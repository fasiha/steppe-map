# -*- encoding: utf-8 -*-

from pyproj import Proj, transform
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import scipy.optimize as opt

import pylab as plt
plt.ion()


def wgs84ToDeg(x, y):
    "Converts lat/lon from epsg:3857 (meters) to degrees"
    P3857 = Proj(init='epsg:3857')  # wgs84 in meters
    P4326 = Proj(init='epsg:4326')  # wgs84 in degrees
    return transform(P3857, P4326, x, y)


def loaddata(fname="gcp.txt", wgs84=False):
    """Load QGIS-style georeference control points (GCP) file.

    The first line will be skipped, and the rest of the lines are assumed to be
    comma-separated numbers.

    Four vectors will be returned:
    - lon, degrees,
    - lat, degrees,
    - x, pixel, and
    - y (pixel).

    If `wgs84` is truthy, lat/lon are treated as meters and converted to degrees.
    """
    arr = np.genfromtxt(fname, skip_header=1, delimiter=',')[:, :4]
    lon, lat, x, y = arr.T
    if wgs84:
        lon, lat = wgs84ToDeg(lon, lat)
    return (lon, lat, x, y)


def remove_polynomial2_2d(t, x):
    """
    For `t = [[t1], [t2]]` and `x = [[x1], [x2]]`, find `A` (2 by 5) and `b` (2 by 1) such that
    
    `t = A @ [[x1], [x2], [x1**2], [x2**2], [x*y]] + b`

    in the least-squares sense. Note that `x` and `t` must have two rows but any number of columns
    (2 by N).

    Returns a 3-tuple:

    0. A version of `x` with the quadratic relationship removed, `xnew`
    1. `A`, a 2 by 5 array
    2. `b`, a 2 by 1 array
    3. `x2t`, a function such that `x2t(x) = t` to machine precision.

    If `t` is exactly a quadratic function of `x`, then `xnew=t` to machine precision.
    """
    maxx = np.max(x)
    y = np.vstack([x / maxx, x**2 / maxx**2, x[0] * x[1] / maxx**2, np.ones_like(x[0])])
    Ab = np.linalg.lstsq(y.T, t.T, rcond=None)[0].T
    A = Ab[:, :-1]
    b = Ab[:, -1:]
    x2t = lambda x: (A @ np.vstack([x/maxx, np.array(x/maxx)**2, np.array(x[0]/maxx) * np.array(x[1]/maxx)]) + b)
    return x2t(x), x2t, A, b


def remove_affine(p, q, q_factor=None, skip_factorization=False):
    """Removes an (unknown) affine transform between two matrixes.

    Given two arrays of the same size, `p` and `q`, finds a matrix `A` and
    column vector `t` such that

    `p = A * q + t`

    in the least-squares sense, and then computes `qnew = A * q + t`. (Notation:
    `matrix + vector` implies the vector is added to each column of the matrix.)

    NB: `p` and the returned `qnew` will be equal if and only if `p` is
    generated from `q` via an affine transform (no noise).

    Returns `(qnew, q_factor, Ahat, that)`. `q_factor` is a matrix factorization
    that can greatly speed up subsequent calls to remove_affine *with the same
    `q`*. If your `q` stays the same for multiple calls, cache `q_factor` and
    pass it in as a keyword argument; `q_factor` won't change from call to call.
    However, if your `q` change from call to call, ignore `q_factor` and pass in
    `skip_factorization=False` to avoid even calculating it. `Ahat` and `that`
    are the estimated values of `A` and `t`.

    NB2: the default `q_factor=None` will trigger computation of the
    factorization unless `skip_factorization=False`. Non-`None` `q_factor` will
    be trusted: no checks will be performed to make sure the given `q_factor` is
    indeed generated by the `q` you pass in. (Example: for `q.shape` of (2, 22),
    the speedup from using `q_factor` is 1.4x with skip_factorization=False, and
    1.3x the case with skip_factorization=True, on a 2009 Mac Book Pro.)

    Implements the algorithm described in H. Spath, "Fitting affine and
    orthogonal transformations between two sets of points" in *Mathematical
    Communications*, vol. 9 (2004), pp. 27--34. http://hrcak.srce.hr/file/1425

    """

    qaug = np.vstack([q, np.ones_like(q[0, :])])
    if q_factor is None:
        Q = np.dot(qaug, qaug.T)

        if skip_factorization:
            sol = la.lstsq(Q, np.dot(qaug, p.T))[0]
            q_factor = None

        else:
            q_factor = scila.cho_factor(Q)
            sol = scila.cho_solve(q_factor, np.dot(qaug, p.T))

    else:
        sol = scila.cho_solve(q_factor, np.dot(qaug, p.T))

    # sol.shape is (n+1, n), for n=p.shape[0]
    Ahat = sol[:-1, :].T  # top square matrix of sol, transposed
    that = sol[-1:, :].T  # bottom row vector of sol, transposed
    qnew = np.dot(Ahat, q) + that
    return (qnew, q_factor, Ahat, that)


def search(lon, lat, x, y, proj, vec2dictfunc, init):
    """
    Given

    - `lon`: an array of longitudes in degrees
    - `lat`: an array of the same size as `lon`, giving latitudes in degrees
    - `x`: horizontal pixel locations, with 0 being leftmost and increasing to the right
    - `y`: vertical pixel locations, with 0 being upmost and DECREASING downwards (all should be
        negative!)
    - `proj`: a string of a projection Proj4 recognizes ('aea', 'wintri', etc.)
    - `vec2dictfunc`: a transformer function that, given a tuple of numbers, creates a dict with
        correct keys for the projection (`make_vector2dictfunc` can make these),
    - `init`: a vector of numbers, one for each projection parameter,
    
    runs several nonlinear least squares methods to find the projection's parameters that bet fit
    the data. Returns a tuple of at least two elements:

    1. a vector of best-fit parameters,
    2. the final error,
    3. whatever else the minimizer might return.
    """
    xy = np.vstack([x, y])

    def minimize(inputvec):
        p = Proj(proj=proj, **vec2dictfunc(inputvec))
        xout, yout = p(lon, lat)
        xout, yout = remove_polynomial2_2d(xy, np.vstack([xout, yout]))[0]
        return np.sum((xy - np.vstack([xout, yout]))**2)

    sols = []
    kws = dict(full_output=True, disp=False, xtol=1e-9, ftol=1e-9, maxiter=10000, maxfun=20000)
    sols.append(opt.fmin(minimize, init, **kws))
    sols.append(opt.fmin_powell(minimize, init, **kws))
    sols.append(opt.fmin_bfgs(minimize, init, full_output=True, disp=False))
    sols.append(opt.fmin_cg(minimize, init, gtol=1e-8, full_output=True, disp=False))
    fixmin = lambda res: [res['x'], res['fun']]
    sols.append(
        fixmin(
            opt.minimize(
                minimize,
                init,
                method="Nelder-Mead",
                tol=1e-9,
                options={
                    "xatol": 1e-9,
                    "fatol": 1e-9,
                    "maxfev": 20000,
                    "maxiter": 10000,
                    'disp': False
                })))

    (idx_x, idx_fx) = (0, 1)
    best = np.argmin(map(lambda x: x[idx_fx], sols))
    return sols[best]


def make_vector2dictfunc(string, delimiter=','):
    """
    Given a string of projection parameters, like `"lat_0,lon_1"`, and a string delimiter,
    returns a function that converts a tuple of numbers to a dict. I.e.,
    ```
    (make_vector2dictfunction('lat_0,lon_1'))(1, 2) # {'lat_0' : 1, 'lon_1' : 2}
    ```
    """
    return lambda x: dict(zip(string.split(delimiter), np.atleast_1d(x)))


def loadshapefile():
    import os.path
    import shapefile
    import pyproj

    coastpath = os.path.join('ne', 'ne_10m_coastline', 'ne_10m_coastline')
    shppath = coastpath + '.shp'
    prjpath = coastpath + '.prj'

    try:
        from osgeo import osr

        # From http://gis.stackexchange.com/questions/17341/
        prjText = open(prjpath, 'r').read()
        srs = osr.SpatialReference()
        if (srs.ImportFromWkt(prjText)):
            print("error importing .prj information from ", prjpath)
            return (None, None)
        inProjection = pyproj.Proj(srs.ExportToProj4())
    except ImportError:
        inProjection = pyproj.Proj('+proj=longlat +ellps=WGS84 +no_defs')

    sf = shapefile.Reader(shppath)
    world = np.vstack(
        [shp.points for (rec, shp) in zip(sf.records(), sf.shapes()) if rec[1] <= 1.0]).T

    return (world, inProjection)


def shape2pixels(inproj, outproj, shape, x2t):
    shape = transform(inproj, outproj, *shape)
    xout, yout = x2t(shape)
    return np.array([xout, yout])


def image_show(x,
               y,
               xout,
               yout,
               imname="TheSteppe.jpg",
               description="",
               shape=None,
               inproj=None,
               outproj=None,
               x2t=None,
               **shapeargs):
    """Load and show an image with control points and estimates.

    Given control points in vectors x and y containing pixels, as well as
    estimates obtained (using some projection), plot both values so they can be
    visually compared.
    """
    try:
        im = plt.imread(imname)
    except IOError:
        print("Couldn't load {}, can't display image!".format(imname))
        return

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(im, interpolation='bicubic')
    imaxis = ax.axis()

    def annot_helper(x, y, c='k', **kwargs):
        ax.annotate(
            "%g,%g" % (x, y),
            xy=(x, y),
            xytext=(x + 30, y + 30),
            size=15,
            arrowprops=dict(facecolor=c, width=1.0),
            color=c,
            **kwargs)

    [annot_helper(xi, yi, 'g') for xi, yi in zip(x, y)]
    [annot_helper(xi, yi, 'r') for xi, yi in zip(xout, yout)]

    plt.title(description + " (Green: control points, red: fit points)")

    if shape is not None:
        (shapex, shapey) = shape2pixels(inproj, outproj, shape, x2t)
        plt.plot(shapex, -shapey, marker='.', markersize=1.0, linestyle='none', **shapeargs)

    ax.axis(imaxis)
    plt.show()

    return ax


def listToParams(params, xformer):
    A0, A1, A2, A3, b0, b1, *rest = params
    A = np.array([A0, A1, A2, A3]).reshape(2, 2) * 1e-6
    b = np.vstack([b0, b1]) * 1e2
    p = Proj(proj='wintri', **xformer(rest))
    return A, b, p


def paramsToInit(A, b, projParams):
    Ainit = A * 1e6
    binit = b.ravel() * 1e-2
    return Ainit.ravel().tolist() + binit.tolist() + np.atleast_1d(projParams).tolist()


def fine(lonsLocs, latsLocs, Ainit, binit, projParamsInit, projParamsXform, sse=True):
    """Fine-tune the estimate using lat-only & lon-only edge ticks. Very custom."""
    pixToLonlat = lambda x, y, A, b, p: p(*np.linalg.solve(A, np.vstack([x, y]) - b), inverse=True)

    def minimize(params):
        A, b, p = listToParams(params, projParamsXform)
        # This is re-factoring `A` for every tick, FIXME
        _, latsHat = pixToLonlat(latsLocs['px'], latsLocs['py'], A, b, p)
        lonsHat, _ = pixToLonlat(lonsLocs['px'], lonsLocs['py'], A, b, p)
        if sse:
            return np.sum((lonsHat - lonsLocs['deg'])**2) + np.sum((latsHat - latsLocs['deg'])**2)
        return np.max(
            [np.max(np.abs(lonsHat - lonsLocs['deg'])),
             np.max(np.abs(latsHat - latsLocs['deg']))])

    init = paramsToInit(Ainit, binit, projParamsInit)
    sols = []
    kws = dict(full_output=True, disp=False, xtol=1e-9, ftol=1e-9, maxiter=10000, maxfun=20000)
    sols.append(opt.fmin(minimize, init, **kws))
    sols.append(opt.fmin_powell(minimize, init, **kws))
    sols.append(opt.fmin_bfgs(minimize, init, full_output=True, disp=False))
    sols.append(opt.fmin_cg(minimize, init, gtol=1e-8, maxiter=10000, full_output=True, disp=False))
    fixmin = lambda res: [res['x'], res['fun']]
    bounds = (np.array(init) + 3.5 * np.vstack([-1, 1.])).T.tolist()
    # bounds = (np.array(init) + [20., 9, 9, 20, 5, 5, 1] * np.vstack([-1, 1.])).T.tolist()
    # bounds = (np.array(init) + [.25, .01, .01, .25, 1, 1, .5] * np.vstack([-1, 1.])).T.tolist()
    if sse == False:
        sols.append(
            fixmin(
                opt.differential_evolution(
                    minimize,
                    bounds,
                    # popsize=20,
                    # maxiter=1000,
                    # mutation=(0.4, 1.1),
                    # recombination=0.5,
                    disp=False,
                    # init='random',
                    polish=True)))
    sols.append(
        fixmin(
            opt.minimize(
                minimize,
                init,
                method="Nelder-Mead",
                tol=1e-9,
                options={
                    "xatol": 1e-9,
                    "fatol": 1e-9,
                    "maxfev": 20000,
                    "maxiter": 10000,
                    'disp': False
                })))

    bestidx = np.argmin([x[1] for x in sols])
    return sols[bestidx], listToParams(sols[bestidx][0], projParamsXform)


def fineLoad(fname='ticks.points'):
    lon, lat, x, y = loaddata(fname)
    latTicks = {'deg': lat[lon == 0], 'px': x[lon == 0], 'py': y[lon == 0]}
    lonTicks = {'deg': lon[lat == 0], 'px': x[lat == 0], 'py': y[lat == 0]}
    return lonTicks, latTicks


def manualinterpolate(im, A, b, p, outLon=None, outLat=None, degPerPix=0.05, fname=None):
    afactor = scila.cho_factor(A)
    pixToLonlat = lambda x, y: p(*scila.cho_solve(afactor, np.vstack([x, y]) - b), inverse=True)

    height, width, _ = im.shape
    xs, ys = np.meshgrid(np.arange(width), -np.arange(height))
    origLon, origLat = pixToLonlat(xs.ravel(), ys.ravel())
    vecToBounds = lambda x: np.array([np.min(x), np.max(x)])
    boundLon = vecToBounds(origLon)
    boundLat = vecToBounds(origLat)

    if outLon is None or outLat is None:
        outLon, outLat = np.meshgrid(
            np.arange(boundLon[0] - 0.5, boundLon[1] + 0.5, degPerPix),
            np.arange(boundLat[0] - 0.5, boundLat[1] + 0.5, degPerPix))

    from scipy.interpolate import griddata
    res = np.dstack([
        griddata(
            np.vstack([origLon, origLat]).T,
            im[:, :, i].ravel(),
            np.vstack([outLon.ravel(), outLat.ravel()]).T,
            method='nearest').reshape(outLon.shape) for i in range(3)
    ])
    if fname:
        plt.imsave(fname=fname, arr=res[::-1, :, :])
        earthRadius = 6378137
        mPerDeg = np.pi / 180 * earthRadius
        print(("gdal_translate -of GTiff -a_ullr {top_left_lon} {top_left_lat} {bottom_right_lon}" +
               " {bottom_right_lat} -a_srs EPSG:32662 out.png output.tif").format(
                   top_left_lon=outLon[0, 0] * mPerDeg,
                   top_left_lat=outLat[-1, -1] * mPerDeg,
                   bottom_right_lon=outLon[-1, -1] * mPerDeg,
                   bottom_right_lat=outLat[0, 0] * mPerDeg))

    return res, outLon, outLat


def myim(x, y, *args, **kwargs):
    def extents(f):
        delta = f[1] - f[0]
        return [f[0] - delta / 2, f[-1] + delta / 2]

    fig, ax = plt.subplots()
    im = ax.imshow(
        *args,
        **kwargs,
        aspect='auto',
        interpolation='none',
        extent=extents(x) + extents(y),
        origin='lower')
    return fig, ax, im


if __name__ == "__main__":
    (shape, shapeproj) = loadshapefile()
    (lon, lat, x, y) = loaddata('gcp29.points', True)

    def searchsolution2xy(proj,
                          parametersString,
                          init,
                          lon=lon,
                          lat=lat,
                          x=x,
                          y=y,
                          plot=True,
                          description="",
                          shape=shape,
                          inproj=shapeproj):
        vec2dictfunc = make_vector2dictfunc(parametersString)
        solutionvec, solutionerr, *_ = search(lon, lat, x, y, proj, vec2dictfunc, init)
        p = Proj(proj=proj, **vec2dictfunc(solutionvec))
        xout, yout = p(lon, lat)
        xy = np.vstack([x, y])
        ((xout, yout), x2t, Ahat, that) = remove_polynomial2_2d(xy, np.vstack([xout, yout]))
        ax = None
        if plot:
            descriptor = "%s%s projection (fit error %.3f)" % (description, proj, solutionerr)
            ax = image_show(
                x,
                -y,
                xout,
                -yout,
                description=descriptor,
                shape=shape,
                inproj=inproj,
                outproj=p,
                x2t=x2t)
        return (xout, yout, p, x2t, Ahat, that, ax)

    pixToLonlat = lambda x, y, A, b, p: p(*np.linalg.solve(A, np.vstack([x, y]) - b), inverse=True)
    ll2pix = lambda lon, lat, A, b, p: A @ np.vstack(p(lon, lat)) + b
    ll2pix2 = lambda lon, lat, x2t, p: x2t(p(lon, lat))
    earthRadius = 6378137
    mPerDeg = np.pi / 180 * earthRadius

    # Some random thing
    # search(lon, lat, x, y, 'wintri', make_vector2dictfunc('lon_0'), [47.])

    searchsolution2xy('aea', "lon_0,lat_0,lat_1,lat_2", [80.0, 50, 40, 60])
    # TWO-parameter Winkel Tripel
    _, _, p, *_ = searchsolution2xy('wintri', "lon_0,lat_1", [47., 0.])
    print("Two-parameter Winkel Tripel SRS: ", p.srs)
    # 1-parameter Winkel Tripel
    _, _, p, *_ = searchsolution2xy('wintri', "lon_0", [47., 0.])
    print("One-parameter Winkel Tripel SRS: ", p.srs)

    ### CUSTOMIZE ME!!! ###
    ### Decide what you want to try to fit
    fit_proj = 'wintri'

    srsParams = 'lon_0,lat_1'
    fit_init = [47., 0]

    _, _, p, x2t, *_ = searchsolution2xy(fit_proj, srsParams, fit_init)

    # Load image
    import pylab as plt
    im = plt.imread('TheSteppe.jpg')
    height, width = im.shape[:2]
    if not (height == 1058 and width == 1600):
        print(
            "Geo-control poins (GCPs) expect a 1600x1058 image but TheSteppe.jpg is not that size")

    def drawImWithGraticules(x2t, p, title, im=im):
        height, width = im.shape[:2]
        myim(np.arange(width), -np.arange(height), im)
        plt.gca().invert_yaxis()
        plt.title(title)

        (shapex, shapey) = shape2pixels(shapeproj, p, shape, x2t)
        plt.plot(shapex, shapey, marker='.', markersize=1.0, linestyle='none')

        for l in range(0, 170, 10):
            plt.plot(*ll2pix2(np.ones(100) * l, np.linspace(0, 70, 100), x2t, p), 'r', lw=.5)
        for l in range(10, 70, 10):
            plt.plot(*ll2pix2(np.linspace(0, 160, 100), np.ones(100) * l, x2t, p), 'r', lw=.5)
        plt.xlim([0, width])
        plt.ylim([-height, 0])

    drawImWithGraticules(x2t, p, 'Optimize GCPs only, blue dots=GCPs')
    s = plt.scatter(x, y, c='b')

if False:
    # Fine estimation with ticks?
    lonTicks, latTicks = fineLoad('ticks.points')

    plt.figure()
    plt.imshow(im)
    plt.scatter(lonTicks['px'], -lonTicks['py'])
    plt.scatter(latTicks['px'], -latTicks['py'])
    annot_helper = lambda ax, x, y, d, c='k', **kwargs: ax.annotate(
        "{}".format(d),
        xy=(x, y),
        xytext=(x - 50, y + 50),
        size=7,
        # arrowprops=dict(facecolor=c, width=1.0, frac=0.5),
        color=c,
        **kwargs)
    for xx, yy, d in zip(lonTicks['px'], -lonTicks['py'], lonTicks['deg']):
        annot_helper(plt.gca(), xx, yy, d, 'g')
    for xx, yy, d in zip(latTicks['px'], -latTicks['py'], latTicks['deg']):
        annot_helper(plt.gca(), xx, yy, d, 'g')

    def solutionToMaxErr(A, b, p):
        e0 = pixToLonlat(lonTicks['px'], lonTicks['py'], A, b, p)[0] - lonTicks['deg']
        e1 = pixToLonlat(latTicks['px'], latTicks['py'], A, b, p)[1] - latTicks['deg']
        return np.max(np.abs(np.hstack([e0, e1])))

    print("Optmized GCPs, worst-case error, in degrees", solutionToMaxErr(Ahat, that, p))

    searchSrs = lambda srs, key: [float(s.split('=')[1]) for s in srs.split(' ') if key in s][0]
    xform = make_vector2dictfunc(srsParams)
    srsToInit = lambda srs, keys: [searchSrs(p.srs, k) for k in keys.split(',')]

    bestsol, (Afine, bfine, pfine) = fine(lonTicks, latTicks, Ahat, that, srsToInit(
        p.srs, srsParams), xform)
    print("Optimized side-ticks, SSE, worst-case error in deg: ",
          solutionToMaxErr(Afine, bfine, pfine))

    bestsol2, (Afine2, bfine2, pfine2) = fine(
        lonTicks, latTicks, Afine, bfine, srsToInit(pfine.srs, srsParams), xform, sse=False)
    print("Optimized SSE+L0, worst-case error in deg ", solutionToMaxErr(Afine2, bfine2, pfine2))

    drawImWithGraticules(Afine2, bfine2, pfine2, 'Optimize side-ticks, SSE+L0')
    drawImWithGraticules(Afine, bfine, pfine, 'Optimize side-ticks, just SSE')

    # Manual interpolation to equirectangular projection. Commented because it doesn't work that great.
    # res, outLon, outLat = manualinterpolate(im, Afine2, bfine2, pfine2, degPerPix=0.05, fname='outfine.png')

    input('All done, hit enter when done')

    fullLat = dict(
        deg=np.hstack([latTicks['deg'], lat]),
        px=np.hstack([latTicks['px'], x]),
        py=np.hstack([latTicks['py'], y]))
    fullLon = dict(
        deg=np.hstack([lonTicks['deg'], lon]),
        px=np.hstack([lonTicks['px'], x]),
        py=np.hstack([lonTicks['py'], y]))

    solfull, (Afull, bfull, pfull) = fine(fullLon, fullLat, Ahat, that, srsToInit(p.srs, srsParams),
                                          xform)
    print("Optimized side-ticks, SSE, worst-case error in deg: ",
          solutionToMaxErr(Afull, bfull, pfull))

    solfull2, (Afull2, bfull2, pfull2) = fine(
        fullLon, fullLat, Afull, bfull, srsToInit(pfull.srs, srsParams), xform, sse=False)
    print("Optimized side-ticks, SSE, worst-case error in deg: ",
          solutionToMaxErr(Afull2, bfull2, pfull2))

    drawImWithGraticules(Afull2, bfull2, pfull2, 'Optimize GCP+side-ticks, SSE+L0')
    drawImWithGraticules(Afull, bfull, pfull, 'Optimize GCP+side-ticks, just SSE')
