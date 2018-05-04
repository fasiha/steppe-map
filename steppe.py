# -*- encoding: utf-8 -*-

from pyproj import Proj
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import scipy.optimize as opt


def loaddata(fname="gcp.txt"):
    """Load QGIS-style georeference control points (GCP) file.

    The first line will be skipped, and the rest of the lines are assumed to be
    comma-separated numbers.

    Four vectors will be returned:
    - lon, degrees,
    - lat, degrees,
    - x, pixel, and
    - y (pixel).

    """
    arr = np.genfromtxt(fname, skip_header=1, delimiter=',')[:, :4]
    lon, lat, x, y = arr.T
    return (lon, lat, x, y)


def L2_norm_error(true, estimated):
    """Compute the normalized L^2 error between two values.

    Compute the L^2 error between two numeric inputs, one of them considered
    "true" (say, x) and the other an "approximation" (say, y), and normalize
    this error using the L^2 norm of x:

    sqrt(|x_1 - y_1|^2 + .. + |x_N - y_N|^2) / sqrt(|x_1|^2 + .. + |x_N|^2).

    The summation is taken over all elements of the inputs: arrays are
    effectively raveled (rasterized).

    To ensure that two numerics that should be "equal" are indeed equal, verify
    that this function's output is near machine precision (see `numpy.spacing`).

    """
    L2 = lambda x: np.sqrt(np.sum(np.abs(x)**2))
    return L2(true - estimated) / L2(true)


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


def remove_linear(x, xp):
    """Find and remove a linear translation between two vectors.

    Given some "good" data in vector x, and an approximation in xp ("p" for
    "prime"), first find the best-fit (in the least-squares sense) slope, m, and
    intercept, b, of

    x = m * xp + b,

    and then return (m*xp + b).

    In linear algebra terms, in compact Matlab/Octave notation, this is (where
    the `(:)` operation is equivalent to `.ravel()` and `\` is the left matrix
    divide):

    ```
    A = [xp(:), ones(size(xp(:)))]
    return A * (A \ x(:))
    ```

    That last line can be replaced with `return A * pinv(A' * A) * A' * x(:)`
    where `'` indicates transpose. Note the hat matrix (or projection matrix) in
    this expression: we're just doing ordinary least squares (OLS).

    This is useful when x and xp are in different units (Celcius versus
    Farenheit, or pixel locations in different images with different origins).
    This function will try to convert xp into the units of x, assuming the units
    are linearly related.

    If the inputs are arrays, the function is called recursively on each axis,
    so that the above operation is done only to vectors. For example, if you
    pass in two 2xN arrays, the output will apply the above operation to the
    first 1xN vector and then the second, finding different slopes/intercepts
    for each row.

    NB: the order of arguments to this function is important! If you have some
    good data x and some scaled approximation xp, you would only expect coherent
    results if you compared the error between x and `remove_linear(x, xp)`.
    Flipping the arguments to remove_linear will give you the wrong answer.

    """
    if x.shape != xp.shape:
        raise ValueError("inputs are not same dimension")

    if x.ndim == 1:
        A = np.vstack([xp, np.ones_like(xp)]).T  # Data matrix for least-squares
        slope_intercept = la.lstsq(A, x)[0]  # lstsq returns other stuff too
        newxp = np.dot(A, slope_intercept)
        return newxp
    else:
        return np.hstack(map(remove_linear, x, xp)).reshape(xp.shape)


def search(lon, lat, x, y, proj, vec2dictfunc, init):
    xy = np.vstack([x, y])

    def minimize(inputvec):
        p = Proj(proj=proj, **vec2dictfunc(inputvec))
        xout, yout = p(lon, lat)
        xout, yout = remove_affine(xy, np.vstack([xout, yout]))[0]
        # return L2_norm_error(np.vstack([x, y]), np.vstack([xout, yout]))
        return np.sum((np.vstack([x, y]) - np.vstack([xout, yout]))**2)

    sols = []
    kws = dict(full_output=True, disp=True, xtol=1e-9, ftol=1e-9, maxiter=10000, maxfun=20000)
    sols.append(opt.fmin(minimize, init, **kws))
    sols.append(opt.fmin_powell(minimize, init, **kws))
    sols.append(opt.fmin_bfgs(minimize, init, full_output=True, disp=True))
    sols.append(opt.fmin_cg(minimize, init, gtol=1e-8, full_output=True, disp=True))
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
                    'disp': True
                })))

    (idx_x, idx_fx) = (0, 1)
    best = np.argmin(map(lambda x: x[idx_fx], sols))
    return sols[best]


def make_vector2dictfunc(string, delimiter=',', initvec=None):
    substrings = string.split(delimiter)
    if initvec and (len(initvec) != len(substrings)):
        print("string=[{}] not split {}-ways".format(string, len(initvec)))
        return None
    return lambda x: dict(zip(substrings, np.atleast_1d(x)))


def grid1dsearch(lon, lat, x, y, proj='moll', plot=True):
    xy = np.vstack([x, y])
    grid = np.arange(-180, 180.0, .25)
    errs = np.empty_like(grid)
    for (idx, l) in enumerate(grid):
        p = Proj(proj=proj, lon_0=l)
        xout, yout = p(lon, lat)
        #xout, yout = remove_linear(xy, np.vstack([xout, yout]))
        xout, yout = remove_affine(xy, np.vstack([xout, yout]))[0]
        errs[idx] = L2_norm_error(np.vstack([x, y]), np.vstack([xout, yout]))

    if plot:
        try:
            import pylab as plt
            plt.ion()
            plt.figure()
            plt.plot(grid, errs)
            title = "%s projection, min=%g at %g degrees" % (proj, np.min(errs),
                                                             grid[np.argmin(errs)])
            plt.title(title)
        except ImportError:
            print("Couldn't plot!")

    bestidx = np.argmin(errs)
    p = Proj(proj=proj, lon_0=grid[bestidx])
    xout, yout = p(lon, lat)
    xout, yout = remove_affine(xy, np.vstack([xout, yout]))[0]

    return (grid, errs, xout, yout)


def proj1d(recompute=True):
    """List of 1D-only projections from [1]. Returns list of strings.

    Ones here but not in pseudocylindrical list: ortho

    In pseudocylindrical but >1D so omitted here: loxim.

    DON'T work with pyproj but are in [1]: hataea, quau, dense, parab.

    [1] ftp://ftp.remotesensing.org/proj/OF90-284.pdf

    """

    if recompute is False:
        return 'sinu,moll,robin,eck1,eck2,eck3,eck4,eck5,eck6,goode,mbtfpp,mbtfpq,mbtfps,putp2,putp5,wink1,boggs,collg,ortho'.split(
            ',')

    pall = 'sinu,moll,robin,eck1,eck2,eck3,eck4,eck5,eck6,goode,hataea,mbtfpp,mbtfpq,mbtfps,putp2,putp5,quau,wink1,boggs,collg,dense,parab,ortho'.split(
        ',')

    (lon, lat, x, y) = loaddata()
    p = []
    for proj in pall:
        # try:
        (grid, errs, _, _) = grid1dsearch(lon, lat, x, y, proj=proj, plot=False)
        print("{}: min={} @ {} degrees".format(proj, np.min(errs), grid[np.argmin(errs)]))
        p.append(proj)
        # except:
        #     print ("{} didn't work".format(proj)
    print(",".join(p))
    return p


def loadshapefile():
    import os.path
    import shapefile
    import pyproj

    # countriespath = os.path.join('ne', 'ne_10m_admin_0_countries', 'ne_10m_admin_0_countries')
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


def shape2pixels(inproj, outproj, shape, Ahat, that):
    import pyproj
    shape = pyproj.transform(inproj, outproj, *shape)

    #countries = [c[3] for c in sf.records()]
    #mongolia = np.array(sf.shape(countries.index('Mongolia')).points).T
    #mongolia = pyproj.transform(inProjection, outproj, *mongolia)

    xout, yout = np.dot(Ahat, shape) + that

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
               Ahat=None,
               that=None,
               **shapeargs):
    """Load and show an image with control points and estimates.

    Given control points in vectors x and y containing pixels, as well as
    estimates obtained (using some projection), plot both values so they can be
    visually compared.

    """
    import pylab as plt
    plt.ion()
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
        (shapex, shapey) = shape2pixels(inproj, outproj, shape, Ahat, that)
        plt.plot(shapex, -shapey, marker='.', markersize=1.0, linestyle='none', **shapeargs)

    ax.axis(imaxis)
    plt.show()

    return ax


def searchsolution2xy(lon,
                      lat,
                      x,
                      y,
                      proj,
                      vec2dictfunc,
                      solvec,
                      plot=True,
                      description="",
                      shape=None,
                      inproj=None):
    xy = np.vstack([x, y])
    solutionvec = solvec[0]
    solutionerr = solvec[1]
    p = Proj(proj=proj, **vec2dictfunc(solutionvec))
    xout, yout = p(lon, lat)
    ((xout, yout), _, Ahat, that) = remove_affine(xy, np.vstack([xout, yout]))

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
            Ahat=Ahat,
            that=that)

    return (xout, yout, p, Ahat, that, ax)


def listToParams(params):
    A0, A1, A2, A3, b0, b1, lon = params
    A = np.array([A0, A1, A2, A3]).reshape(2, 2) * 1e-6
    b = np.vstack([b0, b1]) * 1e2
    p = Proj(proj='wintri', lon_0=lon)
    return A, b, p


def paramsToInit(A, b, pLonInit):
    Ainit = A * 1e6
    binit = b.ravel() * 1e-2
    return Ainit.ravel().tolist() + [binit[0], binit[1], pLonInit]


def fine(lonsLocs, latsLocs, Ainit, binit, pLonInit, sse=True):
    """Fine-tune the estimate using lat-only & lon-only edge ticks"""
    pixToLonlat = lambda x, y, A, b, p: p(*np.linalg.solve(A, np.vstack([x, y]) - b), inverse=True)

    def minimize(params):
        A, b, p = listToParams(params)
        # This is re-factoring `A` for every tick, FIXME
        _, latsHat = pixToLonlat(latsLocs['px'], latsLocs['py'], A, b, p)
        lonsHat, _ = pixToLonlat(lonsLocs['px'], lonsLocs['py'], A, b, p)
        if sse:
            return np.sum((lonsHat - lonsLocs['deg'])**2) + np.sum((latsHat - latsLocs['deg'])**2)
        return np.max(
            [np.max((lonsHat - lonsLocs['deg'])**2),
             np.max((latsHat - latsLocs['deg'])**2)])

    init = paramsToInit(Ainit, binit, pLonInit)
    sols = []
    kws = dict(full_output=True, disp=True, xtol=1e-9, ftol=1e-9, maxiter=10000, maxfun=20000)
    sols.append(opt.fmin(minimize, init, **kws))
    sols.append(opt.fmin_powell(minimize, init, **kws))
    sols.append(opt.fmin_bfgs(minimize, init, full_output=True, disp=True))
    sols.append(opt.fmin_cg(minimize, init, gtol=1e-8, maxiter=10000, full_output=True, disp=True))
    fixmin = lambda res: [res['x'], res['fun']]
    bounds = (np.array(init) + 3.5 * np.vstack([-1, 1.])).T.tolist()
    if sse == False:
        sols.append(
            fixmin(
                opt.differential_evolution(
                    minimize,
                    bounds,
                    popsize=20,
                    mutation=(0.4, 1.1),
                    recombination=0.5,
                    callback=lambda x, convergence: print('ga cb', x, convergence),
                    disp=True,
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
                    'disp': True
                })))

    bestidx = np.argmin([x[1] for x in sols])
    return sols[bestidx], listToParams(sols[bestidx][0])


def fineLoad(fname='ticks.points'):
    lon, lat, x, y = loaddata(fname)
    latTicks = {'deg': lat[lon == 0], 'px': x[lon == 0], 'py': y[lon == 0]}
    lonTicks = {'deg': lon[lat == 0], 'px': x[lat == 0], 'py': y[lat == 0]}
    return lonTicks, latTicks


def manualinterpolate(im, A, b, p, outLon=None, outLat=None, degPerPix=0.05):
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
    return res, outLon, outLat


if __name__ == "__main__":
    (shape, shapeproj) = loadshapefile()
    (lon, lat, x, y) = loaddata()

    fit_description = ''

    # This is how you fit 4 parameters
    fit_proj = 'aea'
    fit_v2dfunc = make_vector2dictfunc("lon_0,lat_0,lat_1,lat_2")
    fit_init = [80.0, 50, 40, 60]
    searchsolution2xy(
        lon,
        lat,
        x,
        y,
        fit_proj,
        fit_v2dfunc,
        search(lon, lat, x, y, fit_proj, fit_v2dfunc, fit_init),
        description=fit_description,
        shape=shape,
        inproj=shapeproj)

    # I believe this is the projection though: 1-parameter Winkel Triple.
    fit_proj = 'wintri'
    fit_v2dfunc = make_vector2dictfunc("lon_0")
    fit_init = [47.0]
    xout, yout, p, Ahat, that, ax = searchsolution2xy(
        lon,
        lat,
        x,
        y,
        fit_proj,
        fit_v2dfunc,
        search(lon, lat, x, y, fit_proj, fit_v2dfunc, fit_init),
        description=fit_description,
        shape=shape,
        inproj=shapeproj)
    print(p.srs)

    recoveredPixels = Ahat @ p(lon, lat) + that
    origPixels = np.vstack([x, y])
    absoluteError = origPixels - recoveredPixels
    relativeError = absoluteError / origPixels

    pixToLonlat0 = lambda x, y: p(*np.linalg.solve(Ahat, np.vstack([x, y]) - that), inverse=True)
    (top_left_lon, top_left_lat) = pixToLonlat0(0, 0)
    import pylab as plt
    im = plt.imread('TheSteppe.jpg')
    height, width = im.shape[:2]
    (bottom_right_lon, bottom_right_lat) = pixToLonlat0(width, -height)
    cmd = ("gdal_translate -of GTiff -a_ullr {top_left_lon} {top_left_lat} {bottom_right_lon}" +
           " {bottom_right_lat} -a_srs SR-ORG:7291 TheSteppe.jpg output.tif").format(
               top_left_lon=top_left_lon[0],
               top_left_lat=top_left_lat[0],
               bottom_right_lon=bottom_right_lon[0],
               bottom_right_lat=bottom_right_lat[0])
    """
    gdal_translate -of GTiff -a_ullr -3.5083634007813402 70.372117747633 131.7449661509387 3.5400212381456004 -a_srs SR-ORG:7291 TheSteppe.jpg output.tif
    """
    # This probably won't work because the Winkel Triple projection isn't widely supported.
    # Fine. We can do interpolations ourselves!
    earthRadius = 6378137
    mPerDeg = np.pi / 180 * earthRadius

    # res, outLon, outLat = manualinterpolate(im, Ahat, that, p, degPerPix=0.05)
    # plt.imsave(fname='out.png', arr=res[::-1, :, :])
    # print(("out.png saved, equirectangular (Plate Carree) projection, with corners: " +
    #        "top_left_lon={top_left_lon} top_left_lat={top_left_lat} " +
    #        "bottom_right_lon={bottom_right_lon} bottom_right_lat={bottom_right_lat} deg").format(
    #            top_left_lon=outLon[0, 0],
    #            top_left_lat=outLat[-1, -1],
    #            bottom_right_lon=outLon[-1, -1],
    #            bottom_right_lat=outLat[0, 0]))
    # cmd = ("gdal_translate -of GTiff -a_ullr {top_left_lon} {top_left_lat} {bottom_right_lon}" +
    #        " {bottom_right_lat} -a_srs EPSG:32662 out.png output.tif").format(
    #            top_left_lon=outLon[0, 0] * mPerDeg,
    #            top_left_lat=outLat[-1, -1] * mPerDeg,
    #            bottom_right_lon=outLon[-1, -1] * mPerDeg,
    #            bottom_right_lat=outLat[0, 0] * mPerDeg)
    # print(cmd)

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
    for x, y, d in zip(lonTicks['px'], -lonTicks['py'], lonTicks['deg']):
        annot_helper(plt.gca(), x, y, d, 'g')
    for x, y, d in zip(latTicks['px'], -latTicks['py'], latTicks['deg']):
        annot_helper(plt.gca(), x, y, d, 'g')

    pixToLonlat = lambda x, y, A, b, p: p(*np.linalg.solve(A, np.vstack([x, y]) - b), inverse=True)
    ll2pix = lambda lon, lat, A, b, p: A @ np.vstack(p(lon, lat)) + b

    def solutionToMaxErr(A, b, p):
        e0 = pixToLonlat(lonTicks['px'], lonTicks['py'], A, b, p)[0] - lonTicks['deg']
        e1 = pixToLonlat(latTicks['px'], latTicks['py'], A, b, p)[1] - latTicks['deg']
        return np.max(np.abs(np.hstack([e0, e1])))

    print(solutionToMaxErr(Ahat, that, p))

    bestsol, (Afine, bfine, pfine) = fine(lonTicks, latTicks, Ahat, that,
                                          float(p.srs.split('+lon_0=')[1]))
    print(solutionToMaxErr(Afine, bfine, pfine))

    bestsol2, (Afine2, bfine2, pfine2) = fine(
        lonTicks, latTicks, Afine, bfine, float(pfine.srs.split('+lon_0=')[1]), sse=False)
    print(solutionToMaxErr(Afine2, bfine2, pfine2))

    a3, b3, p3 = Afine2, bfine2, pfine2
    for iter in range(5):
        a3, b3, p3 = fine(
            lonTicks, latTicks, a3, b3, float(p3.srs.split('+lon_0=')[1]), sse=False)[1]
        print(solutionToMaxErr(a3, b3, p3))

    # resfine, lonfine, latfine = manualinterpolate(im, Afine2, bfine2, pfine2, .05)
    # # resfine, lonfine, latfine = manualinterpolate(im, Afine, bfine, pfine, .05)
    # plt.imsave(fname='outfine.png', arr=resfine[::-1, :, :])
    # cmd = ("gdal_translate -of GTiff -a_ullr {top_left_lon} {top_left_lat} {bottom_right_lon}" +
    #        " {bottom_right_lat} -a_srs EPSG:32662 outfine.png outputfine.tif").format(
    #            top_left_lon=lonfine[0, 0] * mPerDeg,
    #            top_left_lat=latfine[-1, -1] * mPerDeg,
    #            bottom_right_lon=lonfine[-1, -1] * mPerDeg,
    #            bottom_right_lat=latfine[0, 0] * mPerDeg)
    # print(cmd)
    # print(("gdal_translate -of GTiff -a_ullr {top_left_lon} {top_left_lat} {bottom_right_lon}" +
    #        " {bottom_right_lat} -a_srs EPSG:4326 outfine.png outputfine2.tif").format(
    #            top_left_lon=lonfine[0, 0],
    #            top_left_lat=latfine[-1, -1],
    #            bottom_right_lon=lonfine[-1, -1],
    #            bottom_right_lat=latfine[0, 0]))

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

    # myim(lonfine[0], latfine[:, 0], resfine)

    height, width, _ = im.shape

    myim(np.arange(width), -np.arange(height), im)
    plt.gca().invert_yaxis()
    plt.title('fine2')
    for l in range(0, 170, 10):
        plt.plot(*ll2pix(np.ones(100) * l, np.linspace(0, 70, 100), a3, b3, p3))
    for l in range(10, 70, 10):
        plt.plot(*ll2pix(np.linspace(0, 160, 100), np.ones(100) * l, a3, b3, p3))

    myim(np.arange(width), -np.arange(height), im)
    plt.gca().invert_yaxis()
    plt.title('fine')
    for l in range(0, 170, 10):
        plt.plot(*ll2pix(np.ones(100) * l, np.linspace(0, 70, 100), Afine, bfine, pfine))
    for l in range(10, 70, 10):
        plt.plot(*ll2pix(np.linspace(0, 160, 100), np.ones(100) * l, Afine, bfine, pfine))

    myim(np.arange(width), -np.arange(height), im)
    plt.gca().invert_yaxis()
    plt.title('Ahat')
    for l in range(0, 170, 10):
        plt.plot(*ll2pix(np.ones(100) * l, np.linspace(0, 70, 100), Ahat, that, p))
    for l in range(10, 70, 10):
        plt.plot(*ll2pix(np.linspace(0, 160, 100), np.ones(100) * l, Ahat, that, p))
