from pyproj import Proj
import numpy as np
import numpy.linalg as la


p = Proj(proj='moll', lon_0=-90.0)
lon, lat = (-120.108, 34.36116666)
x, y = p(lon, lat)
lon1, lat1 = p(x, y, inverse=True)

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
    arr = np.genfromtxt(fname, skip_header=1, delimiter=',')[:,:4]
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

    To ensure that two numerics that should be "equal" are indeed equal,
    verifying that this function's' output is near machine precision (see
    `numpy.spacing`).

    """
    L2 = lambda x: np.sqrt(np.sum(np.abs(x)**2))
    return L2(true - estimated) / L2(true)

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
        A = np.vstack([xp, np.ones_like(xp)]).T # Data matrix for least-squares
        slope_intercept = la.lstsq(A, x)[0]     # lstsq returns other stuff too
        newxp = np.dot(A, slope_intercept)
        return newxp
    else:
        return np.hstack(map(remove_linear, x, xp)).reshape(xp.shape)

def grid1dsearch(lon, lat, x, y, proj='moll', plot=True):
    xy = np.vstack([x, y])
    grid = np.arange(-180, 180.0, 1.0)
    errs = np.empty_like(grid)
    for (idx, l) in enumerate(grid):
        p = Proj(proj=proj, lon_0=l)
        xout, yout = p(lon, lat)
        xout, yout = remove_linear(xy, np.vstack([xout, yout]))
        errs[idx] = L2_norm_error(np.vstack([x,y]), np.vstack([xout, yout]))
    if plot:
        try:
            import pylab as plt
            plt.ion()
            plt.figure()
            plt.plot(grid, errs)
            title = "%s projection, min=%g at %g degrees" % (proj, np.min(errs),
                                                        grid[np.argmin(errs)])
            plt.title(title)
        except:
            print "Couldn't plot!"
    return (grid, errs, xout, yout)

def proj1d(recompute=True):
    """List of 1D-only projections from [1]. Returns list of strings.

    Ones here but not in pseudocylindrical list: ortho

    In pseudocylindrical but >1D so omitted here: loxim.

    DON'T work with pyproj but are in [1]: hataea, quau, dense, parab.

    [1] ftp://ftp.remotesensing.org/proj/OF90-284.pdf

    """

    if recompute is False:
        return 'sinu,moll,robin,eck1,eck2,eck3,eck4,eck5,eck6,goode,mbtfpp,mbtfpq,mbtfps,putp2,putp5,wink1,boggs,collg,ortho'.split(',')

    pall = 'sinu,moll,robin,eck1,eck2,eck3,eck4,eck5,eck6,goode,hataea,mbtfpp,mbtfpq,mbtfps,putp2,putp5,quau,wink1,boggs,collg,dense,parab,ortho'.split(',')

    (lon, lat, x, y) = loaddata()
    p = []
    for proj in pall:
        try:
            (grid, errs, _, _) = grid1dsearch(lon, lat, x, y,
                                              proj=proj, plot=False)
            print "%s: min=%g @ %g degrees" % (proj, np.min(errs),
                                               grid[np.argmin(errs)])
            p.append(proj)
        except:
            print "%s didn't work" % (proj,)
    print p.join(',')
    return p

def image_show(x, y, xout, yout, imname="small.gif"):
    """Load and show an image with control points and estimates.

    Given control points in vectors x and y containing pixels, as well as
    estimates obtained (using some projection), plot both values so they can be
    visually compared.

    """
    import pylab as plt
    plt.ion()
    im = plt.imread(imname)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(im)

    annot_helper = lambda x, y, **kwargs: ax.annotate("%g,%g"%(x, y),
                                            xy=(x, y), xytext=(x+10, y+10),
                                            size=15,
                                            arrowprops=dict(facecolor='black',
                                                            shrink=0.05),
                                            **kwargs)
    [annot_helper(xi, yi, color='g') for xi,yi in zip(x, y)]
    [annot_helper(xi, yi, color='r') for xi,yi in zip(xout, yout)]

    plt.title("Green: control points. Red: fit points.")



if __name__ == "__main__":
    (lon, lat, x, y) = loaddata()
    (grid, errs, xout, yout) = grid1dsearch(lon, lat, x, y, 'vandg', False)
    image_show(x,-y,xout,-yout)
