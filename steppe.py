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

def grid1dsearch(lon, lat, x, y, proj='moll'):
    xy = np.vstack([x, y])
    grid = np.arange(0, 360.0, 1.0)
    errs = np.empty_like(grid)
    for (idx, l) in enumerate(grid):
        p = Proj(proj=proj, lon_0=l)
        xout, yout = p(lon, lat)
        xout, yout = remove_linear(xy, np.vstack([xout, yout]))
        errs[idx] = L2_norm_error(np.vstack([x,y]), np.vstack([xout, yout]))
    try:
        import pylab as plt
        plt.ion()
        plt.figure()
        plt.plot(grid, errs)
    except:
        print "Couldn't plot!"
    return (grid, errs)

if __name__ == "__main__":
    (lon, lat, x, y) = loaddata()
    (grid, errs) = grid1dsearch(lon, lat, x, y)
