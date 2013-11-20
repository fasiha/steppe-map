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

def remove_affine(p, q, q_factor=None, skip_factorization=False):
    """Removes an (unknown) affine transform between two matrixes.

    Given two arrays of the same size, `p` and `q`, finds a matrix `A` and
    column vector `t`` such that

    `p \approx A * q + t`

    in the least-squares sense, and then computes `qnew = A * q + t`.

    Usage 1 of 3: simple: `qnew = remove_affine(p, q)[0]` returns `qnew`. Note
    how the function returns a tuple and the only the first element is kept.
    Understanding the other return value leads us to...

    Usage 2 of 3: complicated: *if* you plan on calling this function repeatedly
    using the same `q`, you can greatly speeed it up by caching the second
    element of the returned tuple and providing it to this function.

    ```
    # set up p and q
    qnew, q_factor = remove_affine(p, q)

    for i in range(10):
        # change p, keep q the same
        p += 1.0
        qnew, q_factor = remove_affine(p, q, q_factor)
    ```

    q_factor is guaranteed *not* to change after the first call to remove_affine
    so you are free to discard q_factor after the first time it's generated. A
    `q_factor` of `None` automatically causes you to regenerate `q_factor`.

    The first usage example will have the same runtime as the second, because in
    the former, the factorization will be computed, and you're just discarding
    it. There is a third option...

    Usage 3 of 3: simple: `qnew = remove_affine(p, q,
    skip_factorization=True)[0]`. The `skip_factorization` flag will use a
    non-factorizing solver and will return `q_factor` of `None`.

    By skipping factorization and directly solving, #3 will be faster than #2
    the first time you call it. But, assuming `q` doesn't change, subsequent
    calls will be much faster using #2 than #3.

    #1 will be as fast as #2 the first invokation and much slower on subsequent
    ones.

    Suggestions: if your `q` stays the same, use #2. Otherwise, use #3. if you
    don't care, use #1.

    NB: `p` and the returned `qnew` will be equal if and only if `p` and `q` are
    generated as above.

    """

    if q_factor is None:
        qaug = np.vstack([q, np.ones_like(q[0,:])])
        Q = np.dot(qaug, qaug.T)

        if skip_factorization:
            sol = la.lstsq(Q, np.dot(qaug, p.T))
            q_factor = None

        else:
            import scipy.linalg as scila
            q_factor = scila.cho_factor(Q)
            sol = scila.cho_solve(q_factor, np.dot(qaug, p.T))

    else:
        sol = scila.cho_solve(q_factor, p.T)

    Ahat = sol[:-1,:].T
    that = sol[-1,:]
    qnew = np.dot(Ahat, q) + that[:,np.newaxis]
    return (qnew, q_factor)


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
        #xout, yout = remove_linear(xy, np.vstack([xout, yout]))
        xout, yout = remove_affine(xy, np.vstack([xout, yout]))[0]
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
        except ImportError:
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
    (grid, errs, xout, yout) = grid1dsearch(lon, lat, x, y, 'moll', True)
    image_show(x,-y,xout,-yout)
