# -*- encoding: utf-8 -*-

from deproject import *

if __name__ == "__main__":
    plt.ion()

    (shape, shapeproj) = loadshapefile()
    (lon, lat, x, y) = loaddata('gcp29.points', wgs84=True)
    imname = "TheSteppe.jpg"

    def searchsolution2xy(proj,
                          parametersString,
                          init,
                          lon=lon,
                          lat=lat,
                          x=x,
                          y=y,
                          order=2,
                          plot=True,
                          description="",
                          shape=shape,
                          inproj=shapeproj):
        vec2dictfunc = make_vector2dictfunc(parametersString)
        sol, (xout, yout), p, x2t, t2x, rest = search(lon, lat, x, y, proj, vec2dictfunc, init,
                                                      order)
        if plot:
            descriptor = "%s%s projection (%s, %s, fit error %.3f)" % (description, proj,
                                                                       parametersString, 'poly1'
                                                                       if order == 1 else 'poly2',
                                                                       sol[1])
            image_show(
                x,
                y,
                xout,
                yout,
                description=descriptor,
                imname=imname,
                shape=shape,
                inproj=inproj,
                outproj=p,
                x2t=x2t)
        return p, x2t, t2x, rest

    pixToLonlat2 = lambda x, y, t2x, p: p(*t2x([x, y]), inverse=True)
    ll2pix2 = lambda lon, lat, x2t, p: x2t(p(lon, lat))
    earthRadius = 6378137
    mPerDeg = np.pi / 180 * earthRadius

    # Albers equal-area, 4-parameter, poly2 (quadratic)
    p, x2t, t2x, *_ = searchsolution2xy('aea', "lon_0,lat_0,lat_1,lat_2", [80.0, 50, 40, 60])
    plt.savefig('aea.png', dpi=200)

    # TWO-parameter Winkel Tripel, poly1 (affine)
    p, x2t, t2x, *_ = searchsolution2xy('wintri', "lon_0,lat_1", [47., 0.], order=1)
    print("Two-parameter Winkel Tripel, affine, SRS: ", p.srs)
    plt.savefig('wintri-lon_0-lat_1.png', dpi=200)

    # Same as above, but poly2
    p, x2t, t2x, *_ = searchsolution2xy('wintri', "lon_0,lat_1", [47., 0.], order=2)
    plt.savefig('wintri-lon_0-lat_1-quadratic.png', dpi=200)

    # Load image
    import pylab as plt
    im = plt.imread(imname)
    height, width = im.shape[:2]
    if not (height == 1058 and width == 1600):
        print(
            "Geo-control poins (GCPs) expect a 1600x1058 image but TheSteppe.jpg is not that size")

    # Manual interpolation to equirectangular projection.
    print('Equirectangular interpolation started, might take two minutes…', end='')
    # res, outLon, outLat = manualinterpolate(im, t2x, p, degPerPix=0.05, fname='outfine.png')
    print(' done!')
