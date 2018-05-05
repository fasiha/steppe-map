# -*- encoding: utf-8 -*-
"""
For https://pubs.er.usgs.gov/publication/70136641
"""

from deproject import *

if __name__ == "__main__":
    plt.ion()

    imname = 'plate-1-150ppi-preview.png'
    (shape, shapeproj) = loadshapefile()
    (lon, lat, x, y) = loaddata('plate1.points', wgs84=False)

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
    p, x2t, t2x, postfit = searchsolution2xy(
        'aea', "lon_0,lat_0,lat_1,lat_2", [80.0, 50, 40, 30], order=1)
