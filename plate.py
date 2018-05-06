# -*- encoding: utf-8 -*-
"""
For https://pubs.er.usgs.gov/publication/70136641
"""

from deproject import *

if __name__ == "__main__":

    imname = 'plate-1-150ppi-preview.png'
    (shape, shapeproj) = loadshapefile()
    (lon, lat, x, y) = loaddata('plate1.points', wgs84=False)

    import pylab as plt
    plt.ion()
    im = plt.imread(imname)
    height, width = im.shape[:2]

    def searchsolution2xy(proj,
                          parametersString,
                          init,
                          im=im,
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
                im,
                description=descriptor,
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
        'aea', "lon_0,lat_0,lat_1,lat_2", [-80.0, 50, 30, 40], order=1)

    # Can we use GDAL to do the reprojection? Doesn't seem like it >.<
    tl = np.hstack(p(*pixToLonlat2([0], [0], t2x, p)))
    br = np.hstack(p(*pixToLonlat2([width - 1], [-(height - 1)], t2x, p)))
    cmd = ('gdal_translate -a_srs "{srs}" -of GTiff -a_ullr {top_left_lon} {top_left_lat} ' +
           ' {bottom_right_lon} {bottom_right_lat} {imname} image0.tif' +
           ' && gdalwarp -of GTiff -t_srs EPSG:3857 image0.tif image1.tif').format(
               srs=p.srs,
               imname=imname,
               top_left_lon=tl[0],
               top_left_lat=tl[1],
               bottom_right_lon=br[0],
               bottom_right_lat=br[1])
    print(cmd)

    # Manual interpolation always has our back
    res, *_ = manualinterpolate(im, t2x, x2t, p, fname='output.png')
