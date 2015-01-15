The Steppe, or How to fit arbitrary projections to data
=======================================================

![The Steppe, Britannica map](http://i.stack.imgur.com/HkApE.gif)

Introduction
------------

In the current Encyclopedia Britannica's online article on "Steppe, the", there is a fabulous map --- see [http://www.britannica.com/bps/media-view/3658/1/0/0](http://www.britannica.com/bps/media-view/3658/1/0/0) --- that I would love to see broken out of its static fetters and dance with different overlays (e.g., temperature, rainfall, etc.). Choosing to ignore copyright concerns for the sake of narrow educational pursuits, this project is specifically focused on determining the projection used by this map, in order to get as exact georeferencing as possible, i.e., to accurately convert the green pixels to longitude/latitudes.

Hence the subtitle of the project: how can one go from georeferenced control points (GCPs) to a projection's parameters?

Currently, the project leverages Pyproj and allows me to specify

- a projection (e.g., "aea" or the Albers equal-area),
- its parameters (e.g., lon/lat of false origin, and two standard parallels), and
- an initial numeric guess for these parameters

and after running a nonlinear least squares (provided by Scipy), can plot the image with the original and best-fit GCPs, as well as a coastline, courtesy of Natural Earth.

Requirements
------------

Clearly this tiny project stands on the shoulders of giants. Requirements include

- Numpy, for arrays, and Scipy, for function minimization
- Pyproj for all the projections. If it's not in Pyproj, I can't fit it (yet).

To see the original image on top of the GCPs, one needs the original image: it's linked from [http://www.britannica.com/bps/media-view/3658/1/0/0](http://www.britannica.com/bps/media-view/3658/1/0/0) so save it to `TheSteppe.jpg`.

To see the coastline, one needs

- Natural Earth, specifically, the `ne_10m_coastline` physical dataset (the code looks for files `./ne/ne_10m_coastline/ne_10m_coastline.*`)
- pyshp to read Natural Earth shapefiles: [https://code.google.com/p/pyshp/](https://code.google.com/p/pyshp/)
- GDAL's Python module (`osgeo`) to convert Natural Earth's projection descriptor string (WKT) to a Proj4 descriptor.

The project has been tested on MS Windows with Python 2.7 (though I'd like it to be noted that I am not primarily a Windows-user).

Status
------

The parameter fitting aspect of the project is reasonably flexible in fitting any Pyproj-supported projection to be fit with as few or as many unknown parameters.

I have found that the Winkel tripel projection gives excellent accuracy in terms of error between control and fitted points, and in terms of visual correctness:

![Winkel tripel fit of The Steppe map](http://i.stack.imgur.com/Femgi.jpg)

The algorithm found the best fit to be at central longitude of 47 W.

Technical notes
---------------
 §1. Note that most of these projections (see [http://www.remotesensing.org/geotiff/proj_list](http://www.remotesensing.org/geotiff/proj_list)) accept false easting and northing parameters, scalars which are added to all Cartesian locations. While the projection fitting can accommodate these readily, this is unnecessary as the we remove any affine (`a*x + b`) transform between the projection's output (in Cartesian space) and the GCPs' pixel locations using [Späth's algorithm (pdf)](http://hrcak.srce.hr/file/1425). In simpler terms, the projection fitting function will find and remove any rotation and translation that stands between the predicted pixel locations and the actual pixel locations --- treating such factors as unknowns to be estimated from data can possibly be detrimental to fit accuracy, and should not be done.

References
----------

The community at GIS.stackexchange has been very helpful --- see [http://gis.stackexchange.com/questions/43682/](http://gis.stackexchange.com/questions/43682/) --- thank you.
