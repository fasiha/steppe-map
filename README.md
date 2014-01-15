The Steppe, or How to fit arbitrary projections to data
=======================================================

Introduction
------------

In the current Encyclopedia Britannica's online article on "Steppe, the", there is a fabulous map --- see [http://www.britannica.com/bps/media-view/3658/1/0/0](http://www.britannica.com/bps/media-view/3658/1/0/0) --- that I would love to see broken out of its static fetters and dance with different overlays (e.g., temperature, rainfall, etc.). Choosing to ignore copyright concerns for the sake of narrow educational pursuits, this project is specifically focused on determining the projection used by this map, in order to get as exact georeferencing as possible, i.e., to accurately convert the green pixels to longitude/latitudes.

Hence the subtitle of the project. How can one go from georeferenced control points (GCPs) to a projection's paramters.

Currently, the project leverages Pyproj and allows me to specify

- a projection ("aea" or the Albers equal-area),
- its parameters (lon/lat of false origin, and two standard parallels), and
- an initial numeric guess for these parameters

and after running a function minimization (provided by Scipy), can plot the image with the original and best-fit GCPs, as well as a coastline, courtesy of Natural Earth.

Requirements
------------

Clearly this tiny project stands on the shoulders of giants. Requirements include

- Numpy for function minimization
- Pyproj for all the projections. If it's not in Pyproj, I can't fit it (yet).

To see the original image on top of the GCPs, one needs the original image: it's linked from [http://www.britannica.com/bps/media-view/3658/1/0/0](http://www.britannica.com/bps/media-view/3658/1/0/0) so save it to `TheSteppe.jpg`.

To see the coastline, one needs

- Natural Earth, specifically, the `ne_10m_coastline` physical dataset
- pyshp or shapefile to read Natural Earth data: [https://code.google.com/p/pyshp/](https://code.google.com/p/pyshp/)
- GDAL's Python project called `osgeo` to convert Natural Earth's projection descriptor string (WKT) to a Proj4 descriptor.

Status
------

The parameter fitting aspect of the project is reasonably flexible in fitting any Pyproj-supported projection to be fit with as few or as many unknown parameters.

I have found that the Eckert V and Robinson projections give the lowest absolute errors between control points and fitted points, but the resulting projections' coastline can deviate from that of the image, indicating that though close, these aren't this map's projections. I've fitted a couple of other projections (Van Der Grinten and Albers equal-area) as well as numerous pseudocylindrical projections, but none were as good (in terms of error between GCP and fitted points) as these.

Perhaps some kind soul can give me further advice on projections or methods to try, or double-check my GCPs (made carefully in QGIS).

References
----------

The community at GIS.stackexchange has been very helpful to me --- see [http://gis.stackexchange.com/questions/43682/](http://gis.stackexchange.com/questions/43682/) --- and to others, which benefited me.
