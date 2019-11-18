# Generate spherical Voronoi polygons (Thiessen polygons)

This is a Python 2.x script that generates Voronoi polygons on a sphere under the consideration of a curved planetary surface. The script automatically takes care of polar and Date Line intersections. The polygon geometries together with their geodesic area are saved in the ESRI Shapefile format. 

Voronoi polygons can be constructed from either a random set of data points or a list of lat/lon coordinates (Format: [[lon1, lat1],[lon_2, lat_2],[lon_n, lat_n]]). To use a random set of data points, set the variables use_random_set_of_input_points = True and number_of_random_points = 999 to the number of data points you want to generate (20 points minimum). To use an individual list of lat/lon coordinates, modify the input_samples list.

## You need to modify the following variables when using the script (lines 700 ff):

**use_random_set_of_input_points** - True if using randomly distributed points to generate Voronoi polygons, False if using your own point coordinates (20 minimum)

**number_of_random_points** - Number of random points that should be generated, if using randomly distributed points (>= 20 points required)

**input_samples** - Point coordinates to generate Voronoi polygons from, if using own point coordinates (>= 20 points required)

**input_samples_point_shapefile_output_path** - Path to save the input points shapefile: 'C:\path\to\point_file.shp'

**voronoi_polygons_shapefile_output_path** - Path to save the Voronoi polygon shapefile: 'C:\path\to\polygon_file.shp'

**geogr_sr_text** - Geographic coordinate system WKT, defining the reference body the data is located on 

## The script uses the following external libraries

numpy<br>
gdal, ogr, osr (all from the GDAL package)<br>
shapely<br>
scipy<br>

*(Shout out to all the contributors!)* :squirrel:

## Issues

The script uses the polygon vertices that are determined by the SpericalVoronoi function from the scipy library to generate the geodesic Voronoi polygon boundaries. However, I noticed that the SphericalVoronoi function (as of November 2019) has difficulties to calculate the correct polygon vertices when the number of input points are low. This usually happens when at least one polygon is very large (larger than one hemisphere). The likelyhood that we end up with a Voronoi polygon of such size increases, the fewer input points we have. To quantify this, I investigated the rate of failed Voronoi polygon constructions as a function of the number of input points. I did 1000 runs for each randomly distributed set of input points. The error rate reached 0.0% (0 out of 1000 randomly distributed test sets yielded in wrongly determined Voronoi polygons) when the number of input points was at least 20. For this reason, the number of input points should be 20 or more.

## Examples - Voronoi polygons from 1000 randomly located points on a sphere with Earth-like dimensions


Equirectangular projection

Mollweide projection (centered at 0°;0°)

Lambert Azimuthal Equal Area (cetered at the North Pole)




