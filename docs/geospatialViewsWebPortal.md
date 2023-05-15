# Geospatial Views in KNIME Data Apps
In order to use the Geospatial Views in a [KNIME Data App](https://www.knime.com/data-apps) on the KNIME WebPortal 
you need to add the following line to the [knime-server.config file](https://docs.knime.com/latest/webportal_admin_guide/index.html#knime-server-configuration-file):
```
com.knime.server.webportal.csp=default-src 'self'; script-src 'unsafe-inline' 'unsafe-eval' 'self'; style-src 'unsafe-inline' 'self';img-src 'self' 
data: *.arcgisonline.com *.autonavi.com *.basemaps.cartocdn.com *.cloudfront.net/kepler.gl/ *.openrailwaymap.org *.openstreetmap.org 
stamen-tiles-a.a.ssl.fastly.net *.strava.com *.earthdata.nasa.gov; connect-src 'self' *.mapbox.com *.cloudfront.net/kepler.gl/; 
font-src 'self' data:;worker-src blob: 'self';
```
This will allow the views to load the required libraries from the KNIME Server and the background tile images from the different providers such as MapBox or OpenstreetMap.

For further information about the different options see the [WebPortal Admin Guide](https://docs.knime.com/latest/webportal_admin_guide/index.html#knime-server-configuration-file-options-webportal).
