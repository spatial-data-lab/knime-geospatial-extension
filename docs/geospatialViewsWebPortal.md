# Geospatial Views in KNIME Data Apps

The Geospatial Views load visualization libraries from the KNIME Business Hub or KNIME Server during viewing of the data. In addition, dependent on the view and settings they might need to download background tile images from the different providers such as MapBox or OpenstreetMap. By default, these downloads are prevented by the default security settings and need to be enabled as described in the following sections.

## KNIME Business Hub
In order to use the Geospatial Views in a [KNIME Data App](https://www.knime.com/data-apps) on the [KNIME Business Hub](https://www.knime.com/knime-business-hub)
you need to add the following line to the [Content Security Policy for Data Apps](https://docs.knime.com/latest/business_hub_admin_guide/index.html#configure-browser-security):

```
default-src 'self'; script-src 'unsafe-inline' 'unsafe-eval' 'self'; style-src 'unsafe-inline' 'self';img-src 'self' data: *.arcgisonline.com *.autonavi.com *.basemaps.cartocdn.com *.cloudfront.net/kepler.gl/ *.openrailwaymap.org *.openstreetmap.org stamen-tiles-a.a.ssl.fastly.net *.strava.com *.earthdata.nasa.gov; connect-src 'self' *.mapbox.com *.cloudfront.net/kepler.gl/; font-src 'self' data:;
```
For further information see the [Configure Browser Security](https://docs.knime.com/latest/business_hub_admin_guide/index.html#configure-browser-security) section of the [KNIME Business Hub Admin Guide](https://docs.knime.com/latest/business_hub_admin_guide/index.html).

## KNIME Server
In order to use the Geospatial Views in a [KNIME Data App](https://www.knime.com/data-apps) on the KNIME Server 
you need to add the following line to the [knime-server.config file](https://docs.knime.com/latest/webportal_admin_guide/index.html#knime-server-configuration-file):
```
com.knime.server.webportal.csp=default-src 'self'; script-src 'unsafe-inline' 'unsafe-eval' 'self'; style-src 'unsafe-inline' 'self';img-src 'self' data: *.arcgisonline.com *.autonavi.com *.basemaps.cartocdn.com *.cloudfront.net/kepler.gl/ *.openrailwaymap.org *.openstreetmap.org stamen-tiles-a.a.ssl.fastly.net *.strava.com *.earthdata.nasa.gov; connect-src 'self' *.mapbox.com *.cloudfront.net/kepler.gl/; font-src 'self' data:;worker-src blob: 'self';
```
For further information about the different options see the [WebPortal Admin Guide](https://docs.knime.com/latest/webportal_admin_guide/index.html#knime-server-configuration-file-options-webportal).
