# The root category of all Geospatial categories
import knime_extension as knext

# This defines the root Geospatial KNIME category that is displayed in the node repository
category = knext.category(
    path="/",
    level_id="geo",
    name="Geospatial Analytics",
    description="Nodes for Geospatial Analytics",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/Geospatial.png",
)


# The different node files
import nodes.calculation
import nodes.conversion
import nodes.opendata
import nodes.spatialstatistics
import nodes.spatialtool
import nodes.transform
import nodes.visualize
import nodes.geolab
import nodes.spatialmodels
import nodes.io
import nodes.locationanalysis
