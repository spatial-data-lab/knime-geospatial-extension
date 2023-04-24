# The root category of all Geospatial categories
import knime_extension as knext
import util.knime_utils as knut
import sys


# this section is used for things that need to be executed prior any of the nodes is used e.g. global setups such as
# setting lookup directories etc.
def __initialize_pyproj():
    """
    Attach the bundled pyproj_db
    """
    import pyproj
    import os.path as os

    pyproj_path = os.join(knut.get_env_path(), "Library\share\proj")
    pyproj.datadir.set_data_dir(pyproj_path)


__initialize_pyproj()

# Fake the libpysal.examples import to prevent it from downloading data from the internet which might cause problems
# in environments without internet access: https://github.com/spatial-data-lab/knime-geospatial-extension/issues/165
class __LibpysalExamplesModuleMock:
    pass


sys.modules["libpysal.examples"] = __LibpysalExamplesModuleMock()


# This defines the root Geospatial KNIME category that is displayed in the node repository
category = knext.category(
    path="/community",
    level_id="geo",
    name="Geospatial Analytics",
    description="Nodes for Geospatial Analytics",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/Geospatial.png",
)


# The different node files
import nodes.calculation
import nodes.conversion
import nodes.io
import nodes.locationanalysis
import nodes.opendata
import nodes.spatialmodels
import nodes.spatialstatistics
import nodes.spatialtool
import nodes.transform
import nodes.visualize
import nodes.spatialnetwork

# collection of deprecated nodes
import nodes.deprecated


# import nodes.geolab
