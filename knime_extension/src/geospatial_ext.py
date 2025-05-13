# The root category of all Geospatial categories
import nodes.deprecated
import nodes.spatialnetwork
import nodes.visualize
import nodes.transform
import nodes.spatialtool
import nodes.spatialstatistics
import nodes.spatialmodels
import nodes.spatialclustering
import nodes.opendata
import nodes.locationanalysis
import nodes.io
import nodes.conversion
import nodes.calculation
import util.knime_utils as knut
import knime_extension as knext
import sys

# this section is used for things that need to be executed prior the Geopandas lib is imported which happens with the
# import util.knime_utils line

# enforce use of Shapely 2.0 even if PyGEOS is installed which is required by the momepy, segregation and tobler package
import os

os.environ["USE_PYGEOS"] = "0"


# this section is used for things that need to be executed prior any of the nodes is used e.g. global setups such as
# setting lookup directories etc.
def __initialize_pyproj():
    """
    Attach the bundled pyproj_db
    """
    import pyproj
    import os
    import platform
    import sys

    env_path = knut.get_env_path()
    
    # Detect OS and build appropriate path
    if platform.system() == "Windows":
        pyproj_path = os.path.join(env_path, "Library", "share", "proj")
    elif platform.system() == "Darwin":  # macOS
        # Try common macOS paths
        potential_paths = [
            os.path.join(env_path, "share", "proj"),
            os.path.join(env_path, "..", "share", "proj"),
            "/opt/homebrew/share/proj"
        ]
        pyproj_path = None
        for path in potential_paths:
            if os.path.exists(path):
                pyproj_path = path
                break
        
        if pyproj_path is None:
            print("Warning: Could not find PyProj data directory. CRS operations may fail.", file=sys.stderr)
            return
    else:  # Linux and others
        pyproj_path = os.path.join(env_path, "share", "proj")
    
    # Only set the directory if it exists
    if os.path.exists(pyproj_path):
        pyproj.datadir.set_data_dir(pyproj_path)
        print(f"Set PyProj data directory to: {pyproj_path}", file=sys.stderr)
    else:
        print(f"Warning: PyProj data directory not found at {pyproj_path}. CRS operations may fail.", file=sys.stderr)


__initialize_pyproj()


# Fake the libpysal.examples import to prevent it from downloading data from the internet which might cause problems
# in environments without internet access: https://github.com/spatial-data-lab/knime-geospatial-extension/issues/165
# class __LibpysalExamplesModuleMock:
#     pass


# sys.modules["libpysal.examples"] = __LibpysalExamplesModuleMock()


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

# collection of deprecated nodes


# import nodes.geolab
