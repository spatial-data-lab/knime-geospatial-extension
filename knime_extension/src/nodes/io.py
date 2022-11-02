from typing import Callable
from wsgiref.util import shift_path_info
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


__category = knext.category(
    path="/geo",
    level_id="io",
    name="Spatial IO",
    description="Nodes that for reading and writing Geodata.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/IOCategory.png",
    after="",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/IO/"

############################################
# GeoFile Reader
############################################
@knext.node(
    name="GeoFile Reader",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "GeoFileReader.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="Geodata table",
    description="Geodata from the input file path",
)
@knut.geo_node_description(
    short_description="Read single layer GeoFile.",
    description="This node read the Shapefile, zipped  Shapefile or Geojson with geopandas.read_file().",
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
    },
)
class GeoFileReaderNode:
    data_url = knext.StringParameter(
        "Input File Path",
        "The file path for reading data",
        "https://raw.githubusercontent.com/UrbanGISer/Test/main/JsonMap/countries.geojson",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        gdf = gp.read_file(self.data_url)
        gdf = gdf.reset_index(drop=True)
        return knext.Table.from_pandas(gdf)


############################################
# GeoFile Writer
############################################
@knext.node(
    name="GeoFile Writer",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "GeoFileWriter.png",
    category=__category,
    after="GeoPackage Reader",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input portal",
)
@knut.geo_node_description(
    short_description="Write single layer GeoFile.",
    description="This node write the Shapefile,  or Geojson with geopandas.to_file().",
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
    },
)
class GeoFileWriterNode:
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    data_url = knext.StringParameter(
        "Output file path and file name",
        "The file path for writing data ended with .shp or .geo",
        "",
    )

    dataformat = knext.StringParameter(
        "Output File Foramt",
        "The file path for writing data ended with .shp or .geo",
        "Shapefile",
        enum=["Shapefile", "GeoJSON"],
    )

    def configure(self, configure_context, input_schema_1):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        if self.dataformat == "Shapefile":
            fileurl = f"{self.data_url}.shp"
            gdf.to_file(fileurl)
        else:
            fileurl = f"{self.data_url}.geojson"
            gdf.to_file(fileurl, driver="GeoJSON")
        return None


############################################
# GeoPackage Reader
############################################
@knext.node(
    name="GeoPackage Reader",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "GeoPackageReader.png",
    category=__category,
    after="GeoFile Reader",
)
@knext.output_table(
    name="Geodata table",
    description="Geodata from the input file path",
)
@knut.geo_node_description(
    short_description="Read GeoPackage layer.",
    description="This node read the GeoPackage data with geopandas.read_file().",
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
    },
)
class GeoPackageReaderNode:
    data_url = knext.StringParameter(
        "Input File Path",
        "The file path for reading data",
        "https://raw.githubusercontent.com/UrbanGISer/Test/main/JsonMap/countries.gpkg",
    )

    data_layer = knext.StringParameter(
        "Input layer name for readng",
        "The layer name in the GPKG data",
        # TODO we need pre-read layer information
        "countries",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        gdf = gp.read_file(self.data_url, layer=self.data_layer)
        gdf = gdf.reset_index(drop=True)
        return knext.Table.from_pandas(gdf)


############################################
# GeoPackage Writer
############################################
@knext.node(
    name="GeoPackage Writer",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "GeoPackageWriter.png",
    category=__category,
    after="GeoFile Writer",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path",
)
@knut.geo_node_description(
    short_description="Write GeoPackage layer.",
    description="This node write the GeoPackage data with geopandas.to_file().",
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
    },
)
class GeoPackageWriterNode:
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    data_url = knext.StringParameter(
        "Input File Path", "The file path for reading data", ""
    )

    data_layer = knext.StringParameter(
        "Input layer name for writing",
        "The layer name in the GPKG data",
        # TODO we need pre-read layer information and dectect layer conflict
        "",
    )

    def configure(self, configure_context, input_schema_1):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gdf.reset_index(drop=True)
        gdf.to_file(self.data_url, layer=self.data_layer, driver="GPKG")
        return None
