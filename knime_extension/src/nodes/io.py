from typing import Callable
from wsgiref.util import shift_path_info
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import fiona

__category = knext.category(
    path="/community/geo",
    level_id="io",
    name="Spatial IO",
    description="Nodes that for reading and writing Geodata.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/IOCategory.png",
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
    description="Geodata from the input file.",
)
@knut.geo_node_description(
    short_description="Read single layer GeoFile.",
    description="""This node reads a single geospatial file from the path to the file or URL.The supported file 
formats are the popular data types such as [Shapefile (.shp),](https://en.wikipedia.org/wiki/Shapefile)
zipped Shapefiles(.zip) with a single Shapefile, single-layer [Geopackage (.gpkg),](https://www.geopackage.org/) 
or [GeoJSON (.geojson)](https://geojson.org/) files.

Examples of standard local file paths are *C:\\KNIMEworkspace\\test.geojson* for Windows and
*/KNIMEworkspace/test.shp* for Linux. The node can also load resources directly from a web URL, for example to 
load a GeoJSON file from [geojson.xyz](http://geojson.xyz/) you would enter
*http://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_110m_land.geojson*.

**Note:** For larger files the node progress might not change for a time until the file is successfully read.
    """,
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "Read file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.read_file.html",
    },
)
class GeoFileReaderNode:
    data_url = knext.StringParameter(
        "Input file path",
        "The file path for reading data.",
        "",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        exec_context.set_progress(
            0.4, "Reading file (This might take a while without progress changes)"
        )
        gdf = gp.read_file(self.data_url)
        gdf = knut.Turn_all_NA_column_as_str(gdf)
        if "<Row Key>" in gdf.columns:
            gdf = gdf.drop(columns="<Row Key>")
        return knext.Table.from_pandas(gdf)


############################################
# GeoFile Writer
############################################
@knext.node(
    name="GeoFile Writer",
    node_type=knext.NodeType.SINK,
    icon_path=__NODE_ICON_PATH + "GeoFileWriter.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input portal.",
)
@knut.geo_node_description(
    short_description="Write single layer GeoFile.",
    description="""This node writes the data in the format of [Shapefile](https://en.wikipedia.org/wiki/Shapefile) 
    or [GeoJSON](https://geojson.org/).
Examples of standard local file paths are *C:\\KNIMEworkspace\\test.shp* for Windows and
*/KNIMEworkspace/test.geojson* for Linux. 

The file extension e.g. *.shp* or *.geojson* is appended automatically
depending on the selected file format if not specified.

**Note:** Existing files will be overwritten without a warning!""",
    references={
        "Writing Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "To file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html",
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
        "Output file path",
        """The file path for writing data. The file extension e.g. *.shp* or *.geojson* is appended automatically
depending on the selected file format if not specified.""",
        "",
    )

    dataformat = knext.StringParameter(
        "Output file format",
        "The file format to use.",
        "Shapefile",
        enum=["Shapefile", "GeoJSON"],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        exec_context.set_progress(
            0.4, "Writing file (This might take a while without progress changes)"
        )
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        if self.dataformat == "Shapefile":
            fileurl = knut.ensure_file_extension(self.data_url, ".shp")
            gdf.to_file(fileurl)
        else:
            fileurl = knut.ensure_file_extension(self.data_url, ".geojson")
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
    after="",
)
@knext.output_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Geodata Layer",
    description="Layer information from the input file path.",
)
@knut.geo_node_description(
    short_description="Read GeoPackage layer",
    description="""This node reads [Geopackage,](https://www.geopackage.org/) GeoDatabase(GDB) files. 

You can specify the layer to read. If the layer is empty or wrong, the node will read the first layer. 
You can also enter the number of the layer to read starting at 0. The node will output the names of all layers as 
second output table, which can be used to revise the name of the target layer.

Examples of standard local file paths are *C:\\KNIMEworkspace\\test.gpkg* for Windows and
*/KNIMEworkspace/test.gpkg* for Linux. The node can also load resources directly from a web URL.

**Note:** For larger files the node progress might not change for a time until the file is successfully read.
    """,
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "Read file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.read_file.html",
    },
)
class GeoPackageReaderNode:
    data_url = knext.StringParameter(
        "Input file path",
        "The file path for reading data.",
        "",
    )

    data_layer = knext.StringParameter(
        "Input layer name or order for reading",
        "The layer name in the multiple-layer data.",
        "",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        exec_context.set_progress(
            0.4, "Reading file (This might take a while without progress changes)"
        )
        layerlist = fiona.listlayers(self.data_url)
        pnumber = pd.Series(range(0, 100)).astype(str).to_list()
        if self.data_layer in layerlist:
            gdf = gp.read_file(self.data_url, layer=self.data_layer)
        elif self.data_layer in pnumber:
            nlayer = int(self.data_layer)
            gdf = gp.read_file(self.data_url, layer=nlayer)
        else:
            gdf = gp.read_file(self.data_url, layer=0)
        gdf = gdf.reset_index(drop=True)
        listtable = pd.DataFrame({"layerlist": layerlist})
        try:
            test = knext.Table.from_pandas(gdf)
        except:
            gdf = pd.DataFrame(gdf.drop(columns="geometry"))
        return knext.Table.from_pandas(gdf), knext.Table.from_pandas(listtable)


############################################
# GeoPackage Writer
############################################
@knext.node(
    name="GeoPackage Writer",
    node_type=knext.NodeType.SINK,
    icon_path=__NODE_ICON_PATH + "GeoPackageWriter.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knut.geo_node_description(
    short_description="Write GeoPackage layer.",
    description="""This node writes the data as new [Geopackage](https://www.geopackage.org/) file or 
as layer into an existing file.
Examples of standard local file paths are *C:\\KNIMEworkspace\\test.gpkg* for Windows and
*/KNIMEworkspace/test.gpkg* for Linux. 

**Note:** If file and layer already exist, the layer will be overwritten without a warning!
    """,
    references={
        "Writing Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "To file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html",
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
        "Output file path",
        "The file path for saving data.",
        "",
    )

    data_layer = knext.StringParameter(
        "Output layer name for writing",
        "The output layer name in the GeoPackage data.",
        "new",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        exec_context.set_progress(
            0.4, "Writing file (This might take a while without progress changes)"
        )
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gdf.reset_index(drop=True)
        file_name = knut.ensure_file_extension(self.data_url, ".gpkg")
        gdf.to_file(file_name, layer=self.data_layer, driver="GPKG")
        return None
