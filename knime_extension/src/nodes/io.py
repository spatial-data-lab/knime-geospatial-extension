import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut

__category = knext.category(
    path="/community/geo",
    level_id="io",
    name="Spatial IO",
    description="Nodes that read and write spatial data in various formats.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/IOCategory.png",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/IO/"


class Compression(knext.EnumParameterOptions):
    NONE = (
        "None",
        "Does not use any compression at all.",
    )
    BORTLI = (
        "Brotli",
        "Successor to gzip with better compression. For more details see [here.](https://en.wikipedia.org/wiki/Brotli)",
    )
    GZIP = (
        "gzip",
        "Widely used and supported compression format. For more details see [here.](https://en.wikipedia.org/wiki/Gzip)",
    )
    SNAPPY = (
        "Snappy",
        "Compression format aiming for very high speed and reasonable compression. "
        + "For more details see [here.](https://en.wikipedia.org/wiki/Snappy_(compression))",
    )


class ExistingFile(knext.EnumParameterOptions):
    FAIL = (
        "Fail",
        "Will issue an error during the node's execution (to prevent unintentional overwrite).",
    )
    OVERWRITE = (
        "Overwrite",
        "Will replace any existing file.",
    )


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
    description="""This node reads a single geospatial file from the provided local file path or URL. 
    The supported file formats are the popular data types such as [Shapefile (.shp),](https://en.wikipedia.org/wiki/Shapefile)
zipped Shapefiles(.zip) with a single Shapefile, single-layer [Geopackage (.gpkg),](https://www.geopackage.org/) 
[GeoJSON (.geojson),](https://geojson.org/) [GeoParquet,](https://github.com/opengeospatial/geoparquet)
or [MapInfo (.tab)](https://gdal.org/en/latest/drivers/vector/mitab.html) files. 
In addition the node partially supports 
[Keyhole Markup Language (.kml)](https://en.wikipedia.org/wiki/Keyhole_Markup_Language) files or single
entry zipped [.kmz](https://developers.google.com/kml/documentation/kmzarchives) files. 
For more details on the limitations when reading these files see 
[here.](https://gdal.org/drivers/vector/kml.html#kml-reading)

Examples of standard local file paths are *C:\\KNIMEworkspace\\test.geojson* for Windows and
*/KNIMEworkspace/test.shp* for Linux. The node can also load resources directly from a web URL, for example to 
load a GeoJSON file from [geojson.xyz](http://geojson.xyz/) you would enter
*http://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_110m_land.geojson*.

**Note:** For larger files the node progress might not change for a time until the file is successfully read.
    """,
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "Read file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.read_file.html",
        "Read Parquet": "https://geopandas.org/en/stable/docs/reference/api/geopandas.read_parquet.html",
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

        import geopandas as gpd

        def urlread(url: str) -> gpd.GeoDataFrame:
            try:
                gdf = gpd.read_file(url)
                return gdf
            except Exception as e1:
                if url.startswith("http") and url.endswith(".zip"):
                    try:
                        vsizip_url = "/vsizip/vsicurl/" + url
                        gdf = gpd.read_file(vsizip_url)
                        return gdf
                    except Exception as e2:
                        raise RuntimeError(f"Error:{e2}")
                else:
                    raise RuntimeError(f"Error:{e1}")

        if self.data_url.lower().endswith(".kml"):
            import fiona

            fiona.drvsupport.supported_drivers["KML"] = "r"
            gdf = gp.read_file(self.data_url, driver="KML")
        elif self.data_url.lower().endswith(".kmz"):
            import zipfile
            import fiona

            zf = zipfile.ZipFile(self.data_url)
            names = zf.namelist()
            name = None
            for i in range(len(names)):
                if names[i].lower().endswith(".kml"):
                    if name is None:
                        name = names[i]
                    else:
                        raise knext.InvalidParametersError(
                            "Node supports only kmz files with a single kml file"
                        )
            fiona.drvsupport.supported_drivers["KML"] = "r"
            gdf = gp.read_file("/vsizip/" + self.data_url + "/" + name, driver="KML")
        elif (
            self.data_url.lower().endswith(".parquet")
            or self.data_url.lower().endswith(".parquet.br")
            or self.data_url.lower().endswith(".parquet.gz")
            or self.data_url.lower().endswith(".parquet.snappy")
        ):
            gdf = gp.read_parquet(self.data_url)
        else:
            gdf = urlread(self.data_url)

        if "<Row Key>" in gdf.columns:
            gdf = gdf.drop(columns="<Row Key>")
        if "<RowID>" in gdf.columns:
            gdf = gdf.drop(columns="<RowID>")
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
    description="""This node writes the data in the format of [Shapefile](https://en.wikipedia.org/wiki/Shapefile), 
    [GeoJSON](https://geojson.org/), or [GeoParquet](https://github.com/opengeospatial/geoparquet).
Examples of standard local file paths are *C:\\KNIMEworkspace\\test.shp* for Windows and
*/KNIMEworkspace/test.geojson* for Linux. 

The file extension e.g. *.shp*, *.geojson*,  or *.parquet* is appended automatically
depending on the selected file format if not specified.""",
    references={
        "Writing Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "To file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html",
        "To Parquet": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_parquet.html",
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
        """The file path for writing data. The file extension e.g. *.shp*, *.geojson*,  or *.parquet* is appended 
automatically depending on the selected file format if not specified.""",
        "",
    )

    existing_file = knext.EnumParameter(
        "If exists:",
        "Specify the behavior of the node in case the output file already exists.",
        lambda v: (
            ExistingFile.OVERWRITE.name
            if v < knext.Version(1, 2, 0)
            else ExistingFile.FAIL.name
        ),
        enum=ExistingFile,
        since_version="1.2.0",
    )

    dataformat = knext.StringParameter(
        "Output file format",
        "The file format to use.",
        "Shapefile",
        enum=["Shapefile", "GeoJSON", "GeoParquet"],
    )

    parquet_compression = knext.EnumParameter(
        "File compression",
        "The name of the compression to use or none.",
        Compression.NONE.name,
        enum=Compression,
        since_version="1.2.0",
    ).rule(knext.OneOf(dataformat, ["GeoParquet"]), knext.Effect.SHOW)

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
        if "<Row Key>" in gdf.columns:
            gdf = gdf.drop(columns="<Row Key>")
        if "<RowID>" in gdf.columns:
            gdf = gdf.drop(columns="<RowID>")
        if self.dataformat == "Shapefile":
            fileurl = knut.ensure_file_extension(self.data_url, ".shp")
            self.__check_overwrite(fileurl)
            gdf.to_file(fileurl)
        elif self.dataformat == "GeoParquet":
            if self.parquet_compression == Compression.NONE.name:
                file_extension = ".parquet"
                compression = None
            elif self.parquet_compression == Compression.BORTLI.name:
                file_extension = ".parquet.br"
                compression = "brotli"
            elif self.parquet_compression == Compression.GZIP.name:
                file_extension = ".parquet.gz"
                compression = "gzip"
            elif self.parquet_compression == Compression.SNAPPY.name:
                file_extension = ".parquet.snappy"
                compression = "snappy"
            fileurl = knut.ensure_file_extension(self.data_url, file_extension)
            self.__check_overwrite(fileurl)
            gdf.to_parquet(fileurl, compression=compression)
        else:
            fileurl = knut.ensure_file_extension(self.data_url, ".geojson")
            self.__check_overwrite(fileurl)
            gdf.to_file(fileurl, driver="GeoJSON")
        return None

    def __check_overwrite(self, fileurl):
        if self.existing_file == ExistingFile.FAIL.name:
            import os.path

            if os.path.exists(fileurl):
                raise knext.InvalidParametersError(
                    "File already exists and should not be overwritten."
                )


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
        import fiona
        import pandas as pd

        layerlist = fiona.listlayers(self.data_url)
        pnumber = pd.Series(range(0, 100)).astype(str).to_list()
        if self.data_layer in layerlist:
            src = fiona.open(self.data_url, layer=self.data_layer)
        elif self.data_layer in pnumber:
            nlayer = int(self.data_layer)
            src = fiona.open(self.data_url, layer=nlayer)
        else:
            src = fiona.open(self.data_url, layer=0)
        gdf = gp.GeoDataFrame.from_features(src)
        try:
            gdf.crs = src.crs
        except:
            print("Invalid CRS")
        gdf = gdf.reset_index(drop=True)
        if "<Row Key>" in gdf.columns:
            gdf = gdf.drop(columns="<Row Key>")
        if "<RowID>" in gdf.columns:
            gdf = gdf.drop(columns="<RowID>")
        listtable = pd.DataFrame({"layerlist": layerlist})
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
        time_columns = gdf.select_dtypes(
            include=[
                'knime.pandas_type<struct<0:int64,1:int64>, {"value_factory_class":"org.knime.core.data.v2.time.LocalDateTimeValueFactory"}>'
            ]
        ).columns
        if len(time_columns) > 0:
            gdf[time_columns] = gdf[time_columns].astype(str)
        if "<Row Key>" in gdf.columns:
            gdf = gdf.drop(columns="<Row Key>")
        if "<RowID>" in gdf.columns:
            gdf = gdf.drop(columns="<RowID>")
        gdf.to_file(file_name, layer=self.data_layer, driver="GPKG")
        return None
