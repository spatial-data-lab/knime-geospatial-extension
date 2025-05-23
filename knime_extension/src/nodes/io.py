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


def validate_path(path: str) -> None:
    # no path check
    pass


def clean_dataframe(df):
    """
    Cleans the given DataFrame by resetting its index and removing specific columns.

    This function resets the index of the DataFrame, dropping the old index,
    and removes the columns "<Row Key>" and "<RowID>" if they exist in the DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame to be cleaned.

    Returns:
        pandas.DataFrame: A cleaned DataFrame with the index reset and specified columns removed.
    """
    df = df.reset_index(drop=True)
    columns_to_drop = ["<Row Key>", "<RowID>"]
    return df.drop(columns=[col for col in columns_to_drop if col in df.columns])


def check_overwrite(fileurl, existing_file):
    """
    Checks if a file already exists and raises an error if overwriting is not allowed.
    Args:
        fileurl (str): The path to the file to check.
        existing_file (Enum): An enumeration value indicating the overwrite policy.
            It should have a `FAIL` member to signify that overwriting is not allowed.
    Raises:
        knext.InvalidParametersError: If the file exists and the overwrite policy is set to FAIL.
    """
    import os

    if existing_file == ExistingFile.FAIL.name and os.path.exists(fileurl):
        raise knext.InvalidParametersError("File already exists.")


def check_outdir(fileurl):
    """
    Ensures that the directory for the given file path exists. If the directory
    does not exist, it is created.
    Args:
        fileurl (str): The file path for which the directory should be checked
                       and created if necessary.
    Raises:
        OSError: If the directory cannot be created due to an operating system error.
    """
    import os

    output_dir = os.path.dirname(fileurl)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)


class _EncodingOptions(knext.EnumParameterOptions):
    AUTO = (
        "Auto",
        "Automatically detect the encoding from common options",
    )
    UTF8 = (
        "UTF-8",
        "Unicode Transformation Format - 8 bit. Default encoding suitable for most modern GIS data files.",
    )
    GB18030 = (
        "GB18030",
        "Chinese National Standard encoding. More comprehensive than GBK.",
    )
    GBK = (
        "GBK",
        "Chinese internal code specification. Common in Chinese GIS software.",
    )
    GB2312 = (
        "GB2312",
        "Basic Simplified Chinese character encoding.",
    )
    LATIN1 = (
        "ISO-8859-1",
        "Latin-1 encoding. Suitable for Western European languages.",
    )
    WINDOWS1252 = (
        "Windows-1252",
        "Windows Western European encoding. Common in Windows systems.",
    )
    ASCII = (
        "ASCII",
        "Basic ASCII encoding. Only for standard ASCII characters.",
    )

    @classmethod
    def get_default(cls):
        return cls.AUTO


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

The node can load resources directly from a web URL, for example to 
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
    data_url = knext.LocalPathParameter(
        "Input file path",
        "Select the file path or directly enter a remote URL for reading the data.",
        placeholder_text="Select input file path or enter URL...",
        validator=validate_path,
    )

    encoding = knext.EnumParameter(
        label="Encoding",
        description="Select the encoding for reading the data file.",
        default_value=_EncodingOptions.get_default().name,
        enum=_EncodingOptions,
        since_version="1.4.0",
        is_advanced=True,
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        exec_context.set_progress(
            0.4, "Reading file (This might take a while without progress changes)"
        )

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
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf = gp.read_file(self.data_url, engine="pyogrio", on_invalid="ignore")
            else:
                gdf = gp.read_file(
                    self.data_url,
                    encoding=self.encoding,
                    engine="pyogrio",
                    on_invalid="ignore",
                )

        gdf = clean_dataframe(gdf)
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

    data_url = knext.LocalPathParameter(
        "Output file path",
        "Select the file path for saving data.",
        placeholder_text="Select output file path...",
        validator=validate_path,
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
        enum=["Shapefile", "GeoJSON", "GeoParquet", "GML"],
    )

    parquet_compression = knext.EnumParameter(
        "File compression",
        "The name of the compression to use or none.",
        Compression.NONE.name,
        enum=Compression,
        since_version="1.2.0",
    ).rule(knext.OneOf(dataformat, ["GeoParquet"]), knext.Effect.SHOW)

    encoding = knext.EnumParameter(
        label="Encoding",
        description="Select the encoding for saving the data file.",
        default_value=_EncodingOptions.get_default().name,
        enum=_EncodingOptions,
        since_version="1.4.0",
        is_advanced=True,
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

        check_outdir(self.data_url)
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = clean_dataframe(gdf)

        if self.dataformat == "Shapefile":
            fileurl = knut.ensure_file_extension(self.data_url, ".shp")
            check_overwrite(fileurl, self.existing_file)
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf.to_file(fileurl)
            else:
                gdf.to_file(fileurl, encoding=self.encoding)

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
            check_overwrite(fileurl, self.existing_file)
            gdf.to_parquet(fileurl, compression=compression)
        elif self.dataformat == "GeoJSON":
            fileurl = knut.ensure_file_extension(self.data_url, ".geojson")
            check_overwrite(fileurl, self.existing_file)
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf.to_file(fileurl)
            else:
                gdf.to_file(fileurl, driver="GeoJSON", encoding=self.encoding)
        else:
            fileurl = knut.ensure_file_extension(self.data_url, ".gml")
            check_overwrite(fileurl, self.existing_file)
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf.to_file(fileurl)
            else:
                gdf.to_file(fileurl, driver="GML", encoding=self.encoding)
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

The node can load resources directly from a web URL e.g. 
*https://github.com/INSPIRE-MIF/gp-geopackage-encodings/raw/refs/heads/main/examples/GE-gpkg-template.gpkg*.

**Note:** For larger files the node progress might not change for a time until the file is successfully read.
    """,
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
        "Read file": "https://geopandas.org/en/stable/docs/reference/api/geopandas.read_file.html",
    },
)
class GeoPackageReaderNode:
    data_url = knext.LocalPathParameter(
        "Input file path",
        "Select the file path or directly enter a remote URL for reading the data.",
        placeholder_text="Select input file path or enter URL...",
        validator=validate_path,
    )

    data_layer = knext.StringParameter(
        "Input layer name or order for reading",
        "The layer name in the multiple-layer data.",
        "",
    )

    encoding = knext.EnumParameter(
        label="Encoding",
        description="Select the encoding for reading the data file.",
        default_value=_EncodingOptions.get_default().name,
        enum=_EncodingOptions,
        since_version="1.4.0",
        is_advanced=True,
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
        layer = self._get_layer(layerlist)

        if self.encoding == _EncodingOptions.AUTO.name:
            gdf = gp.read_file(
                self.data_url, layer=layer, engine="pyogrio", on_invalid="ignore"
            )
        else:
            gdf = gp.read_file(
                self.data_url,
                layer=layer,
                engine="pyogrio",
                on_invalid="ignore",
                encoding=self.encoding,
            )

        gdf = clean_dataframe(gdf)

        listtable = pd.DataFrame({"layerlist": layerlist})
        return knext.Table.from_pandas(gdf), knext.Table.from_pandas(listtable)

    def _get_layer(self, layerlist):
        if self.data_layer in layerlist:
            return self.data_layer
        elif self.data_layer.isdigit() and 0 <= int(self.data_layer) < 100:
            return int(self.data_layer)
        return 0


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

    data_url = knext.LocalPathParameter(
        "Output file path",
        "Select the file path for saving data.",
        placeholder_text="Select output file path...",
        validator=validate_path,
    )

    data_layer = knext.StringParameter(
        "Output layer name for writing",
        "The output layer name in the GeoPackage data.",
        "new",
    )

    encoding = knext.EnumParameter(
        label="Encoding",
        description="Select the encoding for saving the data file.",
        default_value=_EncodingOptions.get_default().name,
        enum=_EncodingOptions,
        since_version="1.4.0",
        is_advanced=True,
    )

    existing_file = knext.EnumParameter(
        "If exists:",
        "Specify the behavior of the node in case the output file already exists.",
        lambda v: (
            ExistingFile.OVERWRITE.name
            if v < knext.Version(1, 3, 0)
            else ExistingFile.FAIL.name
        ),
        enum=ExistingFile,
        since_version="1.4.0",
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

        check_overwrite(self.data_url, self.existing_file)

        check_outdir(self.data_url)

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

        gdf = clean_dataframe(gdf)

        if self.encoding == _EncodingOptions.AUTO.name:
            gdf.to_file(file_name, layer=self.data_layer, driver="GPKG")
        else:
            gdf.to_file(
                file_name, layer=self.data_layer, driver="GPKG", encoding=self.encoding
            )

        return None
