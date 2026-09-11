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


def check_overwrite(file, existing_file):
    """
    Checks if a file already exists and raises an error if overwriting is not allowed.
    Args:
        file (knext.File): The file to check.
        existing_file (Enum): An enumeration value indicating the overwrite policy.
            It should have a `FAIL` member to signify that overwriting is not allowed.
    Raises:
        knext.InvalidParametersError: If the file exists and the overwrite policy is set to FAIL.
    """
    if existing_file == ExistingFile.FAIL.name and file.exists():
        raise knext.InvalidParametersError("File already exists.")


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
    description="""This node reads a single geospatial file from the selected file or URL. 
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
    data_url = knext.FileSelectionParameter(
        "Input file",
        "Select the file to read the data from or directly enter a remote URL.",
        placeholder_text="Select input file or enter URL...",
        validator=knut.check_file_selected,
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

        file_name = self.data_url.name.lower()
        if knut.is_web_url(self.data_url):
            return self._read_data(self.data_url.path, file_name)
        if file_name.endswith(".shp"):
            with knut.shapefile_to_local(self.data_url) as local_path:
                return self._read_data(str(local_path), file_name)
        with self.data_url.to_local() as local_path:
            return self._read_data(str(local_path), file_name)

    def _read_data(self, data_path: str, file_name: str):
        # the format is taken from the selected file's name because a staged copy of a
        # compressed GeoParquet file keeps only its last suffix
        import geopandas as gpd

        def urlread(url: str) -> gpd.GeoDataFrame:
            # Read the file directly first; modern GDAL/pyogrio can already stream a
            # remote .zip shapefile. If that fails on an online zip, fall back to GDAL's
            # /vsizip/vsicurl/ virtual filesystem, which handles some servers the plain
            # reader trips over.
            try:
                return gpd.read_file(url, engine="pyogrio", on_invalid="ignore")
            except Exception as direct_error:
                if url.lower().startswith("http") and url.lower().endswith(".zip"):
                    try:
                        return gpd.read_file(
                            "/vsizip/vsicurl/" + url,
                            engine="pyogrio",
                            on_invalid="ignore",
                        )
                    except Exception as vsizip_error:
                        raise RuntimeError(
                            f"Could not read {url} directly ({direct_error}) "
                            f"nor via /vsizip/vsicurl/ ({vsizip_error})"
                        )
                raise RuntimeError(f"Could not read {url}: {direct_error}")

        if file_name.endswith(".kml"):
            import fiona

            fiona.drvsupport.supported_drivers["KML"] = "r"
            gdf = gp.read_file(data_path, driver="KML")
        elif file_name.endswith(".kmz"):
            import zipfile
            import fiona

            zf = zipfile.ZipFile(data_path)
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
            gdf = gp.read_file("/vsizip/" + data_path + "/" + name, driver="KML")
        elif file_name.endswith(
            (".parquet", ".parquet.br", ".parquet.gz", ".parquet.snappy")
        ):
            gdf = gp.read_parquet(data_path)

        else:
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf = urlread(data_path)
            else:
                gdf = gp.read_file(
                    data_path,
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

    data_url = knext.FileSelectionParameter(
        "Output file",
        "Select the file to save the data to.",
        placeholder_text="Select output file...",
        is_writer=True,
        validator=knut.check_file_selected,
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

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = clean_dataframe(gdf)

        if self.dataformat == "Shapefile":
            target = knut.file_with_extension(self.data_url, ".shp")
            check_overwrite(target, self.existing_file)
            self._write_staged(target, self._to_file(gdf))

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
            target = knut.file_with_extension(self.data_url, file_extension)
            check_overwrite(target, self.existing_file)
            self._write_staged(
                target,
                lambda local_file: gdf.to_parquet(local_file, compression=compression),
            )
        elif self.dataformat == "GeoJSON":
            target = knut.file_with_extension(self.data_url, ".geojson")
            check_overwrite(target, self.existing_file)
            self._write_staged(target, self._to_file(gdf, driver="GeoJSON"))
        else:
            target = knut.file_with_extension(self.data_url, ".gml")
            check_overwrite(target, self.existing_file)
            self._write_staged(target, self._to_file(gdf, driver="GML"))
        return None

    def _write_staged(self, target, write):
        """
        Writes into a local staging directory and then writes every produced file to the
        target location: an existing target is not downloaded just to be overwritten, and
        the formats that consist of several files (Shapefile, GML) arrive complete.
        """
        import pathlib
        import tempfile

        with tempfile.TemporaryDirectory() as staging:
            local_file = pathlib.Path(staging) / target.name
            write(local_file)
            knut.write_file_set(local_file, target)

    def _to_file(self, gdf, driver=None):
        """Returns a writer that calls ``gdf.to_file`` with the node's encoding setting."""

        def write(local_file):
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf.to_file(local_file)
            elif driver is None:
                gdf.to_file(local_file, encoding=self.encoding)
            else:
                gdf.to_file(local_file, driver=driver, encoding=self.encoding)

        return write


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
    data_url = knext.FileSelectionParameter(
        "Input file",
        "Select the GeoPackage file or GeoDatabase folder to read the data from, or "
        "directly enter a remote URL.",
        placeholder_text="Select input file or enter URL...",
        selection_mode=knext.FileSelectionMode.FILE_OR_FOLDER,
        validator=knut.check_file_selected,
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
        if knut.is_web_url(self.data_url):
            return self._read_data(self.data_url.path)
        with self.data_url.to_local() as local_path:
            return self._read_data(str(local_path))

    def _read_data(self, data_path: str):
        import fiona
        import pandas as pd

        layerlist = fiona.listlayers(data_path)
        layer = self._get_layer(layerlist)

        if self.encoding == _EncodingOptions.AUTO.name:
            gdf = gp.read_file(
                data_path, layer=layer, engine="pyogrio", on_invalid="ignore"
            )
        else:
            gdf = gp.read_file(
                data_path,
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

    data_url = knext.FileSelectionParameter(
        "Output file",
        "Select the file to save the data to.",
        placeholder_text="Select output file...",
        is_writer=True,
        file_extension="gpkg",
        validator=knut.check_file_selected,
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

        target = knut.file_with_extension(self.data_url, ".gpkg")
        check_overwrite(target, self.existing_file)
        target.parent.mkdir()

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gdf.reset_index(drop=True)
        time_columns = gdf.select_dtypes(
            include=[
                'knime.pandas_type<struct<0:int64,1:int64>, {"value_factory_class":"org.knime.core.data.v2.time.LocalDateTimeValueFactory"}>'
            ]
        ).columns
        if len(time_columns) > 0:
            gdf[time_columns] = gdf[time_columns].astype(str)

        gdf = clean_dataframe(gdf)

        # an existing GeoPackage is downloaded first so that the layer is added to it
        with target.to_local(reupload=True) as file_name:
            if self.encoding == _EncodingOptions.AUTO.name:
                gdf.to_file(file_name, layer=self.data_layer, driver="GPKG")
            else:
                gdf.to_file(
                    file_name,
                    layer=self.data_layer,
                    driver="GPKG",
                    encoding=self.encoding,
                )

        return None
