# lingbo
import geopandas as gp
import logging
import knime_extension as knext
import util.knime_utils as knut
import util.projection as kproj

LOGGER = logging.getLogger(__name__)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialTool/"

__category = knext.category(
    path="/community/geo",
    level_id="deprecated",
    name="Deprecated Geospatial Nodes",
    description="Deprecated Geospatial Nodes",
    icon="icons/icon/SpatialToolCategory.png",
)

############################################
# Buffer
############################################


@knext.node(
    name="Buffer",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/SpatialTool/Buffer.png",
    category=__category,
    after="",
    is_deprecated=True,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geodata",
)
@knut.geo_node_description(
    short_description="Generate buffer zone based on a given distance.",
    description="""This node generates polygons representing all points within a given distance of each geometric object 
    based on geopandas.GeoSeries.buffer() with default parameters (resolution=16), which derives from Shapley object.buffer.
    """,
    references={
        "GeoSeries.buffer": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.buffer.html",
        "Shapley object.buffer": "https://shapely.readthedocs.io/en/latest/manual.html#object.buffer",
    },
)
class BufferNode:
    """
    This node aggregate generate buffer zone based on a given distance.
    """

    geo_col = knut.geo_col_parameter()

    bufferdist = knext.DoubleParameter(
        "Buffer distance", "The buffer distance for geometry. ", 1000.0
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = gp.GeoDataFrame(input.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting buffering...")
        gdf[knut.get_unique_column_name("geometry", input.schema)] = gdf.buffer(
            self.bufferdist
        )
        exec_context.set_progress(0.1, "Buffering done")
        LOGGER.debug(
            "Feature geometry " + self.geo_col + " buffered with" + str(self.bufferdist)
        )
        return knext.Table.from_pandas(gdf)


############################################
# Simplify
############################################


@knext.node(
    name="Simplify",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/SpatialTool/Simplify.png",
    category=__category,
    after="",
    is_deprecated=True,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to simplify",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed geo input table",
)
@knut.geo_node_description(
    short_description="Simplify the geometry",
    description="""This node returns a geometry feature containing a simplified representation of each geometry 
    with geopandas.simplify(). The algorithm (Douglas-Peucker) recursively splits the original line into smaller 
    parts and connects these parts’ endpoints by a straight line. Then, it removes all points whose distance 
    to the straight line is smaller than tolerance. It does not move any points and it always preserves endpoints 
    of the original line or polygon.
    """,
    references={
        "GeoSeries.simplify": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.simplify.html",
        "Shapely object.simplify": "http://shapely.readthedocs.io/en/latest/manual.html#object.simplify",
    },
)
class SimplifyNode:
    """
    This node returns a geometry feature containing a simplified representation of each geometry.
    """

    geo_col = knut.geo_col_parameter()

    simplifydist = knext.DoubleParameter(
        label="Simplification tolerance",
        description="""The simplification tolerance distances for geometry.
        All parts of a simplified geometry will be no more than tolerance distance from the original. 
        It has the same units as the coordinate reference system of the GeoSeries. 
        For example, using tolerance=100 in a projected CRS with meters as units means a distance of 100 meters in reality. 
        """,
        default_value=1.0,
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = gp.GeoDataFrame(input.to_pandas(), geometry=self.geo_col)
        gdf[knut.get_unique_column_name("geometry", input.schema)] = (
            gdf.geometry.simplify(self.simplifydist)
        )
        gdf = gdf.reset_index(drop=True)
        exec_context.set_progress(0.1, "Transformation done")
        LOGGER.debug("Feature Simplified")
        return knext.Table.from_pandas(gdf)


############################################
# Nearest Join
############################################


class _JoinModes(knext.EnumParameterOptions):
    INNER = ("Inner", "Retains only matching rows from both input tables.")
    LEFT = (
        "Left",
        "Retains all rows form the left and only matching rows from the right input tables.",
    )
    RIGHT = (
        "Right",
        "Retains all rows form the right and only matching rows from the left input tables.",
    )

    @classmethod
    def get_default(cls):
        return cls.INNER


@knext.node(
    name="Nearest Join",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "NearestJoin.png",
    category=__category,
    after="",
    is_deprecated=True,
)
@knext.input_table(
    name="Left geo table",
    description="Left table with geometry column to join on",
)
@knext.input_table(
    name="Right geo table",
    description="Right table with geometry column to join on",
)
@knext.output_table(
    name="Joined geo table",
    description="Joined geo table",
)
@knut.geo_node_description(
    short_description="Merges the two input tables based on their spatial relationship.",
    description="""This node will merge the left (top) and the right (bottom) table based on the distance between 
    their geometries of the two selected columns to one another. Distance is calculated in CRS units and is returned 
    in the column NearDist. Both layers must be in the same Coordinate Reference System (CRS), otherwise, the CRS of
    right table will be transformed to that of the left table.
    """,
    references={
        "Spatial joins": "https://geopandas.org/en/stable/gallery/spatial_joins.html",
        "Merging data": "https://geopandas.org/en/stable/docs/user_guide/mergingdata.html#nearest-joins",
        "sjoin_nearest": "https://geopandas.org/en/stable/docs/reference/api/geopandas.sjoin_nearest.html",
    },
)
class NearestJoinNode:
    left_geo_col = knext.ColumnParameter(
        "Left geometry column",
        "Select the geometry column from the left (top) input table to join on.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    right_geo_col = knext.ColumnParameter(
        "Right geometry column",
        "Select the geometry column from the right (bottom) input table to join on.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    join_mode = knext.EnumParameter(
        label="Join mode",
        description="The join mode determines of which input table the rows should be retained.",
        default_value=_JoinModes.get_default().name,
        enum=_JoinModes,
    )

    maxdist = knext.DoubleParameter(
        "Maximum distance",
        "Maximum distance within which to query for nearest geometry. Must be greater than 0 ",
        1000.0,
    )

    crs_info = knext.StringParameter(
        label="CRS for distance calculation",
        description=kproj.DEF_CRS_DESCRIPTION,
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        self.left_geo_col = knut.column_exists_or_preset(
            configure_context, self.left_geo_col, left_input_schema, knut.is_geo
        )
        self.right_geo_col = knut.column_exists_or_preset(
            configure_context, self.right_geo_col, right_input_schema, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = knut.load_geo_data_frame(left_input, self.left_geo_col, exec_context)
        right_gdf = knut.load_geo_data_frame(
            right_input, self.right_geo_col, exec_context
        )
        knut.check_canceled(exec_context)
        left_gdf.to_crs(self.crs_info, inplace=True)
        right_gdf.to_crs(self.crs_info, inplace=True)
        gdf = gp.sjoin_nearest(
            left_gdf,
            right_gdf,
            how=self.join_mode.lower(),
            max_distance=self.maxdist,
            distance_col="NearDist",
            lsuffix="1",
            rsuffix="2",
        )
        # reset the index since it might contain duplicates after joining
        gdf.reset_index(drop=True, inplace=True)
        # drop additional index columns if they exist
        gdf.drop(["index_1", "index_2"], axis=1, errors="ignore", inplace=True)
        gdf = gdf[[col for col in gdf.columns if not col.startswith("<RowID>")]]
        return knut.to_table(gdf, exec_context)


############################################
# Multiple Ring Buffer
############################################


@knext.node(
    name="Multiple Ring Buffer",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "MultipleRingBuffer.png",
    category=__category,
    after="",
    is_deprecated=True,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to buffer",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed table by Multiple Ring Buffer",
)
@knut.geo_node_description(
    short_description="This node generate multiple polygons with a series distances of each geometric object.",
    description="""This node generate multiple polygons with a series distances of each geometric object.

**Note:** If the input table contains multiple rows the node first computes the union of all geometries before 
computing the buffers from the union.
    """,
    references={
        "Buffer": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.buffer.html",
    },
)
class MultiRingBufferNode:
    geo_col = knut.geo_col_parameter()

    bufferdist = knext.StringParameter(
        "Serial buffer distances with coma",
        "The buffer distances for geometry ",
        "10,20,30",
    )

    bufferunit = knext.StringParameter(
        label="Serial buffer distances",
        description="The buffer distances for geometry ",
        default_value="Meter",
        enum=[
            "Meter",
            "KiloMeter",
            "Mile",
        ],
    )

    crs_info = knext.StringParameter(
        label="CRS for buffering distance calculation",
        description=kproj.DEF_CRS_DESCRIPTION,
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        gdf.to_crs(self.crs_info, inplace=True)

        from pyproj import CRS  # For CRS Units check

        crsinput = CRS.from_user_input(gdf.crs)
        if crsinput.is_geographic:
            logging.warning("Unit as Degree, Please use Projected CRS")
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting buffering...")
        # transfrom string list to number
        import numpy as np

        bufferlist = np.array(self.bufferdist.split(","), dtype=np.int64)
        if self.bufferunit == "Meter":
            bufferlist = bufferlist
        elif self.bufferunit == "KiloMeter":
            bufferlist = bufferlist * 1000
        else:
            bufferlist = bufferlist * 1609.34
        # sort list
        bufferlist = bufferlist.tolist()
        bufferlist.sort()
        if gdf.shape[0] > 1:
            gdf_union = gdf.unary_union
            gdfunion = gp.GeoDataFrame(geometry=gp.GeoSeries(gdf_union), crs=gdf.crs)
        else:
            gdfunion = gdf
        c1 = gp.GeoDataFrame(geometry=gdfunion.buffer(bufferlist[0]))
        c2 = gp.GeoDataFrame(geometry=gdfunion.buffer(bufferlist[1]))
        gdf0 = gp.overlay(c1, c2, how="union")
        if len(bufferlist) > 2:
            # Construct all other rings by loop
            for i in range(2, len(bufferlist)):
                ci = gp.GeoDataFrame(geometry=gdfunion.buffer(bufferlist[i]))
                gdf0 = gp.overlay(gdf0, ci, how="union")
        # Add ring radius values as a new column
        gdf0["dist"] = bufferlist
        gdf0 = gdf0.reset_index(drop=True)
        exec_context.set_progress(0.1, "Buffering done")
        return knut.to_table(gdf0, exec_context)


############################################
# Euclidean Distance
############################################


@knext.node(
    name="Euclidean Distance",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "EuclideanDistance.png",
    category=__category,
    after="",
    is_deprecated=True,
)
@knext.input_table(
    name="Left geo table",
    description="Left table with geometry column. ",
)
@knext.input_table(
    name="Right geo table",
    description="Right table with geometry column.",
)
@knext.output_table(
    name="Geo table distance",
    description="Euclidean distance between geometry objects.",
)
@knut.geo_node_description(
    short_description="This node will calculate the Euclidean distance between two geometries.",
    description="""This node will calculate the Euclidean distance between two geometries.
    """,
    references={
        "Distance": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.distance.html",
    },
)
class EuclideanDistanceNode:
    left_geo_col = knext.ColumnParameter(
        "Left geometry column",
        "Select the geometry column from the left (top) input table to calculate.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    right_geo_col = knext.ColumnParameter(
        "Right geometry column",
        "Select the geometry column from the right (bottom) input table to calculate.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    crs_info = knext.StringParameter(
        label="CRS for distance calculation",
        description=kproj.DEF_CRS_DESCRIPTION,
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        self.left_geo_col = knut.column_exists_or_preset(
            configure_context, self.left_geo_col, left_input_schema, knut.is_geo
        )
        self.right_geo_col = knut.column_exists_or_preset(
            configure_context, self.right_geo_col, right_input_schema, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(
            right_input.to_pandas(), geometry=self.right_geo_col
        )
        knut.check_canceled(exec_context)
        right_gdf.to_crs(self.crs_info, inplace=True)
        left_gdf.to_crs(self.crs_info, inplace=True)
        # left_gdf['LID'] = range(1,(left_gdf.shape[0]+1))
        # right_gdf['RID'] = range(1,(right_gdf.shape[0]+1))
        mergedf = left_gdf.merge(right_gdf, how="cross")
        mergedf_x = gp.GeoDataFrame(geometry=mergedf["geometry_x"])
        mergedf_y = gp.GeoDataFrame(geometry=mergedf["geometry_y"])
        mergedf["EuDist"] = mergedf_x.distance(mergedf_y, align=False)
        mergedf = mergedf.reset_index(drop=True)
        return knext.Table.from_pandas(mergedf)


############################################
# GeoFile Reader / GeoFile Writer / GeoPackage Reader / GeoPackage Writer
############################################
# these four nodes used a plain local file path (knext.LocalPathParameter) before they
# were switched to knext.FileSelectionParameter, which is not backward compatible with
# workflows that already have a path configured - so the original nodes are kept here
# unchanged and the new file handling was introduced under a new node, see nodes.io


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


@knext.node(
    name="GeoFile Reader",
    node_type=knext.NodeType.SOURCE,
    icon_path="icons/icon/IO/GeoFileReader.png",
    category=__category,
    after="",
    is_deprecated=True,
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
                gdf = urlread(self.data_url)
            else:
                gdf = gp.read_file(
                    self.data_url,
                    encoding=self.encoding,
                    engine="pyogrio",
                    on_invalid="ignore",
                )

        gdf = clean_dataframe(gdf)
        return knext.Table.from_pandas(gdf)


@knext.node(
    name="GeoFile Writer",
    node_type=knext.NodeType.SINK,
    icon_path="icons/icon/IO/GeoFileWriter.png",
    category=__category,
    after="",
    is_deprecated=True,
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


@knext.node(
    name="GeoPackage Reader",
    node_type=knext.NodeType.SOURCE,
    icon_path="icons/icon/IO/GeoPackageReader.png",
    category=__category,
    after="",
    is_deprecated=True,
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


@knext.node(
    name="GeoPackage Writer",
    node_type=knext.NodeType.SINK,
    icon_path="icons/icon/IO/GeoPackageWriter.png",
    category=__category,
    after="",
    is_deprecated=True,
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


############################################
# Spatial Weights
############################################
# the "Get spatial weights matrix from file" option used a plain local file path
# (knext.StringParameter) before it was switched to knext.FileSelectionParameter, which
# is not backward compatible with workflows that already have a path configured - so the
# original node is kept here unchanged and the new file handling was introduced under a
# new node, see nodes.spatialstatistics


@knext.node(
    name="Spatial Weights",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/SpatialStatistics/SpatialWeight.png",
    category=__category,
    after="",
    is_deprecated=True,
)
@knext.input_table(name="Geo table", description="Table with geometry column.")
@knext.output_table(name="Spatial Weights", description="Spatial Weights.")
class spatialWeights:
    """Constructs a contiguity spatial weights matrix from the input data.
    This node constructs a contiguity spatial weights matrix from the input data.
    """

    geo_col = knut.geo_col_parameter(
        description="The name of the geometry column in the input data."
    )

    id_col = knext.ColumnParameter(
        "ID column",
        """Select the column which contains for each observation in the input data a unique ID, it should be an integer column.
        If 'none' is selected, the IDs will be automatically generated from 0 to the number of rows flowing
        the order of the input data.
        The IDs of this column must match with the values of the ID column selected in subsequent ESDA or spatial
        modeling nodes.
        """,
        include_none_column=True,
        column_filter=knut.is_long,
        since_version="1.1.0",
    )

    category = knext.StringParameter(
        "Weights category",
        """ The type of spatial weights to construct. Defaults to 'Queen'. Other options are 'Rook',
        'Binary Distance Band', 'Inverse Distance', 'Lattice', 'K nearest', 'Kernel', and
        'Get spatial weights matrix from file'.

        - `Queen` which will construct a queen contiguity weights matrix, is more robust and more suitable for areal unit data. The queen criterion is somewhat more encompassing and defines
        neighbors as spatial units sharing a common edge or a common vertex.
        - The `Rook` criterion defines neighbors by the existence of a common edge between two spatial units. Therefore, the number of neighbors according to the
        queen criterion will always be at least as large as for the rook criterion.
        - When choosing `K nearest`, select the nearest number 'Nearest k' in the following options. K-nearest are often used for point data.
        - When selecting `Binary Distance Band`, please select
        the distance threshold 'Threshold' in the following options.
        - When selecting `Inverse Distance`, please select
        the distance threshold 'Threshold' and the corresponding power 'Power' in the following options.
        - When 'Your own' is selected, please enter the path of the spatial weights matrix in CSV format in the
        following options.
        - More details about spatial weights, please see the [GeoDa center website](https://geodacenter.github.io/documentation.html).
        """,
        "Queen",
        enum=[
            "Queen",
            "Rook",
            "Binary Distance Band",
            "Inverse Distance",
            "K nearest",
            "Lattice",
            "Kernel",
            "Get spatial weights matrix from file",
        ],
    )
    order = knext.IntParameter(
        "Order",
        """The order of the weight matrix is 1 by default. Users can change the order of the weights, higher order
        weights will treat further units as neighbors.""",
        1,
    ).rule(
        knext.OneOf(
            category,
            [
                "Queen",
                "Rook",
                "Binary Distance Band",
                "Inverse Distance",
                "K nearest",
                "Lattice",
            ],
        ),
        knext.Effect.SHOW,
    )

    Threshold = knext.IntParameter(
        "Threshold for Inverse Distance or Binary Distance Band",
        """The distance threshold for constructing binary distance band and inverse distance weights. Defaults to 1""",
        1,
    ).rule(
        knext.OneOf(category, ["Binary Distance Band", "Inverse Distance"]),
        knext.Effect.SHOW,
    )

    Rows = knext.IntParameter(
        "Rows for Lattice",
        "The number of rows for constructing a lattice spatial weights matrix. Defaults to 5.",
        5,
    ).rule(knext.OneOf(category, ["Lattice"]), knext.Effect.SHOW)

    Columns = knext.IntParameter(
        "Columns for Lattice",
        "The number of columns for constructing a lattice spatial weights matrix. Defaults to 5.",
        5,
    ).rule(knext.OneOf(category, ["Lattice"]), knext.Effect.SHOW)

    Nearest_k = knext.IntParameter(
        "Nearest k",
        "The number of nearest neighbors to use for constructing k-nearest neighbors weights. Defaults to 4.",
        4,
    ).rule(knext.OneOf(category, ["K nearest"]), knext.Effect.SHOW)

    Kernel_K = knext.IntParameter(
        "Kernel K",
        "The number of nearest neighbors to use for determining the bandwidth in kernel weights. Defaults to 12.",
        12,
    ).rule(knext.OneOf(category, ["Kernel"]), knext.Effect.SHOW)

    Kernel_type = knext.StringParameter(
        "Kernel type",
        "The type of kernel to use in constructing kernel weights. Defaults to 'triangular' ",
        "triangular",
        enum=[
            "gaussian",
            "quadratic",
            "quartic",
            "triangular",
            "uniform",
        ],
    ).rule(knext.OneOf(category, ["Kernel"]), knext.Effect.SHOW)

    Kernel_bandwidth = knext.StringParameter(
        "Kernel bandwidth",
        "The type of kernel bandwidth to use in constructing kernel weights. The bandwidth of the kernel. The default is fixed. If adaptive then bandwidth is adaptive across observations.",
        "Fixed",
        enum=[
            "Adaptive",
            "Fixed",
        ],
    ).rule(knext.OneOf(category, ["Kernel"]), knext.Effect.SHOW)

    Your_own_matrix_local_path = knext.StringParameter(
        "Get spatial weights matrix from file",
        """The file path of a user-defined spatial weights matrix in CSV format. Defaults to ''.
        Please enter the path of the spatial weights matrix in CSV format in the following options.
        The weights matrix must be in matrix format and in the order of the samples. """,
        "",
    ).rule(
        knext.OneOf(category, ["Get spatial weights matrix from file"]),
        knext.Effect.SHOW,
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)

        gdf.index = range(len(gdf))
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting projection...")

        import libpysal

        if self.category == "Rook":
            w = libpysal.weights.Rook.from_dataframe(gdf)
            wname = "Rook"
            w.transform = "r"
        if self.category == "Queen":
            w = libpysal.weights.Queen.from_dataframe(gdf)
            wname = "Queen"
            w.transform = "r"
        if self.category == "Inverse Distance":
            w = libpysal.weights.DistanceBand.from_dataframe(
                gdf, self.Threshold, alpha=-1 * self.order, binary=False
            )
            wname = "Inverse Distance"
            w.transform = "r"
        if self.category == "Binary Distance Band":
            import util.projection as kproj

            crs = gdf.crs
            if (crs is not None) and (kproj.is_geographic(crs)):
                gdf = gdf.to_crs("EPSG:3857")
            w = libpysal.weights.DistanceBand.from_dataframe(
                gdf, self.Threshold, binary=True
            )
            wname = "Binary Distance Band"
            w.transform = "r"
        if self.category == "K nearest":
            w = libpysal.weights.KNN.from_dataframe(gdf, k=self.Nearest_k)
            wname = "K nearest"
            w.transform = "r"
        if self.category == "Lattice":
            w = libpysal.weights.lat2W(nrows=self.Rows, ncols=self.Columns, rook=True)
            wname = "Lattice"
            w.transform = "r"
        if self.order != 1:
            w = libpysal.weights.higher_order(w, self.order - 1)
            w.transform = "r"

        import numpy as np

        if self.category == "Get spatial weights matrix from file":
            import pandas as pd
            import numpy as np

            z = pd.read_csv(self.Your_own_matrix_local_path, header=None)
            zz = np.array(z)

            import scipy.sparse

            sparse = scipy.sparse.csr_matrix(zz)

            from libpysal.weights import WSP

            w = WSP(sparse)
            wname = "Get spatial weights matrix from file"

        if self.category == "Kernel":
            bd = False
            if self.Kernel_bandwidth == "Fixed":
                bd = True
            w = libpysal.weights.Kernel.from_dataframe(
                gdf, fixed=bd, k=self.Kernel_K, function=self.Kernel_type
            )
        out = w.to_adjlist(drop_islands=False)

        if "none" not in str(self.id_col).lower():
            # get index id map
            id_map = gdf[self.id_col].to_dict()
            out["focal"] = out["focal"].map(id_map)
            out["neighbor"] = out["neighbor"].map(id_map)
        exec_context.set_progress(
            0.1, "Constructs a contiguity spatial weights matrix done"
        )

        # focal and neighbor should always be int
        out["focal"] = out["focal"].astype(np.int32)
        out["neighbor"] = out["neighbor"].astype(np.int32)

        return knext.Table.from_pandas(out)
