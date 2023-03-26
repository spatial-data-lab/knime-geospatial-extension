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
        gdf[
            knut.get_unique_column_name("geometry", input.schema)
        ] = gdf.geometry.simplify(self.simplifydist)
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
        return knut.to_table(gdf, exec_context)
