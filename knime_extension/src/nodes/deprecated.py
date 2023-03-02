# lingbo
import geopandas as gp
import logging
import knime_extension as knext
import util.knime_utils as knut

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
