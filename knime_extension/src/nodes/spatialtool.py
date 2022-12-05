# lingbo
import logging
from typing import Callable
from typing_extensions import Self
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import numpy as np
from math import radians, cos, sin, asin, sqrt  # For Haversine Distance
from shapely.geometry import Polygon  # For Grid

LOGGER = logging.getLogger(__name__)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialTool/"

__category = knext.category(
    path="/community/geo",
    level_id="spatialtool",
    name="Spatial Manipulation",
    description="Geospatial manipulation nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialToolCategory.png",
    after="calculation",
)


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


############################################
# Buffer
############################################


@knext.node(
    name="Buffer",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Buffer.png",
    category=__category,
    after="",
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
# Dissolve
############################################


@knext.node(
    name="Dissolve",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Dissolve.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column.",
)
@knut.geo_node_description(
    short_description="This node aggregate geometries based on group id (string column) and only keep the two column.",
    description="""This node aggregate geometries based on group id (string column) and only keep the two column with 
    GeoDataFrame.dissolve, which only dissolve geometries here. For grouping attribute values of other than geometry, 
    the [GroupBy node](https://kni.me/n/5stmXk6zY_ORA4bC) in KNIME is a complementary tool.
    """,
    references={
        "GeoDataFrame.dissolve": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.dissolve.html#geopandas.GeoDataFrame.dissolve",
    },
)
class DissolveNode:
    """
    This node aggregate geometries based on a group id (string column) and only keeps the id and geometry column.
    """

    geo_col = knut.geo_col_parameter()

    dissolve_col = knext.ColumnParameter(
        "Dissolve column",
        "Select the dissolve column as group id.",
        # Allow only string columns
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting dissolve...")
        gdf_dissolve = gdf.dissolve(self.dissolve_col, as_index=False)
        gdf = gdf_dissolve[[self.dissolve_col, self.geo_col]]
        exec_context.set_progress(0.1, "Dissolve done")
        LOGGER.debug(
            "Feature geometry " + self.geo_col + "dissolved by " + self.dissolve_col
        )
        return knext.Table.from_pandas(gdf)


############################################
# Spatial Join
############################################


@knext.node(
    name="Spatial Join",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "SpatialJoin.png",
    category=__category,
    after="",
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
    description="""This node will merge the left (top) and the right (bottom) table based on their spatial relationship 
    of the two selected columns to one another.
    Both layers must be in the same Coordinate Reference System (CRS),otherwise,the CRS of right table will be 
    transformed to that of the left table.
    """,
    references={
        "Spatial Joins": "https://geopandas.org/en/stable/gallery/spatial_joins.html",
        "Merging Data": "https://geopandas.org/en/stable/docs/user_guide/mergingdata.html#binary-predicate-joins",
        "sjoin": "https://geopandas.org/en/stable/docs/reference/api/geopandas.sjoin.html",
        "Predicates": "https://shapely.readthedocs.io/en/latest/manual.html#binary-predicates",
    },
)
class SpatialJoinNode:
    class MatchModes(knext.EnumParameterOptions):
        CONTAINS = (
            "Contains",
            """Matches if no points of the right object lies in the exterior of the left object and at least one point 
            of the interior of right object lies in the interior of the left object.""",
        )
        CONTAINS_PROPERLY = (
            "Contains properly",
            "Matches if the right object is completely inside the left object, with no common boundary points.",
        )
        COVERS = (
            "Covers",
            "Matches if every point of the right object is a point on the interior or boundary of the left object.",
        )
        CROSSES = (
            "Crosses",
            """Matches if the interior of the left object intersects the interior of the right object but does not 
            contain it, and the dimension of the intersection is less than the dimension of the one or the other.""",
        )
        INTERSECTS = (
            "Intersects",
            "Matches if the boundary or interior of the left object intersect in any way with those of the right object.",
        )
        OVERLAPS = (
            "Overlaps",
            """Matches if the geometries have more than one but not all points in common, have the same dimension, and 
            the intersection of the interiors of the geometries has the same dimension as the geometries themselves.""",
        )
        TOUCHES = (
            "Touches",
            """Matches if the objects have at least one point in common and their interiors do not intersect with any 
            part of the other.""",
        )
        WITHIN = (
            "Within",
            """Matches if the left object's boundary and interior intersect only with the interior of the right object 
            (not its boundary or exterior).""",
        )

        @classmethod
        def get_default(cls):
            return cls.INTERSECTS

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

    match_mode = knext.EnumParameter(
        label="Match mode",
        description="Defines which predicate is used for matching.",
        default_value=MatchModes.get_default().name,
        enum=MatchModes,
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
        right_gdf.to_crs(left_gdf.crs, inplace=True)
        gdf = left_gdf.sjoin(
            right_gdf, how=self.join_mode.lower(), predicate=self.match_mode.lower()
        )
        # reset the index since it might contain duplicates after joining
        gdf.reset_index(drop=True, inplace=True)
        # drop additional index columns if they exist
        gdf.drop(["index_right", "index_left"], axis=1, errors="ignore", inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# Nearest Join
############################################


@knext.node(
    name="Nearest Join",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "NearestJoin.png",
    category=__category,
    after="",
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
        description=knut.DEF_CRS_DESCRIPTION,
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
        return knext.Table.from_pandas(gdf)


############################################
# Clip
############################################


@knext.node(
    name="Clip",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Clip.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Left geo table",
    description="Left table with geometry column to be clipped",
)
@knext.input_table(
    name="Right geo table",
    description="Right table with geometry column to clip",
)
@knext.output_table(
    name="Joined geo table",
    description="Joined geo table",
)
@knut.geo_node_description(
    short_description="This node will clip target geometries to the mask extent.",
    description="""This node will clip target geometries to the mask extent.
    Both layers must be in the same Coordinate Reference System (CRS),otherwise, the CRS of the right table will be
    transformed to that of the left table. The geometries will be clipped to the full extent of the clip object.
    If there are multiple polygons in the mask geometry column, the geometries in the target geometry column 
    will be clipped to the total boundary of all mask polygons.
    """,
    references={
        "Clip": "https://geopandas.org/en/stable/docs/reference/api/geopandas.clip.html",
    },
)
class ClipNode:
    left_geo_col = knext.ColumnParameter(
        "Target geometry column",
        "Select the geometry column to be clipped to mask.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    right_geo_col = knext.ColumnParameter(
        "Mask  geometry column",
        "Select the geometry column used to clip the target geometry.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
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
        right_gdf.to_crs(left_gdf.crs, inplace=True)
        try:
            gdf_clipnew = gp.clip(left_gdf, right_gdf, keep_geom_type=True)
            return knext.Table.from_pandas(gdf_clipnew)
        except:
            raise ValueError("Improper Mask Geometry")


############################################
# Overlay
############################################


@knext.node(
    name="Overlay",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "SpatialOverlay.png",
    category=__category,
    after="",
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
    short_description="Performs spatial overlay between two geometries.",
    description="""This node will perform spatial overlay between two geometries.
    Currently only supports input tables with uniform geometry types,
    i.e. containing only (Multi)Polygons, or only (Multi)Points, or a combination of (Multi)LineString and 
    LinearRing shapes. Implements several methods that are all effectively subsets of the union.
    """,
    references={
        "Set-Operations with Overlay": "https://geopandas.org/en/stable/docs/user_guide/set_operations.html",
        "overlay": "https://geopandas.org/en/stable/docs/reference/api/geopandas.overlay.html",
    },
)
class OverlayNode:
    class OverlayModes(knext.EnumParameterOptions):
        INTERSECTION = (
            "Intersection",
            "Returns a representation of the intersection of the two geometries.",
        )
        UNION = (
            "Union",
            """Returns a combination of the subset geometry in the left input table not in the right one, and the 
            intersection geometry of the two input tables.""",
        )
        IDENTITY = (
            "Identity",
            "Retains all rows form the right and only matching rows from the left input tables.",
        )
        SYMMETRIC_DIFFERENCE = (
            "Symmetric difference",
            """Returns a combination of the subset geometry in the left input table not in the right one, and the subset 
            geometry in the right one not in the left one.""",
        )
        DIFFERENCE = (
            "Difference",
            "Returns the subset geometry in the left input table not in the right one.",
        )

        @classmethod
        def get_default(cls):
            return cls.INTERSECTION

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

    overlay_mode = knext.EnumParameter(
        label="Overlay mode",
        description="Determines the overlay mode to use.",
        default_value=OverlayModes.get_default().name,
        enum=OverlayModes,
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
        right_gdf.to_crs(left_gdf.crs, inplace=True)
        gdf = gp.overlay(left_gdf, right_gdf, how=self.overlay_mode.lower())
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# Euclidean Distance
############################################


@knext.node(
    name="Euclidean Distance",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "EuclideanDistance.png",
    category=__category,
    after="",
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
        description=knut.DEF_CRS_DESCRIPTION,
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
# Multiple Ring Buffer
############################################


@knext.node(
    name="Multiple Ring Buffer",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "MultipleRingBuffer.png",
    category=__category,
    after="",
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
        description=knut.DEF_CRS_DESCRIPTION,
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gp.GeoDataFrame(geometry=gdf.geometry)
        gdf.to_crs(self.crs_info, inplace=True)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting buffering...")
        # transfrom string list to number
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
        return knext.Table.from_pandas(gdf0)


############################################
# Simplify
############################################


@knext.node(
    name="Simplify",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Simplify.png",
    category=__category,
    after="",
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
# Create Grid Node
############################################


@knext.node(
    name="Create Grid",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "CreateGrid.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table of Create Grid",
)
@knext.output_table(
    name="Output Table",
    description="Output table of Create Grid",
)
class CreateGrid:
    """Create Grid node.
    This node creates an even grid like geospatial object.
    """

    geo_col = knut.geo_col_parameter()

    grid_length = knext.IntParameter(
        "Grid Length",
        "The length in meters of the grid. ",
        default_value=100,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        xmin, ymin, xmax, ymax = gdf.total_bounds
        width = self.grid_length
        height = self.grid_length
        rows = int(np.ceil((ymax - ymin) / height))
        cols = int(np.ceil((xmax - xmin) / width))
        XleftOrigin = xmin
        XrightOrigin = xmin + width
        YtopOrigin = ymax
        YbottomOrigin = ymax - height
        polygons = []
        for i in range(cols):
            Ytop = YtopOrigin
            Ybottom = YbottomOrigin
            for j in range(rows):
                polygons.append(
                    Polygon(
                        [
                            (XleftOrigin, Ytop),
                            (XrightOrigin, Ytop),
                            (XrightOrigin, Ybottom),
                            (XleftOrigin, Ybottom),
                        ]
                    )
                )
                Ytop = Ytop - height
                Ybottom = Ybottom - height
            XleftOrigin = XleftOrigin + width
            XrightOrigin = XrightOrigin + width

        grid = gp.GeoDataFrame({"geometry": polygons}, crs=gdf.crs)

        return knext.Table.from_pandas(grid)


############################################
# Get Geodesic Haversine Distance
############################################
@knext.node(
    name="Haversine Distance",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "HaversineDist.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for haversine distance calculation",
)
@knext.output_table(
    name="Output Table",
    description="Output table containing haversine distance",
)
class HaversineDistGrid:
    """Calculates the Haversine Distance.
    Calculates the [Haversine Distance](https://en.wikipedia.org/wiki/Haversine_formula).
    """

    Lon1 = knext.ColumnParameter(
        "The first longitude column",
        "The column containing the first longitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Lat1 = knext.ColumnParameter(
        "The first latitude column",
        "The column containing the first Latitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Lon2 = knext.ColumnParameter(
        "The second longitude column",
        "The column containing the second  longitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Lat2 = knext.ColumnParameter(
        "The second latitude column",
        "The column containing the second  Latitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema):
        return input_schema.append(knext.Column(knext.double(), name="HDist"))

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        def HaversineDist(x1, y1, x2, y2):
            x1 = radians(x1)
            x2 = radians(x2)
            y1 = radians(y1)
            y2 = radians(y2)
            # Haversine formula
            Dx = x2 - x1
            Dy = y2 - y1
            P = sin(Dy / 2) ** 2 + cos(y1) * cos(y2) * sin(Dx / 2) ** 2
            Q = 2 * asin(sqrt(P))
            # The earth's radius in kilometers.
            EarthR_km = 6371
            # Then we'll compute the outcome.
            return Q * EarthR_km

        df = input_table.to_pandas()
        df["HDist"] = df.apply(
            lambda x: HaversineDist(
                x[self.Lon1], x[self.Lat1], x[self.Lon2], x[self.Lat2]
            ),
            axis=1,
        )
        return knext.Table.from_pandas(df)
