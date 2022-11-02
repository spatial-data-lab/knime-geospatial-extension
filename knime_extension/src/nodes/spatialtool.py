# lingbo
import logging
from typing import Callable
from typing_extensions import Self
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import numpy as np

LOGGER = logging.getLogger(__name__)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialTool/"

__category = knext.category(
    path="/geo",
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class BufferNode:
    """
    This node generate polygons representing all points within a given distance of each geometric object.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    bufferdist = knext.DoubleParameter(
        "Buffer Distance", "The buffer distance for geometry ", 1000.0
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting buffering...")
        gdf["geometry"] = gdf.buffer(self.bufferdist)
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class DissolveNode:
    """
    This node aggregate geometries based on group id (string column) and only keep the two column.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dissolve_col = knext.ColumnParameter(
        "Dissolve column",
        "Select the dissolve column as group id.",
        # Allow only string columns
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
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
        "Right  geometry column",
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
        right_gdf = right_gdf.to_crs(left_gdf.crs)
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
    description="""This node will merge the left (top) and the right (bottom) table based on  the distance between 
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
        "Right  geometry column",
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
        "Maximum Distance",
        "Maximum distance within which to query for nearest geometry. Must be greater than 0 ",
        1000.0,
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
        right_gdf = right_gdf.to_crs(left_gdf.crs)
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
class ClipNode:
    """
    This node will clip target geometries to the mask extent.
    Both layers must be in the same Coordinate Reference System (CRS),otherwise, the CRS of right table will be
    transformed to that of the left table. The gdf will be clipped to the full extent of the clip object.
    If there are multiple polygons in Mask geometry, data from Target geometry will be clipped to the total boundary
    of all polygons in Mask.
    """

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
        right_gdf = right_gdf.to_crs(left_gdf.crs)
        # gdf_clip = gp.clip(left_gdf, right_gdf)
        gdf_clipnew = gp.clip(left_gdf, right_gdf, keep_geom_type=True)
        # gdf_clipnew = gp.GeoDataFrame(geometry=gdf_clip.geometry)
        # =gdf_clip.reset_index(drop=True)
        return knext.Table.from_pandas(gdf_clipnew)


############################################
# Overlay
############################################


@knext.node(
    name="Overlay",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "SpatialOverlay.png",
    category=__category,
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
    Currently only supports data GeoDataFrames with uniform geometry types,
    i.e. containing only (Multi)Polygons, or only (Multi)Points, or a combination of (Multi)LineString and LinearRing shapes.
    Implements several methods that are all effectively subsets of the union.
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
        "Right  geometry column",
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
        right_gdf = right_gdf.to_crs(left_gdf.crs)
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
    name="Geo table Distance",
    description="Euclidean distance between objects",
)
class EuclideanDistanceNode:
    """
    This node will calculate the distance between two geometries.
    """

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
        "Right  geometry column",
        "Select the geometry column from the right (bottom) input table to calculate.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    crs_info = knext.StringParameter(
        label="CRS for distance calculation",
        description="Input the CRS to use",
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
        right_gdf = right_gdf.to_crs(self.crs_info)
        left_gdf = left_gdf.to_crs(self.crs_info)
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to buffer",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed table by Multiple Ring Buffer",
)
class MultiRingBufferNode:
    """
    This node generate multiple polygons with a series  distances of each geometric object.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    bufferdist = knext.StringParameter(
        "Serial Buffer Distances with coma",
        "The buffer distances for geometry ",
        "10,20,30",
    )

    bufferunit = knext.StringParameter(
        label="Serial Buffer Distances",
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
        description="Input the CRS to use",
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
        gdf = gdf.to_crs(self.crs_info)
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
        c1 = gp.GeoDataFrame(geometry=gdf.buffer(bufferlist[0]))
        c2 = gp.GeoDataFrame(geometry=gdf.buffer(bufferlist[1]))
        gdf0 = gp.overlay(c1, c2, how="union")
        if len(bufferlist) > 2:
            # Construct all other rings by loop
            for i in range(2, len(bufferlist)):
                ci = gp.GeoDataFrame(geometry=gdf.buffer(bufferlist[i]))
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to simplify",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class SimplifyNode:
    """
    This node generate the smallest convex Polygon containing all the points in each geometry.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    simplifydist = knext.DoubleParameter(
        label="Simplification tolerance",
        description="The simplification tolerance distances for geometry ",
        default_value="1",
    )

    crs_info = knext.StringParameter(
        label="CRS for simplification ",
        description="Input the CRS to use",
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gdf.to_crs(self.crs_info)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf["geometry"] = gdf.geometry.simplify(self.simplifydist)
        gdf = gdf.reset_index(drop=True)
        exec_context.set_progress(0.1, "Transformation done")
        LOGGER.debug("Feature Simplified")
        return knext.Table.from_pandas(gdf)
