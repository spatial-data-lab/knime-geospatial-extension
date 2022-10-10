import logging
from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import geospatial_types as gt
import numpy as np

LOGGER = logging.getLogger(__name__)


__category = knext.category(
    path="/geo",
    level_id="spatialtool",
    name="Spatial Tool",
    description="Geospatial Tool nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialToolCategory.png",
)

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
        return input_schema_1

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
        return input_schema_1

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
# Bounding Box
############################################

@knext.node(
    name="Bounding Box",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "BoundingBox.png",
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
class BoundingBoxNode:
    """
    This node generate rectangles representing the envelope of each geometry.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf["geometry"] = gdf.envelope
        exec_context.set_progress(0.1, "Transformation done")
        LOGGER.debug("Feature converted to BoundingBox")
        return knext.Table.from_pandas(gdf)


############################################
# Convex Hull
############################################

@knext.node(
    name="Convex Hull",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "ConvexHull.png",
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
class ConvexHullNode:
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

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf["geometry"] = gdf.convex_hull
        exec_context.set_progress(0.1, "Transformation done")
        LOGGER.debug("Feature converted to ConvexHull")
        return knext.Table.from_pandas(gdf)


############################################
# Unary Union
############################################

@knext.node(
    name="Unary Union",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "UnaryUnion.png",
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
class UnaryUnionNode:
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

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf_union = gdf.unary_union
        gdfunion = gp.GeoDataFrame(geometry=gp.GeoSeries(gdf_union), crs=gdf.crs)
        exec_context.set_progress(0.1, "Transformation done")
        LOGGER.debug("Feature converted to ConvexHull")
        return knext.Table.from_pandas(gdfunion)


############################################
# Spatial Join
############################################

@knext.node(
    name="SpatialJoin",
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
class SpatialJoinNode:
    """
    This node will merge the left (top) and the right (bottom) table based on their spatial relationship of the two
    selected columns to one another.
    Both layers must be in the same Coordinate Reference System (CRS),otherwise,the CRS of right table will be tranformed to that of the left table.
    """

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

    join_mode = knext.StringParameter(
        label="Join mode",
        description="""<p>Available join modes are: 
<ul>
<li><b>inner:</b> retains only matching rows from both input tables.</li>
<li><b>left:</b> retains all rows form the left and only matching rows from the right input tables.</li>
<li><b>right:</b> retains all rows form the right and only matching rows from the left input tables.</li>
</ul></p>
""",
        default_value="inner",
        enum=["inner", "left", "right"],
    )

    match_mode = knext.StringParameter(
        label="Match mode",
        description="""<p>Available modes are: 
<ul>
<li><b>contains:</b> Matches if no points of the right object lies in the exterior of the left object and at least 
one point of the interior of right object lies in the interior of the left object.</li>
<li><b>contains_properly:</b> Matches if the right object is completely inside the left object, with no common boundary 
points.</li>
<li><b>covers:</b> Matches if every point of the right object is a point on the interior or boundary of the left 
object.</li>
<li><b>crosses:</b> Matches if the interior of the left object intersects the interior of the right object but does 
not contain it, and the dimension of the intersection is less than the dimension of the one or the other.</li>
<li><b>intersects:</b> Matches if the boundary or interior of the left object intersect in any way with those of the 
right object.</li>
<li><b>overlaps:</b> Matches if the geometries have more than one but not all points in common, have the same dimension, 
and the intersection of the interiors of the geometries has the same dimension as the geometries themselves.</li>
<li><b>touches:</b> Matches if the objects have at least one point in common and their interiors do not intersect with 
any part of the other.</li>
<li><b>within:</b> Matches if the left object's boundary and interior intersect only with the interior of the right 
object (not its boundary or exterior).</li>
</ul></p>
""",
        default_value="intersects",
        enum=[
            "contains",
            "contains_properly",
            "covers",
            "crosses",
            "intersects",
            "overlaps",
            "touches",
            "within",
        ],
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        knut.column_exists(self.left_geo_col, left_input_schema)
        knut.column_exists(self.right_geo_col, left_input_schema)
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(
            right_input.to_pandas(), geometry=self.right_geo_col
        )
        knut.check_canceled(exec_context)
        right_gdf = right_gdf.to_crs(left_gdf.crs)
        gdf = left_gdf.sjoin(right_gdf, how=self.join_mode, predicate=self.match_mode)
        gdf = gdf.reset_index(drop=True).drop(columns='index_right')        
        return knext.Table.from_pandas(gdf)


############################################
# Nearest Join
############################################

@knext.node(
    name="NearestJoin",
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
class NearestJoinNode:
    """
    This node will merge the left (top) and the right (bottom) table based on  the distance between their geometries of the two
    selected columns to one another. Distance is calculated in CRS units and is returned in the column NearDist.
    Both layers must be in the same Coordinate Reference System (CRS),otherwise,the CRS of right table will be tranformed to that of the left table.
    """

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

    join_mode = knext.StringParameter(
        label="Join mode",
        description="""<p>Available join modes are: 
<ul>
<li><b>inner:</b> use intersection of keys from both dfs; retain only left_df geometry column.</li>
<li><b>left:</b> use keys from left_df; retain only left_df geometry column.</li>
<li><b>right:</b> use keys from right_df; retain only right_df geometry column.</li>
</ul></p>
""",
        default_value="inner",
        enum=["inner", "left", "right"],
    )

    maxdist = knext.DoubleParameter(
        "Maximum Distance",
        "Maximum distance within which to query for nearest geometry. Must be greater than 0 ",
        1000.0,
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        knut.column_exists(self.left_geo_col, left_input_schema)
        knut.column_exists(self.right_geo_col, left_input_schema)
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
            how=self.join_mode,
            max_distance=self.maxdist,
            distance_col="NearDist",
            lsuffix="1",
            rsuffix="2",
        )
        gdf = gdf.reset_index(drop=True).drop(columns='index_2')        
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
    Both layers must be in the same Coordinate Reference System (CRS),otherwise, the CRS of right table will be tranformed to that of the left table.
    The gdf will be clipped to the full extent of the clip object.
    If there are multiple polygons in Mask geometry, data from Target geometry will be clipped to the total boundary of all polygons in Mask.
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
        knut.column_exists(self.left_geo_col, left_input_schema)
        knut.column_exists(self.right_geo_col, right_input_schema)
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(
            right_input.to_pandas(), geometry=self.right_geo_col
        )
        knut.check_canceled(exec_context)
        right_gdf = right_gdf.to_crs(left_gdf.crs)
        #gdf_clip = gp.clip(left_gdf, right_gdf)
        gdf_clipnew =gp.clip(left_gdf, right_gdf,keep_geom_type=True)
         #gdf_clipnew = gp.GeoDataFrame(geometry=gdf_clip.geometry)
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
class OverlayNode:
    """
    This node will perform spatial overlay between two geometries.
    Currently only supports data GeoDataFrames with uniform geometry types,
    i.e. containing only (Multi)Polygons, or only (Multi)Points, or a combination of (Multi)LineString and LinearRing shapes.
    Implements several methods that are all effectively subsets of the union.
    """

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

    overlay_mode = knext.StringParameter(
        label="Overlay mode",
        description="""<p>Available overlay modes are: 
<ul>
<li><b>intersection:</b> returns a representation of the intersection of the two geometries</li>
<li><b>union:</b> returns a combiantion of the subset geometry in the left input table not in the right one, and the intersection geometry of the two input tables.</li>
<li><b>identity:</b> retains all rows form the right and only matching rows from the left input tables.</li>
<li><b>symmetric_difference:</b> returns a combiantion of the subset geometry in the left input table not in the right one, and the subset geometry in the right one not in the left one.</li>
<li><b>difference:</b> returns the subset geometry in the left input table not in the right one.</li>
</ul></p>
""",
        default_value="intersection",
        enum=[
            "intersection",
            "union",
            "identity",
            "symmetric_difference",
            "difference",
        ],
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        knut.column_exists(self.left_geo_col, left_input_schema)
        knut.column_exists(self.right_geo_col, left_input_schema)
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(
            right_input.to_pandas(), geometry=self.right_geo_col
        )
        knut.check_canceled(exec_context)
        right_gdf = right_gdf.to_crs(left_gdf.crs)
        gdf = gp.overlay(left_gdf, right_gdf, how=self.overlay_mode)
        gdf = gdf.reset_index(drop=True)
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
        "Select the geometry column from the left (top) input table to calcualte.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    right_geo_col = knext.ColumnParameter(
        "Right  geometry column",
        "Select the geometry column from the right (bottom) input table to calcualte.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    crs_info = knext.StringParameter(
        label="CRS for ditance calculation",
        description="Input the CRS to use",
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        knut.column_exists(self.left_geo_col, left_input_schema)
        knut.column_exists(self.right_geo_col, left_input_schema)
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(right_input.to_pandas(), geometry=self.right_geo_col)
        knut.check_canceled(exec_context)
        right_gdf = right_gdf.to_crs(self.crs_info)
        left_gdf = left_gdf.to_crs(self.crs_info)
        #left_gdf['LID'] = range(1,(left_gdf.shape[0]+1))
        #right_gdf['RID'] = range(1,(right_gdf.shape[0]+1))
        mergedf=left_gdf.merge(right_gdf, how='cross')
        mergedf_x=gp.GeoDataFrame(geometry=mergedf['geometry_x'])
        mergedf_y=gp.GeoDataFrame(geometry=mergedf['geometry_y'])
        mergedf['EuDist']=mergedf_x.distance(mergedf_y, align=False)        
        mergedf= mergedf.reset_index(drop=True)
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
        "Serial Buffer Distances with coma", "The buffer distances for geometry ", "10,20,30"
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
        label="CRS for buffering ditance calculation",
        description="Input the CRS to use",
        default_value="EPSG:3857",
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gp.GeoDataFrame(geometry=gdf.geometry)
        gdf = gdf.to_crs(self.crs_info)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting buffering...")
        # transfrom string list to number 
        bufferlist= np.array(self.bufferdist.split(","), dtype=np.int64)
        if self.bufferunit=="Meter" :
            bufferlist=bufferlist
        elif self.bufferunit=="KiloMeter" :
            bufferlist=bufferlist*1000
        else:
            bufferlist=bufferlist*1609.34 
        # sort list
        bufferlist=bufferlist.tolist() 
        bufferlist.sort()   
        c1 = gp.GeoDataFrame(geometry=gdf.buffer(bufferlist[0]))
        c2 = gp.GeoDataFrame(geometry=gdf.buffer(bufferlist[1]))
        gdf0 = gp.overlay(c1 , c2 , how='union')
        if len(bufferlist)>2:
            # Construct all other rings by loop
            for i in range(2,len(bufferlist)):
                ci = gp.GeoDataFrame(geometry=gdf.buffer(bufferlist[i]));
                gdf0 = gp.overlay(gdf0, ci, how='union')
        # Add ring radius values as a new column    
        gdf0['dist']=bufferlist
        gdf0=gdf0.reset_index(drop=True)
        exec_context.set_progress(0.1, "Buffering done")
        return knext.Table.from_pandas(gdf0)

