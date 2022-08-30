import logging

import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/geo",
    level_id="overlay",
    name="Spatial Overlay",
    description="Spatial overlay nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/overlay.png",
)


@knext.node(name="Spatial Join", node_type=knext.NodeType.MANIPULATOR, icon_path="icons/icon.png", category=category)
@knext.input_table(name="Left geo table", description="Left table with geometry column to join on")
@knext.input_table(name="Right geo table", description="Right table with geometry column to join on")
@knext.output_table(name="Joined geo table", description="Joined geo table")
class SpatialJoinNode:
    """
    This node will merge the left (top) and the right (bottom) table based on their spatial relationship of the two 
    selected columns to one another.
    """
    left_geo_col = knext.ColumnParameter(
        "Left geometry column",
        "Select the geometry column from the left (top) input table to join on.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False
    )

    right_geo_col = knext.ColumnParameter(
        "Right  geometry column",
        "Select the geometry column from the right (bottom) input table to join on.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False
    )

    join_mode = knext.StringParameter(
        label="Join mode",
        description='''<p>Available join modes are: 
<ul>
<li><b>inner:</b> retains only matching rows from both input tables.</li>
<li><b>left:</b> retains all rows form the left and only matching rows from the right input tables.</li>
<li><b>right:</b> retains all rows form the right and only matching rows from the left input tables.</li>
</ul></p>
''',
        default_value="inner",
        enum=["inner", "left", "right"],
    )

    match_mode = knext.StringParameter(
        label="Match mode",
        description='''<p>Available modes are: 
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
''',
        default_value="intersects",
        enum=["contains", "contains_properly", "covers", "crosses", "intersects", "overlaps", "touches", "within"],
    )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        knut.column_exists(self.left_geo_col, left_input_schema)
        knut.column_exists(self.right_geo_col, left_input_schema)
        #TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(right_input.to_pandas(), geometry=self.right_geo_col)
        knut.check_canceled(exec_context)
        gdf = left_gdf.sjoin(right_gdf, how=self.join_mode, predicate=self.match_mode)
        return knext.Table.from_pandas(gdf)
