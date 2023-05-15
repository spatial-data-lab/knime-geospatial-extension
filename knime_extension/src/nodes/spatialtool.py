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
    level_id="spatialtool",
    name="Spatial Manipulation",
    description="Geospatial manipulation nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialToolCategory.png",
    after="calculation",
)


class _JoinModes(knext.EnumParameterOptions):
    INNER = (
        "Inner",
        "Use intersection of index values and retain only the geometry column from the left (top) input table.",
    )
    LEFT = (
        "Left",
        "Use the index and retain the geometry column from the left (top) input table.",
    )
    RIGHT = (
        "Right",
        "Use the index and retain the geometry column from the right (bottom) input table.",
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
    short_description="Generates a buffer based on a given distance.",
    description="""This node generates a buffer with the given distance for each geometric object.""",
    references={
        "GeoSeries.buffer": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.buffer.html",
        "Shapley object.buffer": "https://shapely.readthedocs.io/en/latest/manual.html#object.buffer",
    },
)
class BufferNode2:
    """
    This node aggregate generate buffer zone based on a given distance.
    """

    geo_col = knut.geo_col_parameter()

    distance = kproj.Distance.get_distance_parameter()

    unit = kproj.Distance.get_unit_parameter()

    keep_input_crs = kproj.Distance.get_keep_input_crs_parameter()

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name="Buffer",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.APPEND.name
        self.result_settings.new_column_name = "Buffer"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            knut.TYPE_POLYGON,
        )

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = knut.load_geo_data_frame(input, self.geo_col, exec_context)
        helper = kproj.Distance(self.unit, self.keep_input_crs)
        projected_gdf = helper.pre_processing(exec_context, gdf, False)
        new_distance = helper.convert_input_distance(self.distance)
        knut.check_canceled(exec_context)
        exec_context.set_progress(0.6, "Computing buffer")
        gdf = self.result_settings.get_computed_result_frame(
            exec_context,
            input.schema,
            projected_gdf,
            self.geo_col,
            lambda l: l.buffer(new_distance),
        )
        gdf = helper.post_processing(exec_context, gdf)
        exec_context.set_progress(1, "Buffering done")
        return knut.to_table(gdf, exec_context)


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

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = knut.load_geo_data_frame(input_table, self.geo_col, exec_context)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting dissolve...")
        gdf_dissolve = gdf.dissolve(self.dissolve_col, as_index=False)
        gdf = gdf_dissolve[[self.dissolve_col, self.geo_col]]
        exec_context.set_progress(1, "Dissolve done")
        return knut.to_table(gdf, exec_context)


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
    description="Left (top) table with geometry column to join on",
)
@knext.input_table(
    name="Right geo table",
    description="Right (bottom) table with geometry column to join on",
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
            """Matches if no point of the right object lies in the exterior of the left object and at least one point 
            of the interior of right object lies in the interior of the left object.""",
        )
        CONTAINS_CENTER_OF = (
            "Contains the center of",
            """Matches if the right object's representative point intersects only with the interior of the left object 
            (not its boundary or exterior).""",
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
        HAS_ITS_CENTER_IN = (
            "Has its center in",
            """Matches if the left object's representative point intersects only with the interior of the right object 
            (not its boundary or exterior).""",
        )
        INTERSECTS = (
            "Intersects",
            "Matches if the boundary or interior of the left object intersects in any way with those of the right object.",
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
        left_gdf = knut.load_geo_data_frame(left_input, self.left_geo_col, exec_context)
        right_gdf = knut.load_geo_data_frame(
            right_input, self.right_geo_col, exec_context
        )
        knut.check_canceled(exec_context)
        right_gdf.to_crs(left_gdf.crs, inplace=True)
        if self.match_mode not in [
            self.MatchModes.HAS_ITS_CENTER_IN.name,
            self.MatchModes.CONTAINS_CENTER_OF.name,
        ]:
            gdf = left_gdf.sjoin(
                right_gdf, how=self.join_mode.lower(), predicate=self.match_mode.lower()
            )
        elif self.match_mode == self.MatchModes.HAS_ITS_CENTER_IN.name:
            rep_center = knut.get_unique_column_name("rep_center", left_input.schema)
            left_gdf[rep_center] = left_gdf.representative_point()
            left_gdf_temp = left_gdf.set_geometry(rep_center)
            gdf = gp.sjoin(left_gdf_temp, right_gdf, predicate="within")
            gdf = gdf.set_geometry(self.left_geo_col).drop(columns=[rep_center])
        else:
            rep_center = knut.get_unique_column_name("rep_center", right_input.schema)
            right_gdf[rep_center] = right_gdf.representative_point()
            right_gdf_temp = right_gdf.set_geometry(rep_center).drop(
                columns=[self.right_geo_col]
            )
            gdf = gp.sjoin(left_gdf, right_gdf_temp, predicate="contains")
        # reset the index since it might contain duplicates after joining
        gdf.reset_index(drop=True, inplace=True)
        # drop additional index columns if they exist
        gdf.drop(["index_right", "index_left"], axis=1, errors="ignore", inplace=True)
        return knut.to_table(gdf, exec_context)


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
class NearestJoinNode2:
    left_geo_column = knext.ColumnParameter(
        "Left geometry column",
        "Select the geometry column from the left (top) input table to join on.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    right_geo_column = knext.ColumnParameter(
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
        description="The join mode specifies the type of join that will occur and which geometry column is "
        + "retained in the output table.",
        default_value=_JoinModes.get_default().name,
        enum=_JoinModes,
    )

    distance = kproj.Distance.get_distance_parameter(
        "Maximum distance",
        "Maximum distance within which to query for nearest geometry. Must be greater than 0 ",
        1000.0,
    )

    unit = kproj.Distance.get_unit_parameter()

    keep_input_crs = kproj.Distance.get_keep_input_crs_parameter(
        label="Keep CRS from left input table",
        description="If checked the CRS of the left input table is retained even if a re-projection was necessary "
        + "for the selected distance unit.",
    )

    # left_include_columns = knext.MultiColumnParameter(
    #     "Left columns",
    #     "Select columns which should be included in the join result from the left (top) input table.",
    #     port_index=0,
    # )

    # right_include_columns = knext.MultiColumnParameter(
    #     "Right columns",
    #     "Select columns which should be included in the join result from the right (bottom) input table.",
    #     port_index=1,
    # )

    def configure(self, configure_context, left_input_schema, right_input_schema):
        self.left_geo_col = knut.column_exists_or_preset(
            configure_context, self.left_geo_column, left_input_schema, knut.is_geo
        )
        self.right_geo_col = knut.column_exists_or_preset(
            configure_context, self.right_geo_column, right_input_schema, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        # keep only selected columns
        # left_include_columns = {self.left_geo_column}
        # if self.left_include_columns is not None:
        #     left_include_columns.update(self.left_include_columns)
        # filtered_left = left_input[list(left_include_columns)]
        # right_include_columns = {self.right_geo_column}
        # if self.right_include_columns is not None:
        #     right_include_columns.update(self.right_include_columns)
        # filtered_right = right_input[list(right_include_columns)]

        # left_gdf = knut.load_geo_data_frame(
        #     filtered_left, self.left_geo_column, exec_context
        # )
        # right_gdf = knut.load_geo_data_frame(
        #     filtered_right, self.right_geo_column, exec_context
        # )
        left_gdf = knut.load_geo_data_frame(
            left_input, self.left_geo_column, exec_context
        )
        right_gdf = knut.load_geo_data_frame(
            right_input, self.right_geo_column, exec_context
        )
        knut.check_canceled(exec_context)
        distance_helper = kproj.Distance(self.unit, self.keep_input_crs)
        distance_helper.pre_processing(exec_context, right_gdf, True)
        # process left last to keep its CRS as described in the node description
        distance_helper.pre_processing(exec_context, left_gdf, True)

        # left_include_columns.update(right_include_columns)
        # distance_col_name = knut.get_unique_name("Distance", left_include_columns)
        col_names = set(left_input.column_names)
        col_names.update(right_input.column_names)
        distance_col_name = knut.get_unique_name(
            kproj.DEFAULT_DISTANCE_COLUMN_NAME, col_names
        )

        left_suffix = "left"
        right_suffix = "right"

        gdf = gp.sjoin_nearest(
            left_gdf,
            right_gdf,
            how=self.join_mode.lower(),
            max_distance=distance_helper.convert_input_distance(self.distance),
            distance_col=distance_col_name,
            lsuffix=left_suffix,
            rsuffix=right_suffix,
        )

        # convert returned distance to selected distance unit
        knut.check_canceled(exec_context)
        exec_context.set_progress(
            0.8, f"Convert {distance_col_name} column to selected unit"
        )
        gdf[distance_col_name] = (
            gdf[distance_col_name] / distance_helper.get_distance_factor()
        )

        # reset the index since it might contain duplicates after joining
        gdf.reset_index(drop=True, inplace=True)
        # drop additional index columns if they exist
        gdf.drop(
            ["index_" + left_suffix, "index_" + right_suffix],
            axis=1,
            errors="ignore",
            inplace=True,
        )

        distance_helper.post_processing(exec_context, gdf, True)

        return knut.to_table(gdf, exec_context)


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
    description="Left (top) table with geometry column to be clipped",
)
@knext.input_table(
    name="Right geo table",
    description="Right (bottom) table with geometry column to clip",
)
@knext.output_table(
    name="Clipped geo table",
    description="Clipped geo table",
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
        "Clip vector data with GeoPandas": "https://geopandas.org/en/stable/gallery/plot_clip.html",
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

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Clipped",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.REPLACE.name
        self.result_settings.new_column_name = "Clipped"

    def configure(self, configure_context, left_input_schema, right_input_schema):
        self.left_geo_col = knut.column_exists_or_preset(
            configure_context, self.left_geo_col, left_input_schema, knut.is_geo
        )
        self.right_geo_col = knut.column_exists_or_preset(
            configure_context, self.right_geo_col, right_input_schema, knut.is_geo
        )

        return self.result_settings.get_result_schema(
            configure_context,
            left_input_schema,
            self.left_geo_col,
            knut.TYPE_GEO,
        )

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        left_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.left_geo_col)
        right_gdf = gp.GeoDataFrame(
            right_input.to_pandas(), geometry=self.right_geo_col
        )
        knut.check_canceled(exec_context)
        # ensure that both dataframe use the same projection
        right_gdf.to_crs(left_gdf.crs, inplace=True)
        try:
            gdf_clip = gp.clip(left_gdf, right_gdf, keep_geom_type=True)
        except:
            raise ValueError("Improper Mask Geometry")

        if self.result_settings.mode == knut.ResultSettingsMode.APPEND.name:
            left_gdf[self.result_settings.new_column_name] = gdf_clip[self.left_geo_col]
            gdf_clip = left_gdf
        return knut.to_table(gdf_clip, exec_context)


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
    description="Left (top) table with geometry column to join on",
)
@knext.input_table(
    name="Right geo table",
    description="Right (bottom) table with geometry column to join on",
)
@knext.output_table(
    name="Overlayed geo table",
    description="Overlayed geo table",
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
        "Overlay": "https://geopandas.org/en/stable/docs/reference/api/geopandas.overlay.html",
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

    keep_geom_type = knext.BoolParameter(
        label="Return only geometries of the same geometry type",
        description="""If selected, the node returns only geometries of the same geometry type as the geometry column 
        from the left (top) input table, otherwise it returns all resulting geometries.""",
        default_value=lambda v: True if v < knext.Version(1, 1, 0) else False,
        since_version="1.1.0",
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
        right_gdf.to_crs(left_gdf.crs, inplace=True)
        gdf = gp.overlay(
            left_gdf,
            right_gdf,
            how=self.overlay_mode.lower(),
            keep_geom_type=self.keep_geom_type,
        )
        gdf.reset_index(drop=True, inplace=True)
        return knut.to_table(gdf, exec_context)


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
    name="Origin geo table",
    description="Origin table with geometry column. ",
)
@knext.input_table(
    name="Destination geo table",
    description="Destination table with geometry column.",
)
@knext.output_table(
    name="Euclidean distance table",
    description="Euclidean distance between all origin and destination geometry objects.",
)
@knut.geo_node_description(
    short_description="This node will calculate the Euclidean distance between each origin and destination pair.",
    description="This node will calculate the Euclidean distance in the selected unit between each origin and destination pair.",
    references={
        "Distance": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.distance.html",
    },
)
class EuclideanDistanceNode2:
    o_geo_col = knext.ColumnParameter(
        "Origin geometry column",
        "Select the geometry column from the origin table to calculate.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    o_id_col = knext.ColumnParameter(
        "Origin ID column",
        """Select the column which contains for each origin a unique ID. The selected column will be returned
        in the result table and can be used to link back to the original data.""",
        port_index=0,
        include_row_key=False,
        include_none_column=False,
    )

    d_geo_col = knext.ColumnParameter(
        "Destination geometry column",
        "Select the geometry column from the destination table to calculate.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    d_id_col = knext.ColumnParameter(
        "Destination ID column",
        """Select the column which contains for each destination a unique ID. The selected column will be returned
        in the result table and can be used to link back to the original data.""",
        port_index=1,
        include_row_key=False,
        include_none_column=False,
    )

    unit = kproj.Distance.get_unit_parameter()

    __COL_ORIGIN = "Origin ID"
    __COL_DESTINATION = "Destination ID"
    __COL_DISTANCE = "Distance"

    def configure(self, configure_context, o_schema, d_schema):
        self.o_geo_col = knut.column_exists_or_preset(
            configure_context, self.o_geo_col, o_schema, knut.is_geo
        )
        self.d_geo_col = knut.column_exists_or_preset(
            configure_context, self.d_geo_col, d_schema, knut.is_geo
        )
        knut.column_exists(self.o_id_col, o_schema)
        o_id_type = o_schema[self.o_id_col].ktype
        knut.column_exists(self.d_id_col, d_schema)
        d_id_type = o_schema[self.d_id_col].ktype
        return knext.Schema.from_columns(
            [
                knext.Column(o_id_type, self.__COL_ORIGIN),
                knext.Column(d_id_type, self.__COL_DESTINATION),
                knext.Column(knext.double(), self.__COL_DISTANCE),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, o_input, d_input):
        o_gdf = gp.GeoDataFrame(
            o_input.to_pandas()[[self.o_geo_col, self.o_id_col]],
            geometry=self.o_geo_col,
        )
        d_gdf = gp.GeoDataFrame(
            d_input.to_pandas()[[self.d_geo_col, self.d_id_col]],
            geometry=self.d_geo_col,
        )
        # rename the id columns to origin and destination and the geometry to geometry
        o_gdf.rename(columns={self.o_id_col: self.__COL_ORIGIN}, inplace=True)
        o_gdf.rename_geometry("geometry", inplace=True)
        d_gdf.rename(columns={self.d_id_col: self.__COL_DESTINATION}, inplace=True)
        d_gdf.rename_geometry("geometry", inplace=True)

        helper = kproj.Distance(self.unit, True)
        helper.pre_processing(exec_context, o_gdf, True)
        helper.pre_processing(exec_context, d_gdf, True)

        # this could happen if the user selects to keep the input projection
        if o_gdf.crs != d_gdf.crs:
            d_gdf.to_crs(o_gdf.crs, inplace=True)

        merged_df = o_gdf.merge(d_gdf, how="cross")
        merged_df_x = gp.GeoDataFrame(geometry=merged_df["geometry_x"])
        merged_df_y = gp.GeoDataFrame(geometry=merged_df["geometry_y"])
        # compute and adjust distance if necessary
        merged_df[self.__COL_DISTANCE] = (
            merged_df_x.distance(merged_df_y, align=False)
            / helper.get_distance_factor()
        )
        merged_df = merged_df[
            [self.__COL_ORIGIN, self.__COL_DESTINATION, self.__COL_DISTANCE]
        ].reset_index(drop=True)

        # computes only the upper part of the distance matrix
        # data = []
        # idx = 1
        # length = len(o_gdf)
        # for o_id, o_geo in zip(o_gdf[self.__COL_ORIGIN], o_gdf[self.o_geo_col]):
        #     knut.check_canceled(exec_context)
        #     exec_context.set_progress(
        #         0.9 * idx / float(length), f"Processing origin row {idx} of {length}"
        #     )
        #     for d_id, d_geo in zip(
        #         d_gdf[self.__COL_DESTINATION], d_gdf[self.d_geo_col]
        #     ):
        #         data.append(
        #             (o_id, d_id, o_geo.distance(d_geo) / helper.get_distance_factor())
        #         )
        #     idx = idx + 1
        # import pandas as pd

        # merged_df = pd.DataFrame(
        #     data,
        #     columns=(self.__COL_ORIGIN, self.__COL_DESTINATION, self.__COL_DISTANCE),
        # )
        return knut.to_table(merged_df, exec_context)


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
    description="Table with geometry column to compute the buffers for.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed table with multiple ring buffers.",
)
@knut.geo_node_description(
    short_description="This node generates a buffer for each of the given distances for the input geometric object.",
    description="""This node generates a buffer for each of the given distances for the input geometric 
    object. If the input table contains multiple rows the node first computes the union of all geometries before 
    computing the buffers from the union.
    """,
    references={
        "GeoSeries.buffer": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.buffer.html",
        "Shapley object.buffer": "https://shapely.readthedocs.io/en/latest/manual.html#object.buffer",
    },
)
class MultipleRingBufferNode:
    geo_col = knut.geo_col_parameter()

    distance = knext.StringParameter(
        "Buffer distances (comma separated)",
        "Comma separated list of buffer distances in the selected distance unit.",
        "10,20,30",
        validator=kproj.string_distances_parser,
    )

    unit = kproj.Distance.get_unit_parameter()

    keep_input_crs = kproj.Distance.get_keep_input_crs_parameter()

    __BUFFER_COL_NAME = "Buffer"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )

        return knext.Schema(
            [knut.TYPE_GEO, knext.double()],
            [self.__BUFFER_COL_NAME, kproj.DEFAULT_DISTANCE_COLUMN_NAME],
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # parse and sort the buffer distances
        distances = kproj.string_distances_parser(self.distance)
        distances.sort()
        dist_length = len(distances)

        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        # compute union of all shape if the input table has multiple rows
        if gdf.shape[0] > 1:
            gdf_union = gdf.unary_union
            gdf_union = gp.GeoDataFrame(geometry=gp.GeoSeries(gdf_union), crs=gdf.crs)
        else:
            gdf_union = gdf

        # initialize distance helper to adjust distances and project if necessary
        distance_helper = kproj.Distance(self.unit, self.keep_input_crs)
        projected_gdf = distance_helper.pre_processing(exec_context, gdf_union, False)

        exec_context.set_progress(0.4, f"Processing distance 1 of {dist_length}")
        knut.check_canceled(exec_context)
        # Compute the buffers and append them to the result table
        result_gdf = gp.GeoDataFrame(
            geometry=projected_gdf.buffer(
                distance_helper.convert_input_distance(distances[0])
            )
        )
        if dist_length > 1:
            for i in range(1, dist_length):
                exec_context.set_progress(
                    (i + 1) * 1.0 / dist_length,
                    f"Processing distance {i + 1} of {dist_length}",
                )
                knut.check_canceled(exec_context)
                gdf_plus = gp.GeoDataFrame(
                    geometry=projected_gdf.buffer(
                        distance_helper.convert_input_distance(distances[i])
                    )
                )
                result_gdf = gp.overlay(
                    result_gdf, gdf_plus, how="union", keep_geom_type=False
                )

        # Append original distances as new column
        result_gdf[kproj.DEFAULT_DISTANCE_COLUMN_NAME] = distances
        result_gdf.rename_geometry(self.__BUFFER_COL_NAME, inplace=True)
        result_gdf = result_gdf.reset_index(drop=True)
        distance_helper.post_processing(exec_context, result_gdf, in_place=True)
        return knut.to_table(result_gdf, exec_context)


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
    description="Table with geometries to simplify.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Input table with simplified geometries.",
)
@knut.geo_node_description(
    short_description="Simplify the geometry",
    description="""This node returns for each input geometry a simplified representation. 
    The applied [algorithm](https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm) 
    recursively splits the original lines into smaller parts and connects the endpoints of each part by a straight line. 
    Then, it removes all points whose distance to the straight line is smaller than the provided tolerance distance. 
    It does not move any points and it always preserves endpoints of the original line or polygon.
    """,
    references={
        "GeoSeries.simplify": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.simplify.html",
        "Shapely object.simplify": "http://shapely.readthedocs.io/en/latest/manual.html#object.simplify",
    },
)
class SimplifyNode2:
    """
    This node returns a geometry feature containing a simplified representation of each geometry.
    """

    geo_col = knut.geo_col_parameter()

    distance = kproj.Distance.get_distance_parameter(
        label="Tolerance distance",
        description="The tolerance distances in the selected distance unit. "
        + "All parts of a simplified geometry will be no more than tolerance distance from their original position.",
        default_value=1.0,
    )

    unit = kproj.Distance.get_unit_parameter()

    keep_input_crs = kproj.Distance.get_keep_input_crs_parameter()

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name="Simplified",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.APPEND.name
        self.result_settings.new_column_name = "Simplified"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema_1,
            self.geo_col,
            knut.TYPE_GEO,
        )

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = knut.load_geo_data_frame(input, self.geo_col, exec_context)
        helper = kproj.Distance(self.unit, self.keep_input_crs)
        helper.pre_processing(exec_context, gdf, True)
        new_distance = helper.convert_input_distance(self.distance)

        result_gdf = self.result_settings.get_computed_result_frame(
            exec_context,
            input.schema,
            gdf,
            self.geo_col,
            lambda l: l.simplify(new_distance),
        )
        helper.post_processing(exec_context, result_gdf, True)
        return knut.to_table(result_gdf, exec_context)


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

    _COL_ID = "Grid ID"
    _COL_GEOMETRY = "geometry"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return knext.Schema(
            [knut.TYPE_POLYGON, knext.int64()], [self._COL_GEOMETRY, self._COL_ID]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = knut.load_geo_data_frame(input_table, self.geo_col, exec_context)

        xmin, ymin, xmax, ymax = gdf.total_bounds
        width = self.grid_length
        height = self.grid_length
        import numpy as np

        rows = int(np.ceil((ymax - ymin) / height))
        cols = int(np.ceil((xmax - xmin) / width))
        XleftOrigin = xmin
        XrightOrigin = xmin + width
        YtopOrigin = ymax
        YbottomOrigin = ymax - height
        polygons = []

        from shapely.geometry import Polygon  # For Grid

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

        grid = gp.GeoDataFrame({self._COL_GEOMETRY: polygons}, crs=gdf.crs)
        grid[self._COL_ID] = list(range(1, grid.shape[0] + 1))
        return knut.to_table(grid, exec_context)


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
        "The column containing the first latitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Lon2 = knext.ColumnParameter(
        "The second longitude column",
        "The column containing the second longitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Lat2 = knext.ColumnParameter(
        "The second latitude column",
        "The column containing the second latitude coordinates. ",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema):
        return input_schema.append(knext.Column(knext.double(), name="HDist"))

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        def HaversineDist(x1, y1, x2, y2):
            from math import radians, cos, sin, asin, sqrt  # For Haversine Distance

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
        result_col = knut.get_unique_column_name("HDist", input_table.schema)
        df[result_col] = df.apply(
            lambda x: HaversineDist(
                x[self.Lon1], x[self.Lat1], x[self.Lon2], x[self.Lat2]
            ),
            axis=1,
        )
        return knext.Table.from_pandas(df)


############################################
# Create Voronoi
############################################


@knext.node(
    name="Voronoi (Thiessen) Polygon",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Voronoi.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input Point Table",
    description="Input point data for Voronoi (Thiessen) polygons",
)
@knext.input_table(
    name="Reference Boundary Table",
    description="Input boundary data for Voronoi (Thiessen) polygons",
)
@knext.output_table(
    name="Voronoi (Thiessen) polygons",
    description="Output table of Voronoi (Thiessen) polygons",
)
class CreateVoronoi:
    """Create Voronoi (Thiessen) polygons.
    This node creates [Voronoi (Thiessen) polygons](https://en.wikipedia.org/wiki/Voronoi_diagram) from the input
    point data according to the reference boundary. The input data for the reference boundary should be a
    Polygon or MultiPolygon.

    The buffer distance with the provided unit is used to create dummy points that define a virtual boundary around the
    given reference boundary that controls the output of the Voronoi polygons. If the final Voronoi polygons are
    smaller than the given reference boundary, you might want to increase the buffer distance. For an illustration
    of the buffer distance see
    [here.](https://raw.githubusercontent.com/spatial-data-lab/knime-geospatial-extension/main/docs/imgs/voronoiDistance.png)
    """

    point_geo_col = knext.ColumnParameter(
        "Point geometry column",
        "Select the geometry column with the points for the Voronoi (Thiessen) polygons.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

    boundary_geo_col = knext.ColumnParameter(
        "Boundary geometry column",
        """Select the geometry column with the Polygon or MultiPolygon that defines the reference boundary for the 
        Voronoi (Thiessen) polygons.""",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    distance = kproj.Distance.get_distance_parameter(
        label="Buffer distance",
        description="""The buffer distance defines the distance of the dummy points that define the virtual boundary 
        from the given reference boundary.
              
        If the buffer distance is too small, the resulting polygons may not fill the entire reference 
        boundary, resulting in gaps or missing parts within it. This is because the buffer distance 
        determines how far away from the reference boundary the Voronoi polygons will be generated, and if the 
        distance is too small, some of the polygons may not intersect with the reference boundary or may only 
        partially intersect with it. For an illustration of the buffer distance see 
        [here.](https://raw.githubusercontent.com/spatial-data-lab/knime-geospatial-extension/main/docs/imgs/voronoiDistance.png)""",
        default_value=100,
        min_distance=1,
    )

    unit = kproj.Distance.get_unit_parameter()

    keep_input_crs = kproj.Distance.get_keep_input_crs_parameter()

    _COL_GEOMETRY = "Geometry"
    _COL_ID = "Thiessen ID"

    def configure(self, configure_context, input_schema1, input_schema2):
        self.point_geo_col = knut.column_exists_or_preset(
            configure_context, self.point_geo_col, input_schema1, knut.is_geo
        )
        self.boundary_geo_col = knut.column_exists_or_preset(
            configure_context, self.boundary_geo_col, input_schema2, knut.is_geo_polygon
        )
        return knext.Schema(
            [knut.TYPE_POLYGON, knext.int64()], [self._COL_GEOMETRY, self._COL_ID]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table1, input_table2):
        import pandas as pd
        from shapely.geometry import Point, Polygon
        import numpy as np
        from scipy.spatial import Voronoi

        origin_point = knut.load_geo_data_frame(
            input_table1, self.point_geo_col, exec_context
        )
        boundary = knut.load_geo_data_frame(
            input_table2, self.boundary_geo_col, exec_context
        )

        helper = kproj.Distance(self.unit, self.keep_input_crs)
        helper.pre_processing(exec_context, boundary, True)
        helper.pre_processing(exec_context, origin_point, True)
        dist_buffer = helper.convert_input_distance(self.distance)

        # Convert the GeoDataFrame to a numpy array
        points = np.array([[pt.x, pt.y] for pt in origin_point.geometry])

        # get the envelope of the union of the points with a buffer of 1
        knut.check_canceled(exec_context)
        exec_context.set_progress(0.3, "Computing virtual boundary")
        boundary_union = boundary.unary_union
        envelope = boundary_union.envelope.buffer(dist_buffer).envelope

        # get the coordinates of the envelope
        xmin, ymin, xmax, ymax = envelope.bounds
        dummy_points = np.array(
            [[xmin, ymin], [xmin, ymax], [xmax, ymax], [xmax, ymin]]
        )

        # combine the sample points and the dummy points
        points = np.concatenate((points, dummy_points))

        knut.check_canceled(exec_context)
        exec_context.set_progress(0.4, "Computing Voronoi regions")
        # compute Voronoi tessellation
        vor = Voronoi(points)
        knut.check_canceled(exec_context)
        exec_context.set_progress(0.8, "Converting Voronoi regions to table")
        # extract the vertices of each Voronoi region
        polygons = []
        for region in vor.regions:
            if -1 not in region:
                polygon = Polygon([vor.vertices[i] for i in region])
                polygons.append(polygon)

        knut.check_canceled(exec_context)
        # convert the list of Polygon objects to a GeoSeries object
        geo_series = gp.GeoSeries(polygons)
        # create GeoDataFrame
        gdf = pd.DataFrame({self._COL_GEOMETRY: geo_series})
        gdf = gp.GeoDataFrame(gdf, geometry=self._COL_GEOMETRY, crs=origin_point.crs)

        gdf = gdf[~gdf.geometry.is_empty]
        extent = gp.GeoDataFrame(
            geometry=gp.GeoSeries(boundary_union), crs=origin_point.crs
        )
        gdf = gp.overlay(gdf, extent, how="intersection")
        gdf.rename_geometry(self._COL_GEOMETRY, inplace=True)

        # project back if necessary
        helper.post_processing(exec_context, gdf, True)
        # append region id column
        gdf[self._COL_ID] = range(1, (gdf.shape[0] + 1))
        return knut.to_table(gdf, exec_context)
