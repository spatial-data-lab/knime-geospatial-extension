# lingbo
import geopandas as gp
import logging
import knime_extension as knext
import util.knime_utils as knut
import util.projection as kproj

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/community/geo",
    level_id="transform",
    name="Spatial Transformation",
    description="Nodes that transform, decompose, and generate new geometric entities from single geometric objects.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/TransformationCategory.png",
    after="spatialtool",
)

############################################
# CRS Transformer
############################################


@knext.node(
    name="Projection",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/Projection.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed geo output table.",
)
@knut.geo_node_description(
    short_description="Projection Transformation",
    description="""This node transforms the 
    [Coordinate reference system (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) of the selected 
    geometry column to the entered new coordinate reference system. The node will transform the points in all 
    objects individually. It has no notion of projecting entire geometries. All segments joining points are assumed 
    to be lines in the current projection, not geodesics. Objects crossing the dateline (or other projection boundary) 
    will have undesirable behavior.
    """,
    references={
        "Map projection (Wikipedia)": "https://en.wikipedia.org/wiki/Map_projection",
        "Comparison of map projection": "https://map-projections.net/index.php",
        "Collection of common map projections and their properties": "https://www.icsm.gov.au/sites/default/files/projections.pdf",
        "Projection wizard that helps to find a good projection": "https://projectionwizard.org/",
        "Coordinate Reference System (CRS) EPSG:4326": "https://epsg.io/4326",
        "geopandas.GeoSeries.to_crs()": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.to_crs.html",
    },
)
class CrsTransformerNode:
    """
    This node projects the data from its original CRS to the entered CRS.
    """

    geo_col = knut.geo_col_parameter()

    new_crs = knext.StringParameter(
        "New CRS", kproj.DEF_CRS_DESCRIPTION, kproj.DEFAULT_CRS
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Projected",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.REPLACE.name
        self.result_settings.new_column_name = "Projected"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        # use the data type of the selected column as result type
        result_type = input_schema[self.geo_col].ktype
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            result_type,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = knut.load_geo_data_frame(input_table, self.geo_col, exec_context)
        if self.result_settings.mode == knut.ResultSettingsMode.APPEND.name:
            result_col = knut.get_unique_column_name(
                self.result_settings.new_column_name, input_table.schema
            )
            gdf[result_col] = gdf[self.geo_col]
            gdf.set_geometry(result_col, inplace=True)
        gdf.to_crs(self.new_crs, inplace=True)
        return knut.to_table(gdf, exec_context)


############################################
# Geometry To Point
############################################


@knext.node(
    name="Geometry To Point",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/FeatureToPoint.png",
    category=category,
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
    short_description="Returns a GeoSeries of points representing each geometry.",
    description="""This node returns a GeoSeries of points representing each geometry.
    There are two types of points, centroids and representative points. 
    The latter is guaranteed to be within each geometry.
    """,
    references={
        "GeoSeries.representative_point": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.representative_point.html",
        "GeoSeries.centroid": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.centroid.html",
    },
)
class GeometryToPointNode:
    """
    This node returns a GeoSeries of points representing each geometry.
    """

    geo_col = knut.geo_col_parameter()

    pointtype = knext.StringParameter(
        "Point Type Selection",
        "The point type to choose from.",
        "centroid",
        enum=["centroid", "representative_point"],
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Point",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.REPLACE.name
        self.result_settings.new_column_name = "Point"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            knut.TYPE_POINT,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        if self.pointtype == "centroid":
            func = lambda l: l.centroid
        else:
            func = lambda l: l.representative_point()

        return self.result_settings.get_computed_result_table(
            exec_context, input_1, self.geo_col, func
        )


############################################
# Multipart to Singlepart
############################################


@knext.node(
    name="Multipart To Singlepart",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/Explode.png",
    category=category,
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
    short_description="Explode multi-part geometries into multiple single geometries.",
    description="""This node explodes multi-part geometries into multiple single geometries.
    Each row containing a multi-part geometry will be split into multiple rows with single geometries,
    thereby increasing the number rows in the output table.
    """,
    references={
        "GeoDataFrame.explode": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.explode.html",
    },
)
class ExplodeNode:
    """
    This node dismantles the multiparts into single parts.
    """

    geo_col = knut.geo_col_parameter()

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting explosion...")
        exploded = gdf.explode(ignore_index=True)
        # gdf[self.geo_col] = exploded.geometry
        exec_context.set_progress(0.1, "Explosion done")
        LOGGER.debug("Feature geometry " + self.geo_col + "exploded")
        return knext.Table.from_pandas(exploded)


############################################
# Polygon To Line
############################################


@knext.node(
    name="Polygon To Line",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/PolygonToLine.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column. ",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column. ",
)
@knut.geo_node_description(
    short_description="Returns the boundaries of each polygon.",
    description="""This node return the boundaries of each polygon with geopandas.GeoSeries.boundary,
    which Returns a GeoSeries of lower dimensional objects representing each geometry’s set-theoretic boundary.
    """,
    references={
        "GeoSeries.boundary": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.boundary.html",
    },
)
class PolygonToLineNode:
    """
    This node generate lines from the boundaries of polygons.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.boolean_or(knut.is_geo_polygon, knut.is_geo_multi_polygon),
        include_row_key=False,
        include_none_column=False,
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Line",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.REPLACE.name
        self.result_settings.new_column_name = "Line"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context,
            self.geo_col,
            input_schema_1,
            knut.boolean_or(knut.is_geo_polygon, knut.is_geo_multi_polygon),
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema_1,
            self.geo_col,
            knut.TYPE_LINE,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # extract the boundary for each geometry
        return self.result_settings.get_computed_result_table(
            exec_context, input_table, self.geo_col, lambda l: l.boundary
        )


############################################
# Points To Line
############################################


@knext.node(
    name="Points To Line",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/PointToLine.png",
    category=category,
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
    short_description="This node generate lines from points according to group id and serial label.",
    description="""This node generate lines from points according to group id and serial label. The result
    table contains one LineString object per group value. Each constructed LineString object represents one or 
    more connected linear splines between the given points of a group ordered by the serial column. 
    Repeated points in the ordered sequence are allowed, but may incur performance penalties and should be avoided. 
    A LineString may cross itself. A LineString has zero area and non-zero length.
    """,
    references={
        "Shapely.LineStrings": "https://shapely.readthedocs.io/en/stable/manual.html#linestrings",
    },
)
class PointsToLineNode:
    """
    This node generate lines from points according to group id and serial label.
    """

    geo_col = knut.geo_point_col_parameter()

    group_col = knext.ColumnParameter(
        "Group column",
        "Select the group column (string) as group id for points.",
        # Allow only string columns
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    seiral_col = knext.ColumnParameter(
        "Serial column",
        "Select the serial column (numeric) for each group .",
        # Allow only string columns
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo_point
        )
        self.group_col = knut.column_exists_or_preset(
            configure_context, self.group_col, input_schema, knut.is_int_or_string
        )
        self.seiral_col = knut.column_exists_or_preset(
            configure_context, self.seiral_col, input_schema, knut.is_numeric
        )
        return gNone

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = gp.GeoDataFrame(input.to_pandas(), geometry=self.geo_col)
        if self.geo_col != "geometry":
            gdf.rename_geometry("geometry", inplace=True)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting grouping...")
        from shapely.geometry import MultiPoint, LineString

        line_gdf = (
            gdf.sort_values(by=[self.seiral_col])
            .groupby([self.group_col], as_index=False)["geometry"]
            .apply(lambda x: LineString(x.tolist()))
        )
        line_gdf = gp.GeoDataFrame(line_gdf, geometry="geometry", crs=gdf.crs)
        exec_context.set_progress(0.1, "PointsToLine done")
        return knext.Table.from_pandas(line_gdf)


############################################
# Line To MultiPoint
############################################


@knext.node(
    name="Line To MultiPoint",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/LinePolygonToPoints.png",
    category=category,
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
    short_description="This node generate points from the lines.",
    description="""This node generate points from the lines.
    The list of coordinates that describe a geometry are represented as the CoordinateSequence object in Shapely 
    which is the dependence of GeoPandas. 
    """,
    references={
        "Coordinate sequences": "https://shapely.readthedocs.io/en/stable/manual.html",
    },
)
class GeometryToMultiPointNode:
    """
    This node generate points from the lines.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_line,
        include_row_key=False,
        include_none_column=False,
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Multipoint",
    )

    def __init__(self):
        # set twice as workaround until fixed in KNIME framework
        self.result_settings.mode = knut.ResultSettingsMode.REPLACE.name
        self.result_settings.new_column_name = "Multipoint"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo_line
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            knut.TYPE_MULTI_POINT,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # extract coordinates of each geometry into a new MultiPoint geometry
        from shapely.geometry import MultiPoint

        return self.result_settings.get_computed_result_table(
            exec_context, input_table, self.geo_col, lambda l: MultiPoint(l.coords)
        )


############################################
# Create points in polygon
############################################


@knext.node(
    name="Create Random Points",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/RandomPoint.png",
    category=category,
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
    short_description="This node generate random points in polygons.",
    description="""This node generate random points in polygons.
    The Create Random Points node enables you to generate random points inside polygons based on a numerical value column and an ID column.
    The numerical value column is used to determine the number of points to be generated inside each polygon.
    Additionally, you will need to provide an ID column that will be used to identify each polygon.
    The node will create a new MultiPoint geometry that includes a random set of points for each polygon,which can be exploded into Points by the node Multipart To Singlepart.
    """,
    references={
        "Spatial Relationships": "https://shapely.readthedocs.io/en/stable/manual.html",
    },
)
class RandomPointNode:
    """
    This node generate points from the polygons.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    num_col = knext.ColumnParameter(
        "Number column",
        "Select the integer column for the number of points.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_int,
        include_row_key=False,
        include_none_column=False,
    )
    id_col = knext.ColumnParameter(
        "Unique ID column",
        "Select the column as polygon ID for points.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gdf[(gdf[self.num_col] >= 1) & (gdf.area > 0)].reset_index(drop=True)

        ###--Old Code Start--###
        # import random
        # from shapely.geometry import box, Point, MultiPoint
        # import pandas as pd
        # import numpy as np
        # import math
        # from itertools import product

        # def off(n, x):
        #     return np.array([random.uniform(-0.5 * x, 0.5 * x) for _ in range(n)])

        # def offset_sample(pt, n, grid_size_x, grid_size_y):
        #     indices = np.random.choice(len(pt), size=n)
        #     x = pt.x[indices] + off(n, grid_size_x)
        #     y = pt.y[indices] + off(n, grid_size_y)
        #     pointf = gp.GeoSeries(gp.points_from_xy(x, y)).unary_union
        #     return pointf

        # def offset_point(pt, n, grid_size_x, grid_size_y, polygon):
        #     points0 = []
        #     loop = math.ceil(n / len(pt)) + 1
        #     for i in range(loop):
        #         x = pt.x + off(len(pt), grid_size_x)
        #         y = pt.y + off(len(pt), grid_size_y)
        #         pointx = gp.GeoSeries(gp.points_from_xy(x, y))
        #         points0.extend(pointx)
        #     points1 = gp.GeoSeries(points0).unary_union.intersection(polygon)
        #     points1 = gp.GeoSeries(points1).explode(index_parts=True).values
        #     pointf = offset_sample(points1, n, grid_size_x, grid_size_y)
        #     return pointf

        # # define a function to generate grid cells for a single row
        # def generate_points(row, num_col):
        #     polygon = row.geometry
        #     n = int(row[num_col])
        #     bbox = polygon.bounds
        #     ncol = 10 if n < 30 else 30
        #     x_grid = np.linspace(bbox[0], bbox[2], ncol + 1)
        #     y_grid = np.linspace(bbox[1], bbox[3], ncol + 1)
        #     grid_size_x = (bbox[2] - bbox[0]) / ncol
        #     grid_size_y = (bbox[3] - bbox[1]) / ncol
        #     grid = pd.DataFrame(list(product(x_grid, y_grid)), columns=["x", "y"])
        #     pointinter = (
        #         gp.points_from_xy(grid.x, grid.y).unary_union().intersection(polygon)
        #     )
        #     points = gp.GeoSeries(pointinter).explode(index_parts=True).values
        #     if len(points) >= n:
        #         return offset_sample(points, n, grid_size_x, grid_size_y)
        #     else:
        #         return offset_point(points, n, grid_size_x, grid_size_y, polygon)

        # exec_context.set_progress(0.2, "Geo data frame loaded. Starting explosion...")

        # Generate all representative_point for point =1
        # gdf1 = gdf[gdf[self.num_col] == 1]
        # if gdf1.shape[0] > 0:
        #     gdf1 = gp.GeoDataFrame(
        #         gdf1[self.id_col], geometry=gdf1.representative_point(), crs=gdf.crs
        #     )
        # else:
        #     gdf1 = gp.GeoDataFrame()
        # exec_context.set_progress(0.3, "Geo data frame loaded. Starting explosion...")
        # gdf2 = gdf[gdf[self.num_col] > 1]
        # if gdf2.shape[0] > 0:
        #     points_df = gdf2.apply(
        #         lambda row: generate_points(row, self.num_col), axis=1
        #     )
        #     gdf3 = gp.GeoDataFrame(gdf2[self.id_col], geometry=points_df, crs=gdf.crs)
        # else:
        #     gdf3 = gp.GeoDataFrame()
        # dfs = pd.concat([gdf1, gdf3], axis=0)
        # dfs.sort_values(by=[self.id_col], inplace=True)
        # dfs.reset_index(drop=True, inplace=True)
        ####--Old Code End--###
        gdf2 = gdf[[self.id_col, self.num_col]]
        gdf2["geometry"] = gdf[self.geo_col].sample_points(size=gdf[self.num_col])
        return knext.Table.from_pandas(gdf2)
