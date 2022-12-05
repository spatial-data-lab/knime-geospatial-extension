# lingbo
import logging
from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
from shapely.geometry import Point, MultiPoint, LineString

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/community/geo",
    level_id="transform",
    name="Spatial Transformation",
    description="Geospatial transformation nodes",
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
    description="""This node transforms the Coordinate Reference System (CRS) of the geometry column  with the input 
    parameter by geopandas.to_crs(). This method will transform all points in all objects. 
    It has no notion or projecting entire geometries. 
    All segments joining points are assumed to be lines in the current projection, not geodesics. 
    Objects crossing the dateline (or other projection boundary) will have undesirable behavior.
    """,
    references={
        "geopandas.GeoSeries.to_crs()": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.to_crs.html",
        "Coordinate Reference System (CRS) EPSG:4326": "https://epsg.io/4326",
    },
)
class CrsTransformerNode:
    """
    This node projects the data from its original CRS to the entered CRS.
    """

    geo_col = knut.geo_col_parameter()

    new_crs = knext.StringParameter(
        "New CRS", knut.DEF_CRS_DESCRIPTION, knut.DEFAULT_CRS
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = knut.load_geo_data_frame(input_table, self.geo_col, exec_context)
        gdf.to_crs(self.new_crs, inplace=True)
        crs = gdf.crs
        LOGGER.debug("CRS converted to " + self.new_crs)
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
        "Selection",
        "The point type to choose from.",
        "centroid",
        enum=["centroid", "representative_point"],
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        if self.pointtype == "centroid":
            gdf["point"] = gdf.centroid
        else:
            gdf["point"] = gdf.representative_point()
        gdf = gdf.set_geometry("point")
        gdf = gdf.drop(columns=self.geo_col)
        gdf = gdf.rename(columns={"point": self.geo_col})
        exec_context.set_progress(0.1, "Transformation done")
        LOGGER.debug("Feature converted to " + self.pointtype)
        return knext.Table.from_pandas(gdf)


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
        column_filter=knut.is_geo_polygon,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo_polygon
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting explosion...")
        gdf[self.geo_col] = gdf.boundary
        exec_context.set_progress(0.1, "PolygonToLine done")
        LOGGER.debug("Polygon feature " + self.geo_col + "transformed to line")
        return knext.Table.from_pandas(gdf)


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
    description="""This node generate lines from points according to group id and serial label by Shapely.LineString().
    The constructed LineString object represents one or more connected linear splines between the points. 
    Repeated points in the ordered sequence are allowed, but may incur performance penalties and should be avoided. 
    A LineString may cross itself.A LineString has zero area and non-zero length.
    """,
    references={
        "Shapely.LineStrings": "https://shapely.readthedocs.io/en/stable/manual.html",
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
        column_filter=knut.is_string,
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
            configure_context, self.group_col, input_schema, knut.is_string
        )
        self.seiral_col = knut.column_exists_or_preset(
            configure_context, self.seiral_col, input_schema, knut.is_numeric
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = gp.GeoDataFrame(input.to_pandas(), geometry=self.geo_col)
        gdf = gdf.rename(columns={self.geo_col: "geometry"})
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting explosion...")
        line_gdf = (
            gdf.sort_values(by=[self.seiral_col])
            .groupby([self.group_col], as_index=False)["geometry"]
            .apply(lambda x: LineString(x.tolist()))
        )
        line_gdf = gp.GeoDataFrame(line_gdf, geometry="geometry", crs=gdf.crs)
        exec_context.set_progress(0.1, "PolygonToLine done")
        LOGGER.debug(
            "Point feature "
            + self.geo_col
            + "transformed to line by group column"
            + self.group_col
            + "according to the order of"
            + self.seiral_col
        )
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

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo_line
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting explosion...")
        gdf["points"] = gdf.apply(lambda l: l["geometry"].coords, axis=1)
        gdf["geometry"] = gdf["points"].apply(lambda l: MultiPoint(l))
        gdf = gdf.drop(columns="points")
        exec_context.set_progress(0.1, "LineToMultiPoint done")
        LOGGER.debug("Line feature " + self.geo_col + "transformed to MultiPoint")
        return knext.Table.from_pandas(gdf)
