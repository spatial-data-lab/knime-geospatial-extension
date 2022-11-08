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
    path="/geo",
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class CrsTransformerNode:
    """
    This node projects the data from its original CRS to the entered CRS.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "EPSG:4326")

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting projection...")
        gdf = gdf.to_crs(self.new_crs)
        exec_context.set_progress(0.1, "Projection done")
        LOGGER.debug("CRS converted to " + self.new_crs)
        return knext.Table.from_pandas(gdf)


############################################
# Geometry To Point
############################################


@knext.node(
    name="Geometry To Point",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/FeatureToPoint.png",
    category=category,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class GeometryToPointNode:
    """
    This node generate centroids of or respresentative points inside polygons,lines.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    pointtype = knext.StringParameter(
        "Selection",
        "The point type to choose from.",
        "centroid",
        enum=["centroid", "representative_point"],
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class ExplodeNode:
    """
    This node dismantles the multiparts into single parts.
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
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
        knut.column_exists(self.geo_col, input_schema_1)
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class PointsToLineNode:
    """
    This node generate lines from points according to group id and serial label.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

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

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
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
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
class GeometryToMultiPointNode:
    """
    This node generate line from the boundary of polygons.
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
