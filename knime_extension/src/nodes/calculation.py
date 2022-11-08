from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


__category = knext.category(
    path="/geo",
    level_id="calculation",
    name="Spatial Calculation",
    description="Nodes that calculate properties for given geometric objects.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/CalculationCategory.png",
    after="io",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/GeometryCalculation/"

############################################
# Simple transformation helper class
############################################
class _SingleCalculator:
    """
    Helper class for transformations that append a single result column to an input table with geo column.
    """

    def __init__(
        self,
        func: Callable[[gp.GeoDataFrame], gp.GeoDataFrame],
        col_name: str,
        col_type: knext.KnimeType = knext.double(),
    ):
        self.func = func
        self.col_name = col_name
        self.coltype = col_type

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description=f"Input table with a geospatial column to compute the '{self.col_name}' for.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name=f"Geo table with {self.col_name}",
                description=f"""Geo input table with additional '{self.col_name}' column containing the
                {self.col_name} of each input geometry in the units of the CRS.
                """,
            )
        ]

    def configure(self, configure_context, input_schema):
        # geo_col needs to be defined in the child class
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return input_schema.append(knext.Column(self.coltype, self.col_name))

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column.
        # The geo_col needs to be defined in the child class
        gdf = knut.load_geo_data_frame(input_table, self.geo_col, exec_context)
        gdf[self.col_name] = self.func(
            gdf
        )  # execute the function and append the result
        return knut.to_table(gdf, exec_context)


############################################
# Area node
############################################
@knext.node(
    name="Area",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Area.png",
    category=__category,
)
@knut.geo_node_description(
    short_description="Calculates the area of geometric objects.",
    description="""This node will add a new column named 'area' to the input data with the calculated area of
    each geometry in the units of the CRS.  
    Area may be invalid for a geographic CRS using degrees as units. Use the Projection node to project 
    geometries to a planar CRS before using this node.
    This operation is planar, i.e. the potential third dimension is not taken into account.
    """,
    references={
        "Area": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.area.html"
    },
)
class AreaNode(_SingleCalculator):

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the area."
    )

    def __init__(self):
        super().__init__(lambda gdf: gdf.area, "area")


############################################
# Length node
############################################
@knext.node(
    name="Length",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Length.png",
    category=__category,
)
@knut.geo_node_description(
    short_description="Calculates the length of geometric objects.",
    description="""This node will add a new column named 'length' to the input data
    with the calculated length of each geometry in the units of the CRS.  
    Length may be invalid for a geographic CRS using degrees as units. Use the Projection node to project 
    geometries to a planar CRS before using this node.
    This operation is planar, i.e. the potential third dimension is not taken into account.
    """,
    references={
        "Length": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.length.html"
    },
)
class LengthNode(_SingleCalculator):

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the length."
    )

    def __init__(self):
        super().__init__(lambda gdf: gdf.length, "length")


############################################
# Bounds node
############################################
@knext.node(
    name="Bounds",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Bounds.png",
    category=__category,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to compute the bounds per geometric object",
)
@knext.output_table(
    name="Geo table with bounds",
    description="Input geo table with the bounds per geometric object",
)
@knut.geo_node_description(
    short_description="Computes the bound for geometric objects.",
    description="This node computes the bounds for each geometry object. The bounds are appended to the input table as x/y coordinate columns.",
    references={
        "Bounds": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.bounds.html"
    },
)
class BoundsNode:
    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the bounds for."
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1.append(
            [
                knext.Column(knext.double(), "minx"),
                knext.Column(knext.double(), "miny"),
                knext.Column(knext.double(), "maxx"),
                knext.Column(knext.double(), "maxy"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        bounds = gdf.bounds
        gdf = pd.concat([gdf, bounds], axis=1, copy=False)
        return knut.to_table(gdf, exec_context)


############################################
# Total bounds node
############################################
@knext.node(
    name="Total Bounds",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "TotalBounds.png",
    category=__category,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to compute the total bounds for.",
)
@knext.output_table(
    name="Table with the total bounds",
    description="Table with the total bounds for all geometric input objects",
)
@knut.geo_node_description(
    short_description="Computes the total bounds for all given geometric objects.",
    description="This node computes the total bounds for all given geometric input objects.",
    references={
        "Total bounds": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.total_bounds.html"
    },
)
class TotalBoundsNode:

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the total bounds for."
    )

    col_min_x = "minx"
    col_min_y = "miny"
    col_max_x = "maxx"
    col_max_y = "maxy"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return knext.Schema.from_columns(
            [
                knext.Column(knext.double(), self.col_min_x),
                knext.Column(knext.double(), self.col_min_y),
                knext.Column(knext.double(), self.col_max_x),
                knext.Column(knext.double(), self.col_max_y),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        bounds = gdf.total_bounds
        result = pd.DataFrame(
            [bounds],
            columns=[self.col_min_x, self.col_min_y, self.col_max_x, self.col_max_y],
        )
        return knut.to_table(result, exec_context)


############################################
# Coordinates node
############################################
@knext.node(
    name="Coordinates XYZ",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "XYZcoordinates.png",
    category=__category,
)
@knext.input_table(
    name="Geo table",
    description="Table with a point geometry column to extract the coordinates from.",
)
@knext.output_table(
    name="Table with coordinates",
    description="Table with the coordinates for all point geometries",
)
@knut.geo_node_description(
    short_description="Extracts the XYZ coordinates.",
    description="This node extracts the XYZ coordinates for all given point objects."
    + "The coordinates are appended to the input table as xyz coordinate columns.",
    references={
        "X coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.x.html",
        "Y coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.y.html",
        "Z coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.z.html",
    },
)
class CoordinatesNode:

    geo_col = knut.geo_point_col_parameter(
        description="Select the point geometry column to extract the coordinates from."
    )

    add_z = knext.BoolParameter(
        "Extract Z coordinate",
        "Select this option to extract the Z coordinate. If not selected only the X and Y coordinate are extracted",
        False,
    )

    col_x = "x"
    col_y = "y"
    col_z = "z"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo_point
        )
        result = input_schema_1.append(
            [
                knext.Column(knext.double(), self.col_x),
                knext.Column(knext.double(), self.col_y),
            ]
        )
        if self.add_z:
            result = result.append(knext.Column(knext.double(), self.col_z))
        return result

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        gs = gdf.loc[:, self.geo_col]
        gdf[self.col_x] = gs.x
        gdf[self.col_y] = gs.y
        if self.add_z:
            gdf[self.col_z] = gs.z
        return knut.to_table(gdf, exec_context)


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
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf["geometry"] = gdf.envelope
        exec_context.set_progress(0.1, "Transformation done")
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
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf["geometry"] = gdf.convex_hull
        exec_context.set_progress(0.1, "Transformation done")
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
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf_union = gdf.unary_union
        gdfunion = gp.GeoDataFrame(geometry=gp.GeoSeries(gdf_union), crs=gdf.crs)
        exec_context.set_progress(0.1, "Transformation done")
        return knext.Table.from_pandas(gdfunion)
