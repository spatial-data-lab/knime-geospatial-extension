from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


category = knext.category(
    path="/geo",
    level_id="calculation",
    name="Geometry Calculation",
    description="Nodes that calculate properties for given geometric objects.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/CalculationCategory.png",
)


class _SingelCalculator:
    """
    Helper class for transformations that append a single result column to an input table with geo column.
    """

    def __init__(
        self,
        geo_col: knext.Column,
        func: Callable[[gp.GeoDataFrame], gp.GeoDataFrame],
        colname: str,
        coltype: knext.KnimeType = knext.double(),
    ):
        self.geocol = geo_col
        self.func = func
        self.colname = colname
        self.coltype = coltype

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description=f"Input table with a geospatial column to compute the '{self.colname}' for.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name=f"Geo table with {self.colname}",
                description="""Geo input table with additional '{self.colname}' column containing the "
                {self.colname} of each input geometry in the units of the CRS.  
                
                Area may be invalid for a geographic CRS using degrees as units; use the Projection node to project 
                geometries to a planar CRS before using this node.
                This operation is planar, i.e. the potential third dimension is not taken into account.
                """,
            )
        ]

    def configure(self, configure_context, input_schema):
        knut.column_exists(self.geo_col, input_schema)
        return input_schema.append(knext.Column(self.coltype, self.colname))

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)
        gdf[self.colname] = self.func(gdf)  # execute the function and append the result
        return knext.Table.from_pandas(gdf)


# Node definition based on _SingelManipulator works without the class but would have no node description:
# @knext.node(name="Compute Area", node_type=knext.NodeType.MANIPULATOR, icon_path="icons/icon.png", category=category)
# def compute_area_node():
#     __doc__="This node will add a new column named 'area' to the input data, which calculate the area of each "
#               +"geometry in the units of the CRS."
#     return _SingelManipulator("area", lambda gdf: gdf.area)


@knext.node(
    name="Area",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryCalculation/LengthArea.png",
    category=category,
)
class ComputeAreaNode(_SingelCalculator):
    """
    This node will add a new column named 'area' to the input data, which calculate the area of each geometry in
    the units of the CRS.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to compute the area.",
        column_filter=knut.is_geo,  # Allow all GeoValue compatible columns
        include_row_key=False,
        include_none_column=False,
    )

    def __init__(self):
        super().__init__(self.geo_col, lambda gdf: gdf.area, "area")


@knext.node(
    name="Length",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryCalculation/LengthArea.png",
    category=category,
)
class ComputeLengthNode(_SingelCalculator):
    """
    This node will add a new column named 'length' to the input data, which calculate the length of each geometry in
    the units of the CRS.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to compute the length.",
        column_filter=knut.is_geo,  # Allow all GeoValue compatible columns
        include_row_key=False,
        include_none_column=False,
    )

    def __init__(self):
        super().__init__(self.geo_col, lambda gdf: gdf.length, "length")


@knext.node(
    name="Bounds",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryCalculation/Bounds.png",
    category=category,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to compute the bounds per geometric object",
)
@knext.output_table(
    name="Geo table with bounds",
    description="Input geo table with the bounds per geometric object",
)
class BoundsNode:
    """
    This node computes the bounds for each geometry object.
    The bounds are appended to the input table as x/y coordinate columns.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to compute the bounds for.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return input_schema_1.append(
            [
                knext.Column(knext.double(), "minx"),
                knext.Column(knext.double(), "miny"),
                knext.Column(knext.double(), "maxx"),
                knext.Column(knext.double(), "maxy"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting computation...")
        bounds = gdf.bounds
        gdf = pd.concat([gdf, bounds], axis=1, copy=False)
        exec_context.set_progress(1.0, "Computation done")
        return knext.Table.from_pandas(gdf)


@knext.node(
    name="Total Bounds",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryCalculation/TotalBounds.png",
    category=category,
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to compute the total bounds for.",
)
@knext.output_table(
    name="Table with the total bounds",
    description="Table with the total bounds for all geometric input objects",
)
class TotalBoundsNode:
    """
    This node computes the total bounds for all given geometric input objects.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to compute the total bounds for.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return knext.Schema.from_columns(
            [
                knext.Column(knext.double(), "minx"),
                knext.Column(knext.double(), "miny"),
                knext.Column(knext.double(), "maxx"),
                knext.Column(knext.double(), "maxy"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting computation...")
        bounds = gdf.total_bounds
        exec_context.set_progress(1.0, "Computation done")
        return knext.Table.from_pandas(
            pd.DataFrame([bounds], columns=["minx", "miny", "maxx", "maxy"])
        )


@knext.node(
    name="Coordinates (XYZ)",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryCalculation/XYZcoordinates.png",
    category=category,
)
@knext.input_table(
    name="Geo table",
    description="Table with a point geometry column to extract the coordinates from.",
)
@knext.output_table(
    name="Table with coordinates",
    description="Table with the coordinates for all point geometries",
)
class CoordinatesNode:
    """
    This node extracts the coordinates for each given point geometry.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the point geometry column to extract the coordinates from.",
        # Allow only GeoPointValue compatible columns
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return input_schema_1.append(
            [
                knext.Column(knext.double(), "x"),
                knext.Column(knext.double(), "y"),
                knext.Column(knext.double(), "z"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting computation...")
        gs = gdf.loc[:, self.geo_col]
        gdf["x"] = gs.x
        gdf["y"] = gs.y
        gdf["z"] = gs.z
        exec_context.set_progress(1.0, "Computation done")
        return knext.Table.from_pandas(gdf)
