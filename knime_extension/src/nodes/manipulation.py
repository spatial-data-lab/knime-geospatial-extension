
from typing import Callable
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


category = knext.category(
    path="/geo",
    level_id="manipulation",
    name="Geometric Manipulation",
    description="Geometric manipulation nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/edit.png",
)


class _SingelManipulator:
    """
    Helper class for transformations that append a single result column to an input table with geo column.
    """
    def __init__(self, geo_col: knext.Column, func: Callable[[gp.GeoDataFrame], gp.GeoDataFrame],
                 colname: str, coltype: knext.KnimeType = knext.double()):
        self.geocol = geo_col
        self.func = func
        self.colname = colname
        self.coltype = coltype

        #standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description=f"Input table with a geospatial column to compute the '{self.colname}' for.",
            )
        ]
        #standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name=f"Geo table with {self.colname}",
                description=f"Geo input table with additional '{self.colname}' column, which calculate the "
                            + f"{self.colname} of each geometry in the units of the CRS.",
            )
        ]

    def configure(self, configure_context, input_schema):
        knut.column_exists(self.geo_col, input_schema)
        return input_schema.append(knext.Column(self.coltype, self.colname))

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)
        gdf[self.colname] = self.func(gdf)  #execute the function and append the result
        return knext.Table.from_pandas(gdf)

# Node definition based on _SingelManipulator works without the class but would have no node description:
# @knext.node(name="Compute Area", node_type=knext.NodeType.MANIPULATOR, icon_path="icons/icon.png", category=category)
# def compute_area_node():
#     __doc__="This node will add a new column named 'area' to the input data, which calculate the area of each "
#               +"geometry in the units of the CRS."
#     return _SingelManipulator("area", lambda gdf: gdf.area)

@knext.node(name="Compute Area", node_type=knext.NodeType.MANIPULATOR, icon_path="icons/icon.png", category=category)
class ComputeAreaNode(_SingelManipulator):
    """
    This node will add a new column named 'area' to the input data, which calculate the area of each geometry in 
    the units of the CRS.
    """
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to compute the area.",
        column_filter=knut.is_geo,  # Allow all GeoValue compatible columns
        include_row_key=False,
        include_none_column=False
    )

    def __init__(self):
        super().__init__(self.geo_col, lambda gdf: gdf.area, "area")


@knext.node(name="Compute Length", node_type=knext.NodeType.MANIPULATOR, icon_path="icons/icon.png", category=category)
class ComputeLengthNode(_SingelManipulator):
    """
    This node will add a new column named 'length' to the input data, which calculate the length of each geometry in 
    the units of the CRS.
    """
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to compute the length.",
        column_filter=knut.is_geo,  # Allow all GeoValue compatible columns
        include_row_key=False,
        include_none_column=False
    )

    def __init__(self):
        super().__init__(self.geo_col, lambda gdf: gdf.length, "length")
