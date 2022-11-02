from re import S
from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import pyproj as pyp
from shapely.geometry import shape


__category = knext.category(
    path="/geo",
    level_id="conversion",
    name="Spatial Conversion",
    description="Nodes that perform conversions from/to geometric objects.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/ConversionCategory.png",
    after="transform",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/GeometryConversion/"

########################################################################################
# Conversion helper methods
########################################################################################

DEFAULT_CRS = "epsg:4326"
"""Default coordinate reference system."""

__DEF_CRS_DESCRIPTION = """[Coordinate reference system (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) column
        Supports the following input types:
        
        - PROJ string
        - JSON string with PROJ parameters
        - CRS WKT string
        - An authority string [i.e. 'epsg:4326']
        - An EPSG integer code [i.e. 4326]
        - A tuple of ('auth_name': 'auth_code') [i.e ('epsg', '4326')]
        """


def crs_input_parameter(
    label: str = "CRS",
    description: str = __DEF_CRS_DESCRIPTION,
) -> knext.StringParameter:
    """Returns a CRS (coordinate reference system) string input parameter."""
    return knext.StringParameter(
        label=label,
        description=description,
        default_value=DEFAULT_CRS,
    )


def validate_crs(crs: str) -> None:
    """Validates the input crs and throws an appropriate exception if it is invalid."""
    try:
        parse_crs(crs)
    except Exception as error:
        raise knext.InvalidParametersError(str(error))


def parse_crs(crs: str) -> pyp.CRS:
    """Parses the input crs into a CRS object and throws an exception if it is invalid."""
    return pyp.CRS.from_user_input(crs)


########################################################################################
# Converter that converts a single input column to a geometry column
########################################################################################
class _ToGeoConverter:
    """
    Helper class for conversion that creates a new geometry column from a single input column.
    """

    DEF_GEO_COL_NAME = "geometry"

    def __init__(
        self,
        converter: Callable[[pd.DataFrame, knext.Column], gp.GeoSeries],
        column_type_check: Callable[[knext.Column], bool] = None,
    ):
        self.converter = converter
        self.column_type_check = column_type_check

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Data table",
                description=f"Input table to extract the geometry column from.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description="Input table with additional geometry column in the units of the provided CRS.",
            )
        ]

    def configure(self, configure_context, input_schema):
        # input_column and crs need to be defined in the child class
        self.input_column = knut.column_exists_or_preset(
            configure_context, self.input_column, input_schema, self.column_type_check
        )
        # check that the output column does NOT exist to prevent overwriting
        knut.fail_if_column_exists(
            self.DEF_GEO_COL_NAME,
            input_schema,
            f"Output column '{self.DEF_GEO_COL_NAME}'  exists in input table",
        )
        validate_crs(self.crs)

        # we can not return a schema since we do not know the geometric type of the column without parsing all rows

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        df = input_table.to_pandas()
        try:
            # input_column and crs need to be defined in the child class
            df[self.DEF_GEO_COL_NAME] = self.converter(df, self.input_column)
            gdf = gp.GeoDataFrame(df, geometry=self.DEF_GEO_COL_NAME, crs=self.crs)
        except Exception as error:
            raise RuntimeError(
                f"Error converting '{self.input_column}' column to geometry: {str(error)}"
            )
        # execute the function and append the result
        return knut.to_table(gdf, exec_context)


############################################
# WKT to Geometry
############################################
@knext.node(
    name="WKT to Geometry",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "WKTtoGeo.png",
    category=__category,
)
@knut.geo_node_description(
    short_description="Converts the input Well-known-text (WKT) column to a geometry column.",
    description="""This node converts the selected 
    [Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) input column to 
    a geometry column in the units of the provided CRS.
    """,
    references={
        "From WKT": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.from_wkt.html",
        "CRS parser from pyproj": "https://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input",
    },
)
class WKTtoGeoNode(_ToGeoConverter):

    input_column = knext.ColumnParameter(
        label="WKT Column",
        description="[Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column to convert",
        column_filter=knut.is_string,
        include_row_key=False,
        include_none_column=False,
    )

    crs = crs_input_parameter()

    def __init__(self):
        super().__init__(
            lambda df, col: gp.GeoSeries.from_wkt(df[col].astype("str")), knut.is_string
        )


############################################
# WKB to Geometry
############################################
# binary data cells not supported yet in the KNIME Python Extension
# @knext.node(
#     name="WKB to Geometry",
#     node_type=knext.NodeType.MANIPULATOR,
#     icon_path=__NODE_ICON_PATH + "WKBtoGeo.png",
#     category=__category,
# )
# @knut.geo_node_description(
#     short_description="Converts the input Well-known-binary (WKB) column to a geometry column.",
#     description="""This node converts the selected
#     [Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) input column to
#     a geometry column in the units of the provided CRS.
#     """,
#     references={
#         "From WKB": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.from_wkb.html",
#         "CRS parser from pyproj": "https://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input",
#     },
# )
# class WKBtoGeoNode(_ToGeoConverter):

#     input_column = knext.ColumnParameter(
#         label="WKB Column",
#         description="[Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column to convert",
#         column_filter=knut.is_binary,
#         include_row_key=False,
#         include_none_column=False,
#     )

#     crs = crs_input_parameter()

#     def __init__(self):
#         super().__init__(lambda df, col: gp.GeoSeries.from_wkb(df[col]), knut.is_binary)


############################################
# WKT to Geometry
############################################
@knext.node(
    name="GeoJSON to Geometry",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoJSONtoGeo.png",
    category=__category,
)
@knut.geo_node_description(
    short_description="Converts the input GeoJSON column to a geometry column.",
    description="""This node converts the selected [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON) input column to 
    a geometry column in the units of the provided CRS.
    """,
    references={
        "From GeoJSON": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.from_wkt.html",
        "CRS parser from pyproj": "https://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input",
    },
)
class GeoJSONtoGeoNode(_ToGeoConverter):

    input_column = knext.ColumnParameter(
        label="GeoJSON Column",
        description="[GeoJSON](https://en.wikipedia.org/wiki/GeoJSON) column to convert",
        column_filter=knut.is_string,
        include_row_key=False,
        include_none_column=False,
    )

    crs = crs_input_parameter()

    def __init__(self):
        super().__init__(
            lambda df, col: df,
            knut.is_string,
        )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        df = input_table.to_pandas()
        try:
            # input_column and crs need to be defined in the child class
            mygdf = gp.read_file(df[self.input_column][0], driver="GeoJSON")

            gdf = gp.GeoDataFrame(df, geometry=self.DEF_GEO_COL_NAME, crs=self.crs)
        except Exception as error:
            raise RuntimeError(
                f"Error converting '{self.input_column}' column to geometry: {str(error)}"
            )
        # execute the function and append the result
        return knut.to_table(gdf, exec_context)


############################################
# Lat/Lon to Geometry
############################################
@knext.node(
    name="Lat/Lon to Geometry",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "LatLongToGeo.png",
    category=__category,
)
@knext.input_table(
    name="Table with latitude and longitude",
    description="Table with columns for latitude and longitude to convert to a geometric point.",
)
@knext.output_table(
    name="Geo table",
    description="Input table with additional geometry column in the units of the provided CRS.",
)
@knut.geo_node_description(
    short_description="Converts the given latitude and longitude column to a geometric point column.",
    description="This node converts the given latitude and longitude into a geometric point."
    + "The geometric point is appended to the input table as geometry column.",
    references={
        "Points from XY": "https://geopandas.org/en/stable/docs/reference/api/geopandas.points_from_xy.html",
    },
)
class LatLongToGeoNode:

    lat_col = knext.ColumnParameter(
        label="Latitude column",
        description="Please select the latitude column",
        column_filter=knut.is_numeric,
    )

    lon_col = knext.ColumnParameter(
        label="Longitude column",
        description="Please select the longitude column",
        column_filter=knut.is_numeric,
    )

    crs = crs_input_parameter()

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.lat_col, input_schema_1, knut.is_numeric)
        knut.column_exists(self.lon_col, input_schema_1, knut.is_numeric)
        knut.fail_if_column_exists(
            _ToGeoConverter.DEF_GEO_COL_NAME,
            input_schema_1,
            f"Output column '{_ToGeoConverter.DEF_GEO_COL_NAME}'  exists in input table",
        )
        validate_crs(self.crs)

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # https://geopandas.org/en/stable/docs/user_guide/projections.html#the-axis-order-of-a-crs
        # In GeoPandas the coordinates are always stored as (x, y), and thus as (lon, lat) order, regardless of the CRS
        df = input_table.to_pandas()
        try:
            # input_column and crs need to be defined in the child class
            df[_ToGeoConverter.DEF_GEO_COL_NAME] = gp.points_from_xy(
                df[self.lon_col], df[self.lat_col]
            )
            gdf = gp.GeoDataFrame(
                df, geometry=_ToGeoConverter.DEF_GEO_COL_NAME, crs=self.crs
            )
        except Exception as error:
            raise RuntimeError(
                f"Error converting '{self.lat_col}' and '{self.lon_col}' column to geometry: {str(error)}"
            )
        # execute the function and append the result
        return knut.to_table(gdf, exec_context)


########################################################################################
# Converter that converts a geometry column to another column
########################################################################################
class _FromGeoConverter:
    """
    Helper class for conversion that creates a new column from a geometry column.
    """

    def __init__(
        self,
        column_name: str,
        column_type: knext.KnimeType,
        converter: Callable[[gp.GeoSeries], gp.GeoSeries],
    ):
        self.column_name = column_name
        self.column_type = column_type
        self.converter = converter

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Data table",
                description=f"Input table with geometry column to extract the {self.column_name} column from.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description=f"Input table with additional {self.column_name} column extracted from the specified geometry column.",
            )
        ]

    def configure(self, configure_context, input_schema):
        # geo_column needs to be defined in the child class
        self.geo_column = knut.column_exists_or_preset(
            configure_context, self.geo_column, input_schema, knut.is_geo
        )
        # check that the output column does NOT exist to prevent overwriting
        knut.fail_if_column_exists(
            self.column_name,
            input_schema,
            f"Output column '{self.column_name}'  exists in input table",
        )
        return input_schema.append(knext.Column(self.column_type, self.column_name))

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        gdf = knut.load_geo_data_frame(input_table, self.geo_column, exec_context)
        # execute the function and append the result
        try:
            result = self.converter(gdf[self.geo_column])
            gdf[self.column_name] = result
        except Exception as error:
            raise knext.InvalidParametersError(
                f"Error converting to {self.column_name}: {str(error)}"
            )
        return knut.to_table(gdf, exec_context)


############################################
# Geometry to WKT
############################################
@knext.node(
    name="Geometry to WKT",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoToWKT.png",
    category=__category,
)
@knut.geo_node_description(
    short_description="Converts the input geometry column to a Well-known-text (WKT) column.",
    description="""This node converts the selected geometry column into a
    [Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry).
    """,
    references={
        "To WKT": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_wkt.html"
    },
)
class GeoToWKTNode(_FromGeoConverter):

    geo_column = knut.geo_col_parameter(
        description="Geometry column to convert to [Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry)"
    )

    def __init__(self):
        super().__init__(
            "WKT",
            knext.string(),
            lambda gs: gs.to_wkt(),
        )


############################################
# Geometry to WKB
############################################
# binary data cells not supported yet in the KNIME Python Extension
# @knext.node(
#     name="Geometry to WKB",
#     node_type=knext.NodeType.MANIPULATOR,
#     icon_path=__NODE_ICON_PATH + "GeoToWKB.png",
#     category=__category,
# )
# @knut.geo_node_description(
#     short_description="Converts the input geometry column to a Well-known-binary (WKB) column.",
#     description="""This node converts the selected geometry column into a
#     [Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column.
#     """,
#     references={
#         "To WKB": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_wkb.html"
#     },
# )
# class GeoToWKBNode(_FromGeoConverter):

#     geo_column = knut.geo_col_parameter(
#         description="Geometry column to convert to [Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry)"
#     )

#     def __init__(self):
#         super().__init__(
#             "WKB",
#             knext.blob(),
#             lambda gs: gs.to_wkb(),
#         )


############################################
# Geometry to GeoJSON
############################################
@knext.node(
    name="Geometry to GeoJSON",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoToGeoJSON.png",
    category=__category,
)
@knut.geo_node_description(
    short_description="Converts the input geometry column to a Well-known-binary (WKB) column.",
    description="""This node converts the selected geometry column into a 
    [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON) column.
    """,
    references={
        "To JSON": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.to_json.html"
    },
)
class GeoToGeoJSONNode(_FromGeoConverter):

    geo_column = knut.geo_col_parameter(
        description="Geometry column to convert to [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON)"
    )

    def __init__(self):
        super().__init__(
            "GeoJSON",
            # TODO: JSON
            knext.string(),
            lambda gs: gs.to_json(),
        )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        gdf = knut.load_geo_data_frame(input_table, self.geo_column, exec_context)
        # execute the function and append the result
        try:
            gdf[self.column_name] = gdf[self.geo_column].apply(
                lambda e: gp.GeoSeries(e).to_json()
            )
        except Exception as error:
            raise knext.InvalidParametersError(
                f"Error converting to {self.column_name}: {str(error)}"
            )
        return knut.to_table(gdf, exec_context)


############################################
# Geometry to Lat/Lon
############################################
@knext.node(
    name="Geometry to Lat/Long",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoToLatLon.png",
    category=__category,
)
@knext.input_table(
    name="Geo table",
    description="Table with a point geometry column to extract the latitude and longitude from.",
)
@knext.output_table(
    name="Table with latitude and longitude",
    description="Table with the latitude and longitude for all point geometries",
)
@knut.geo_node_description(
    short_description="Extracts the latitude and longitude.",
    description="This node extracts the latitude and longitude for all given point objects."
    + "The longitude and latitude columns are appended to the input table as lat and long column.",
    references={
        "X coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.x.html",
        "Y coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.y.html",
    },
)
class GeoToLatLongNode:

    geo_col = knut.geo_point_col_parameter(
        description="Select the point geometry column to extract the latitude and longitude from."
    )

    col_lat = "lat"
    col_lon = "lon"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo_point
        )
        knut.fail_if_column_exists(
            self.col_lat,
            input_schema_1,
            f"Output column '{self.col_lat}'  exists in input table",
        )
        knut.fail_if_column_exists(
            self.col_lon,
            input_schema_1,
            f"Output column '{self.col_lon}'  exists in input table",
        )
        return input_schema_1.append(
            [
                knext.Column(knext.double(), self.col_lat),
                knext.Column(knext.double(), self.col_lon),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        gs = gdf.loc[:, self.geo_col]
        # https://geopandas.org/en/stable/docs/user_guide/projections.html#the-axis-order-of-a-crs
        # In GeoPandas the coordinates are always stored as (x, y), and thus as (lon, lat) order, regardless of the CRS
        gdf[self.col_lon] = gs.x
        gdf[self.col_lat] = gs.y
        return knut.to_table(gdf, exec_context)
