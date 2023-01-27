import logging
from typing import Callable
from typing import List

import geopandas as gp
import knime.types.geospatial as gt
import knime_extension as knext

from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPolygon
from shapely.geometry import GeometryCollection


LOGGER = logging.getLogger(__name__)

DEFAULT_CRS = "epsg:4326"
"""Default coordinate reference system."""

DEF_CRS_DESCRIPTION = """Enter the 
        [Coordinate reference system (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) to use.

        Common [EPSG codes](https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset) that can be universally 
        used for mapping coordinates everywhere in the world are 
        [epsg:4326 (WGS 84, Unit: degree)](https://epsg.io/4326) (Latitude/longitude coordinate system based 
        on the Earth's center of mass;  Used by the Global Positioning System among others) and 
        [epsg:3857 (Unit: meter)](https://epsg.io/3857) (Web Mercator projection used by many web-based mapping tools,
        including Google Maps and OpenStreetMap.). 
        
        There are EPSG codes for specific regions that provide a higher accuracy in these regions, such as 
        [epsg:4269 (NAD83, Unit: degree)](https://epsg.io/4269) 
        and [epsg:26918 (NAD83 18N, Unit: meter)](https://epsg.io/26918) for North America, 
        and [epsg:4490 (CGCS2000, Unit: degree)](https://epsg.io/4490) 
        and [epsg:4479 (CGCS2000, Unit: meter)](https://epsg.io/4479) for China.
        
        The input field supports the following input types:
        
        - An authority string [i.e. 'epsg:4326']
        - An EPSG integer code [i.e. 4326]
        - A tuple of ('auth_name': 'auth_code') [i.e ('epsg', '4326')]
        - [CRS WKT string](https://www.ogc.org/standards/wkt-crs)
        - [PROJ string](https://proj.org/usage/quickstart.html)
        - JSON string with [PROJ parameters](https://proj.org/specifications/projjson.html)
        """


############################################
# Geometry value types
############################################

TYPE_GEO = knext.logical(gt.GeoValue)

TYPE_POINT = knext.logical(Point)
TYPE_LINE = knext.logical(LineString)
TYPE_POLYGON = knext.logical(Polygon)

TYPE_MULTI_POINT = knext.logical(MultiPoint)
TYPE_MULTI_LINE = knext.logical(MultiLineString)
TYPE_MULTI_POLYGON = knext.logical(MultiPolygon)

TYPE_GEO_COLLECTION = knext.logical(GeometryCollection)


############################################
# Column selection helper
############################################

__DEF_GEO_COL_LABEL = "Geometry column"
__DEF_GEO_COL_DESC = "Select the geometry column to use"

__CELL_TYPE_GEO = "org.knime.geospatial.core.data.cell.Geo"
__CELL_TYPE_POINT = "GeoPointCell"
__CELL_TYPE_LINE = "GeoLineCell"
__CELL_TYPE_POLYGON = "GeoPolygonCell"
__CELL_TYPE_COLLECTION = "GeoCollectionCell"
__CELL_TYPE_MULTI_POINT = "GeoMultiPointCell"
__CELL_TYPE_MULTI_LINE = "GeoMultiLineCell"
__CELL_TYPE_MULTI_POLYGON = "GeoMultiPolygonCell"


def geo_point_col_parameter(
    label: str = __DEF_GEO_COL_LABEL,
    description: str = __DEF_GEO_COL_DESC,
) -> knext.ColumnParameter:
    """
    Returns a column selection parameter that only supports GeoPoint columns.
    @return: Column selection parameter for point columns
    """
    return typed_geo_col_parameter(label, description, is_geo_point)


def geo_col_parameter(
    label: str = __DEF_GEO_COL_LABEL,
    description: str = __DEF_GEO_COL_DESC,
) -> knext.ColumnParameter:
    """
    Returns a column selection parameter that supports Geo columns.
    @return: Column selection parameter for all Geo columns
    """
    return typed_geo_col_parameter(label, description, is_geo)


def typed_geo_col_parameter(
    label: str = __DEF_GEO_COL_LABEL,
    description: str = __DEF_GEO_COL_DESC,
    type_filter: Callable[[knext.Column], bool] = None,
) -> knext.ColumnParameter:
    """
    Returns a column selection parameter that allows the user to select all columns that are compatible with
    the provided type filter.
    @return: Column selection parameter
    """
    return knext.ColumnParameter(
        label=label,
        description=description,
        column_filter=type_filter,
        include_row_key=False,
        include_none_column=False,
    )

def negate(function):
    """
    Negates the incoming function e.g. negate(is_numeric) can be used in a column parameter to allow the user 
    to select from all none numeric columns. 
    @return: the negated input function e.g. if the input function returns true this function returns false
    """
    def new_function(*args, **kwargs):
       return not function(*args, **kwargs)
    return new_function

def is_numeric(column: knext.Column) -> bool:
    """
    Checks if column is numeric e.g. int, long or double.
    @return: True if Column is numeric
    """
    return (
        column.ktype == knext.double()
        or column.ktype == knext.int32()
        or column.ktype == knext.int64()
    )


def is_string(column: knext.Column) -> bool:
    """
    Checks if column is string
    @return: True if Column is string
    """
    return column.ktype == knext.string()


def is_boolean(column: knext.Column) -> bool:
    """
    Checks if column is boolean
    @return: True if Column is boolean
    """
    return column.ktype == knext.boolean()


def is_numeric_or_string(column: knext.Column) -> bool:
    """
    Checks if column is numeric or string
    @return: True if Column is numeric or string
    """
    return column.ktype in [
        knext.double(),
        knext.int32(),
        knext.int64(),
        knext.string(),
    ]


def is_binary(column: knext.Column) -> bool:
    """
    Checks if column is binary
    @return: True if Column is binary
    """
    return column.ktype == knext.blob


def is_date(column: knext.Column) -> bool:
    """
    Checks if column is compatible to the GeoValue interface and thus returns true for all geospatial types such as:
    GeoPointCell, GeoLineCell, GeoPolygonCell, GeoMultiPointCell, GeoMultiLineCell, GeoMultiPolygonCell, ...
    @return: True if Column Type is GeoValue compatible
    """
    return __is_type_x(column, "org.knime.core.data.v2.time.LocalDateValueFactory")


def is_geo(column: knext.Column) -> bool:
    """
    Checks if column is compatible to the GeoValue interface and thus returns true for all geospatial types such as:
    GeoPointCell, GeoLineCell, GeoPolygonCell, GeoMultiPointCell, GeoMultiLineCell, GeoMultiPolygonCell, ...
    @return: True if Column Type is GeoValue compatible
    """
    return __is_type_x(column, __CELL_TYPE_GEO)


def is_geo_point(column: knext.Column) -> bool:
    """
    Checks if column is a GeoPointCell.
    @return: True if Column Type is a GeoPointCell
    """
    return __is_type_x(column, __CELL_TYPE_POINT)


def is_geo_line(column: knext.Column) -> bool:
    """
    Checks if column is a GeoLineCell.
    @return: True if Column Type is a GeoLineCell
    """
    return __is_type_x(column, __CELL_TYPE_LINE)


def is_geo_polygon(column: knext.Column) -> bool:
    """
    Checks if column is a GeoPolygonCell.
    @return: True if Column Type is a GeoPolygonCell
    """
    return __is_type_x(column, __CELL_TYPE_POLYGON)


def is_geo_collection(column: knext.Column) -> bool:
    """
    Checks if column is a GeoCollectionCell.
    @return: True if Column Type is a GeoCollectionCell
    """
    return __is_type_x(column, __CELL_TYPE_COLLECTION)


def is_geo_multi_point(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMultiPointCell.
    @return: True if Column Type is a GeoMultiPointCell
    """
    return __is_type_x(column, __CELL_TYPE_MULTI_POINT)


def is_geo_multi_line(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMultiLineCell.
    @return: True if Column Type is a GeoMultiLineCell
    """
    return __is_type_x(column, __CELL_TYPE_MULTI_LINE)


def is_geo_multi_polygon(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMultiPolygonCell.
    @return: True if Column Type is a GeoMultiPolygonCell
    """
    return __is_type_x(column, __CELL_TYPE_MULTI_POLYGON)


def __is_type_x(column: knext.Column, type: str) -> bool:
    """
    Checks if column contains the given type whereas type can be :
    GeoPointCell, GeoLineCell, GeoPolygonCell, GeoMultiPointCell, GeoMultiLineCell, GeoMultiPolygonCell, ...
    @return: True if Column Type is a GeoLogical Point
    """
    return (
        isinstance(column.ktype, knext.LogicalType)
        and type in column.ktype.logical_type
    )


############################################
# GeoPandas node class decorator
############################################
def geo_node_description(short_description: str, description: str, references: dict):
    """This decorator takes the provided information and generates a standardized node description
    for nodes that are based on GeoPandas functionality."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the [GeoPandas](https://geopandas.org/) project and uses the following related information and function"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        node_factory.__doc__ = s
        return node_factory

    return set_description


def census_node_description(short_description: str, description: str, references: dict):
    """This decorator takes the provided information and generates a standardized node description
    for nodes that are based on GeoPandas functionality."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the data from [US Census](https://www.census.gov/) and [FIPS code](https://www.census.gov/library/reference/code-lists/ansi.html) and here are related data sources and references"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        node_factory.__doc__ = s
        return node_factory

    return set_description


def osm_node_description(short_description: str, description: str, references: dict):
    """This decorator takes the provided information and generates a standardized node description
    for nodes that are based on GeoPandas functionality."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the project [OSMnx](https://github.com/gboeing/osmnx)  and here are related tools and references"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        node_factory.__doc__ = s
        return node_factory

    return set_description


def pd_node_description(short_description: str, description: str, references: dict):
    """This decorator takes the provided information and generates a standardized node description
    for nodes that are based on GeoPandas functionality."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the package [pandas](https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html)  and here are related tools and references"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        node_factory.__doc__ = s
        return node_factory

    return set_description


def pulp_node_description(short_description: str, description: str, references: dict):
    """This decorator takes the provided information and generates a standardized node description
    for nodes that are based on GeoPandas functionality."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the package [PuLP](https://coin-or.github.io/pulp/)  and here are related tools and references"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        node_factory.__doc__ = s
        return node_factory

    return set_description


############################################
# GeoPandas helper
############################################
def load_geo_data_frame(
    input_table: knext.Table,
    column: knext.Column,
    exec_context: knext.ExecutionContext = None,
    load_msg: str = "Loading Geo data frame...",
    done_msg: str = "Geo data frame loaded. Start computation...",
) -> gp.GeoDataFrame:
    """Creates a GeoDataFrame from the given input table using the provided column as geo column."""
    if exec_context:
        exec_context.set_progress(0.0, load_msg)
    gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=column)
    if exec_context:
        exec_context.set_progress(0.3, done_msg)
    return gdf


def to_table(
    gdf: gp.GeoDataFrame,
    exec_context: knext.ExecutionContext = None,
    done_msg: str = "Computation done",
) -> knext.Table:
    """Returns a KNIME table representing the given GeoDataFrame."""
    if exec_context:
        exec_context.set_progress(1.0, done_msg)
    return knext.Table.from_pandas(gdf)


############################################
# General helper
############################################

__DEF_PLEASE_SELECT_COLUMN = "Please select a column"


def column_exists_or_preset(
    context: knext.ConfigurationContext,
    column: str,
    schema: knext.Schema,
    func: Callable[[knext.Column], bool] = None,
    none_msg: str = "No compatible column found in input table",
) -> str:
    """
    Checks that the given column is not None and exists in the given schema. If none is selected it returns the
    first column that is compatible with the provided function. If none is compatible it throws an exception.
    """
    if column is None:
        for c in schema:
            if func(c):
                context.set_warning(f"Preset column to: {c.name}")
                return c.name
        raise knext.InvalidParametersError(none_msg)
    __check_col_and_type(column, schema, func)
    return column


def geo_column_exists(
    column: str,
    schema: knext.Schema,
    none_msg: str = __DEF_PLEASE_SELECT_COLUMN,
) -> None:
    """
    Checks that the given column is not None, of type geo and exists in the given schema otherwise it throws an exception
    """
    return column_exists(column, schema, is_geo, none_msg)


def column_exists(
    column: str,
    schema: knext.Schema,
    func: Callable[[knext.Column], bool] = None,
    none_msg: str = __DEF_PLEASE_SELECT_COLUMN,
) -> None:
    """
    Checks that the given column is not None and exists in the given schema otherwise it throws an exception
    """
    if column is None:
        raise knext.InvalidParametersError(none_msg)
    __check_col_and_type(column, schema, func)


def __check_col_and_type(
    column: str,
    schema: knext.Schema,
    check_type: Callable[[knext.Column], bool] = None,
) -> None:
    """
    Checks that the given column exists in the given schema and that it matches the given type_check function.
    """
    # Check that the column exists in the schema and that it has a compatible type
    try:
        existing_column = schema[column]
        if check_type is not None and not check_type(existing_column):
            raise knext.InvalidParametersError(
                f"Column '{str(column)}' has incompatible data type"
            )
    except IndexError:
        raise knext.InvalidParametersError(
            f"Column '{str(column)}' not available in input table"
        )


def columns_exist(
    columns: List[str],
    schema: knext.Schema,
    func: Callable[[knext.Column], bool] = lambda c: True,
    none_msg: str = __DEF_PLEASE_SELECT_COLUMN,
) -> None:
    """
    Checks that the given columns are not None and exist in the given schema otherwise it throws an exception
    """
    for col in columns:
        column_exists(col, schema, func, none_msg)


def fail_if_column_exists(
    column_name: str, input_schema: knext.Schema, msg: str = None
):
    """Checks that the given column name does not exists in the input schema.
    Can be used to check that a column is not accidentally overwritten."""
    if column_name in input_schema.column_names:
        if msg is None:
            msg = f"Column '{column_name}' exists"
        raise knext.InvalidParametersError(msg)


def get_unique_column_name(column_name: str, input_schema: knext.Schema) -> str:
    """Checks if the column name exists in the given schema and if so appends a number to it to make it unique.
    The unique name if returned or the original if it was already unique."""
    if column_name is None:
        raise knext.InvalidParametersError("Column name must not be None")
    uniquifier = 1
    result = column_name
    while result in input_schema.column_names:
        result = column_name + f"(#{uniquifier})"
        uniquifier += 1
    return result


def check_canceled(exec_context: knext.ExecutionContext) -> None:
    """
    Checks if the user has canceled the execution and if so throws a RuntimeException
    """
    if exec_context.is_canceled():
        raise RuntimeError("Execution canceled")


def ensure_file_extension(file_name: str, file_extension: str) -> str:
    """
    Checks if the given file_name ends with the given file_extension and if not appends it to the returned file_name.
    """
    if not file_name:
        raise knext.InvalidParametersError("Please enter a valid file name")
    if file_name.lower().endswith(file_extension):
        return file_name
    return file_name + file_extension


def Turn_all_NA_column_as_str(gdf) -> None:
    """
    Checks if the GeoDataFrame has columns of all NAs, otherwise it cannot be write out and makes the node crush.
    """
    # Transform the NA columns to string
    NotNacol = list(gdf.dropna(axis=1, how="all").columns)
    Nacol = gdf.loc[:, ~gdf.columns.isin(NotNacol)].columns.tolist()
    if len(Nacol) > 0:
        gdf[Nacol] = gdf[Nacol].astype(str)
    gdf = gdf.reset_index(drop=True)
    return gdf

def get_env_path():
    """
    Returns the path to the used Python environment e.g. the Python packages of this environment can be found in
    <RETURNED_VAL>\Lib\site-packages.
    """
    import sys
    import os.path as os

    # path to the Python executable in the used Python environment
    exec_path = sys.executable
    env_path = os.dirname(exec_path)
    return env_path