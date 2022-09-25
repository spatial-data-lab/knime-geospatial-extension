import logging
from typing import Callable
from typing import List

import geopandas as gp
import geospatial_types as gt
import knime_extension as knext

LOGGER = logging.getLogger(__name__)


############################################
# Column selection helper
############################################

DEF_GEO_COL_LABEL = "Geometry column"
DEF_GEO_COL_DESC = "Select the Geometry column to use"


def geo_point_col_parameter(
    label: str = DEF_GEO_COL_LABEL,
    description: str = DEF_GEO_COL_DESC,
):
    """
    Returns a column selection parameter that only supports GeoPoint columns.
    @return: Column selection parameter for point columns
    """
    return typed_geo_col_parameter(label, description, is_geo_point)


def geo_col_parameter(
    label: str = DEF_GEO_COL_LABEL,
    description: str = DEF_GEO_COL_DESC,
):
    """
    Returns a column selection parameter that supports Geo columns.
    @return: Column selection parameter for all Geo columns
    """
    return typed_geo_col_parameter(label, description, is_geo)


def typed_geo_col_parameter(
    label: str = DEF_GEO_COL_LABEL,
    description: str = DEF_GEO_COL_DESC,
    type_filter: Callable[[knext.Column], bool] = None,
):
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
    return is_numeric(column) or is_string(column)

def is_geo(column: knext.Column) -> bool:
    """
    Checks if column is compatible to the GeoValue interface and thus returns true for all geospatial types such as:
    GeoPointCell, GeoLineCell, GeoPolygonCell, GeoMultiPointCell, GeoMultiLineCell, GeoMultiPolygonCell, ...
    @return: True if Column Type is GeoValue compatible
    """
    return (
        isinstance(column.ktype, knext.LogicalType)
        and column.ktype.value_type == knext.logical(gt.GeoValue).value_type
    )


def is_geo_point(column: knext.Column) -> bool:
    """
    Checks if column is a GeoPointCell.
    @return: True if Column Type is a GeoPointCell
    """
    return __is_type_x(column, "GeoPointCell")


def is_geo_line(column: knext.Column) -> bool:
    """
    Checks if column is a GeoLineCell.
    @return: True if Column Type is a GeoLineCell
    """
    return __is_type_x(column, "GeoLineCell")


def is_geo_polygon(column: knext.Column) -> bool:
    """
    Checks if column is a GeoPolygonCell.
    @return: True if Column Type is a GeoPolygonCell
    """
    return __is_type_x(column, "GeoPolygonCell")


def is_geo_collection(column: knext.Column) -> bool:
    """
    Checks if column is a GeoCollectionCell.
    @return: True if Column Type is a GeoCollectionCell
    """
    return __is_type_x(column, "GeoCollectionCell")


def is_geo_multi_point(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMultiPointCell.
    @return: True if Column Type is a GeoMultiPointCell
    """
    return __is_type_x(column, "GeoMultiPointCell")


def is_geo_multi_line(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMultiLineCell.
    @return: True if Column Type is a GeoMultiLineCell
    """
    return __is_type_x(column, "GeoMultiLineCell")


def is_geo_multi_polygon(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMultiPolygonCell.
    @return: True if Column Type is a GeoMultiPolygonCell
    """
    return __is_type_x(column, "GeoMultiPolygonCell")


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
        s += "The node is based on the [GeoPandas](https://geopandas.org/) project and uses the following function"
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
def column_exists(
    col_name: str,
    schema: knext.Schema,
    func: Callable[[knext.Column], bool] = lambda c: True,
    none_msg: str = "Please select a column",
) -> None:
    """
    Checks that the given column is not None and exists in the given schema otherwise it throws an exception
    """
    if col_name is None:
        raise knext.InvalidParametersError(none_msg)
    # Check that the column exists in the schema and that it has a compatible type
    try:
        column = schema[col_name]
        if not func(column):
            raise knext.InvalidParametersError(
                f"Column '{col_name}' has incompatible data type"
            )
    except IndexError:
        raise knext.InvalidParametersError(
            f"Column '{col_name}' not available in input table"
        )


def columns_exist(
    columns: List[str],
    schema: knext.Schema,
    func: Callable[[knext.Column], bool] = lambda c: True,
    none_msg: str = "Please select a column",
) -> None:
    """
    Checks that the given columns are not None and exist in the given schema otherwise it throws an exception
    """
    for col in columns:
        column_exists(col, schema, func, none_msg)


def check_canceled(exec_context: knext.ExecutionContext) -> None:
    """
    Checks if the user has canceled the execution and if so throws a RuntimeException
    """
    if exec_context.is_canceled():
        raise RuntimeError("Execution canceled")
