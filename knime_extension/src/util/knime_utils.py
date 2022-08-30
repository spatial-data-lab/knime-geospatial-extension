import logging
from typing import List

import geospatial_types as gt
import knime_extension as knext

LOGGER = logging.getLogger(__name__)

############################################
# Column selection helper
############################################
def is_numeric(column: knext.Column) -> bool:
    """
    Checks if column is numeric e.g. int, long or double.
    @return: True if Column is numeric
    """
    return column.ktype == knext.double()

def is_string(column: knext.Column) -> bool:
    """
    Checks if column is string 
    @return: True if Column is string
    """
    return column.ktype == knext.string()


def is_geo(column: knext.Column) -> bool:
    """
    Checks if column is compatible to the GeoValue interface and thus returns true for allgeospatial types such as:
    GeoPointCell, GeoLineCell, GeoPolygonCell, GeoMultiPointCell, GeoMultiLineCell, GeoMultiPolygonCell, ...
    @return: True if Column Type is GeoValue compatible
    """
    return isinstance(column.ktype, knext.LogicalType) and column.ktype.value_type == knext.logical(gt.GeoValue).value_type


def is_geo_point(column: knext.Column) -> bool:
    """
    Checks if column is a GeoPointCell.
    @return: True if Column Type is a GeoPointCell
    """
    return __is_geo_x(column, "GeoPointCell")


def is_geo_line(column: knext.Column) -> bool:
    """
    Checks if column is a GeoLineCell.
    @return: True if Column Type is a GeoLineCell
    """
    return __is_geo_x(column, "GeoLineCell")


def is_geo_polygon(column: knext.Column) -> bool:
    """
    Checks if column is a GeoPolygonCell.
    @return: True if Column Type is a GeoPolygonCell
    """
    return __is_geo_x(column, "GeoPolygonCell")


def is_geo_collection(column: knext.Column) -> bool:
    """
    Checks if column is a GeoCollectionCell.
    @return: True if Column Type is a GeoCollectionCell
    """
    return __is_geo_x(column, "GeoCollectionCell")


def is_geo_multi_point(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMulitPointCell.
    @return: True if Column Type is a GeoMultiPointCell
    """
    return __is_geo_x(column, "GeoMultiPointCell")


def is_geo_multi_line(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMulitLineCell.
    @return: True if Column Type is a GeoMultiLineCell
    """
    return __is_geo_x(column, "GeoMultiLineCell")


def is_geo_multi_polygon(column: knext.Column) -> bool:
    """
    Checks if column is a GeoMulitPolygonCell.
    @return: True if Column Type is a GeoMultiPolygonCell
    """
    return __is_geo_x(column, "GeoMultiPolygonCell")


def __is_geo_x(column: knext.Column, type: str) -> bool:
    """
    Checks if column contains the given type whereas type can be :
    GeoPointCell, GeoLineCell, GeoPolygonCell, GeoMultiPointCell, GeoMultiLineCell, GeoMultiPolygonCell, ...
    @return: True if Column Type is a GeoLogical Point
    """
    return isinstance(column.ktype, knext.LogicalType) and type in column.ktype.logical_type


############################################
# General helper
############################################
def column_exists(column: knext.Column, schema: knext.Schema, none_msg : str = "Please select a column") -> None:
    """
    Checks that the given column is not None and exists in the given schema otherwise it throws an exception
    """
    if column is None:
        raise knext.InvalidParametersError(none_msg)
    #Check that the column exists in the schema     
    if column not in schema.column_names:
        raise knext.InvalidParametersError(f"Column '{column}' not available in input table")

def columns_exist(columns: List[knext.Column], schema: knext.Schema, none_msg : str = "Please select a column") -> None:
    """
    Checks that the given columns are not None and exist in the given schema otherwise it throws an exception
    """
    for col in columns:
        column_exists(col, schema, none_msg)


def check_canceled(exec_context: knext.ExecutionContext) -> None:
    """
    Checks if the user has canceled the execution and if so throws a RuntimeException
    """
    if exec_context.is_canceled():
        raise RuntimeError("Execution canceled")