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
__DEF_GEO_COL_DESC = "Select the geometry column to use."

__CELL_TYPE_GEO = "org.knime.geospatial.core.data.cell.Geo"
__CELL_TYPE_POINT = "GeoPointCell"
__CELL_TYPE_LINE = "GeoLineCell"
__CELL_TYPE_POLYGON = "GeoPolygonCell"
__CELL_TYPE_COLLECTION = "GeoCollectionCell"
__CELL_TYPE_MULTI_POINT = "GeoMultiPointCell"
__CELL_TYPE_MULTI_LINE = "GeoMultiLineCell"
__CELL_TYPE_MULTI_POLYGON = "GeoMultiPolygonCell"

# The default request header to use in all nodes that perform a web request
WEB_REQUEST_HEADER = {"User-Agent": "KNIME-Geospatial/1.1"}


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


def boolean_or(*functions):
    """
    Return True if any of the given functions returns True
    @return: True if any of the functions returns True
    """

    def new_function(*args, **kwargs):
        return any(f(*args, **kwargs) for f in functions)

    return new_function


def boolean_and(*functions):
    """
    Return True if all of the given functions return True
    @return: True if all of the functions return True
    """

    def new_function(*args, **kwargs):
        return all(f(*args, **kwargs) for f in functions)

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


def is_int(column: knext.Column) -> bool:
    """
    Checks if column is integer.
    @return: True if Column is integer
    """
    return column.ktype in [
        knext.int32(),
        knext.int64(),
    ]


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
    return boolean_or(is_numeric, is_string)(column)


def is_int_or_string(column: knext.Column) -> bool:
    """
    Checks if column is int or string
    @return: True if Column is numeric or string
    """
    return column.ktype in [
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
    for nodes that are based on US Census data."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the data from [US Census](https://www.census.gov/) and "
        s += "[FIPS code](https://www.census.gov/library/reference/code-lists/ansi.html) "
        s += "and uses the following related information and function"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        s += "\n\n##Note\n"
        s += "This node uses the Census Bureau Data API but is not endorsed or certified by the Census Bureau. "
        s += "For the terms of service click [here.](https://www.census.gov/data/developers/about/terms-of-service.html)"
        node_factory.__doc__ = s
        return node_factory

    return set_description


def osm_node_description(short_description: str, description: str, references: dict):
    """This decorator takes the provided information and generates a standardized node description
    for nodes that are based on OpenStreetMap data. It also places the proper copyright notice which
    is necessary."""

    def set_description(node_factory):
        s = f"{short_description}\n"
        s += f"{description}\n\n"
        # s += "___\n\n"  # separator line between description and general part
        s += "The node is based on the [OpenStreetMap project](https://www.openstreetmap.org/about) "
        s += "and uses the following related information and function"
        if references is not None:
            if len(references) > 1:
                s += "s"
            s += ":"
            s += "\n\n"
            for key in references:
                s += f"- [{key}]({references[key]})\n"
        s += "\n\n##Note\n"
        s += (
            "Data copyright by [OpenStreetMap](https://www.openstreetmap.org/copyright)"
        )
        s += "[(ODbl)](https://opendatacommons.org/licenses/odbl/index.html) and provided under "
        s += "[CC-BY-SA.](https://creativecommons.org/licenses/by-sa/2.0/)"
        s += "To report a problem and contribute to OpenStreetMap click [here.](https://www.openstreetmap.org/fixthemap)"
        s += "Please note the OpenStreetMap licence and attribution guidelines as described "
        s += "[here.](https://wiki.osmfoundation.org/wiki/Licence/Attribution_Guidelines)"
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
        exec_context.set_progress(0.1, done_msg)
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
    return get_unique_name(column_name, input_schema.column_names)
    # if column_name is None:
    #     raise knext.InvalidParametersError("Column name must not be None")
    # uniquifier = 1
    # result = column_name
    # while result in input_schema.column_names:
    #     result = column_name + f"(#{uniquifier})"
    #     uniquifier += 1
    # return result


def get_unique_name(column_name: str, existing_col_names) -> str:
    """Checks if the column name exists in the given schema and if so appends a number to it to make it unique.
    The unique name if returned or the original if it was already unique."""
    if column_name is None:
        raise knext.InvalidParametersError("Column name must not be None")
    uniquifier = 1
    result = column_name
    while result in existing_col_names:
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


class ResultSettingsMode(knext.EnumParameterOptions):
    REPLACE = (
        "Replace",
        "Replace the selected input column with the result.",
    )
    APPEND = (
        "Append",
        "Append a new column with the name provided below.",
    )

    @classmethod
    def get_default(cls):
        return cls.REPLACE


@knext.parameter_group(label="Output", since_version="1.1.0")
class ResultSettings:
    """
    Group of settings that define the format of the result table.
    """

    mode = knext.EnumParameter(
        label="Output column",
        description="Choose where to place the result column:",
        default_value=ResultSettingsMode.get_default().name,
        enum=ResultSettingsMode,
    )

    new_column_name = knext.StringParameter(
        "New column name",
        "The name of the new column that is appended if 'Append' is selected.",
        default_value="geometry",
    )

    def __init__(self, mode=ResultSettingsMode.get_default().name, new_name="geometry"):
        self.mode = mode
        self.new_column_name = new_name

    def get_result_schema(
        self,
        configure_context: knext.ConfigurationContext,
        schema: knext.Schema,
        selected_col: knext.Column,
        result_type,
    ) -> knext.Schema:
        """
        Either replaces the selected column or appends a new column to the end.
        """
        if self.mode == ResultSettingsMode.REPLACE.name:
            col_names = schema.column_names
            i = 0
            while i < len(col_names):
                if col_names[i] == selected_col:
                    result_schema = schema.remove(i)
                    return result_schema.insert(
                        knext.Column(result_type, selected_col), i
                    )
                i += 1
            raise knext.InvalidParametersError(
                f"Selected column '{selected_col}' not found"
            )
        # make sure the appended column is unique
        result_col = get_unique_column_name(self.new_column_name, schema)
        return schema.append(knext.Column(result_type, result_col))

    def get_result_table(
        self,
        exec_context: knext.ExecutionContext,
        gdf: gp.GeoDataFrame,
        selected_col: knext.Column,
        result_col: str,
    ) -> gp.GeoDataFrame:
        """
        Assumes that the result_col and the select_col are part of the input data frame.
        The (altered) input data frame is returned.
        """
        if self.mode == ResultSettingsMode.REPLACE.name:
            check_canceled(exec_context)
            exec_context.set_progress(0.9, "Replace input column with result column")
            gdf[selected_col] = gdf[result_col]
            gdf.drop(result_col, axis=1, inplace=True)
            gdf.rename(columns={self.new_column_name: selected_col}, inplace=True)
            return gdf
        return gdf

    def get_computed_result_table(
        self,
        exec_context: knext.ExecutionContext,
        input_table: knext.Table,
        selected_col: knext.Column,
        func: Callable,
    ) -> knext.Table:
        """
        Uses the given function to either append a new column or replace the existing column to the given input table and
        returns the result as a table depending on the user chosen settings.
        """
        gdf = load_geo_data_frame(input_table, selected_col, exec_context)
        gdf = self.get_computed_result_frame(
            exec_context, input_table.schema, gdf, selected_col, func
        )
        return to_table(gdf, exec_context)

    def get_computed_result_frame(
        self,
        exec_context: knext.ExecutionContext,
        schema: knext.Schema,
        gdf: gp.GeoDataFrame,
        selected_col: knext.Column,
        func: Callable,
    ) -> knext.Table:
        """
        Uses the given function to either append a new column or replace the existing column to the given input
        GeoDataFrame and returns the result as a table depending on the user chosen settings.
        """
        result_col = selected_col
        if self.mode == ResultSettingsMode.APPEND.name:
            result_col = get_unique_column_name(self.new_column_name, schema)
        gdf[result_col] = gdf.apply(lambda l: func(l[selected_col]), axis=1)
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


# def re_order_weight_rows(gdf, adjust_list, id_col):
#     """
#     Reorder the spatial weight according to the id_col.
#     """
#     # Reorder the rows based on the weight column
#     import libpysal

#     gdf.index = range(len(gdf))
#     w_ref = libpysal.weights.Rook.from_dataframe(gdf)
#     id_map = gdf[id_col].to_dict()
#     w_ref.transform = "r"
#     adjust_list_ref = w_ref.to_adjlist()
#     for k, row in adjust_list_ref.iterrows():
#         focal = id_map[row["focal"]]
#         neighbor = id_map[row["neighbor"]]
#         row["weight"] = adjust_list[
#             (adjust_list["focal"] == focal) & (adjust_list["neighbor"] == neighbor)
#         ]["weight"]

#     return adjust_list_ref
