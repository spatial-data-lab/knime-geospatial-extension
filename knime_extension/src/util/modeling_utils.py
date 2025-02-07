import knime_extension as knext
import util.knime_utils as knut


def re_order_weight_rows(gdf, adjust_list, id_col):
    """
    Reorder the spatial weight according to the id_col.
    """
    # Reorder the rows based on the weight column
    import libpysal

    gdf.index = range(len(gdf))
    w_ref = libpysal.weights.Rook.from_dataframe(gdf)
    id_map = gdf[id_col].to_dict()
    w_ref.transform = "r"
    adjust_list_ref = w_ref.to_adjlist()
    for k, row in adjust_list_ref.iterrows():
        focal = id_map[row["focal"]]
        neighbor = id_map[row["neighbor"]]
        row["weight"] = adjust_list[
            (adjust_list["focal"] == focal) & (adjust_list["neighbor"] == neighbor)
        ]["weight"]

    return adjust_list_ref


# get w from adjust list
def get_w_from_adjust_list(adj_list, gdf, id_col):
    """
    Create a W object from an adjacency list.
    """

    from libpysal.weights import W

    adjust_list = adj_list.to_pandas()
    if "none" not in str(id_col).lower():
        gdf.index = range(len(gdf))
        id_map = dict(zip(gdf[id_col], gdf.index))
        adjust_list["focal"] = adjust_list["focal"].map(id_map)
        adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)
    w = W.from_adjlist(adjust_list)

    return w, gdf


def get_id_col_parameter(
    label: str = "ID column",
    description: str = """Select the column which contains for each observation in the input data a unique ID, it should be an integer column.
    The IDs must match with the values of the 
    [Spatial Weights node](https://hub.knime.com/center%20for%20geographic%20analysis%20at%20harvard%20university/extensions/sdl.harvard.features.geospatial/latest/org.knime.python3.nodes.extension.ExtensionNodeSetFactory$DynamicExtensionNodeFactory:4d710eae/)
    ID column.
    If 'none' is selected, the IDs will be automatically generated from 0 to the number of rows following the order of 
    the first input table.
    """,
    since_version: str = "1.1.0",
):
    """
    Returns the unique ID column. It should always keep the same as the ID column in the spatial weights matrix node.
    The selected column should contain unique IDs for each observation in the input data.
    """
    return knext.ColumnParameter(
        label=label,
        description=description,
        include_none_column=True,
        column_filter=knut.is_long,
        since_version=since_version,
    )


def get_dependent_and_independent_variables():
    """
    Returns the dependent and independent variables.
    """
    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the model.",
        column_filter=knut.is_numeric,
    )

    return (dependent_variable, independent_variables)
