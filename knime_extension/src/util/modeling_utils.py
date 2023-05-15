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
