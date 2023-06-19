import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


__category = knext.category(
    path="/community/geo",
    level_id="spatialstatistic",
    name="Exploratory Spatial Data Analysis",
    description="Spatial Statistic Nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialStatisticsCategory.png",
    after="viz",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialStatistics/"

__global_statistics_output_table_description = """
The output table contains the following columns:
the local statistic value for the input geometry,
the p-value for the local statistic value,
the z-score for the local statistic value.
"""
__global_statistics_interactive_view_description = """
The interactive view shows the density plot of the statistic values for permuted samples.
The red line is the value of Moran’s I
The blue line is the expected value under normality assumption
"""

__local_statistics_output_table_description = """
The output table contains the original input table with the following additional columns:
the local statistic value,
`p-value` is the p-value for the local statistic value,
`z-score` is the z-score for the local statistic value.
"""

__spots = """
`spots` are the values that indicate quadrant location 1 HH, 2 LH, 3 LL, 4 HL, 
`spots_type` has the values of HH (High-High), LH (Low-High), LL (Low-Low),
HL (High-Low), Not Significant (the p-value is greater than the significance level).
"""


def _var_col_exists_or_preset(
    context: knext.ConfigurationContext,
    field,
    schema: knext.Schema,
) -> str:
    return knut.column_exists_or_preset(
        context,
        field,
        schema,
        knut.is_numeric,
        "No compatible variable column found in input table",
    )


############################################
# Spatial Weights
############################################
@knext.node(
    name="Spatial Weights",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "SpatialWeight.png",
    category=__category,
    after="",
)
@knext.input_table(name="Geo table", description="Table with geometry column.")
@knext.output_table(name="Spatial Weights", description="Spatial Weights.")
class spatialWeights:
    """This node constructs a contiguity spatial weights matrix from the input data.
    This node constructs a contiguity spatial weights matrix from the input data.
    """

    geo_col = knut.geo_col_parameter(
        description="The name of the geometry column in the input data."
    )

    id_col = knext.ColumnParameter(
        "ID column",
        """Select the column which contains for each observation in the input data a unique ID, it should be an integer column.
        If 'none' is selected, the IDs will be automatically generated from 0 to the number of rows flowing 
        the order of the input data.
        The IDs of this column must match with the values of the ID column selected in subsequent ESDA or spatial 
        modeling nodes.
        """,
        include_none_column=True,
        column_filter=knut.is_int,
        since_version="1.1.0",
    )

    category = knext.StringParameter(
        "Weights category",
        """ The type of spatial weights to construct. Defaults to 'Queen'. Other options are 'Rook', 
        'Binary Distance Band', 'Inverse Distance', 'Lattice', 'K nearest', 'Kernel', and 
        'Get spatial weights matrix from file'.

        - `Queen` which will construct a queen contiguity weights matrix, is more robust and more suitable for areal unit data. The queen criterion is somewhat more encompassing and defines 
        neighbors as spatial units sharing a common edge or a common vertex. 
        - The `Rook` criterion defines neighbors by the existence of a common edge between two spatial units. Therefore, the number of neighbors according to the 
        queen criterion will always be at least as large as for the rook criterion. 
        - When choosing `K nearest`, select the nearest number 'Nearest k' in the following options. K-nearest are often used for point data. 
        - When selecting `Binary Distance Band`, please select 
        the distance threshold 'Threshold' in the following options. 
        - When selecting `Inverse Distance`, please select 
        the distance threshold 'Threshold' and the corresponding power 'Power' in the following options. 
        - When 'Your own' is selected, please enter the path of the spatial weights matrix in CSV format in the 
        following options. 
        - More details about spatial weights, please see the [GeoDa center website](https://geodacenter.github.io/documentation.html).
        """,
        "Queen",
        enum=[
            "Queen",
            "Rook",
            "Binary Distance Band",
            "Inverse Distance",
            "Lattice",
            "K nearest",
            "Lattice",
            "Kernel",
            "Get spatial weights matrix from file",
        ],
    )
    order = knext.IntParameter(
        "Order for Queen or Rook",
        """The order of the weight matrix is 1 by default. Users can change the order of the weights, higher order 
        weights will treat further units as neighbors.""",
        1,
    )

    Threshold = knext.IntParameter(
        "Threshold for Inverse Distance or Binary Distance Band",
        """The distance threshold for constructing binary distance band and inverse distance weights. Defaults to 1""",
        1,
    )

    Power = knext.IntParameter(
        "Power for Inverse Distance",
        """The power for constructing inverse distance weights. Defaults to 1.""",
        1,
    )

    Rows = knext.IntParameter(
        "Rows for Lattice",
        "The number of rows for constructing a lattice spatial weights matrix. Defaults to 5.",
        5,
    )
    Columns = knext.IntParameter(
        "Columns for Lattice",
        "The number of columns for constructing a lattice spatial weights matrix. Defaults to 5.",
        5,
    )

    # k = knext.IntParameter(
    #     "K for K nearest or Kernel",
    #     "k is the number of the nearest neighbor.",
    #     4,
    # )

    Nearest_k = knext.IntParameter(
        "Nearest k",
        "The number of nearest neighbors to use for constructing k-nearest neighbors weights. Defaults to 4.",
        4,
    )

    Kernel_K = knext.IntParameter(
        "Kernel K",
        "The number of nearest neighbors to use for determining the bandwidth in kernel weights. Defaults to 12.",
        12,
    )

    Kernel_type = knext.StringParameter(
        "Kernel type",
        "The type of kernel to use in constructing kernel weights. Defaults to 'triangular' ",
        "triangular",
        enum=[
            "gaussian",
            "quadratic",
            "quartic",
            "triangular",
            "uniform",
        ],
    )

    Kernel_bandwidth = knext.StringParameter(
        "Kernel bandwidth",
        "The type of kernel bandwidth to use in constructing kernel weights. The bandwidth of the kernel. The default is fixed. If adaptive then bandwidth is adaptive across observations.",
        "Fixed",
        enum=[
            "Adaptive",
            "Fixed",
        ],
    )

    Your_own_matrix_local_path = knext.StringParameter(
        "Get spatial weights matrix from file",
        """The file path of a user-defined spatial weights matrix in CSV format. Defaults to ''.
        Please enter the path of the spatial weights matrix in CSV format in the following options. 
        The weights matrix must be in matrix format and in the order of the samples. """,
        "",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)

        gdf.index = range(len(gdf))
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting projection...")

        import libpysal

        if self.category == "Rook":
            w = libpysal.weights.Rook.from_dataframe(gdf)
            wname = "Rook"
            w.transform = "r"
        if self.category == "Queen":
            w = libpysal.weights.Queen.from_dataframe(gdf)
            wname = "Queen"
            w.transform = "r"
        if self.category == "Inverse Distance":
            w = libpysal.weights.DistanceBand.from_dataframe(
                gdf, self.Threshold, alpha=-1 * self.order, binary=False
            )
            wname = "Inverse Distance"
            w.transform = "r"
        if self.category == "Binary Distance Band":

            import util.projection as kproj 

            crs = gdf.crs
            if (crs is not None) and (kproj.is_geographic(crs)):
                gdf = gdf.to_crs("EPSG:3857")
            w = libpysal.weights.DistanceBand.from_dataframe(
                gdf, self.Threshold, binary=True
            )
            wname = "Binary Distance Band"
            w.transform = "r"
        if self.category == "K nearest":
            # self.Nearest_k =  knext.IntParameter("Nearest k", "K-nearest are often used for point data. k is the number of the nearest neighbor", 4)
            w = libpysal.weights.KNN.from_dataframe(gdf, k=self.Nearest_k)
            wname = "K nearest"
            w.transform = "r"
        if self.category == "Lattice":
            w = libpysal.weights.lat2W(nrows=self.Rows, ncols=self.Columns, rook=True)
            wname = "Lattice"
            w.transform = "r"
        if self.order != 1:
            w = libpysal.weights.higher_order(w, self.order - 1)
            w.transform = "r"

        import numpy as np

        if self.category == "Get spatial weights matrix from file":
            import pandas as pd

            z = pd.read_csv(self.Your_own_matrix_local_path, header=None)
            zz = np.array(z)

            import scipy.sparse

            sparse = scipy.sparse.csr_matrix(zz)

            from libpysal.weights import WSP

            w = WSP(sparse)
            wname = "Get spatial weights matrix from file"

        if self.category == "Kernel":
            bd = False
            if self.Kernel_bandwidth == "Fixed":
                bd = True
            w = libpysal.weights.Kernel.from_dataframe(
                gdf, fixed=bd, k=self.Kernel_K, function=self.Kernel_type
            )
        # path = flow_variables["knime.workspace"] + os.sep +"spatialweights.gal"
        # w.to_file(path)
        # flow_variables["weights"] =path
        out = w.to_adjlist(drop_islands=False)

        if "none" not in str(self.id_col).lower():
            # get index id map
            id_map = gdf[self.id_col].to_dict()
            out["focal"] = out["focal"].map(id_map)
            out["neighbor"] = out["neighbor"].map(id_map)
        # gdf.to_crs(self.new_crs, inplace=True)
        exec_context.set_progress(
            0.1, "Constructs a contiguity spatial weights matrix done"
        )

        # focal and neighbor should always be int
        out["focal"] = out["focal"].astype(int)
        out["neighbor"] = out["neighbor"].astype(int)

        # focal and neighbor should always be int
        out["focal"] = out["focal"].astype(int)
        out["neighbor"] = out["neighbor"].astype(int)

        return knext.Table.from_pandas(out)


@knext.parameter_group(label="Variable Setting")
class VariableSetting:
    """
    Select the variable you want to use for the analysis.
    """

    Field_col = knext.ColumnParameter(
        "Variable column",
        "The variable column you want to use for the analysis",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )


@knext.parameter_group(label="ID Setting", since_version="1.1.0")
class IDSetting:
    """
    The unique ID column. The values need to match the values from the ID column selected in the
    [Spatial Weights](https://hub.knime.com/center%20for%20geographic%20analysis%20at%20harvard%20university/extensions/sdl.harvard.features.geospatial/latest/org.knime.python3.nodes.extension.ExtensionNodeSetFactory$DynamicExtensionNodeFactory:4d710eae/) node.
    The selected column must contain unique IDs for each observation in the input data of type integer.
    """

    Field_col = knext.ColumnParameter(
        "ID column",
        """The selected column should contain unique IDs for each observation in the input data and should be of 
        type integer. The values need to match the values from the ID column selected in the
        [Spatial Weights](https://hub.knime.com/center%20for%20geographic%20analysis%20at%20harvard%20university/extensions/sdl.harvard.features.geospatial/latest/org.knime.python3.nodes.extension.ExtensionNodeSetFactory$DynamicExtensionNodeFactory:4d710eae/) node.
        If you selected 'none' in the Spatial Weights node select it here as well.""",
        column_filter=knut.is_int,
        include_none_column=True,
    )


############################################
# Global Moran's I node
############################################


@knext.node(
    name="Global Moran's I",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GlobalMoran.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Global Moran's I.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Global Moran's I.",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Global Moran's I. "
    + __global_statistics_output_table_description,
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Moran's I",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="Output View",
    description="Output view of Global Moran's I. "
    + __global_statistics_interactive_view_description,
)
class GlobalMoransI:
    """Global Moran's I.
    Moran’s I Global Autocorrelation Statistic.
    """

    # input parameters
    geo_col = knut.geo_col_parameter(
        description="The column containing the geometry to use for the spatial weights matrix."
    )

    id_col_setting = IDSetting()

    variable_setting = VariableSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        self.variable_setting.Field_col = _var_col_exists_or_preset(
            configure_context, self.variable_setting.Field_col, input_schema_1
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col_setting.Field_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col_setting.Field_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        y = gdf[self.variable_setting.Field_col]

        import numpy as np

        np.random.seed(12345)

        import esda

        mi = esda.moran.Moran(y, w)
        result = {"Moran's I": mi.I, "p-value": mi.p_norm, "z-score": mi.z_norm}

        import pandas as pd

        out = pd.DataFrame(result, index=[0])

        import seaborn as sbn

        ax = sbn.kdeplot(mi.sim, shade=True)

        import matplotlib.pyplot as plt

        plt.vlines(mi.I, 0, 1, color="r")
        plt.vlines(mi.EI, 0, 1)
        plt.xlabel("Moran's I")

        return knext.Table.from_pandas(out), knext.view_matplotlib(ax.get_figure())


############################################
# Local Moran's I node
############################################


@knext.node(
    name="Local Moran's I",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "LocalMoran.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Local Moran's I",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Local Moran's I",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Local Moran's I. "
    + __local_statistics_output_table_description
    + __spots,
)
@knext.output_view(
    name="Output View",
    description="The scatter plot of Local Moran's I",
)
class LocalMoransI:
    """Local Moran's I.
    Local Moran's I Statistics.
    """

    # input parameters
    geo_col = knut.geo_col_parameter(
        description="The column containing the geometry to use for local Moran's I.",
    )

    id_col_setting = IDSetting()

    variable_setting = VariableSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        self.variable_setting.Field_col = _var_col_exists_or_preset(
            configure_context, self.variable_setting.Field_col, input_schema_1
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col_setting.Field_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col_setting.Field_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        y = gdf[self.variable_setting.Field_col]

        import numpy as np

        np.random.seed(12345)

        import esda

        li = esda.moran.Moran_Local(y, w)

        # gdf.loc[:,"spots_type"] = gdf["spots_type"].fillna("Not Significant")

        gdf.loc[:, "Moran's I"] = li.Is
        gdf.loc[:, "p-value"] = li.p_sim
        gdf.loc[:, "z-score"] = li.z_sim
        gdf.loc[:, "spots"] = li.q

        gdf.loc[:, "spots_type"] = gdf["spots"].replace(
            {1: "HH", 2: "LH", 3: "LL", 4: "HL"}
        )
        gdf.loc[gdf["p-value"] > 0.05, "spots_type"] = "Not Significant"
        # out = pd.merge(gdf, out, left_index=True, right_index=True)

        import pysal.lib as lps

        lag_index = lps.weights.lag_spatial(w, gdf[self.variable_setting.Field_col])
        index_v = gdf[self.variable_setting.Field_col]

        import numpy as np
        import matplotlib.pyplot as plt

        b, a = np.polyfit(index_v, lag_index, 1)
        f, ax = plt.subplots(1, figsize=(9, 9))

        plt.plot(index_v, lag_index, ".", color="firebrick")

        # dashed vert at mean of the index_v
        plt.vlines(index_v.mean(), lag_index.min(), lag_index.max(), linestyle="--")
        # dashed horizontal at mean of lagged index_v
        plt.hlines(lag_index.mean(), index_v.min(), index_v.max(), linestyle="--")

        # red line of best fit using global I as slope
        plt.plot(index_v, a + b * index_v, "r")
        plt.title("Moran Scatterplot")
        plt.ylabel("Spatial Lag of %s" % self.variable_setting.Field_col)
        plt.xlabel("%s" % self.variable_setting.Field_col)

        # out.drop(columns=[self.geo_col], inplace=True)
        # out.reset_index(inplace=True)

        # return knext.Table.from_pandas(out)
        return knext.Table.from_pandas(gdf), knext.view_matplotlib(f)


############################################
# Global  Geary’s C node
############################################


@knext.node(
    name="Global Geary’s C",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GlobalGeary.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Global Geary’s C",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Global Geary’s C",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Global Geary’s C. "
    + __global_statistics_output_table_description,
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Geary’s C",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="Output View",
    description="Output view of Global Geary’s C. "
    + __global_statistics_interactive_view_description,
)
class GlobalGearysC:
    """Global Geary’s C.
    Global Geary C Autocorrelation statistic.
    """

    # input parameters
    geo_col = knut.geo_col_parameter(
        description="The column containing the geometry to use for global Geary’s C.",
    )

    id_col_setting = IDSetting()

    variable_setting = VariableSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        self.variable_setting.Field_col = _var_col_exists_or_preset(
            configure_context, self.variable_setting.Field_col, input_schema_1
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col_setting.Field_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col_setting.Field_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        y = gdf[self.variable_setting.Field_col]

        import numpy as np

        np.random.seed(12345)

        import esda

        gc = esda.geary.Geary(y, w)

        import pandas as pd

        out = pd.DataFrame(
            {"Geary's C": [gc.C], "p-value": [gc.p_sim], "z-score": [gc.z_sim]}
        )
        out.reset_index(inplace=True)

        import seaborn as sbn
        import matplotlib.pyplot as plt

        ax = sbn.kdeplot(gc.sim, shade=True)
        plt.vlines(gc.C, 0, 1, color="r")
        plt.vlines(gc.EC, 0, 1)
        plt.xlabel("Geary's C")

        return knext.Table.from_pandas(out), knext.view_matplotlib(ax.get_figure())


############################################
# Global Getis-Ord node
############################################


@knext.node(
    name="Global Getis-Ord G",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GlobalGetis.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Global Getis-Ord",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Global Getis-Ord",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Global Getis-Ord. "
    + __global_statistics_output_table_description,
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Getis-Ord",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="Output View",
    description="Output view of Global Getis-Ord. "
    + __global_statistics_interactive_view_description,
)
class GlobalGetisOrd:
    """Global Getis-Ord G.
    Global G Autocorrelation Statistic.
    """

    # input parameters
    geo_col = knut.geo_col_parameter(
        description="The column containing the geometry to use for global Getis-Ord.",
    )

    id_col_setting = IDSetting()

    variable_setting = VariableSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        self.variable_setting.Field_col = _var_col_exists_or_preset(
            configure_context, self.variable_setting.Field_col, input_schema_1
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col_setting.Field_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col_setting.Field_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        y = gdf[self.variable_setting.Field_col]

        import numpy as np

        np.random.seed(12345)

        import esda

        go = esda.getisord.G(y, w)

        import pandas as pd

        out = pd.DataFrame(
            {"Getis-Ord G": [go.G], "p-value": [go.p_sim], "z-score": [go.z_sim]}
        )
        out.reset_index(inplace=True)

        import seaborn as sbn
        import matplotlib.pyplot as plt

        ax = sbn.kdeplot(go.sim, shade=True)
        plt.vlines(go.G, 0, 1, color="r")
        plt.vlines(go.EG, 0, 1)
        plt.xlabel("Getis-Ord G")

        return knext.Table.from_pandas(out), knext.view_matplotlib(ax.get_figure())


############################################
# Local Getis-Ord node
############################################


@knext.node(
    name="Local Getis-Ord G",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "LocalGeary.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Local Getis-Ord",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Local Getis-Ord",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Local Getis-Ord. "
    + __local_statistics_output_table_description
    + "`Standardized Gs` is the standardization of Gs."
    + __spots,
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Local Getis-Ord",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="Output View",
    description="The scatter plot Local Getis-Ord",
)
class LocalGetisOrd:
    """Local Getis-Ord.
    Local Getis-Ord G Statistics.
    """

    # input parameters

    geo_col = knut.geo_col_parameter(
        description="The column containing the geometry to use for local Getis-Ord.",
    )

    id_col_setting = IDSetting()

    variable_setting = VariableSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        self.variable_setting.Field_col = _var_col_exists_or_preset(
            configure_context, self.variable_setting.Field_col, input_schema_1
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col_setting.Field_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col_setting.Field_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        y = gdf[self.variable_setting.Field_col]

        import numpy as np

        np.random.seed(12345)

        import esda

        lo = esda.getisord.G_Local(y, w)

        gdf.loc[:, "Local Getis-Ord G"] = lo.Gs
        gdf.loc[:, "p-value"] = lo.p_sim
        gdf.loc[:, "z-score"] = lo.z_sim

        gdf.loc[:, "Standardized Gs"] = lo.Zs

        gdf.loc[gdf["Standardized Gs"] > 0, "spots_type"] = "HH"
        gdf.loc[gdf["Standardized Gs"] < 0, "spots_type"] = "LL "
        gdf.loc[gdf["Standardized Gs"] > 0, "spots"] = 1
        gdf.loc[gdf["Standardized Gs"] < 0, "spots"] = 3
        gdf.loc[gdf["p-value"] > 0.05, "spots_type"] = "Not Significant"

        import pysal.lib as lps

        lag_index = lps.weights.lag_spatial(w, gdf[self.variable_setting.Field_col])
        index_v = gdf[self.variable_setting.Field_col]
        b, a = np.polyfit(index_v, lag_index, 1)

        import matplotlib.pyplot as plt

        f, ax = plt.subplots(1, figsize=(9, 9))

        plt.plot(index_v, lag_index, ".", color="firebrick")
        plt.vlines(index_v.mean(), lag_index.min(), lag_index.max(), linestyle="--")
        plt.hlines(lag_index.mean(), index_v.min(), index_v.max(), linestyle="--")

        plt.plot(index_v, a + b * index_v, "r")
        plt.xlabel(self.variable_setting.Field_col)
        plt.ylabel("Spatial Lag of " + self.variable_setting.Field_col)
        plt.title("Local Getis-Ord G Scatterplot")

        return knext.Table.from_pandas(gdf), knext.view_matplotlib(f)
