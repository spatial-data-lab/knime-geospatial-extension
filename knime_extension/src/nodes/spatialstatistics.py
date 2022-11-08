import knime_extension as knext
import util.knime_utils as knut
import pandas as pd
import geopandas as gp
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW
import numpy as np
import libpysal
import scipy.sparse
from libpysal.weights import WSP
import esda
import pysal.lib as lps
import pickle
import seaborn as sbn
import matplotlib.pyplot as plt
from libpysal.weights import W
import spreg
from io import StringIO
import sys


__category = knext.category(
    path="/geo",
    level_id="spatialstatistic",
    name="Exploratory Spatial Data Analysis",
    description="Spatial Statistic Nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialStatisticsCategory.png",
    after="viz",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialStatistics/"

############################################
# Spatial Weights
############################################
@knext.node(
    name="Spatial Weights",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "SpatialWeight.png",
    category=__category,
)
@knext.input_table(name="Geo table", description="Table with geometry column")
@knext.output_table(name="Spatial Weights", description="Spatial Weights")
class spatialWeights:
    """

    This node constructs a contiguity spatial weights matrix from the input data.
    """

    # FIXME:
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    category = knext.StringParameter(
        "Weights category",
        "The default value is ‘Queen’ which will construct a queen contiguity weights matrix. Queen weights is more robust and more suitable for areal unit data. The queen criterion is somewhat more encompassing and defines neighbors as spatial units sharing a common edge or a common vertex. The rook criterion defines neighbors by the existence of a common edge between two spatial units. Therefore, the number of neighbors according to the queen criterion will always be at least as large as for the rook criterion. When choosing k-nearest, select the nearest number 'Nearest k' in the following options. When selecting Binary Distance Band, please select the distance threshold 'Threshold' in the following options. When selecting Inverse Distance, please select the distance threshold 'Threshold' and the corresponding power 'Power' in the following options. When 'Your own' is selected, please enter the path of the spatial weights matrix in CSV format in the following options. More details about spatial weights, please see the GeoDa center website: https://geodacenter.github.io/documentation.html",
        "Queen",
        enum=[
            "Queen",
            "Rook",
            "Inverse Distance",
            "Binary Distance Band",
            "K nearest",
            "Lattice",
            "Kernel",
            "Your own",
        ],
    )
    order = knext.IntParameter(
        "Order",
        "The order of the weight matrix is 1 by default. Users can change the order of the weights, higher order weights will treat further units as neighbors.",
        1,
    )

    Nearest_k = knext.IntParameter(
        "Nearest k",
        "K-nearest are often used for point data. k is the number of the nearest neighbor",
        4,
    )

    Threshold = knext.IntParameter(
        "Threshold",
        "Distance band weights are often used for point data. The weights within the threshold are 1 and otherwise 0. Inverse distance weights are often used for point data. The weights within the threshold are distance^-power, and otherwise 0. The distance is Euclidean distance.",
        1,
    )
    Power = knext.IntParameter(
        "Power",
        "Distance band weights are often used for point data. The weights within the threshold are 1 and otherwise 0. Inverse distance weights are often used for point data. The weights within the threshold are distance^-power, and otherwise 0. The distance is Euclidean distance.",
        1,
    )
    Rows = knext.IntParameter(
        "Rows", "Please choose your rows and colunns of your lattice.", 5
    )
    Columns = knext.IntParameter(
        "Columns", "Please choose your rows and columns of your lattice.", 5
    )
    Your_own_matrix_local_path = knext.StringParameter(
        "Your own matrix local path",
        "Please enter the path of the spatial weights matrix in CSV format in the following options. The weights matrix must be in matrix format and in the order of the samples.",
        "",
    )
    # FIXME:
    Kernel_type = knext.StringParameter(
        "Kernel type",
        " ",
        "triangular",
        enum=["triangular", "uniform", "quadratic", "quartic", "gaussian"],
    )
    Kernel_K = knext.IntParameter("Kernel K", " ", 12)
    # FIXME:
    Kernel_bandwidth = knext.StringParameter(
        "Kernel bandwidth", " ", "Fixed", enum=["Fixed", "Adaptive"]
    )

    # new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "3857")
    # new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "3857")
    # new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "3857")
    # new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "3857")
    # new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "3857")
    # new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "3857")
    # li = knext.ColumnParameter()

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting projection...")

        # my functions here TODO:

        if self.category == "Rook":
            w = libpysal.weights.Rook.from_dataframe(gdf)
            wname = "Rook"
            w.transform = "r"
            # FIXME:
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

        if self.category == "Your own":
            z = pd.read_csv(self.Your_own_matrix_local_path, header=None)
            zz = np.array(z)
            sparse = scipy.sparse.csr_matrix(zz)
            w = WSP(sparse)
            wname = "Your own"

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
        out = w.to_adjlist()

        # gdf=gdf.to_crs(self.new_crs)
        exec_context.set_progress(
            0.1, "Constructs a contiguity spatial weights matrix done"
        )

        return knext.Table.from_pandas(out)


############################################
# Global Moran's I node
############################################
#


@knext.node(
    name="Global Moran's I",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GlobalMoran.png",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Global Moran's I",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Global Moran's I",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Global Moran's I",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Moran's I",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="output view",
    description="Output view of Global Moran's I",
)
class GlobalMoransI:
    """
    Global Moran's I
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for the spatial weights matrix.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col = knext.ColumnParameter(
        "Field column",
        "The column containing the field to use for the calculation of Global Moran's I.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)

        y = gdf[self.Field_col]
        np.random.seed(12345)
        mi = esda.moran.Moran(y, w)
        result = {"Moran's I": mi.I, "p-value": mi.p_norm, "z-score": mi.z_norm}
        out = pd.DataFrame(result, index=[0])

        ax = sbn.kdeplot(mi.sim, shade=True)
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
    description="Output table results of Local Moran's I",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Moran's I",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="output view",
    description="Output view of Local Moran's I",
)
class LocalMoransI:
    """
    Local Moran's I
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for local Moran's I.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col = knext.ColumnParameter(
        "Field column",
        "The column containing the field to use for the calculation of Local Moran's I.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)

        y = gdf[self.Field_col]
        np.random.seed(12345)
        li = esda.moran.Moran_Local(y, w)

        # gdf.loc[:,"spots_type"] = gdf["spots_type"].fillna("Not Significant")

        gdf.loc[:, "Moran's I"] = li.Is
        gdf.loc[:, "p-value"] = li.p_sim
        gdf.loc[:, "z-score"] = li.z_sim
        gdf.loc[:, "spots"] = li.q

        gdf.loc[:, "spots_type"] = gdf["spots"].replace(
            {1: "HH", 2: "LL", 3: "LH", 4: "HL"}
        )
        gdf.loc[gdf["p-value"] < 0.05, "spots_type"] = "Not Significant"
        # out = pd.merge(gdf, out, left_index=True, right_index=True)

        lag_index = lps.weights.lag_spatial(w, gdf[self.Field_col])
        index_v = gdf[self.Field_col]
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
        plt.ylabel("Spatial Lag of %s" % self.Field_col)
        plt.xlabel("%s" % self.Field_col)

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
    description="Output table results of Global Geary’s C",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Geary’s C",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="output view",
    description="Output view of Global Geary’s C",
)
class GlobalGearysC:
    """
    Global Geary’s C
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for global Geary’s C.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col = knext.ColumnParameter(
        "Field column",
        "The column containing the field to use for the calculation of Global Geary’s C.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)

        y = gdf[self.Field_col]
        np.random.seed(12345)
        gc = esda.geary.Geary(y, w)

        out = pd.DataFrame(
            {"Geary's C": [gc.C], "p-value": [gc.p_sim], "z-score": [gc.z_sim]}
        )
        out.reset_index(inplace=True)

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
    description="Output table results of Global Getis-Ord",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Global Getis-Ord",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="output view",
    description="Output view of Global Getis-Ord",
)
class GlobalGetisOrd:
    """
    Global Getis-Ord
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for global Getis-Ord.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col = knext.ColumnParameter(
        "Field column",
        "The column containing the field to use for the calculation of Global Getis-Ord.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)

        y = gdf[self.Field_col]
        np.random.seed(12345)
        go = esda.getisord.G(y, w)

        out = pd.DataFrame(
            {"Getis-Ord G": [go.G], "p-value": [go.p_sim], "z-score": [go.z_sim]}
        )
        out.reset_index(inplace=True)

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
    description="Output table results of Local Getis-Ord",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Local Getis-Ord",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="output view",
    description="Output view of Local Getis-Ord",
)
class LocalGetisOrd:
    """
    Local Getis-Ord
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for local Getis-Ord.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col = knext.ColumnParameter(
        "Field column",
        "The column containing the field to use for the calculation of Local Getis-Ord.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)

        y = gdf[self.Field_col]
        np.random.seed(12345)
        lo = esda.getisord.G_Local(y, w)

        gdf.loc[:, "Local Getis-Ord G"] = lo.Gs
        gdf.loc[:, "p-value"] = lo.p_sim
        gdf.loc[:, "z-score"] = lo.z_sim

        lag_index = lps.weights.lag_spatial(w, gdf[self.Field_col])
        index_v = gdf[self.Field_col]
        b, a = np.polyfit(index_v, lag_index, 1)
        f, ax = plt.subplots(1, figsize=(9, 9))

        plt.plot(index_v, lag_index, ".", color="firebrick")
        plt.vlines(index_v.mean(), lag_index.min(), lag_index.max(), linestyle="--")
        plt.hlines(lag_index.mean(), index_v.min(), index_v.max(), linestyle="--")

        plt.plot(index_v, a + b * index_v, "r")
        plt.xlabel(self.Field_col)
        plt.ylabel("Spatial Lag of " + self.Field_col)
        plt.title("Local Getis-Ord G Scatterplot")

        return knext.Table.from_pandas(gdf), knext.view_matplotlib(f)
