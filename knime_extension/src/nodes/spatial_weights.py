import logging
import knime_extension as knext
import util.knime_utils as knut
import pandas as pd
import geopandas as gp

import numpy as np
import libpysal
import scipy.sparse
from libpysal.weights import WSP


LOGGER = logging.getLogger(__name__)

category = knext.category(
    path="/geo",
    level_id="spatialstatistic",
    name="Spatial Statistic",
    description="Spatial Statistic Nodes",
    icon="icons/icon.png",
)


@knext.node(
    name="Spatial Weights",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon.png",
    category=category,
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
        knut.column_exists(self.geo_col, input_schema_1)
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
