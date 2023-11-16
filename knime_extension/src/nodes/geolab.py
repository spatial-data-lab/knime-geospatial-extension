# lingbo
from typing import Callable
from wsgiref.util import shift_path_info
from jmespath import search
import pandas as pd
import geopandas as gp
import knime_extension as knext
from sympy import content
import util.knime_utils as knut
import requests
import io
import numpy as np
from shapely.geometry import Polygon
import requests  # for OSRM
import json  # for OSRM
import urllib.request as urllib2  # for Google Drive

__category = knext.category(
    path="/community/geo",
    level_id="geolab",
    name="Spatial Data Lab",
    description="Nodes that for testing and future exploration.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/GeolabCategroy.png",
    after="opendataset",
)


############################################
# Bivariate Global Moran’s I node
############################################
# FIXME: add test workflow
@knext.node(
    name="Bivariate Global Moran’s I",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "BivariateGlobal.png",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Bivariate Global Moran’s I",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Bivariate Global Moran’s I",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Bivariate Global Moran’s I",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Bivariate Global Moran’s I",
#     id="pysal.esda.moran.Moran",
# )
@knext.output_view(
    name="output view",
    description="Output view of Bivariate Global Moran’s I",
)
class BivariateGlobalMoran:
    """
    Bivariate Global Moran’s I
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for Bivariate Global Moran’s I.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col1 = knext.ColumnParameter(
        "Field column 1",
        "The column containing the field to use for the calculation of Bivariate Global Moran’s I.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    Field_col2 = knext.ColumnParameter(
        "Field column 2",
        "The column containing the field to use for the calculation of Bivariate Global Moran’s I.",
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

        x = gdf[self.Field_col1]
        y = gdf[self.Field_col2]
        np.random.seed(12345)
        bi = esda.moran.Moran_BV(x, y, w)

        out = pd.DataFrame(
            {
                "Bivariate Moran’s I": [bi.I],
                "p-value": [bi.p_sim],
                "z-score": [bi.z_sim],
            }
        )
        out.reset_index(inplace=True)

        ax = sbn.kdeplot(bi.sim, shade=True)
        plt.vlines(bi.I, 0, 1, color="r")
        # plt.vlines(bi.EI, 0, 1)
        plt.xlabel("Bivariate Moran’s I")

        return knext.Table.from_pandas(out), knext.view_matplotlib(ax.get_figure())


############################################
# Bivariate Local Moran Statistics
############################################
# FIXME: add test workflow
@knext.node(
    name="Bivariate Local Moran Statistics",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "BivariateLocal.png",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Bivariate Local Moran Statistics",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Bivariate Local Moran Statistics",
)
@knext.output_table(
    name="Output Table",
    description="Output table results of Bivariate Local Moran Statistics",
)
# @knext.output_binary(
#     name="output model",
#     description="Output model of Bivariate Local Moran Statistics",
#     id="pysal.esda.moran.Moran_Local",
# )
@knext.output_view(
    name="output view",
    description="Output view of Bivariate Local Moran Statistics",
)
class BivariateLocalMoran:
    """
    Bivariate Local Moran Statistics
    """

    # input parameters
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "The column containing the geometry to use for Bivariate Local Moran Statistics.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    Field_col1 = knext.ColumnParameter(
        "Field column 1",
        "The column containing the field to use for the calculation of Bivariate Local Moran Statistics.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    Field_col2 = knext.ColumnParameter(
        "Field column 2",
        "The column containing the field to use for the calculation of Bivariate Local Moran Statistics.",
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

        x = gdf[self.Field_col1]
        y = gdf[self.Field_col2]
        np.random.seed(12345)
        bi = esda.moran.Moran_Local_BV(x, y, w)

        gdf.loc[:, "Bivariate Local Moran’s I"] = bi.Is
        gdf.loc[:, "p-value"] = bi.p_sim
        gdf.loc[:, "z-score"] = bi.z_sim
        gdf.loc[:, "spots"] = bi.q

        gdf.loc[:, "spots_type"] = gdf["spots"].replace(
            {1: "HH", 2: "LL", 3: "LH", 4: "HL"}
        )
        gdf.loc[gdf["p-value"] < 0.05, "spots_type"] = "Not Significant"

        lag_index = lps.weights.lag_spatial(w, gdf[self.Field_col1])
        index_v = gdf[self.Field_col1]
        b, a = np.polyfit(index_v, lag_index, 1)
        f, ax = plt.subplots(1, figsize=(9, 9))

        plt.plot(index_v, lag_index, ".", color="firebrick")
        plt.vlines(index_v.mean(), lag_index.min(), lag_index.max(), linestyle="--")
        plt.hlines(lag_index.mean(), index_v.min(), index_v.max(), linestyle="--")

        plt.plot(index_v, a + b * index_v, "r")
        plt.title("Moran Scatterplot")
        plt.ylabel("Spatial Lag of %s" % self.Field_col1)
        plt.xlabel("%s" % self.Field_col1)

        return knext.Table.from_pandas(gdf), knext.view_matplotlib(f)
