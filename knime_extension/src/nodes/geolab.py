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


# The following two models are rarely used and complex, so I comment them out.

# ###########################################
# spatial GM_Endog_Error_Het node
# ###########################################
# @knext.node(
#     name="Spatial GM Endog Error Het",
#     node_type=knext.NodeType.LEARNER,
#     category=__category,
#     icon_path=__NODE_ICON_PATH + "GMendogErrHet.png",
# )
# @knext.input_table(
#     name="Input Table",
#     description="Input Table with dependent and independent variables for calculation of the spatial GM Endog Error Het model.",
# )
# @knext.input_table(
#     name="Spatial Weights",
#     description="Input Table with spatial weights for calculation of the spatial GM Endog Error Het model.",
# )
# @knext.output_table(
#     name="Output Table",
#     description="Description of the spatial GM Endog Error Het model.",
# )
# @knext.output_table(
#     name="Variable and Coefficient Table",
#     description="Variable and Coefficient Table of the spatial GM Endog Error Het model.",
# )
# @knext.output_view(
#     name="Model summary view",
#     description="Model summary view of the spatial GM Endog Error Het model.",
# )
# class SpatialGM_Endog_Error_Het:
#     (
#         """
#     Spatial GM Endog Error Het Model.
#     GMM method for a spatial error model with heteroskedasticity and endogenous variables, with results and diagnostics; based on [ADKP10], following [Ans11]. More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Endog_Error_Het.html#spreg.GM_Endog_Error_Het). Please refer the following papers for more details.

#     - Irani Arraiz, David M. Drukker, Harry H. Kelejian, and Ingmar R. Prucha. A spatial Cliff-Ord-type model with heteroskedastic innovations: Small and large sample results. Journal of Regional Science, 50(2):592–614, 2010. doi:10.1111/j.1467-9787.2009.00618.x.
#     - Luc Anselin. GMM estimation of spatial error autocorrelation with and without heteroskedasticity. Technical Report, GeoDa Center for Geospatial Analysis and Computation, 2011.
#     """
#         + SHOULD_NOT_CONTAIN_MISSING_VALUES
#     )

#     geo_col = knut.geo_col_parameter()

#     id_col = mut.get_id_col_parameter()

#     dependent_variable = knext.ColumnParameter(
#         "Dependent variable",
#         "The column containing the dependent variable to use for the calculation of the spatial GM Endog Error Het model.",
#         column_filter=knut.is_numeric,
#         include_none_column=False,
#     )

#     independent_variables = knext.MultiColumnParameter(
#         "Independent variables",
#         "The columns containing the independent variables to use for the calculation of the spatial GM Endog Error Het model.",
#         column_filter=knut.is_numeric,
#     )

#     yend = knext.ColumnParameter(
#         "Endogenous variable",
#         "The column containing the endogenous variable to use for the calculation of the spatial GM Endog Error Het model.",
#         column_filter=knut.is_numeric,
#         include_none_column=False,
#     )

#     q = knext.ColumnParameter(
#         "External exogenous variable",
#         "The column containing the external exogenous variable to use for the calculation of the spatial GM Endog Error Het model.",
#         column_filter=knut.is_numeric,
#         include_none_column=False,
#     )

#     def configure(self, configure_context, input_schema, input_schema_2):
#         self.geo_col = knut.column_exists_or_preset(
#             configure_context, self.geo_col, input_schema, knut.is_geo
#         )

#         return None

#     def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
#         gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
#         adjust_list = input_2.to_pandas()

#         if "none" not in str(self.id_col).lower():
#             gdf.index = range(len(gdf))
#             id_map = dict(zip(gdf[self.id_col], gdf.index))
#             adjust_list["focal"] = adjust_list["focal"].map(id_map)
#             adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

#         import pandas as pd
#         import spreg
#         from libpysal.weights import W

#         w = W.from_adjlist(adjust_list)
#         # Prepare Georgia dataset inputs
#         X = gdf[self.independent_variables].values
#         y = gdf[self.dependent_variable].values
#         yend = gdf[self.yend].values
#         q = gdf[self.q].values

#         model = spreg.GM_Endog_Error_Het(y=y, x=X, w=w, yend=yend, q=q)

#         results = pd.DataFrame(
#             [model.name_x, model.betas, model.std_err, model.z_stat]
#         ).T
#         results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
#         results = results.dropna()
#         results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
#             lambda x: x[0]
#         )
#         results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
#             lambda x: x[1]
#         )
#         results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
#             lambda x: x[0]
#         )
#         # #
#         results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
#             lambda x: round(x, 7)
#         )
#         results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
#             lambda x: round(x, 7)
#         )
#         results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
#             lambda x: round(x, 7)
#         )
#         results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
#             lambda x: round(x, 7)
#         )

#         result2 = pd.DataFrame(
#             {
#                 "Pseudo R-squared ": model.pr2,
#                 "Number of Observations": model.n,
#                 "Number of Variables": model.k,
#             },
#             index=[0],
#         )
#         result2 = result2.round(7)

#         html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

#         return (
#             knext.Table.from_pandas(result2),
#             knext.Table.from_pandas(results),
#             knext.view_html(html),
#         )


# ############################################
# # spatial GM_Endog_Error_Hom node
# ############################################
# # FIXME: add another two parameters
# @knext.node(
#     name="Spatial GM Endog Error_Hom",
#     node_type=knext.NodeType.LEARNER,
#     category=__category,
#     icon_path=__NODE_ICON_PATH + "GMendogErrHom.png",
# )
# @knext.input_table(
#     name="Input Table",
#     description="Input Table with dependent and independent variables for calculation of the spatial GM Endog Error_Hom model.",
# )
# @knext.input_table(
#     name="Spatial Weights",
#     description="Input Table with spatial weights for calculation of the spatial GM Endog Error_Hom model.",
# )
# @knext.output_table(
#     name="Output Table",
#     description="Description of the spatial GM Endog Error_Hom model.",
# )
# @knext.output_table(
#     name="Variable and Coefficient Table",
#     description="Variable and Coefficient Table of the spatial GM Endog Error_Hom model.",
# )
# @knext.output_view(
#     name="Model summary view",
#     description="Model summary view of the spatial GM Endog Error_Hom model.",
# )
# class SpatialGM_Endog_Error_Hom:
#     (
#         """
#     Spatial GM Endog Error Hom Model.
#     GMM method for a spatial error model with homoskedasticity and endogenous variables, with results and diagnostics; based on Drukker et al. (2013) [DEP13], following Anselin (2011) [Ans11]. More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Endog_Error_Hom.html#spreg.GM_Endog_Error_Hom). Please refer the following papers for more details.

#     - %s
#     - %s
#     """
#         % (mrs.model_references["DEP13"], mrs.model_references["Ans11"])
#         + SHOULD_NOT_CONTAIN_MISSING_VALUES
#     )

#     geo_col = knext.ColumnParameter(
#         "Geometry Column",
#         "The column containing the geometry of the input table.",
#         column_filter=knut.is_geo,
#         include_row_key=False,
#         include_none_column=False,
#     )

#     dependent_variable = knext.ColumnParameter(
#         "Dependent variable",
#         "The column containing the dependent variable to use for the calculation of the spatial GM Endog Error Hom model.",
#         column_filter=knut.is_numeric,
#         include_none_column=False,
#     )

#     independent_variables = knext.MultiColumnParameter(
#         "Independent variables",
#         "The columns containing the independent variables to use for the calculation of the spatial GM Endog Error Hom model.",
#         column_filter=knut.is_numeric,
#     )

#     yend = knext.ColumnParameter(
#         "Endogenous variable",
#         "The column containing the endogenous variable to use for the calculation of the spatial GM Endog Error Hom model.",
#         column_filter=knut.is_numeric,
#         include_none_column=False,
#     )

#     q = knext.ColumnParameter(
#         "External exogenous variable",
#         "The column containing the external exogenous variable to use for the calculation of the spatial GM Endog Error Hom model.",
#         column_filter=knut.is_numeric,
#         include_none_column=False,
#     )

#     def configure(self, configure_context, input_schema, input_schema_2):
#         self.geo_col = knut.column_exists_or_preset(
#             configure_context, self.geo_col, input_schema, knut.is_geo
#         )

#         return None

#     def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
#         gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
#         adjust_list = input_2.to_pandas()

#         import pandas as pd
#         import spreg
#         from libpysal.weights import W

#         w = W.from_adjlist(adjust_list)

#         if "none" not in str(self.id_col).lower():
#             gdf.index = range(len(gdf))
#             id_map = dict(zip(gdf[self.id_col], gdf.index))
#             adjust_list["focal"] = adjust_list["focal"].map(id_map)
#             adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

#         # Prepare Georgia dataset inputs
#         X = gdf[self.independent_variables].values
#         y = gdf[self.dependent_variable].values
#         yend = gdf[self.yend].values
#         q = gdf[self.q].values

#         model = spreg.GM_Endog_Error_Hom(y=y, x=X, w=w, yend=yend, q=q)

#         results = pd.DataFrame(
#             [model.name_x, model.betas, model.std_err, model.z_stat]
#         ).T
#         results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
#         results = results.dropna()
#         results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
#             lambda x: x[0]
#         )
#         results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
#             lambda x: x[1]
#         )
#         results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
#             lambda x: x[0]
#         )
#         # #
#         results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
#             lambda x: round(x, 7)
#         )
#         results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
#             lambda x: round(x, 7)
#         )
#         results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
#             lambda x: round(x, 7)
#         )
#         results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
#             lambda x: round(x, 7)
#         )

#         result2 = pd.DataFrame(
#             {
#                 "Pseudo R-squared ": model.pr2,
#                 "Number of Observations": model.n,
#                 "Number of Variables": model.k,
#             },
#             index=[0],
#         )
#         result2 = result2.round(7)

#         html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

#         return (
#             knext.Table.from_pandas(result2),
#             knext.Table.from_pandas(results),
#             knext.view_html(html),
#         )

