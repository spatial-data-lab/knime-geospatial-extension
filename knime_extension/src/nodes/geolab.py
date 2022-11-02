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
from pandarallel import pandarallel
from pyDataverse.api import NativeApi, DataAccessApi
from pyDataverse.models import Dataverse
import io
import numpy as np
from shapely.geometry import Polygon

__category = knext.category(
    path="/geo",
    level_id="geolab",
    name="Spatial Data Lab",
    description="Nodes that for testing and future exploration.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/GeolabCategroy.png",
    after="opendataset",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Geolab/"

############################################
# spatial GM_Error node
############################################


@knext.node(
    name="Spatial GM_Error",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMErr.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Error model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Error model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Error model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Error model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Error model.",
)
class SpatialGM_Error:
    """
    Spatial GM_Error
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Error model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Error model.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        model = spreg.GM_Error(y, X, w)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        # results =  results.dropna()

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Error_Het node
############################################


@knext.node(
    name="Spatial GM_Error_Het",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMErrHet.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Error_Het model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Error_Het model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Error_Het model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Error_Het model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Error_Het model.",
)
class SpatialGM_Error_Het:
    """
    Spatial GM_Error_Het
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Error_Het model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Error_Het model.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        model = spreg.GM_Error_Het(y, X, w)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        # results =  results.dropna()

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Error_Hom node
############################################


@knext.node(
    name="Spatial GM_Error_Hom",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMErrHom.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Error_Hom model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Error_Hom model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Error_Hom model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Error_Hom model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Error_Hom model.",
)
class SpatialGM_Error_Hom:
    """
    Spatial GM_Error_Hom
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Error_Hom model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Error_Hom model.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        model = spreg.GM_Error_Hom(y, X, w)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        # results =  results.dropna()

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Combo node
############################################


@knext.node(
    name="Spatial GM_Combo",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMcombo.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Combo model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Combo model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Combo model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Combo model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Combo model.",
)
class SpatialGM_Combo:
    """
    Spatial GM_Combo
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Combo model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Combo model.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        model = spreg.GM_Combo(y=y, x=X, w=w)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Spatial Pseudo R-squared ": model.pr2_e,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Combo_Het node
############################################


@knext.node(
    name="Spatial GM_Combo_Het",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMcomboHet.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Combo_Het model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Combo_Het model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Combo_Het model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Combo_Het model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Combo_Het model.",
)
class SpatialGM_Combo_Het:
    """
    Spatial GM_Combo_Het
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Combo_Het model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Combo_Het model.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        model = spreg.GM_Combo_Het(y=y, x=X, w=w)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Spatial Pseudo R-squared ": model.pr2_e,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Combo_Hom node
############################################


@knext.node(
    name="Spatial GM_Combo_Hom",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMcomboHom.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Combo_Hom model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Combo_Hom model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Combo_Hom model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Combo_Hom model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Combo_Hom model.",
)
class SpatialGM_Combo_Hom:
    """
    Spatial GM_Combo_Hom
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Combo_Hom model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Combo_Hom model.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        model = spreg.GM_Combo_Hom(y=y, x=X, w=w)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Endog_Error node
############################################
# FIXME: add another two parameters
@knext.node(
    name="Spatial GM_Endog_Error",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMendogErr.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Endog_Error model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Endog_Error model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Endog_Error model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Endog_Error model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Endog_Error model.",
)
class SpatialGM_Endog_Error:
    """
    Spatial GM_Endog_Error
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Endog_Error model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Endog_Error model.",
        column_filter=knut.is_numeric,
    )

    yend = knext.ColumnParameter(
        "Endogenous variable",
        "The column containing the endogenous variable to use for the calculation of the spatial GM_Endog_Error model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    q = knext.ColumnParameter(
        "External exogenous variable",
        "The column containing the external exogenous variable to use for the calculation of the spatial GM_Endog_Error model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values
        yend = gdf[self.yend].values
        q = gdf[self.q].values
        model = spreg.GM_Endog_Error(y=y, x=X, w=w, yend=yend, q=q)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Endog_Error_Het node
############################################
# FIXME: add another two parameters
@knext.node(
    name="Spatial GM_Endog_Error_Het",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMendogErrHet.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Endog_Error_Het model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Endog_Error_Het model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Endog_Error_Het model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Endog_Error_Het model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Endog_Error_Het model.",
)
class SpatialGM_Endog_Error_Het:
    """
    Spatial GM_Endog_Error_Het
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Endog_Error_Het model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Endog_Error_Het model.",
        column_filter=knut.is_numeric,
    )

    yend = knext.ColumnParameter(
        "Endogenous variable",
        "The column containing the endogenous variable to use for the calculation of the spatial GM_Endog_Error_Het model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    q = knext.ColumnParameter(
        "External exogenous variable",
        "The column containing the external exogenous variable to use for the calculation of the spatial GM_Endog_Error_Het model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values
        yend = gdf[self.yend].values
        q = gdf[self.q].values

        model = spreg.GM_Endog_Error_Het(y=y, x=X, w=w, yend=yend, q=q)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
        )


############################################
# spatial GM_Endog_Error_Hom node
############################################
# FIXME: add another two parameters
@knext.node(
    name="Spatial GM_Endog_Error_Hom",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMendogErrHom.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM_Endog_Error_Hom model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM_Endog_Error_Hom model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM_Endog_Error_Hom model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM_Endog_Error_Hom model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM_Endog_Error_Hom model.",
)
class SpatialGM_Endog_Error_Hom:
    """
    Spatial GM_Endog_Error_Hom
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the input table.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial GM_Endog_Error_Hom model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial GM_Endog_Error_Hom model.",
        column_filter=knut.is_numeric,
    )

    yend = knext.ColumnParameter(
        "Endogenous variable",
        "The column containing the endogenous variable to use for the calculation of the spatial GM_Endog_Error_Hom model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    q = knext.ColumnParameter(
        "External exogenous variable",
        "The column containing the external exogenous variable to use for the calculation of the spatial GM_Endog_Error_Hom model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()
        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values
        yend = gdf[self.yend].values
        q = gdf[self.q].values

        model = spreg.GM_Endog_Error_Hom(y=y, x=X, w=w, yend=yend, q=q)

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results = results.dropna()
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        # #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Number of Observations": model.n,
                "Number of Variables": model.k,
            },
            index=[0],
        )
        result2 = result2.round(7)

        html = """<p><pre>%s</pre>""" % model.summary.replace("\n", "<br/>")

        return (
            knext.Table.from_pandas(result2),
            knext.Table.from_pandas(results),
            knext.view_html(html),
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


############################################
# Harvard DataVerse Search Node
############################################


@knext.node(
    name="Harvard DataVerse Search",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "DvSearch.png",
    category=__category,
)
@knext.output_table(
    name="Output Table",
    description="Output table of Harvard DataVerse Search",
)
class HarvardDataVerseSearch:
    """
    Harvard DataVerse Search
    """

    # input parameters
    search_term = knext.StringParameter(
        "Search Term",
        "The search term or terms. Using “title:data” will search only the “title” field. “*” can be used as a wildcard either alone or adjacent to a term (i.e. “bird*”). ",
        default_value="mobility",
    )

    search_type = knext.StringParameter(
        "Search Type",
        "Can be either “Dataverse”, “dataset”, or “file”. ",
        enum=["dataset", "file", "dataverse", "all"],
        default_value="dataset",
    )

    def configure(self, configure_context):

        return None

    def execute(self, exec_context: knext.ExecutionContext):
        # get the search term
        pandarallel.initialize(progress_bar=False, nb_workers=15)

        search_term = self.search_term

        start = 0
        type_ = self.search_type
        per_page = 1000
        url = (
            "https://dataverse.harvard.edu/api/search?q=%s&start=%d&type=%s&per_page=%d"
            % (search_term, start, type_, per_page)
        )
        r = requests.get(url)
        data = r.json()
        pages = data["data"]["total_count"] // per_page + 1

        # Get the data all pages from the API and save it to a dataframe

        urls = []
        for start in range(pages):
            temp = {}
            url = (
                "https://dataverse.harvard.edu/api/search?q=%s&start=%d&type=%s&per_page=%d"
                % (search_term, start, type_, per_page)
            )
            temp["url"] = url
            temp["query"] = search_term
            urls.append(temp)
        urls = pd.DataFrame(urls)

        def get_data(url):
            import requests
            import pandas as pd

            r = requests.get(url)
            data = r.json()
            return pd.DataFrame(data["data"]["items"])

        urls["data"] = urls["url"].parallel_apply(get_data)
        df = pd.concat(urls["data"].values, ignore_index=True)

        def float_list(x):
            try:
                f = ";".join([str(i) for i in x])
                return f
            except:
                return x

        df[
            [
                "subjects",
                "contacts",
                "authors",
                "keywords",
                "producers",
                "relatedMaterial",
                "geographicCoverage",
                "dataSources",
            ]
        ] = df[
            [
                "subjects",
                "contacts",
                "authors",
                "keywords",
                "producers",
                "relatedMaterial",
                "geographicCoverage",
                "dataSources",
            ]
        ].applymap(
            float_list
        )
        df.drop("publications", axis=1, inplace=True)

        # return the results as a table
        return knext.Table.from_pandas(df)


############################################
# Harvard DataVerse Query Data Files Source Node
############################################


@knext.node(
    name="Harvard DataVerse GlobalDOI Search",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "DvDOIsource.png",
    category=__category,
)
@knext.output_table(
    name="Output Table",
    description="Output table of Harvard DataVerse Query Data Files",
)
class HarvardDataVerseQueryDataFilesSource:
    """
    Harvard DataVerse Query Data Files
    """

    # input parameters
    global_doi = knext.StringParameter(
        "Global DOI",
        "The global DOI of the dataset. ",
        default_value="doi:10.7910/DVN/ZAKKCE",
    )

    def configure(self, configure_context):

        return None

    def execute(self, exec_context: knext.ExecutionContext):

        base_url = "https://dataverse.harvard.edu/"
        global_doi = self.global_doi
        api = NativeApi(base_url)
        data_api = DataAccessApi(base_url)
        dataset = api.get_dataset(global_doi)
        files_list = dataset.json()["data"]["latestVersion"]["files"]
        df = pd.json_normalize(files_list)
        df = df.fillna(method="bfill")

        return knext.Table.from_pandas(df)


############################################
# Harvard DataVerse Query Data Files Node
############################################


@knext.node(
    name="Harvard DataVerse GlobalDOI Link ",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "DvGlobalDOIlink.png",
    category=__category,
)
@knext.input_table(
    name="Input Table",
    description="Input table of Harvard DataVerse Query Data Files",
)
@knext.output_table(
    name="Output Table",
    description="Output table of Harvard DataVerse Query Data Files",
)
class HarvardDataVerseQueryDataFiles:
    """
    Harvard DataVerse Query Data Files
    """

    # input parameters
    global_doi_column = knext.ColumnParameter(
        "Global DOI Column",
        "The column containing the global DOI of the dataset. ",
    )

    def configure(self, configure_context, input_schema):

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        base_url = "https://dataverse.harvard.edu/"
        global_doi = input_table.to_pandas()[self.global_doi_column].values[0]
        api = NativeApi(base_url)
        data_api = DataAccessApi(base_url)
        dataset = api.get_dataset(global_doi)
        files_list = dataset.json()["data"]["latestVersion"]["files"]
        df = pd.json_normalize(files_list)
        df = df.fillna(method="bfill")

        return knext.Table.from_pandas(df)


############################################
# Harvard DataVerse Read Data File Node
############################################


@knext.node(
    name="Harvard DataVerse DataID Reader",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "DvFileReader.png",
    category=__category,
)
@knext.input_table(
    name="Input Table",
    description="Input table of Harvard DataVerse Read Data File",
)
@knext.output_table(
    name="Output Table",
    description="Output table of Harvard DataVerse Read Data File",
)
class HarvardDataVerseReadDataFile:
    """
    Harvard DataVerse Read Data File
    """

    # input parameters
    dataFile_id_column = knext.ColumnParameter(
        "DataFile ID Column",
        "The column containing the DataFile ID of the dataset. ",
    )

    is_geo = knext.BoolParameter(
        "Is Geo",
        "Is the file a geo file?",
        default_value=False,
    )

    def configure(self, configure_context, input_schema):

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        base_url = "https://dataverse.harvard.edu/"
        api = NativeApi(base_url)
        data_api = DataAccessApi(base_url)

        file_id = input_table.to_pandas()[self.dataFile_id_column].values[0]
        response = data_api.get_datafile(file_id)
        content_ = io.BytesIO(response.content)

        if self.is_geo:
            df = gp.read_file(content_)
            df.reset_index(inplace=True, drop=True)

        else:
            df = pd.read_csv(content_, encoding="utf8", sep="\t")
        return knext.Table.from_pandas(df)


############################################
# Harvard DataVerse Replace Data File Node
############################################


@knext.node(
    name="Harvard DataVerse Data File Replacer",
    node_type=knext.NodeType.SINK,
    icon_path=__NODE_ICON_PATH + "DvReplace.png",
    category=__category,
)
@knext.input_table(
    name="Input Table",
    description="Input table of Harvard DataVerse Replace Data File",
)
class HarvardDataVerseReplaceDataFile:
    """
    Harvard DataVerse Replace Data File
    """

    # input parameters
    dataFile_id_column = knext.ColumnParameter(
        "DataFile ID Column",
        "The column containing the DataFile ID of the dataset. ",
    )

    dataFile_name_column = knext.ColumnParameter(
        "DataFile Name Column",
        "The column containing the DataFile Name of the dataset. ",
    )

    upload_file_path = knext.StringParameter(
        "Upload File Path", "The path to the file to be uploaded. "
    )

    API_TOKEN = knext.StringParameter(
        "API Token", "The API Token for the Harvard DataVerse. "
    )

    def configure(self, configure_context, input_schema):

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        base_url = "https://dataverse.harvard.edu/"
        api = NativeApi(base_url, self.API_TOKEN)
        json_str = """{"description":"My description.","categories":["Data"],"forceReplace":true}"""

        def check_need_replace(x, file_list):
            a = False
            if x.endswith(".tab"):
                if x.split(".")[0] in file_list:
                    a = True
            return a

        df = input_table.to_pandas()
        file_list = [
            file_name.split(".")[0] for file_name in os.listdir(self.upload_file_path)
        ]
        need_upload = df[
            df[self.dataFile_name_column].map(
                lambda x: check_need_replace(x, file_list)
            )
        ]

        for k, v in need_upload.iterrows():
            remote_file_id = v["dataFile.id"]
            local_file_path = input_path + v["dataFile.filename"].replace(
                ".tab", ".csv"
            )
            for i in range(60):
                test = api.replace_datafile(
                    remote_file_id, local_file_path, json_str, False
                )
                if json.loads(test.content)["status"] != "ERROR":
                    break
                if (
                    json.loads(test.content)["message"]
                    == "Error! You may not replace a file with a file that has duplicate content."
                ):
                    break
                # print(json.loads(test.content))
                time.sleep(1)
                print("retry %d" % (i + 1))

        return None


############################################
# Harvard DataVerse Publish Node
############################################


@knext.node(
    name="Harvard DataVerse Publish",
    node_type=knext.NodeType.SINK,
    icon_path=__NODE_ICON_PATH + "DvPublish.png",
    category=__category,
)
class HarvardDataVersePublish:
    """
    Harvard DataVerse Publish
    """

    # input parameters
    dataset_doi = knext.StringParameter(
        "Dataset DOI",
        "The DOI of the dataset to be published. ",
    )

    API_TOKEN = knext.StringParameter(
        "API Token", "The API Token for the Harvard DataVerse. "
    )

    def configure(self, configure_context):

        return None

    def execute(self, exec_context: knext.ExecutionContext):

        base_url = "https://dataverse.harvard.edu/"
        api = NativeApi(base_url, self.API_TOKEN)
        api.publish_dataset(pid=self.dataset_doi, release_type="major")
        return None


############################################
# Create Grid Node
############################################


@knext.node(
    name="Create Grid",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "CreateGrid.png",
    category=__category,
)
@knext.input_table(
    name="Input Table",
    description="Input table of Create Grid",
)
@knext.output_table(
    name="Output Table",
    description="Output table of Create Grid",
)
class CreateGrid:
    """
    Create Grid
    """

    geo_col = knext.ColumnParameter(
        "Geometry Column",
        "The column containing the geometry of the dataset. ",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    grid_length = knext.IntParameter(
        "Grid Length",
        "The length in meters of the grid. ",
        default_value=100,
    )

    def configure(self, configure_context, input_schema):

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        xmin, ymin, xmax, ymax = gdf.total_bounds
        width = self.grid_length
        height = self.grid_length
        rows = int(np.ceil((ymax - ymin) / height))
        cols = int(np.ceil((xmax - xmin) / width))
        XleftOrigin = xmin
        XrightOrigin = xmin + width
        YtopOrigin = ymax
        YbottomOrigin = ymax - height
        polygons = []
        for i in range(cols):
            Ytop = YtopOrigin
            Ybottom = YbottomOrigin
            for j in range(rows):
                polygons.append(
                    Polygon(
                        [
                            (XleftOrigin, Ytop),
                            (XrightOrigin, Ytop),
                            (XrightOrigin, Ybottom),
                            (XleftOrigin, Ybottom),
                        ]
                    )
                )
                Ytop = Ytop - height
                Ybottom = Ybottom - height
            XleftOrigin = XleftOrigin + width
            XrightOrigin = XrightOrigin + width

        grid = gp.GeoDataFrame({"geometry": polygons}, crs=gdf.crs)

        return knext.Table.from_pandas(grid)
