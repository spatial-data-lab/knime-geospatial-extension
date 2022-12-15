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
import  requests # for OSRM
import json # for OSRM

__category = knext.category(
    path="/community/geo",
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
# OSRM
############################################
@knext.node(
    name="OSRM Drive Matrix",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "BivariateLocal.png",
)
@knext.input_table(
    name="Input Table",
    description="Input table with startx,y,endx,y",
)

@knext.output_table(
    name="Output Table",
    description="Output table with travel Cost",
)


class OSRMDriveMatrix:
    """
    OSRM Distance Matrix
    """

    # input parameters
    StartX = knext.ColumnParameter(
        "Start X-Longitude",
        "The column containing the value for longitude of start point.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    StartY = knext.ColumnParameter(
        "Start Y-Latitude",
        "The column containing the value for latitude of start point.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    EndX = knext.ColumnParameter(
        "End X-Longitude",
        "The column containing the value for longitude of end point.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    EndY = knext.ColumnParameter(
        "End Y-Latitude",
        "The column containing the value for latitude  of end point.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1):
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        df = gp.GeoDataFrame(input_1.to_pandas())
        df = df.reset_index(drop=True)
        # set minimun set
        slnum=100
        nlength=df.shape[0]
        nloop=nlength//slnum
        ntail=nlength%slnum
        df['duration']=0.0
        df['distance']=0.0
        osrm_route_service="http://router.project-osrm.org/route/v1/driving/"
        if nloop>=1:
            for i in range(nloop):
                ns=slnum*i
                ne=ns+slnum-1
                dfs=df.copy().loc[ns:ne]
                dfs=dfs.rename(columns={self.StartX: "StartX", self.StartY: "StartY",self.EndX:"EndX", self.EndY: "EndY"})
                dfs=dfs[['StartX','StartY','EndX','EndY']]
                dfs=dfs.astype(str)
                dfs["period"] = dfs["StartX"] +","+ dfs["StartY"]+";"+dfs["EndX"] +","+ dfs["EndY"]
                Querylist = [';'.join(dfs['period'])]
                address = osrm_route_service + Querylist[0]
                try:
                    r = requests.get(address, params={'continue_straight':'false'}, timeout=None)
                    data = json.loads(r.text)
                    if data['code']=='Ok':                        
                        dfr = pd.DataFrame(data['routes'][0]['legs'])[['duration','distance']].iloc[::2]
                        df.loc[ns:ne,"duration"]=dfr.duration.to_list()
                        df.loc[ns:ne,"distance"]=dfr.distance.to_list()
                    else:
                        print("error from:{} to :{}".format(ns,ne))
                except:
                    print("error from:{} to :{}".format(ns,ne)) 
        if ntail>0:
            ns=slnum*nloop
            ne=ns+ntail-1
            dfs=df.copy().loc[ns:ne]
            dfs=dfs.rename(columns={self.StartX: "StartX", self.StartY: "StartY",self.EndX:"EndX", self.EndY: "EndY"})
            dfs=dfs[['StartX','StartY','EndX','EndY']]
            dfs=dfs.astype(str)
            dfs["period"] = dfs["StartX"] +","+ dfs["StartY"]+";"+dfs["EndX"] +","+ dfs["EndY"]
            Querylist = [';'.join(dfs['period'])]
            address = osrm_route_service + Querylist[0]
            try:
                r = requests.get(address, params={'continue_straight':'false'}, timeout=None)
                data = json.loads(r.text)
                if data['code']=='Ok':                    
                    dfr = pd.DataFrame(data['routes'][0]['legs'])[['duration','distance']].iloc[::2]
                    df.loc[ns:ne,"duration"]=dfr.duration.to_list()
                    df.loc[ns:ne,"distance"]=dfr.distance.to_list()
                else:
                    print("error from:{} to :{}".format(ns,ne))
            except:
                print("error from:{} to :{}".format(ns,ne)) 
            else:
                print("error from:{} to :{}".format(ns,ne))
      
        return knext.Table.from_pandas(df)