import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import util.modeling_utils as mut


__category = knext.category(
    path="/community/geo",
    level_id="spatialmodels",
    name="Spatial Modelling",
    description="Spatial Models Nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialModelCategory.png",
    after="spatialstatistic",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialModel/"


def get_id_col_parameter(
    label: str = "ID column",
    description: str = """Select the column which contains for each observation in the input data a unique ID, it should be an integer column.
    The IDs must match with the values of the 
    [Spatial Weights node](https://hub.knime.com/center%20for%20geographic%20analysis%20at%20harvard%20university/extensions/sdl.harvard.features.geospatial/latest/org.knime.python3.nodes.extension.ExtensionNodeSetFactory$DynamicExtensionNodeFactory:4d710eae/)
    ID column.
    If 'none' is selected, the IDs will be automatically generated from 0 to the number of rows flowing the order of 
    the first input table.
    """,
):
    """
    Returns the unique ID column. It should always keep the same as the ID column in the spatial weights matrix node.
    The selected column should contain unique IDs for each observation in the input data.
    """
    return knext.ColumnParameter(
        label=label,
        description=description,
        include_none_column=True,
        column_filter=knut.is_int,
        since_version="1.1.0",
    )


############################################
# spatial 2SlS node
############################################
@knext.node(
    name="2SLS with Spatial Test",
    node_type=knext.NodeType.LEARNER,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "2SLS.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for calculation of Spatial 2SlS",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights table for calculation of Spatial 2SlS",
)
@knext.output_table(
    name="Model Description Table",
    description="Description of Spatial 2SlS, including Pseudo R-squared, Spatial Pseudo R-squared, Number of Observations, and Number of Variables",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of Spatial 2SlS",
)
@knext.output_view(
    name="Model Summary View",
    description="Model Summary View of Spatial 2SlS",
)
class Spatial2SLSModel:
    """Spatial two stage least squares (S2SLS) with results and diagnostics.
    Spatial two stage least squares (S2SLS) with results and diagnostics. More details can be found in the following reference, Luc Anselin. Spatial Econometrics: Methods and Models. Kluwer. Dordrecht, 1988.
    """

    # input parameters
    geo_col = knut.geo_col_parameter()

    id_col = get_id_col_parameter()

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of Spatial 2SlS.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of Spatial 2SlS.",
        column_filter=knut.is_numeric,
    )

    Orders_of_W = knext.IntParameter(
        "Orders of W",
        """Orders of W to include as instruments for the spatially lagged dependent variable. For example, w_lags=1, 
        then instruments are WX; if w_lags=2, then WX, WWX; and so on.""",
        default_value=1,
    )

    Spatial_Diagnostics = knext.BoolParameter(
        "Spatial Diagnostics",
        "If selected, the node computes the Anselin-Kelejian test",
        default_value=False,
    )

    robust = knext.StringParameter(
        "Robust",
        """If ‘white’, then a White consistent estimator of the variance-covariance matrix is given. If ‘hac’, 
        then a HAC consistent estimator of the variance-covariance matrix is given. Set to None for default.""",
        enum=["white", "hac", "none"],
        default_value="none",
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)
        y = gdf[self.dependent_variable].values
        x = gdf[self.independent_variables].values

        kws = {
            "y": y,
            "x": x,
            "w": w,
            "w_lags": self.Orders_of_W,
            "spat_diag": self.Spatial_Diagnostics,
            "name_y": self.dependent_variable,
            "name_x": self.independent_variables,
            "name_w": "Spatial Weights",
            "name_ds": "Input Table",
        }
        if "none" not in str(self.robust).lower():
            kws["robust"] = self.robust

        import spreg

        model = spreg.GM_Lag(**kws)

        # model = spreg.GM_Lag(y, x, w=w,w_lags=self.Orders_of_W, robust= self.robust,
        # name_y=self.dependent_variable, name_x=self.independent_variables, name_w="Spatial Weights", name_ds="Input Table")

        import pandas as pd

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "z-Statistic"]
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "z-Statistic"] = results.loc[:, "z-Statistic"].map(
            lambda x: x[0]
        )
        #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "z-Statistic"] = results.loc[:, "z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        results = results.dropna()

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Spatial Pseudo R-squared": model.pr2_e,
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


############################################################################################################
# Spatial Lag Panel Model with Fixed Effects node
############################################################################################################


@knext.node(
    name="Spatial Lag Panel Model",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "SpatialLag.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input Table for Spatial Lag Panel Model with Fixed Effects",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights for Spatial Lag Panel Model with Fixed Effects",
)
@knext.output_table(
    name="Model Description Table",
    description="Model Description Table for Spatial Lag Panel Model with Fixed Effects",
)
@knext.output_table(
    name="Model Coefficients Table",
    description="Model Coefficients Table for Spatial Lag Panel Model with Fixed Effects",
)
@knext.output_view(
    name="Model Summary",
    description="Model Summary for Spatial Lag Panel Model with Fixed Effects",
)
class SpatialLagPanelModelwithFixedEffects:
    """Spatial Lag Panel Model with Fixed Effects.
    Spatial Lag Panel Model with Fixed Effects. ML estimation of the fixed effects spatial lag model with all results and diagnostics. More details can be found at J. Paul Elhorst. Specification and estimation of spatial panel data models. International Regional Science Review, 26(3):244–268, 2003. doi:10.1177/0160017603253791.
    """

    geo_col = knut.geo_col_parameter()

    id_col = get_id_col_parameter()

    dependent_variable = knext.MultiColumnParameter(
        "Dependent variables",
        "The column containing the dependent variables to use for the calculation of Spatial Lag Panel Model with Fixed Effects.",
        column_filter=knut.is_numeric,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of Spatial Lag Panel Model with Fixed Effects.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        x = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        kws = {
            "y": y,
            "x": x,
            "w": w,
            "name_y": self.dependent_variable,
            "name_x": self.independent_variables,
            "name_w": "Spatial Weights",
            "name_ds": "Input Table",
        }

        import spreg

        model = spreg.Panel_FE_Lag(**kws)

        import pandas as pd

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "z-Statistic"]
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "z-Statistic"] = results.loc[:, "z-Statistic"].map(
            lambda x: x[0]
        )
        #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "z-Statistic"] = results.loc[:, "z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        results = results.dropna()

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Spatial Pseudo R-squared": model.pr2_e,
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


############################################################################################################
# Spatial Error Panel Model with Fixed Effects node
############################################################################################################


@knext.node(
    name="Spatial Error Panel Model",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "SpatialError.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input Table for Spatial Error Panel Model with Fixed Effects",
)
@knext.input_table(
    name="Spatial Weights",
    description="Spatial Weights for Spatial Error Panel Model with Fixed Effects",
)
@knext.output_table(
    name="Model Description Table",
    description="Model Description Table for Spatial Error Panel Model with Fixed Effects",
)
@knext.output_table(
    name="Model Coefficients Table",
    description="Model Coefficients Table for Spatial Error Panel Model with Fixed Effects",
)
@knext.output_view(
    name="Model Summary",
    description="Model Summary for Spatial Error Panel Model with Fixed Effects",
)
class SpatialErrorPanelModelwithFixedEffects:
    """Spatial Error Panel Model with Fixed Effects node.
    Spatial Error Panel Model with Fixed Effects node. ML estimation of the fixed effects spatial error model with all results and diagnostics. More details can be found at J. Paul Elhorst. Specification and estimation of spatial panel data models. International Regional Science Review, 26(3):244–268, 2003. doi:10.1177/0160017603253791.
    """

    geo_col = knut.geo_col_parameter()

    id_col = get_id_col_parameter()

    dependent_variable = knext.MultiColumnParameter(
        "Dependent variables",
        "The column containing the dependent variables to use for the calculation of Spatial Error Panel Model with Fixed Effects.",
        column_filter=knut.is_numeric,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of Spatial Error Panel Model with Fixed Effects.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        adjust_list = input_2.to_pandas()

        if "none" not in str(self.id_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)

        x = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        kws = {
            "y": y,
            "x": x,
            "w": w,
            "name_y": self.dependent_variable,
            "name_x": self.independent_variables,
            "name_w": "Spatial Weights",
            "name_ds": "Input Table",
        }

        import spreg

        model = spreg.Panel_FE_Error(**kws)

        import pandas as pd

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "z-Statistic"]
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "z-Statistic"] = results.loc[:, "z-Statistic"].map(
            lambda x: x[0]
        )
        #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "z-Statistic"] = results.loc[:, "z-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        results = results.dropna()

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


############################################################################################################
# Geographically Weighted Regression (GWR) node
############################################################################################################


@knext.node(
    name="GWR Model",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GWR.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input Table for Geographically Weighted Regression",
)
# @knext.output_table(
#     name="Model Description Table",
#     description="Model Description Table for Geographically Weighted Regression",
# )
@knext.output_table(
    name="Model Coefficients Table",
    description="Model Coefficients Table for Geographically Weighted Regression",
)
@knext.output_binary(
    name="Model",
    description="Model for Geographically Weighted Regression",
    id="mgwr.gwr.GWR",
)
@knext.output_view(
    name="Model Summary",
    description="Model Summary for Geographically Weighted Regression",
)
class GeographicallyWeightedRegression:
    """Geographically Weighted Regression node.
    Performs Geographically Weighted Regression (GWR), a local form of linear regression used to model spatially varying relationships. Can currently estimate Gaussian, Poisson, and logistic models(built on a GLM framework).
    More details can be found at [here](https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-statistics-toolbox/geographically-weighted-regression.htm).
    """

    geo_col = knut.geo_col_parameter()

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of Geographically Weighted Regression.",
        column_filter=knut.is_numeric,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of Geographically Weighted Regression.",
        column_filter=knut.is_numeric,
    )

    search_method = knext.StringParameter(
        "Search method",
        "Bw search method: ‘golden’, ‘interval’",
        default_value="golden",
        enum=["golden", "interval"],
    )

    bandwith_min = knext.IntParameter(
        "Bandwith min",
        "Min value used in bandwidth search",
        default_value=2,
    )

    # bandwith_max = knext.IntParameter(
    #     "Bandwith Max",
    #     "max value used in bandwidth search",
    #     default=200,
    # )

    kernel = knext.StringParameter(
        "Kernel",
        "Type of kernel function used to weight observations; available options: ‘gaussian’ ‘bisquare’ ‘exponential’",
        default_value="bisquare",
        enum=["gaussian", "bisquare", "exponential"],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        # Prepare Georgia dataset inputs
        g_y = gdf[self.dependent_variable].values.reshape((-1, 1))
        g_X = gdf[self.independent_variables].values
        u = gdf.centroid.x
        v = gdf.centroid.y
        g_coords = list(zip(u, v))
        # g_X = (g_X - g_X.mean(axis=0)) / g_X.std(axis=0)
        g_y = g_y.reshape((-1, 1))
        # g_y = (g_y - g_y.mean(axis=0)) / g_y.std(axis=0)

        # Calibrate GWR model
        from mgwr.sel_bw import Sel_BW

        gwr_selector = Sel_BW(g_coords, g_y, g_X)
        gwr_bw = gwr_selector.search(bw_min=self.bandwith_min)
        # gwr_bw = gwr_selector.search(bw_min=self.bandwith_min,search_method=self.search_method)
        # print(gwr_bw)
        from mgwr.gwr import GWR

        gwr_model = GWR(g_coords, g_y, g_X, gwr_bw, fixed=False).fit()

        gdf.loc[:, "predy"] = gwr_model.predy
        gdf.loc[:, "resid"] = gwr_model.resid_response.reshape(-1, 1)
        gdf.loc[
            :,
            ["Intercept_beta"]
            + ["%s_beta" % item for item in self.independent_variables],
        ] = gwr_model.params
        gdf.loc[
            :, ["Intercept_t"] + ["%s_t" % item for item in self.independent_variables]
        ] = gwr_model.filter_tvals()
        intervals = gwr_model.get_bws_intervals(gwr_selector)

        import numpy as np

        intervals = np.asarray(intervals)
        if gwr_bw.shape == ():
            gdf.loc[:1, "bw"] = gwr_bw
            gdf.loc[:1, ["bw_lower", "bw_upper"]] = intervals
        else:
            gdf.loc[: (gwr_bw.shape[0]), "bw"] = gwr_bw
            gdf.loc[: (intervals.shape[0]), ["bw_lower", "bw_upper"]] = intervals
        # gdf.loc[:,"localR2"] = results.localR2
        # gdf.drop(columns=["<Row Key>"], inplace=True, axis=1)
        gdf.reset_index(drop=True, inplace=True)

        from io import StringIO
        import sys

        buffer = StringIO()
        sys.stdout = buffer
        gwr_model.summary()
        summary = buffer.getvalue()
        sys.stdout = sys.__stdout__

        html = """<p><pre>%s</pre>""" % summary.replace("\n", "<br/>")

        import pickle

        model_string = pickle.dumps(gwr_model)

        return knext.Table.from_pandas(gdf), model_string, knext.view_html(html)


#############################################################################################################
# Geographically Weighted Regression Predictor Node
#############################################################################################################
@knext.node(
    name="GWR Predictor",
    node_type=knext.NodeType.PREDICTOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GWRp.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input table for Geographically Weighted Regression Predictor",
)
@knext.input_binary(
    name="Model",
    description="Model for Geographically Weighted Regression Predictor",
    id="mgwr.gwr.GWR",
)
@knext.output_table(
    name="Output Table",
    description="Output table with predictions for Geographically Weighted Regression Predictor",
)
class GeographicallyWeightedRegressionPredictor:
    """Geographically Weighted Regression Predictor node.
    Geographically Weighted Regression Predictor. It will predict the dependent variable using the model and the input table.
    """

    geo_col = knut.geo_col_parameter()

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of Geographically Weighted Regression.",
        column_filter=knut.is_numeric,
    )

    def configure(self, configure_context, input_schema, input_binary_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, model):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        # Prepare Georgia dataset inputs
        g_X = gdf[self.independent_variables].values
        u = gdf["geometry"].x
        v = gdf["geometry"].y

        import numpy as np

        g_coords = np.array(list(zip(u, v)))
        # g_X = (g_X - g_X.mean(axis=0)) / g_X.std(axis=0)

        import pickle

        gwr_model = pickle.loads(model)

        gdf.loc[:, "predy"] = gwr_model.model.predict(g_coords, g_X).predictions
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


##############################################################################################################
# Multiscale Geographically Weighted Regression (MGWR)
##############################################################################################################


@knext.node(
    name="MGWR Model",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "MGWR.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="The input table containing the data to use for the calculation of Multiscale Geographically Weighted Regression.",
)
@knext.output_table(
    name="Model Coefficients Table",
    description="The output table containing the model coefficients for Multiscale Geographically Weighted Regression.",
)
# @knext.output_binary(
#     name="Model",
#     description="The output binary containing the model for Multiscale Geographically Weighted Regression.",
#     id = "mgwr.gwr.MGWR",
# )
@knext.output_view(
    name="Model Summary",
    description="Model Summary for Multiscale Geographically Weighted Regression",
)
class MultiscaleGeographicallyWeightedRegression:
    """Multiscale Geographically Weighted Regression node.
    Multiscale Geographically Weighted Regression estimation. 
    More details can be found at 
    1. A. Stewart Fotheringham, Wenbai Yang, and Wei Kang. Multiscale geographically weighted regression (mgwr). Annals of the American Association of Geographers, 107(6):1247–1265, 2017. URL: http://dx.doi.org/10.1080/24694452.2017.1352480, arXiv:http://dx.doi.org/10.1080/24694452.2017.1352480, doi:10.1080/24694452.2017.1352480. and Hanchen Yu, Alexander Stewart Fotheringham, Ziqi Li, Taylor Oshan, Wei Kang, and Levi John Wolf. Inference in multiscale geographically weighted regression. Geographical Analysis, 2019. URL: https://onlinelibrary.wiley.com/doi/abs/10.1111/gean.12189, arXiv:https://onlinelibrary.wiley.com/doi/pdf/10.1111/gean.12189, doi:10.1111/gean.12189.
    2. https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/how-multiscale-geographically-weighted-regression-mgwr-works.htm
    """

    geo_col = knut.geo_col_parameter()

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of Multiscale Geographically Weighted Regression.",
        column_filter=knut.is_numeric,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of Multiscale Geographically Weighted Regression.",
        column_filter=knut.is_numeric,
    )

    search_method = knext.StringParameter(
        "Search method",
        """Bw search method: ‘golden’, ‘interval’. Golden Search— Determines either the number of neighbors or distance band for each 
        explanatory variable using the Golden Search algorithm. This method searches multiple combinations 
        of values for each explanatory variable between a specified minimum and maximum value. Intervals— Determines the number of neighbors or distance band for each 
        explanatory variable by incrementing the number of neighbors or distance band from a minimum value.""",
        default_value="golden",
        enum=["golden", "interval"],
    )

    bandwith_min = knext.IntParameter(
        "Bandwith min",
        "Min value used in bandwidth search",
        default_value=2,
    )

    # bandwith_max = knext.IntParameter(
    #     "Bandwith Max",
    #     "max value used in bandwidth search",
    #     default=200,
    # )

    kernel = knext.StringParameter(
        "Kernel",
        "type of kernel function used to weight observations; available options: ‘gaussian’ ‘bisquare’ ‘exponential’",
        default_value="bisquare",
        enum=["gaussian", "bisquare", "exponential"],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        # Prepare Georgia dataset inputs
        g_y = gdf[self.dependent_variable].values.reshape((-1, 1))
        g_X = gdf[self.independent_variables].values
        u = gdf["geometry"].centroid.x
        v = gdf["geometry"].centroid.y
        g_coords = list(zip(u, v))
        # g_X = (g_X - g_X.mean(axis=0)) / g_X.std(axis=0)
        g_y = g_y.reshape((-1, 1))
        # g_y = (g_y - g_y.mean(axis=0)) / g_y.std(axis=0)

        # Calibrate MGWR model
        from mgwr.sel_bw import Sel_BW

        mgwr_selector = Sel_BW(g_coords, g_y, g_X, multi=True)
        mgwr_bw = mgwr_selector.search(multi_bw_min=[self.bandwith_min])
        from mgwr.gwr import MGWR

        mgwr_model = MGWR(
            g_coords,
            g_y,
            g_X,
            mgwr_selector,
            fixed=False,
            sigma2_v1=True,
            kernel="bisquare",
        ).fit()

        gdf.loc[:, "predy"] = mgwr_model.predy
        gdf.loc[:, "resid"] = mgwr_model.resid_response.reshape(-1, 1)
        gdf.loc[
            :,
            ["Intercept_beta"]
            + ["%s_beta" % item for item in self.independent_variables],
        ] = mgwr_model.params
        gdf.loc[
            :, ["Intercept_t"] + ["%s_t" % item for item in self.independent_variables]
        ] = mgwr_model.filter_tvals()
        intervals = mgwr_model.get_bws_intervals(mgwr_selector)

        import numpy as np

        intervals = np.asarray(intervals)
        if mgwr_bw.shape == ():
            gdf.loc[:1, "bw"] = mgwr_bw
            gdf.loc[:1, ["bw_lower", "bw_upper"]] = intervals
        else:
            gdf.loc[: (mgwr_bw.shape[0]), "bw"] = mgwr_bw.reshape(-1, 1)
            gdf.loc[: (intervals.shape[0]), ["bw_lower", "bw_upper"]] = intervals
        # gdf.loc[:,"localR2"] = results.localR2
        # gdf.drop(columns=["<Row Key>"], inplace=True, axis=1)
        gdf.reset_index(drop=True, inplace=True)

        from io import StringIO
        import sys

        buffer = StringIO()
        sys.stdout = buffer
        mgwr_model.summary()
        summary = buffer.getvalue()
        sys.stdout = sys.__stdout__

        html = """<p><pre>%s</pre>""" % summary.replace("\n", "<br/>")

        import pickle

        model_string = pickle.dumps(mgwr_model)

        return knext.Table.from_pandas(gdf), knext.view_html(html)
        # return knext.Table.from_pandas(gdf),model_string,knext.view_html(html)


##############################################################################################################
# Multiscale Geographically Weighted Regression Predictor Node
# Notice that the prediction of MGWR is not yet implemented in the mgwr package
##############################################################################################################

# @knext.node(
#     name="Multiscale Geographically Weighted Regression Predictor",
#     node_type=knext.NodeType.PREDICTOR,
#     category=__category,
#     icon_path=__NODE_ICON_PATH + "SpatialWeight.png",
#     after=""
# )
# @knext.input_table(
#     name="Input Table",
#     description="The input table containing the data to use for the prediction of Multiscale Geographically Weighted Regression.",
# )
# @knext.input_binary(
#     name="Model",
#     description="The model to use for the prediction of Multiscale Geographically Weighted Regression.",
#     id="mgwr.gwr.MGWR",
# )
# @knext.output_table(
#     name="Output Table",
#     description="The output table containing the prediction of Multiscale Geographically Weighted Regression.",
# )
# class MGWRPredictorNode:
#     """
#     Multiscale Geographically Weighted Regression Predictor Node
#     """

#     geo_col = knext.ColumnParameter(
#         "Geometry Column",
#         "The column containing the geometry of the input table.",
#         column_filter=knut.is_geo,
#         include_row_key=False,
#         include_none_column=False,
#     )

# # do we need this?
#     independent_variables = knext.MultiColumnParameter(
#         "Independent variables",
#         "The columns containing the independent variables to use for the prediction of Multiscale Geographically Weighted Regression.",
#         column_filter=knut.is_numeric
#     )

#     def configure(self, configure_context, input_schema,input_schema_2):
#         self.geo_col = knut.column_exists_or_preset(configure_context, self.geo_col, input_schema, knut.is_geo)
#         return None

#     def execute(self, exec_context:knext.ExecutionContext, input_1, model):

#             gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
#             #Prepare Georgia dataset inputs
#             g_X = gdf[self.independent_variables].values
#             u = gdf['geometry'].x
#             v = gdf['geometry'].y
#             g_coords = np.array(list(zip(u,v)))
#             # g_X = (g_X - g_X.mean(axis=0)) / g_X.std(axis=0)

#             #Calibrate MGWR model
#             mgwr_model = pickle.loads(model)
#             gdf.loc[:,'predy'] = mgwr_model.model.predict(g_coords, g_X).predictions
#             gdf.reset_index(drop=True, inplace=True)

#             return knext.Table.from_pandas(gdf)


############################################
# spatial OLS node
############################################


@knext.node(
    name="OLS with Spatial Test",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "SpatialOLS.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial OLS model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial OLS model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial OLS model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial OLS model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial OLS model.",
)
class SpatialOLS:
    """Spatial OLS node.
    Ordinary least squares with results and diagnostics. More information can be found at
    [here](https://spreg.readthedocs.io/en/latest/generated/spreg.OLS.html)
    """

    geo_col = knut.geo_col_parameter()

    id_col = get_id_col_parameter()

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial OLS model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial OLS model.",
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

        if "none" not in str(self.id_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        import spreg

        model = spreg.OLS(
            y,
            X,
            w=w,
            spat_diag=True,
            moran=True,
            white_test=True,
            name_y=self.dependent_variable,
            name_x=self.independent_variables,
        )

        import pandas as pd

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.t_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "t-Statistic"]
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "t-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "t-Statistic"] = results.loc[:, "t-Statistic"].map(
            lambda x: x[0]
        )
        #
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Std.Error"] = results.loc[:, "Std.Error"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "t-Statistic"] = results.loc[:, "t-Statistic"].map(
            lambda x: round(x, 7)
        )
        results.loc[:, "Probability"] = results.loc[:, "Probability"].map(
            lambda x: round(x, 7)
        )
        results = results.dropna()

        result2 = pd.DataFrame(
            {
                "R squared": model.r2,
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
# spatial ML_Lag node
############################################


@knext.node(
    name="Spatial Lag Model",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "MLLag.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial ML_Lag model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial ML_Lag model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial ML_Lag model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial ML_Lag model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial ML_Lag model.",
)
class SpatialML_Lag:
    """Spatial ML_Lag.
    ML estimation of the spatial lag model with all results and diagnostics. More details can be found at Luc Anselin. Spatial Econometrics: Methods and Models. Kluwer, Dordrecht, 1988.
    """

    geo_col = knut.geo_col_parameter()

    id_col = get_id_col_parameter()

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial ML_Lag model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial ML_Lag model.",
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

        if "none" not in str(self.id_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        import spreg

        model = spreg.ML_Lag(
            y,
            X,
            w=w,
            method="ord",
            name_x=self.independent_variables,
            name_y=self.dependent_variable,
        )

        import pandas as pd

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        #
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
        results = results.dropna()

        result2 = pd.DataFrame(
            {
                "Pseudo R-squared ": model.pr2,
                "Spatial Pseudo R-squared": model.pr2_e,
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
# spatial ML_Error node
############################################


@knext.node(
    name="Spatial Error Model",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "MLErr.png",
    after="",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial ML_Error model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial ML_Error model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial ML_Error model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial ML_Error model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial ML_Error model.",
)
class SpatialML_Error:
    """Spatial ML_Error.
    ML estimation of the spatial error model with all results and diagnostics. More details can be found at Luc Anselin. Spatial Econometrics: Methods and Models. Kluwer, Dordrecht, 1988.
    """

    geo_col = knut.geo_col_parameter()

    id_col = get_id_col_parameter()

    dependent_variable = knext.ColumnParameter(
        "Dependent variable",
        "The column containing the dependent variable to use for the calculation of the spatial ML_Error model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    independent_variables = knext.MultiColumnParameter(
        "Independent variables",
        "The columns containing the independent variables to use for the calculation of the spatial ML_Error model.",
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

        if "none" not in str(self.id_col).lower():
            gdf.index = range(len(gdf))
            id_map = dict(zip(gdf[self.id_col], gdf.index))
            adjust_list["focal"] = adjust_list["focal"].map(id_map)
            adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)

        from libpysal.weights import W

        w = W.from_adjlist(adjust_list)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        import spreg

        model = spreg.ML_Error(
            y,
            X,
            w=w,
            method="ord",
            name_x=self.independent_variables,
            name_y=self.dependent_variable,
        )

        import pandas as pd

        results = pd.DataFrame(
            [model.name_x, model.betas, model.std_err, model.z_stat]
        ).T
        results.columns = ["Variable", "Coefficient", "Std.Error", "Z-Statistic"]
        results.loc[:, "Coefficient"] = results.loc[:, "Coefficient"].map(
            lambda x: x[0]
        )
        results.loc[:, "Probability"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[1]
        )
        results.loc[:, "Z-Statistic"] = results.loc[:, "Z-Statistic"].map(
            lambda x: x[0]
        )
        #
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
        results = results.dropna()

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
