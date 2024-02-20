import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import util.modeling_utils as mut
import util.model_references as mrs

__category = knext.category(
    path="/community/geo",
    level_id="spatialmodels",
    name="Spatial Modelling",
    description="Nodes that conduct spatial regression, panel modelling, and geographically weighted regression analyses.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialModelCategory.png",
    after="spatialstatistic",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/SpatialModel/"


@knext.parameter_group("Advanced Settings", since_version="1.2.0")
class AdvancedLsSetting:
    """Advanced Setting"""

    vm = knext.BoolParameter(
        "VM",
        """If True, include variance-covariance matrix in summary results. Default set to False.""",
        default_value=False,
        is_advanced=True,
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
    Spatial two stage least squares (S2SLS) with results and diagnostics. More details can be found in the following reference, Luc Anselin.
    Spatial Econometrics: Methods and Models. Kluwer. Dordrecht, 1988.

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    # input parameters
    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

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
        "If selected, the node computes the Anselin-Kelejian test.",
        default_value=False,
    )

    robust = knext.StringParameter(
        "Robust",
        """If ‘white’, then a White consistent estimator of the variance-covariance matrix is given. If ‘hac’, 
        then a HAC consistent estimator of the variance-covariance matrix is given. Set to None for default.""",
        enum=["white", "hac", "none"],
        default_value="none",
    )

    advanced_settings = AdvancedLsSetting()

    sig2n_k = knext.BoolParameter(
        "Sig2n k",
        """If True, then use n-k to estimate sigma^2. If False, use n. Default set to True.""",
        default_value=False,
        is_advanced=True,
        since_version="1.2.0",
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
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
            "sig2n_k": self.sig2n_k,
            "vm": self.advanced_settings.vm,
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


@knext.parameter_group("Advanced Settings", since_version="1.2.0")
class AdvancedLagErrorSetting:
    """Advanced Setting"""

    epsilon = knext.DoubleParameter(
        "Epsilon",
        """tolerance criterion in mimimize_scalar function and inverse_product""",
        default_value=1e-07,
        is_advanced=True,
    )

    vm = knext.BoolParameter(
        "VM",
        """if True, include variance-covariance matrix in summary results.""",
        default_value=False,
        is_advanced=True,
    )

    # name_y = knext.StringParameter(
    #     "Name Y",
    #     """Name of dependent variable for use in output, if None use the column name of dependent variable""",
    #     default_value="None",
    #     is_advanced=True,
    # )

    # name_x = knext.StringParameter(
    #     "Name X",
    #     """Names of independent variables for use in output, if None use the column names of independent variables""",
    #     default_value="None",
    #     is_advanced=True,
    # )

    # name_w = knext.StringParameter(
    #     "Name W",
    #     """Name of weights matrix for use in output""",
    #     default_value="None",
    #     is_advanced=True,
    # )

    # name_ds = knext.StringParameter(
    #     "Name DS",
    #     """Name of dataset for use in output""",
    #     default_value="None",
    #     is_advanced=True,
    # )


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
    Spatial Lag Panel Model with Fixed Effects. ML estimation of the fixed effects spatial lag model with all results and diagnostics. More details can be found at J. Paul Elhorst.
    Specification and estimation of spatial panel data models. International Regional Science Review, 26(3):244–268, 2003. doi:10.1177/0160017603253791.

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    advanced_settings = AdvancedLagErrorSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
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
            "epsilon": self.advanced_settings.epsilon,
            "vm": self.advanced_settings.vm,
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
    Spatial Error Panel Model with Fixed Effects node. ML estimation of the fixed effects spatial error model with all results and diagnostics.
    More details can be found at J. Paul Elhorst. Specification and estimation of spatial panel data models. International Regional Science Review, 26(3):244–268, 2003. doi:10.1177/0160017603253791.

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    advanced_settings = AdvancedLagErrorSetting()

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
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
            "epsilon": self.advanced_settings.epsilon,
            "vm": self.advanced_settings.vm,
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


@knext.parameter_group("Advanced Settings", since_version="1.2.0")
class AdvancedGWRSetting:
    """Advanced Setting"""

    sigma2_v1 = knext.BoolParameter(
        "Sigma2 v1",
        """specify form of corrected denominator of sigma squared to use for model diagnostics; 
        Acceptable options are:
        ‘True’: n-tr(S) (default) ‘False’: n-2(tr(S)+tr(S’S))""",
        default_value=True,
        is_advanced=True,
    )

    constant = knext.BoolParameter(
        "Constant",
        """True to include intercept (default) in model and False to exclude intercept.""",
        default_value=True,
        is_advanced=True,
    )

    spherical = knext.BoolParameter(
        "Spherical",
        """True for spherical coordinates (long-lat), False for projected coordinates (default).""",
        default_value=False,
        is_advanced=True,
    )

    hat_matrix = knext.BoolParameter(
        "Hat Matrix",
        """True to store full n by n hat matrix, False to not store full hat matrix to minimize memory footprint (default).""",
        default_value=False,
        is_advanced=True,
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
    Performs Geographically Weighted Regression (GWR), a local form of linear regression used to model spatially varying relationships. Can currently estimate Gaussian, Poisson, and logistic models (built on a GLM framework).
    More details can be found at [here](https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-statistics-toolbox/geographically-weighted-regression.htm).

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    fixed = knext.BoolParameter(
        "Fixed bandwidth",
        """True for distance-based kernel function and False for adaptive (nearest neighbor) kernel function (default)""",
        default_value=False,
        since_version="1.2.0",
    )

    kernel = knext.StringParameter(
        "Kernel",
        "Type of kernel function used to weight observations; available options: ‘gaussian’, ‘bisquare’, ‘exponential’.",
        default_value="bisquare",
        enum=["gaussian", "bisquare", "exponential"],
    )

    # use_bindwidth_search = knext.BoolParameter(
    #     "Use Bindwidth Search",
    #     "If True, then use bandwidth search. If False, use bandwidth value.",
    #     default_value=True,
    #     since_version="1.2.0",
    # )

    search_method = knext.StringParameter(
        "Search method",
        "Bw search method: ‘golden’, ‘interval’.",
        default_value="golden",
        enum=["golden", "interval"],
    )
    # .rule(knext.OneOf(use_bindwidth_search, [True]), knext.Effect.SHOW)

    bandwith_min = knext.IntParameter(
        "Bandwidth min",
        "Min value used in bandwidth search.",
        default_value=2,
    )
    # .rule(knext.OneOf(use_bindwidth_search, [True]), knext.Effect.SHOW)

    # bandwith_max = knext.IntParameter(
    #     "Bandwidth Max",
    #     "max value used in bandwidth search",
    #     default_value=200,
    #     since_version="1.2.0",
    # ).rule(knext.OneOf(use_bindwidth_search, [True]), knext.Effect.SHOW)

    # interval = knext.IntParameter(
    #     "Interval",
    #     "Interval used in bandwidth search",
    #     default_value=1,
    #     since_version="1.2.0",
    # ).rule(knext.OneOf(use_bindwidth_search, [True]), knext.Effect.SHOW)

    # criterion = knext.StringParameter(
    #     "Criterion",
    #     "Criterion used in bandwidth search: ‘AICc’, ‘AIC’, ‘BIC’, ‘CV’",
    #     default_value="AICc",
    #     enum=["AICc", "AIC", "BIC", "CV"],
    #     since_version="1.2.0",
    # ).rule(knext.OneOf(use_bindwidth_search, [True]), knext.Effect.SHOW)

    # bw = knext.IntParameter(
    #     "Bandwidth",
    #     "bandwidth value consisting of either a distance or N nearest neighbors",
    #     default_value=100,
    #     since_version="1.2.0",
    #     # is_advanced=True,
    # ).rule(knext.OneOf(use_bindwidth_search, [False]), knext.Effect.SHOW)

    # advanced_settings = AdvancedGWRSetting()

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

        from mgwr.gwr import GWR

        # kws = {
        #     "coords": g_coords,
        #     "y": g_y,
        #     "X": g_X,
        #     "fixed": self.fixed,
        #     "kernel": self.kernel,
        #     "constant": self.advanced_settings.constant,
        #     "sigma2_v1": self.advanced_settings.sigma2_v1,
        #     "spherical": self.advanced_settings.spherical,
        #     "hat_matrix": self.advanced_settings.hat_matrix,
        # }

        # # Calibrate GWR model
        # if self.use_bindwidth_search:
        #     from mgwr.sel_bw import Sel_BW

        # gwr_selector = Sel_BW(
        #     g_coords, g_y, g_X, kernel=self.kernel, fixed=self.fixed
        # )

        #     gwr_bw = gwr_selector.search(
        #         bw_min=self.bandwith_min,
        #         bw_max=self.bandwith_max,
        #         interval=self.interval,
        #         criterion=self.criterion,
        #         search_method=self.search_method,
        #     )

        #     kws["bw"] = gwr_bw
        # else:
        #     kws["bw"] = self.bw
        # kws["fixed"] = self.fixed

        from mgwr.sel_bw import Sel_BW

        gwr_selector = Sel_BW(g_coords, g_y, g_X)
        gwr_bw = gwr_selector.search(bw_min=self.bandwith_min)
        # print(gwr_bw)

        # gwr_model = GWR(**kws).fit()
        gwr_model = GWR(g_coords, g_y, g_X, gwr_bw, fixed=self.fixed).fit()

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

        # if self.use_bindwidth_search:
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

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    independent_variables = mut.get_dependent_and_independent_variables()[-1]

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

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    fixed = knext.BoolParameter(
        "Fixed bandwidth",
        """True for distance-based kernel function and False for adaptive (nearest neighbor) kernel function (default).""",
        default_value=False,
        since_version="1.2.0",
    )

    kernel = knext.StringParameter(
        "Kernel",
        "Type of kernel function used to weight observations; available options: ‘gaussian’, ‘bisquare’, ‘exponential’.",
        default_value="bisquare",
        enum=["gaussian", "bisquare", "exponential"],
        # since_version="1.2.0",
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
        "Bandwidth min",
        "Min value used in bandwidth search.",
        default_value=2,
    )

    bandwith_max = knext.IntParameter(
        "Bandwidth Max",
        "Max value used in bandwidth search.",
        default_value=200,
        since_version="1.2.0",
    )

    interval = knext.IntParameter(
        "Interval",
        "Interval used in bandwidth search.",
        default_value=1,
        since_version="1.2.0",
    )

    criterion = knext.StringParameter(
        "Criterion",
        "Criterion used in bandwidth search: ‘AICc’, ‘AIC’, ‘BIC’, ‘CV’.",
        default_value="AICc",
        enum=["AICc", "AIC", "BIC", "CV"],
        since_version="1.2.0",
    )

    advanced_settings = AdvancedGWRSetting()

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

        mgwr_selector = Sel_BW(
            g_coords, g_y, g_X, multi=True, kernel=self.kernel, fixed=self.fixed
        )
        mgwr_bw = mgwr_selector.search(
            multi_bw_min=[self.bandwith_min],
            #    multi_bw_max=[self.bandwith_max],
            #    multi_bw_interval=[self.interval],
            #    criterion=self.criterion,
            #    search_method=self.search_method
        )

        from mgwr.gwr import MGWR

        mgwr_model = MGWR(
            g_coords,
            g_y,
            g_X,
            mgwr_selector,
            fixed=self.fixed,
            sigma2_v1=self.advanced_settings.sigma2_v1,
            kernel=self.kernel,
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

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    Spatial_Diagnostics = knext.BoolParameter(
        "Spatial Diagnostics",
        "If selected, the node computes the Anselin-Kelejian test",
        default_value=True,
        since_version="1.2.0",
    )

    robust = knext.StringParameter(
        "Robust",
        """If ‘white’, then a White consistent estimator of the variance-covariance matrix is given. If ‘hac’, 
        then a HAC consistent estimator of the variance-covariance matrix is given. Set to None for default.""",
        enum=["white", "hac", "none"],
        default_value="none",
        since_version="1.2.0",
    )

    advanced_settings = AdvancedLsSetting()

    sig2n_k = knext.BoolParameter(
        "Sig2n k",
        """If True, then use n-k to estimate sigma^2. If False, use n. Default set to True.""",
        default_value=True,
        is_advanced=True,
        since_version="1.2.0",
    )

    moran = knext.BoolParameter(
        "Moran",
        """If True, compute Moran’s I on the residuals. Note: requires spat_diag=True.
        Default set to True.""",
        default_value=True,
        is_advanced=True,
        since_version="1.2.0",
    )

    white_test = knext.BoolParameter(
        "White Test",
        """If True, compute White’s specification robust test (requires nonspat_diag=True).
        Default set to True.""",
        default_value=True,
        is_advanced=True,
        since_version="1.2.0",
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        import spreg

        kws = {
            "y": y,
            "x": X,
            "w": w,
            "spat_diag": self.Spatial_Diagnostics,
            "name_y": self.dependent_variable,
            "name_x": self.independent_variables,
            "name_w": "Spatial Weights",
            "name_ds": "Input Table",
            "sig2n_k": self.sig2n_k,
            "vm": self.advanced_settings.vm,
            "moran": self.moran,
            "white_test": self.white_test,
        }
        if "none" not in str(self.robust).lower():
            kws["robust"] = self.robust

        model = spreg.OLS(**kws)
        # model = spreg.OLS(
        #     y,
        #     X,
        #     w=w,
        #     spat_diag=True,
        #     moran=True,
        #     white_test=True,
        #     name_y=self.dependent_variable,
        #     name_x=self.independent_variables,
        # )

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

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    advanced_settings = AdvancedLagErrorSetting()

    method = knext.StringParameter(
        "Method",
        """if ‘full’, brute force calculation (full matrix expressions) if ‘ord’, 
        Ord eigenvalue method if ‘LU’, LU sparse matrix decomposition""",
        default_value="ord",
        enum=["full", "ord", "LU"],
        since_version="1.2.0",
        is_advanced=True,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        import spreg

        model = spreg.ML_Lag(
            y,
            X,
            w=w,
            method=self.method,
            name_x=self.independent_variables,
            name_y=self.dependent_variable,
            # name_w=self.advanced_settings.name_w,
            # name_ds=self.advanced_settings.name_ds,
            epsilon=self.advanced_settings.epsilon,
            vm=self.advanced_settings.vm,
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

    **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them.
    """

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    advanced_settings = AdvancedLagErrorSetting()

    method = knext.StringParameter(
        "Method",
        """if ‘full’, brute force calculation (full matrix expressions) if ‘ord’, 
        Ord eigenvalue method if ‘LU’, LU sparse matrix decomposition""",
        default_value="ord",
        enum=["full", "ord", "LU"],
        since_version="1.2.0",
        is_advanced=True,
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        # Prepare Georgia dataset inputs
        X = gdf[self.independent_variables].values
        y = gdf[self.dependent_variable].values

        import spreg

        model = spreg.ML_Error(
            y,
            X,
            w=w,
            method=self.method,
            name_x=self.independent_variables,
            name_y=self.dependent_variable,
            # name_w=self.advanced_settings.name_w,
            # name_ds=self.advanced_settings.name_ds,
            epsilon=self.advanced_settings.epsilon,
            vm=self.advanced_settings.vm,
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


# New added from Geolab in 2023-11-20

SHOULD_NOT_CONTAIN_MISSING_VALUES = """
   **Note:** The input table should not contain missing values. You can use the
    [Missing Value](https://hub.knime.com/knime/extensions/org.knime.features.base/latest/org.knime.base.node.preproc.pmml.missingval.compute.MissingValueHandlerNodeFactory/) node to replace them."""

# Root path for all node icons in this file
############################################
# spatial GM_Error node
############################################


@knext.node(
    name="Spatial GM Error",
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
    description="Input Table with spatial weights for calculation of the spatial GM Error model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Error model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Error model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Error model.",
)
class SpatialGM_Error:
    (
        """
    Spatial GM Error Model.
    GMM method for a spatial error model, with results and diagnostics; based on Kelejian and Prucha (1998, 1999) [KP98] [KP99]. More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Error.html#spreg.GM_Error). Please refer the following papers for more details.
    - %s
    - %s
    """
        % (mrs.model_references["KP98"], mrs.model_references["KP99"])
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        import pandas as pd
        import spreg

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
    name="Spatial GM Error Het",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMErrHet.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM Error Het model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM Error Het model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Error Het model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Error Het model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Error Het model.",
)
class SpatialGM_Error_Het:
    (
        """
    Spatial GM Error Het.
    GMM method for a spatial error model with heteroskedasticity, with results and diagnostics; based on Kelejian and Prucha (1998). More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Error_Het.html#spreg.GM_Error_Het). Please refer the following papers for more details.
    - Irani Arraiz, David M. Drukker, Harry H. Kelejian, and Ingmar R. Prucha. A spatial Cliff-Ord-type model with heteroskedastic innovations: Small and large sample results. Journal of Regional Science, 50(2):592–614, 2010. doi:10.1111/j.1467-9787.2009.00618.x
    - Luc Anselin. GMM estimation of spatial error autocorrelation with and without heteroskedasticity. Technical Report, GeoDa Center for Geospatial Analysis and Computation, 2011.
    """
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)

        import pandas as pd
        import spreg

        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
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
    name="Spatial GM Error Hom",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMErrHom.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM Error Hom model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM Error Hom model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Error Hom model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Error Hom model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Error Hom model.",
)
class SpatialGM_Error_Hom:
    (
        """
    Spatial GM Error Hom Model.
    GMM method for a spatial error model with homoskedasticity, with results and diagnostics; based on Drukker et al. (2013),  following Anselin (2011). More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Error_Hom.html#spreg.GM_Error_Hom). Please refer the following papers for more details.

    - David M Drukker, Peter Egger, and Ingmar R Prucha. On two-step estimation of a spatial autoregressive model with autoregressive disturbances and endogenous regressors. Econometric Reviews, 32(5-6):686–733, 2013.
    - Luc Anselin. GMM estimation of spatial error autocorrelation with and without heteroskedasticity. Technical Report, GeoDa Center for Geospatial Analysis and Computation, 2011.
    """
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        import pandas as pd
        import spreg

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
    name="Spatial GM Combo",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMcombo.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM Combo model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM Combo model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Combo model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Combo model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Combo model.",
)
class SpatialGM_Combo:
    (
        """
    Spatial GM Combo Model.
    GMM method for a spatial lag and error model with endogenous variables, with results and diagnostics; based on Kelejian and Prucha (1998, 1999). More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Combo.html#spreg.GM_Combo). Please refer the following papers for more details.

    - Harry H Kelejian and Ingmar R Prucha. A generalized spatial two-stage least squares procedure for estimating a spatial autoregressive model with autoregressive disturbances. J. Real Estate Fin. Econ., 17(1):99–121, 1998.
    - H H Kelejian and I R Prucha. A generalized moments estimator for the autoregressive parameter in a spatial model. Int. Econ. Rev., 40:509–534, 1999.
    """
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        import pandas as pd
        import spreg

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
    name="Spatial GM Combo Het",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMcomboHet.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM Combo Het model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM Combo Het model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Combo Het model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Combo Het model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Combo Het model.",
)
class SpatialGM_Combo_Het:
    (
        """
    Spatial GM Combo Het Model.
    GMM method for a spatial lag and error model with heteroskedasticity and endogenous variables, with results and diagnostics; More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Combo_Het.html#spreg.GM_Combo_Het). Please refer the following papers for more details.

    - Irani Arraiz, David M. Drukker, Harry H. Kelejian, and Ingmar R. Prucha. A spatial Cliff-Ord-type model with heteroskedastic innovations: Small and large sample results. Journal of Regional Science, 50(2):592–614, 2010. doi:10.1111/j.1467-9787.2009.00618.x.
    - Luc Anselin. GMM estimation of spatial error autocorrelation with and without heteroskedasticity. Technical Report, GeoDa Center for Geospatial Analysis and Computation, 2011.
    """
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        import pandas as pd
        import spreg

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
    name="Spatial GM Combo Hom",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMcomboHom.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM Combo Hom model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM Combo Hom model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Combo Hom model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Combo Hom model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Combo Hom model.",
)
class SpatialGM_Combo_Hom:
    (
        """
    Spatial GM Combo Hom Model.
    GMM method for a spatial lag and error model with homoskedasticity and endogenous variables, with results and diagnostics; based on Drukker et al. (2013) [DEP13], following Anselin (2011) [Ans11]. More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Combo_Hom.html#spreg.GM_Combo_Hom). Please refer the following papers for more details.
    - David M Drukker, Peter Egger, and Ingmar R Prucha. On two-step estimation of a spatial autoregressive model with autoregressive disturbances and endogenous regressors. Econometric Reviews, 32(5-6):686–733, 2013.
    - Luc Anselin. GMM estimation of spatial error autocorrelation with and without heteroskedasticity. Technical Report, GeoDa Center for Geospatial Analysis and Computation, 2011.
    """
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    def configure(self, configure_context, input_schema, input_schema_2):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        import pandas as pd
        import spreg

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
    name="Spatial GM Endog Error",
    node_type=knext.NodeType.LEARNER,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GMendogErr.png",
)
@knext.input_table(
    name="Input Table",
    description="Input Table with dependent and independent variables for calculation of the spatial GM Endog Error model.",
)
@knext.input_table(
    name="Spatial Weights",
    description="Input Table with spatial weights for calculation of the spatial GM Endog Error model.",
)
@knext.output_table(
    name="Output Table",
    description="Description of the spatial GM Endog Error model.",
)
@knext.output_table(
    name="Variable and Coefficient Table",
    description="Variable and Coefficient Table of the spatial GM Endog Error model.",
)
@knext.output_view(
    name="Model summary view",
    description="Model summary view of the spatial GM Endog Error model.",
)
class SpatialGM_Endog_Error:
    (
        """
    Spatial GM Endog Error Model.
    GMM method for a spatial error model with endogenous variables, with results and diagnostics; based on Kelejian and Prucha (1998, 1999) [KP98] [KP99]. More information can be found at [here](https://spreg.readthedocs.io/en/latest/generated/spreg.GM_Endog_Error.html#spreg.GM_Endog_Error). Please refer the following papers for more details.

    - Harry H Kelejian and Ingmar R Prucha. A generalized spatial two-stage least squares procedure for estimating a spatial autoregressive model with autoregressive disturbances. J. Real Estate Fin. Econ., 17(1):99–121, 1998.
    - H H Kelejian and I R Prucha. A generalized moments estimator for the autoregressive parameter in a spatial model. Int. Econ. Rev., 40:509–534, 1999.
    """
        + SHOULD_NOT_CONTAIN_MISSING_VALUES
    )

    geo_col = knut.geo_col_parameter()

    id_col = mut.get_id_col_parameter()

    (dependent_variable, independent_variables) = (
        mut.get_dependent_and_independent_variables()
    )

    yend = knext.ColumnParameter(
        "Endogenous variable",
        "The column containing the endogenous variable to use for the calculation of the spatial GM Endog Error model.",
        column_filter=knut.is_numeric,
        include_none_column=False,
    )

    q = knext.ColumnParameter(
        "External exogenous variable",
        "The column containing the external exogenous variable to use for the calculation of the spatial GM Endog Error model.",
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
        w, gdf = mut.get_w_from_adjust_list(input_2, gdf, self.id_col)
        import pandas as pd
        import spreg

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
