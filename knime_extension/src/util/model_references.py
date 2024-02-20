model_references = {
    "Aka74": "Hirotugu Akaike. A new look at the statistical model identification. IEEE transactions on automatic control, 19(6):716–723, 1974.",
    "Ans88": "Luc Anselin. Spatial Econometrics: Methods and Models. Kluwer, Dordrecht, 1988.",
    "Ans11": "Luc Anselin. GMM estimation of spatial error autocorrelation with and without heteroskedasticity. Technical Report, GeoDa Center for Geospatial Analysis and Computation, 2011.",
    "ABFY96": "Luc Anselin, Anil K Bera, Raymond Florax, and Mann J Yoon. Simple diagnostic tests for spatial dependence. Regional science and urban economics, 26(1):77–104, 1996.",
    "AK97": "Luc Anselin and Harry H Kelejian. Testing for spatial error autocorrelation in the presence of endogenous regressors. International Regional Science Review, 20(1-2):153–182, 1997.",
    "ADKP10": "Irani Arraiz, David M. Drukker, Harry H. Kelejian, and Ingmar R. Prucha. A spatial Cliff-Ord-type model with heteroskedastic innovations: Small and large sample results. Journal of Regional Science, 50(2):592–614, 2010. doi:10.1111/j.1467-9787.2009.00618.x.",
    "BKW05": "David A Belsley, Edwin Kuh, and Roy E Welsch. Regression diagnostics: Identifying influential data and sources of collinearity. Volume 571. John Wiley & Sons, 2005.",
    "BP79": "Trevor S Breusch and Adrian R Pagan. A simple test for heteroscedasticity and random coefficient variation. Econometrica: Journal of the Econometric Society, pages 1287–1294, 1979.",
    "DEP13": "David M Drukker, Peter Egger, and Ingmar R Prucha. On two-step estimation of a spatial autoregressive model with autoregressive disturbances and endogenous regressors. Econometric Reviews, 32(5-6):686–733, 2013.",
    "DPR13": "David M. Drukker, Ingmar R. Prucha, and Rafal Raciborski. A command for estimating spatial-autoregressive models with spatial-autoregressive disturbances and additional endogenous variables. The Stata Journal, 13(2):287–301, 2013. URL: https://journals.sagepub.com/doi/abs/10.1177/1536867X1301300203.",
    "Gre03": "William H Greene. Econometric analysis. Pearson Education India, 2003.",
    "JB80": "Carlos M Jarque and Anil K Bera. Efficient tests for normality, homoscedasticity and serial independence of regression residuals. Economics letters, 6(3):255–259, 1980.",
    "KP99": "H H Kelejian and I R Prucha. A generalized moments estimator for the autoregressive parameter in a spatial model. Int. Econ. Rev., 40:509–534, 1999.",
    "KP98": "Harry H Kelejian and Ingmar R Prucha. A generalized spatial two-stage least squares procedure for estimating a spatial autoregressive model with autoregressive disturbances. J. Real Estate Fin. Econ., 17(1):99–121, 1998.",
    "KBJ82": "Roger Koenker and Gilbert Bassett Jr. Robust tests for heteroscedasticity based on regression quantiles. Econometrica: Journal of the Econometric Society, pages 43–61, 1982.",
    "S+78": "Gideon Schwarz and others. Estimating the dimension of a model. The annals of statistics, 6(2):461–464, 1978.",
    "Whi80": "Halbert White. A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. Econometrica: Journal of the Econometric Society, pages 817–838, 1980.",
}

# get w from adjust list

def get_w_from_adjust_list(adj_list, gdf, id_col):
    """
    Create a W object from an adjacency list.
    """

    from libpysal.weights import W

    adjust_list = adj_list.to_pandas()
    if "none" not in str(id_col).lower():
        gdf.index = range(len(gdf))
        id_map = dict(zip(gdf[id_col], gdf.index))
        adjust_list["focal"] = adjust_list["focal"].map(id_map)
        adjust_list["neighbor"] = adjust_list["neighbor"].map(id_map)
    w = W.from_adjlist(adjust_list)

    return w, gdf