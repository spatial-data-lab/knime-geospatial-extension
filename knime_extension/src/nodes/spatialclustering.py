import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut

__category = knext.category(
    path="/community/geo",
    level_id="spatialclustering",
    name="Spatial Clustering",
    description="Spatial Clustering (Regionalization).",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialclusteringCategory.png",
    after="opendataset",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Spatialclustering/"

# def  geodataframe to PPP
def gdf2ppp(gdf):
    from pointpats import PointPattern

    ngdf = gp.GeoDataFrame(geometry=gdf.geometry, crs=gdf.crs)
    ngdf["coords"] = ngdf.apply(lambda x: list(x.geometry.coords[0]), axis=1)
    ppts = ngdf.coords.to_list()
    pp = PointPattern(ppts)
    return pp


############################################
# MeanCenter and StandarDistance
############################################
@knext.node(
    name="Mean Center",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "MeanCenter.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with Mean center and standard distance .",
)
class MeanCenterNode:
    """
    Mean center and standard distance.

    This node use the Pysal package [pointpats](https://pointpats.readthedocs.io/en/master/api.html)
    to  measure the compactness of a spatial distribution of features around its mean center.
    Standard distance (or standard distance deviation) is usually represented as a circle where the radius of the circle is the standard distance.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

    weight_col = knext.ColumnParameter(
        "Weight column",
        "Select the weighting column for calculatethe meann center.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=True,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pointpats

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        pp = gdf2ppp(gdf)
        if (
            (self.weight_col is None)
            or (self.weight_col == "None")
            or (self.weight_col == "<none>")
        ):
            mc = pointpats.mean_center(pp.points)
        else:
            weights = gdf[self.weight_col]
            mc = pointpats.weighted_mean_center(pp.points, weights)
        mcgdf = gp.GeoDataFrame(
            {"x": [mc[0]], "y": [mc[1]]},
            geometry=gp.GeoSeries.from_xy(x=[mc[0]], y=[mc[1]], crs=gdf.crs),
        )
        # Standard distance
        mcgdf["std_distance"] = pointpats.std_distance(pp.points)
        mcgdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(mcgdf)


############################################
# Ellipse
############################################
@knext.node(
    name="Standard Deviational Ellipse",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "StandardEllipse.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with Standard deviational ellips.",
)
class StandardEllipseNode:
    """
    Standard deviational ellipse.

    This node use the Pysal package [pointpats](https://pointpats.readthedocs.io/en/master/api.html)
    to  calculate parameters of standard deviational ellipse for a point pattern, which is a common way of measuring
    the trend for a set of points or areas. These measures define the axes of an ellipse (or ellipsoid) encompassing the distribution of features.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pointpats
        import numpy as np
        from matplotlib.patches import Ellipse
        from shapely.geometry import Polygon

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        pp = gdf2ppp(gdf)
        sx, sy, theta = pointpats.ellipse(pp.points)
        theta_degree = np.degrees(theta)
        e = Ellipse(
            xy=pointpats.mean_center(pp.points),
            width=sx * 2,
            height=sy * 2,
            angle=-theta_degree,
        )
        # angle is rotation in degrees (anti-clockwise)
        # get the vertices from the ellipse object
        vertices = e.get_verts()
        ellipse = Polygon(vertices)
        gdf = gp.GeoDataFrame(geometry=gp.GeoSeries(ellipse), crs=gdf.crs)
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# SKATER
############################################
@knext.node(
    name="SKATER",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "skater.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class SKATERNode:
    """
    The Spatial C(K)luster Analysis by Tree Edge Removal(SKATER).

    This node use the Pysal package [pointpats](https://pointpats.readthedocs.io/en/master/api.html)
    to  calculate parameters of standard deviational ellipse for a point pattern, which is a common way of measuring
    The Spatial C(K)luster Analysis by Tree Edge Removal(SKATER) algorithm introduced by Assuncao et al. (2006)
    is based on the optimal pruning of a minimum spanning tree that reflects the contiguity structure among the observations.
    It provides an optimized algorithm to prune to tree into several clusters that their values of selected variables are as similar as possible.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    bound_col = knext.ColumnParameter(
        "Bound column for minibound",
        "Select the bound column for clusters with minibound.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    attribute_list = knext.StringParameter(
        "Attribute columns for Clustering",
        "Input the column names with Semicolon.",
        "",
    )
    cluster_k = knext.IntParameter(
        "Number of clusters",
        "The Number of user-defined clusters.",
        default_value=4,
    )
    minibound = knext.DoubleParameter(
        "Minimum sum value for bound variable in output clusters",
        "A minimum value that the sum value of bounding variable in each cluster should be greater than.",
    )
    weight_mode = knext.StringParameter(
        "Spatial Weight model",
        "Input spatial weight mode.",
        "Queen",
        enum=["Queen", "Rook"],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):

        import pygeoda

        k = self.cluster_k
        controlVar = self.bound_col
        m_bound = self.minibound
        attributelist = self.bound_col.split(";")
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        geodadf = pygeoda.open(gdf)
        data = geodadf[attributelist]
        b_vals = geodadf.GetRealCol(controlVar)
        if self.weight_mode == "Queen":
            w = pygeoda.queen_weights(geodadf)
        else:
            w = pygeoda.rook_weights(geodadf)
        skatercluster = pygeoda.skater(
            k,
            w,
            data,
            distance_method="euclidean",
            bound_variable=b_vals,
            min_bound=m_bound,
        )
        gdf["clusterid"] = skatercluster["Clusters"]
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# REDCAP
############################################
@knext.node(
    name="REDCAP",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "redcap.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class REDCAPNode:
    """
    REDCAP (Regionalization with dynamically constrained agglomerative clustering and partitioning).

    It is developed by D. Guo (2008). Like SKATER, REDCAP starts from building a spanning tree with
    4 different ways (single-linkage, average-linkage, complete-linkage and wards-linkage).
    Then,REDCAP provides 2 different ways (first‐order and full-order constraining) to prune
    the tree to find clusters. The first-order approach with a minimum spanning tree is exactly
    the same with SKATER. In GeoDa and pygeoda, the following methods are provided:

     - First-order and Single-linkage
     - Full-order and Complete-linkage
     - Full-order and Average-linkage
     - Full-order and Single-linkage
     - Full-order and Wards-linkage
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    bound_col = knext.ColumnParameter(
        "Bound column for minibound",
        "Select the bound column for clusters with minibound.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    attribute_list = knext.StringParameter(
        "Attribute columns for Clustering",
        "Input the column names with Semicolon.",
        "",
    )
    cluster_k = knext.IntParameter(
        "Number of clusters",
        "The Number of user-defined clusters.",
        default_value=4,
    )
    minibound = knext.DoubleParameter(
        "Minimum sum value for bound variable in output clusters",
        "A minimum value that the sum value of bounding variable in each cluster should be greater than.",
    )
    weight_mode = knext.StringParameter(
        "Spatial Weight model",
        "Input spatial weight mode.",
        "Queen",
        enum=["Queen", "Rook"],
    )
    link_mode = knext.StringParameter(
        "Linkage model",
        "Input linkage mode.",
        "fullorder-completelinkage",
        enum=[
            "firsorder-singlelinkage",
            "fullorder-singlelinkage",
            "fullorder-completelinkage",
            "fullorder-averagelinkage",
            "fullorder-wardlinkage",
        ],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):

        import pygeoda

        k = self.cluster_k
        controlVar = self.bound_col
        m_bound = self.minibound
        attributelist = self.bound_col.split(";")
        linkage = self.link_mode
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        geodadf = pygeoda.open(gdf)
        data = geodadf[attributelist]
        b_vals = geodadf.GetRealCol(controlVar)

        if self.weight_mode == "Queen":
            w = pygeoda.queen_weights(geodadf)
        else:
            w = pygeoda.rook_weights(geodadf)
        final_cluster = pygeoda.redcap(
            k,
            w,
            data,
            linkage,
            bound_variable=b_vals,
            min_bound=m_bound,
        )
        gdf["clusterid"] = final_cluster["Clusters"]
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# SCHC
############################################
@knext.node(
    name="SCHC",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "schc.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class SCHCNode:
    """
    SCHC (Spatially constrained hierarchical clustering).

    It is a special form of constrained clustering, where the constraint is based on contiguity
     (common borders). The method builds up the clusters using agglomerative hierarchical clustering methods:
     single linkage, complete linkage, average linkage and Ward’s method (a special form of centroid linkage).
     Meanwhile, it also maintains the spatial contiguity when merging two clusters. The method builds up the
     clusters using agglomerative hierarchical clustering methods:

     - single
     - complete
     - average
     - ward
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    bound_col = knext.ColumnParameter(
        "Bound column for minibound",
        "Select the bound column for clusters with minibound.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    attribute_list = knext.StringParameter(
        "Attribute columns for Clustering",
        "Input the column names with Semicolon.",
        "",
    )
    cluster_k = knext.IntParameter(
        "Number of clusters",
        "The Number of user-defined clusters.",
        default_value=4,
    )
    minibound = knext.DoubleParameter(
        "Minimum sum value for bound variable in output clusters",
        "A minimum value that the sum value of bounding variable in each cluster should be greater than.",
    )
    weight_mode = knext.StringParameter(
        "Spatial Weight model",
        "Input spatial weight mode.",
        "Queen",
        enum=["Queen", "Rook"],
    )
    link_mode = knext.StringParameter(
        "Linkage model",
        "Input linkage mode.",
        "ward",
        enum=[
            "single",
            "complete",
            "average",
            "ward",
        ],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):

        import pygeoda

        k = self.cluster_k
        controlVar = self.bound_col
        m_bound = self.minibound
        attributelist = self.bound_col.split(";")
        linkage = self.link_mode
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        geodadf = pygeoda.open(gdf)
        data = geodadf[attributelist]
        b_vals = geodadf.GetRealCol(controlVar)

        if self.weight_mode == "Queen":
            w = pygeoda.queen_weights(geodadf)
        else:
            w = pygeoda.rook_weights(geodadf)
        final_cluster = pygeoda.schc(
            k,
            w,
            data,
            linkage,
            bound_variable=b_vals,
            min_bound=m_bound,
        )
        gdf["clusterid"] = final_cluster["Clusters"]
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# MaxP Greedy
############################################
@knext.node(
    name="MaxP-Greedy",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "maxp.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class MaxPgreedyNode:
    """
    A greedy algorithm to solve the max-p-region problem.

    The so-called max-p regions model (outlined in Duque, Anselin, and Rey 2012)
    uses a different approach and considers the regionalization problem as an application of
    integer programming. In addition, the number of regions is determined endogenously.

    The algorithm itself consists of a search process that starts with an initial feasible
    solution and iteratively improves upon it while maintaining contiguity among the elements of each cluster.

    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    bound_col = knext.ColumnParameter(
        "Bound column for minibound",
        "Select the bound column for clusters with minibound.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    attribute_list = knext.StringParameter(
        "Attribute columns for Clustering",
        "Input the column names with Semicolon.",
        "",
    )
    minibound = knext.DoubleParameter(
        "Minimum sum value for bound variable in output clusters",
        "A minimum value that the sum value of bounding variable in each cluster should be greater than.",
    )
    weight_mode = knext.StringParameter(
        "Spatial Weight model",
        "Input spatial weight mode.",
        "Queen",
        enum=["Queen", "Rook"],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):

        import pygeoda

        controlVar = self.bound_col
        m_bound = self.minibound
        attributelist = self.bound_col.split(";")

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        geodadf = pygeoda.open(gdf)
        data = geodadf[attributelist]
        b_vals = geodadf.GetRealCol(controlVar)

        if self.weight_mode == "Queen":
            w = pygeoda.queen_weights(geodadf)
        else:
            w = pygeoda.rook_weights(geodadf)
        final_cluster = pygeoda.maxp_greedy(
            w,
            data,
            bound_variable=b_vals,
            min_bound=m_bound,
        )
        gdf["clusterid"] = final_cluster["Clusters"]
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# AZP Greedy
############################################
@knext.node(
    name="AZP-Greedy",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "azp.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class AZPgreedyNode:
    """
    A greedy algorithm to solve the max-p-region problem.

    The so-called max-p regions model (outlined in Duque, Anselin, and Rey 2012)
    uses a different approach and considers the regionalization problem as an application of
    integer programming. In addition, the number of regions is determined endogenously.

    The algorithm itself consists of a search process that starts with an initial feasible
    solution and iteratively improves upon it while maintaining contiguity among the elements of each cluster.

    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    bound_col = knext.ColumnParameter(
        "Bound column for minibound",
        "Select the bound column for clusters with minibound.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    attribute_list = knext.StringParameter(
        "Attribute columns for Clustering",
        "Input the column names with Semicolon.",
        "",
    )
    cluster_k = knext.IntParameter(
        "Number of clusters",
        "The Number of user-defined clusters.",
        default_value=4,
    )
    minibound = knext.DoubleParameter(
        "Minimum sum value for bound variable in output clusters",
        "A minimum value that the sum value of bounding variable in each cluster should be greater than.",
    )
    weight_mode = knext.StringParameter(
        "Spatial Weight model",
        "Input spatial weight mode.",
        "Queen",
        enum=["Queen", "Rook"],
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):

        import pygeoda
        import numpy as np

        k = self.cluster_k
        controlVar = self.bound_col
        m_bound = self.minibound
        attributelist = self.bound_col.split(";")

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        geodadf = pygeoda.open(gdf)
        data = geodadf[attributelist]
        b_vals = np.array(geodadf.GetRealCol(controlVar))

        if self.weight_mode == "Queen":
            w = pygeoda.queen_weights(geodadf)
        else:
            w = pygeoda.rook_weights(geodadf)
        final_cluster = pygeoda.azp_greedy(
            k,
            w,
            data,
            bound_variable=b_vals,
            min_bound=m_bound,
        )
        gdf["clusterid"] = final_cluster["Clusters"]
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# Peano Curve
############################################
@knext.node(
    name="Peano Curve",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "peanocurve.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with peano curve spatal order.",
)
class PeanoCurveNode:
    """
    Spatial Order by Peano Curve.

    Peano curves are space-filling curves first introduced by Italian mathematician Giuseppe Peano (1890).
    A variation and more complicated form is called Hilbert curves because Hilbert (1891) visualized the
    spacefilling idea described in Peano curves and later referred to it as “topological monsters”(Bartholdi and Platzman 1988).

    Peano and Hilbert curves have been used to find all-nearest-neighbors (Chen and Chang 2011) and spatial ordering of
    geographic data (Guo and Gahegan 2006). Conceptually, Peano curves use algorithms to assign spatial orders to
    points in 2D space and map the points onto one-dimensional (1D) space.

    Following the spatial order,  a point in 2D space can be mapped onto the 1D line underneath,
    and the connected line in 2D space (on the right) is the Peano curve. Following the spatial orders along the 1D line,
    spatial clustering can be achieved by classification with many methods.

    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    grid_k = knext.IntParameter(
        "Binary-digit scale",
        "Grid scale defined for spatial filling.",
        default_value=32,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import math
        import pandas as pd
        import numpy as np

        gdf0 = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gp.GeoDataFrame(geometry=gdf0.geometry, crs=gdf0.crs)
        # Scale Coordinates
        gdf["ctrX"] = gdf.centroid.x
        gdf["ctrY"] = gdf.centroid.y
        Xmin, Ymin, Xmax, Ymax = gdf.total_bounds
        dx = Xmax - Xmin
        dy = Ymax - Ymin

        if dx >= dy:
            offsetx = 0.0
            offsety = (1.0 - dy / dx) / 2.0
            scale = dx
        else:
            offsetx = (1.0 - dy / dx) / 2.0
            offsety = 0.0
            scale = dy

        gdf["unitx"] = (gdf["ctrX"] - Xmin) / scale + offsetx
        gdf["unity"] = (gdf["ctrY"] - Ymin) / scale + offsety

        def GetFractionalPart(dbl):
            return dbl - math.floor(dbl)

        # return the peano curve coordinate for a given x,y value
        def Peano(x, y, k):
            if k == 0 or (x == 1 and y == 1):
                return 0.5
            if x <= 0.5:
                if y <= 0.5:
                    quad = 0
                else:
                    quad = 1
            elif y <= 0.5:
                quad = 3
            else:
                quad = 2
            subpos = Peano(2 * abs(x - 0.5), 2 * abs(y - 0.5), k - 1)

            if quad == 1 or quad == 3:
                subpos = 1 - subpos

            return GetFractionalPart((quad + subpos - 0.5) / 4.0)

        gdf["theid"] = np.arange(gdf.shape[0]).tolist()
        gdf["peanoorder"] = 0
        for i in range(gdf.shape[0]):
            x = gdf.loc[(gdf.theid == i), "unitx"].item()
            y = gdf.loc[(gdf.theid == i), "unity"].item()
            gdf.loc[(gdf.theid == i), "peanoorder"] = Peano(x, y, self.grid_k)
        gdf0["peanoorder"] = gdf.peanoorder
        gdf0.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf0)


############################################
#  MSSC
############################################
@knext.node(
    name=" MSSC Initialization",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "mssc.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the input file path.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class MSSCNode:
    """
    Modified scale-space clustering(MSSC) Initialization.

    This node performs the initial MSSC clustering based on the input data.
    The MSSC model (Mu and Wang 2008) follows the values of spatial order along the Peano curve with breaking points that are defined by a threshold population size.
    Through iterations in programming, each cluster satisfies the criterias of ascending spatial order and aggregation volumn minimum constraints.
    MSSC(Mu and Wang 2008) is developed based on scale-space theory, an earlier algorithm and applications of the theory in remote sensing and GIS.
    Using analogies of solid melting and viewing images, scale-space theory treats “scale”—corresponding to temperature in solid melting or distance in viewing
    images—as a parameter in describing the processes and phenomena. With the increase of scale (as temperature in
    the melting algorithm), a piece of metal will melt into liquid but not evenly,
    showing a clustering pattern; with the increase of scale (as distance in the blurring algorithm),
    the same image can reveal different levels of generalizations and details, or different cluster centers.

    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the column containing the geometry data for spatial clustering.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    order_col = knext.ColumnParameter(
        "Weighted order column",
        "Select the order column for MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    constraint_list = knext.StringParameter(
        "Constraints columns names for clustering",
        "Input the column names with Semicolon.",
        "",
    )
    const_capacity_list = knext.StringParameter(
        "Minimum constraint values for clustering",
        "Input the capacity list with Semicolon.",
        "",
    )
    _clusterid = "clusterid"
    _isolated = "isolate"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        result = input_schema.append(
            [
                knext.Column(
                    knext.int64(),
                    knut.get_unique_column_name(self._isolated, input_schema),
                ),
                knext.Column(
                    knext.int64(),
                    knut.get_unique_column_name(self._clusterid, input_schema),
                ),
            ]
        )
        return result

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import libpysal

        ConstraintList = self.constraint_list.split(";")
        CapacityList = [float(x) for x in self.const_capacity_list.split(";")]

        exec_context.flow_variables["constraint_variable"] = self.constraint_list
        exec_context.flow_variables["constraint_capacity"] = self.const_capacity_list

        newlist = ConstraintList + [self.order_col] + [self.geo_col]

        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)

        df = gdf[newlist].copy()
        df["OriginalID"] = list(range(df.shape[0]))
        df["TheID"] = df[self.order_col].rank().astype(int) - 1
        df = df.set_index("TheID", drop=False).sort_index().rename_axis(None)

        # Create spatial weight matrix
        wq = libpysal.weights.Rook.from_dataframe(df)
        w = wq.neighbors

        df["isolate"] = 0
        df["included"] = 0
        tmp_df = df

        class_value = 1
        count = 0
        classValueDict = {class_value: []}

        tot_rec = df.shape[0]
        roundCount = 0
        newRoundCount = tot_rec
        orphanIteration = 0

        acc_capacity = []
        for i in range(len(ConstraintList)):
            acc_capacity.append(0)

        # Major round of clustering starts

        while tmp_df.shape[0] > 0:
            tmp_df = tmp_df[tmp_df["included"] == 0]
            # print('dfFC:',dfFC.shape[0])
            for indexs, row in tmp_df.iterrows():
                index = int(row["TheID"])
                count += 1
                roundCount += 1
                cur_capacity = []
                for i in range(len(ConstraintList)):
                    cur_capacity.append(row[ConstraintList[i]])

                if count == 1:
                    tmp_df.loc[indexs, "included"] = 1
                    df.loc[indexs, "included"] = 1

                if count > 1:
                    if len(classValueDict[class_value]) == 0:
                        tmp_df.loc[indexs, "included"] = 1
                        df.loc[indexs, "included"] = 1
                    else:
                        for NID in w[index]:
                            if NID in classValueDict[class_value]:
                                tmp_df.loc[indexs, "included"] = 1
                                df.loc[indexs, "included"] = 1
                                break

                if tmp_df.loc[indexs, "included"] == 1:
                    for i in range(len(ConstraintList)):
                        acc_capacity[i] += cur_capacity[i]
                    classValueDict[class_value].append(index)

                    satisfyAll = 1
                    for i in range(len(ConstraintList)):
                        if float(acc_capacity[i]) >= float(CapacityList[i]):
                            satisfyAll = satisfyAll * 1
                        else:
                            satisfyAll = satisfyAll * 0

                    if satisfyAll == 1 and count < tot_rec:
                        for i in range(len(ConstraintList)):
                            acc_capacity[i] = 0

                        class_value += 1
                        classValueDict[class_value] = []

                    roundCount = 0
                    tmp_df = tmp_df[tmp_df["included"] == 0]
                    # tmp_df = dfFC.sort_values('OrdVal', ascending=True)
                    newRoundCount = tmp_df.shape[0]
                    break

                elif tmp_df.loc[indexs, "included"] == 0:
                    count -= 1

                if roundCount == newRoundCount:
                    dforphan = tmp_df[tmp_df["included"] == 0]
                    orphanCount = dforphan.shape[0]
                    if orphanCount > 0:
                        orphanIteration += 1
                        if len(classValueDict[class_value]) > 0:
                            for i in range(len(ConstraintList)):
                                acc_capacity[i] = 0
                            class_value += 1
                            classValueDict[class_value] = []

                        tmp_df = tmp_df[tmp_df["included"] == 0]
                        newRoundCount = tmp_df.shape[0]
                        roundCount = 0
                        break

        df["SubClusID"] = 0
        for k in classValueDict.keys():
            df.loc[(df.TheID.isin(classValueDict[k])), "SubClusID"] = k
        df = df.set_index("OriginalID", drop=False).sort_index()

        gdf[self._isolated] = df.isolate.tolist()
        gdf[self._clusterid] = df.SubClusID.tolist()
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
#  MSSC Refiner
############################################
@knext.node(
    name=" MSSC Refiner",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "msscmerge.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the MSSC clustering.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with modified cluster tag.",
)
class MSSCmodifierNode:
    """
    Modified scale-space clustering(MSSC) Refiner.

    This node takes the clusters generated by the MSSC Initialization node and applies additional
    steps to refine those clusters that do not meet the certain criteria.
    The method combinies spatial weight matrix with spatial orders to address the jumping problem— disconnected members in a cluster.

    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    clusterid_col = knext.ColumnParameter(
        "Cluster id column from MSSC",
        "Select the cluster id column generated by MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    isolate_col = knext.ColumnParameter(
        "isolate column from MSSC",
        "Select the column 'isoloate'generated by MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    constraint_list = knext.StringParameter(
        "Constraints columns names for clustering",
        "Input the column names with Semicolon.",
        "",
    )
    const_capacity_list = knext.StringParameter(
        "Minimum constraint values for clustering",
        "Input the capacity list with Semicolon.",
        "",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return input_schema

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import libpysal
        import pandas as pd

        ConstraintList = self.constraint_list.split(";")
        CapacityList = [float(x) for x in self.const_capacity_list.split(";")]

        SubClusID = self.clusterid_col
        TheID = "TheID"
        isolateid = self.isolate_col
        keycolum = ConstraintList + [SubClusID] + [isolateid] + [self.geo_col]
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        df = gdf[keycolum].copy()

        df = df.rename(columns={isolateid: "isolate", SubClusID: "SubClusID"})
        df["TheID"] = list(range(df.shape[0]))
        df = df.set_index("TheID", drop=False).rename_axis(None)

        keycolum1 = ["SubClusID"] + ConstraintList
        sum_tmp = df[keycolum1].groupby("SubClusID").sum().reset_index()
        classValueDict = df.groupby("SubClusID")["TheID"].apply(list).to_dict()

        # Create spatial weight matrix
        wq = libpysal.weights.Rook.from_dataframe(df)
        w = wq.neighbors

        df_list = []
        # Loop through the constraint and capacity lists
        for j in range(len(ConstraintList)):
            df_filtered = sum_tmp[sum_tmp[ConstraintList[j]] < CapacityList[j]]
            df_list.append(df_filtered)
        df2 = pd.concat(df_list).drop_duplicates()

        if df2.shape[0] > 0:
            clusAdjList = []
            isoClus = []
            rowCount = 0
            clusAdjDict = {}
            class_value = max(df2["SubClusID"])
            for index, rowLocal in df2.iterrows():
                rowCount += 1
                FoundIt = False
                foundCount = 0
                currentClus = rowLocal["SubClusID"]
                theCluster = currentClus

                capA = []
                for j in range(len(ConstraintList)):
                    capA.append(rowLocal[ConstraintList[j]])

                swmRows = []
                for ids in classValueDict[currentClus]:
                    swmRows.extend(w[ids])
                swmRows = list(set(swmRows))

                matching_indices = set()
                for swmRow in swmRows:
                    for i in range(class_value, 0, -1):
                        if i != currentClus:
                            if swmRow in classValueDict[i]:
                                matching_indices.add(i)

                for i in matching_indices:
                    foundRow = sum_tmp[sum_tmp.SubClusID == i]
                    capB = []
                    for j in range(len(ConstraintList)):
                        capB.append(foundRow[ConstraintList[j]])
                    newcap = []
                    for j in range(len(ConstraintList)):
                        newcap.append(capA[j] + capB[j])

                    newSatisfyAll = 1
                    for j in range(len(ConstraintList)):
                        if float(newcap[j]) >= float(CapacityList[j]):
                            newSatisfyAll = newSatisfyAll * 1
                        else:
                            newSatisfyAll = newSatisfyAll * 0

                    # if newcap1 >= capacity1 and newcap2 >= capacity2:
                    if newSatisfyAll == 1:
                        foundCount += 1
                        if foundCount == 1:
                            mincap = []
                            for k in range(len(ConstraintList)):
                                mincap.append(newcap[k])
                            FoundIt = True
                            theCluster = i
                            break
                        # elif foundCount > 1 and newcap1 < mincap1 and newcap2 < mincap2:
                        elif foundCount > 1:
                            minSatisfyAll = 1
                            for j in range(len(ConstraintList)):
                                if float(newcap[j]) < float(mincap[j]):
                                    minSatisfyAll = minSatisfyAll * 1
                                else:
                                    minSatisfyAll = minSatisfyAll * 0

                            if minSatisfyAll == 1:
                                for k in range(len(ConstraintList)):
                                    mincap[k] = newcap[k]
                                FoundIt = True
                                theCluster = i
                                break

                if FoundIt == True:
                    clusAdjDict[currentClus] = theCluster
                    df.loc[(df.SubClusID == currentClus), "SubClusID"] = theCluster
                    df.loc[(df.SubClusID == currentClus), "isolate"] = 0
                else:
                    df.loc[(df.SubClusID == currentClus), "isolate"] = 1
                    isoClus.append(currentClus)

                clusAdjList.append(FoundIt)

            clusAdjust = 0
            listInd = 0

            for indexs, rowLocal in df2.iterrows():
                if len(clusAdjList) > 0:
                    if clusAdjList[listInd] == True:
                        index = rowLocal["SubClusID"]
                        df.loc[(df.SubClusID > (index - clusAdjust)), "SubClusID"] -= 1
                        clusAdjust += 1
                        listInd += 1

        # df = df.set_index("TheID", drop=False).sort_index()
        gdf[SubClusID] = df.SubClusID.tolist()
        gdf[isolateid] = df.isolate.tolist()

        # gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
#  MSSC Isolation Tackler
############################################
@knext.node(
    name="Isolation Tackler",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "msscisolate.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata from the clustering.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with modified cluster tag.",
)
class MSSCisolationNode:
    """
    Modified scale-space clustering(MSSC).

    The node forces a cluster membership for unclaimed or loose-end units, especially useful for the enclave units in Mixed level regionalization model.

    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column for Geodata.",
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    clusterid_col = knext.ColumnParameter(
        "Cluster id column from MSSC",
        "Select the cluster id column generated by MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    isolate_col = knext.ColumnParameter(
        "isolate column from MSSC",
        "Select the column 'isoloate'generated by MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    constraint_list = knext.StringParameter(
        "Constraints columns names for clustering",
        "Input the column names with Semicolon.",
        "",
    )
    const_capacity_list = knext.StringParameter(
        "Minimum constraint values for clustering",
        "Input the capacity list with Semicolon.",
        "",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return input_schema

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import libpysal
        import pandas as pd

        ConstraintList = self.constraint_list.split(";")
        CapacityList = [float(x) for x in self.const_capacity_list.split(";")]

        SubClusID = self.clusterid_col
        isolateid = self.isolate_col
        keycolum = ConstraintList + [SubClusID] + [isolateid] + [self.geo_col]
        gdf0 = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf = gdf0[keycolum].copy()

        gdf = gdf.rename(
            columns={
                isolateid: "isolate",
                SubClusID: "SubClusID",
                self.geo_col: "geometry",
            }
        )
        gdf["FinalClus"] = gdf["SubClusID"]

        # Create dictionary of Constraint list
        mixStatFlds = {}
        for i in ConstraintList:
            mixStatFlds[i] = "sum"
        mixStatFlds["SubClusID"] = "first"
        mixStatFlds["isolate"] = "min"

        # Dissolve gdf based on  Constraint list
        tmpMixedClusFC = (
            gdf.dissolve(by="FinalClus", aggfunc=mixStatFlds)
            .reset_index()
            .rename({"index": "FinalClus"})
        )

        def nearest(
            row,
            geom_union,
            df1,
            df2,
            geom1_col="geometry",
            geom2_col="geometry",
            src_column=None,
        ):
            """Find the nearest point and return the corresponding value from specified column."""
            from shapely.ops import nearest_points

            # Find the geometry that is closest
            nearest = df2[geom2_col] == nearest_points(row[geom1_col], geom_union)[1]
            # Get the corresponding value from df2 (matching is based on the geometry)
            value = df2[nearest][src_column].to_numpy()[0]
            return value

        def MergeIsolated(tmpMixedClusFC, theClus):
            # get isolate layer out
            lyr5 = tmpMixedClusFC[tmpMixedClusFC.index == theClus]
            theClusID = lyr5.FinalClus.to_list()[0]
            # get satified layer out
            lyr6 = tmpMixedClusFC[tmpMixedClusFC[isolateid] == 0]
            # get weight
            wq = libpysal.weights.Rook.from_dataframe(tmpMixedClusFC)
            w_mix = wq.neighbors

            lyr6select = lyr6[lyr6.index.isin(w_mix[theClus])]
            count = lyr6select.shape[0]
            if count > 0:
                # get the smaller group
                minClus = lyr6select[ConstraintList[0]].idxmin()
                theNewClusID = minClus
            else:
                lyr6["centroid"] = lyr6.centroid
                lyr5["centroid"] = lyr5.centroid
                unary_union = lyr6.centroid.unary_union
                lyr5["nearest_id"] = lyr5.apply(
                    nearest,
                    geom_union=unary_union,
                    df1=lyr5,
                    df2=lyr6,
                    geom1_col="centroid",
                    geom2_col="centroid",
                    src_column="SubClusID",
                    axis=1,
                )
                theNewClusID = lyr5["nearest_id"].values[0]
            return theClusID, theNewClusID

        lyrIso = tmpMixedClusFC[tmpMixedClusFC["isolate"] > 0]
        fldIso = "isolate"
        isoCount1 = lyrIso.shape[0]
        isoClusters = []

        for index, row in lyrIso.iterrows():
            isoClusters.append(index)
        for i in range(isoCount1):
            theClusID, theNewClusID = MergeIsolated(tmpMixedClusFC, isoClusters[i])
            gdf.loc[(gdf["SubClusID"] == theClusID), fldIso] = 0
            gdf.loc[(gdf["SubClusID"] == theClusID), "SubClusID"] = theNewClusID

        gdf0[SubClusID] = gdf.SubClusID.tolist()
        gdf0[isolateid] = gdf.isolate.tolist()

        gdf0.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf0)
