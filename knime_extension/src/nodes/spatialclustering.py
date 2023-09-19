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
__NODE_ICON_PATH = "icons/icon/SpatialClustering/"
_CLUSTER_ID = "Cluster ID"
_ISOLATED = "Isolate"


def validate_k(
    value: int, min_val: int = 1, msg: str = "Number of clusters must be greater 1."
) -> int:
    """
    Checks if the cluster k value is larger than 1.
    """
    if value <= min_val:
        raise knext.InvalidParametersError(msg)
    return int(value)


# def  geodataframe to PPP
def gdf2ppp(gdf):
    from pointpats import PointPattern

    ngdf = gp.GeoDataFrame(geometry=gdf.geometry, crs=gdf.crs)
    ngdf["coords"] = ngdf.apply(lambda x: list(x.geometry.coords[0]), axis=1)
    ppts = ngdf.coords.to_list()
    pp = PointPattern(ppts)
    return pp


def get_cluster_k():
    return knext.IntParameter(
        "Number of clusters",
        "The number of user-defined clusters.",
        default_value=4,
        validator=validate_k,
    )


def get_seed():
    return knext.IntParameter(
        "Seed",
        "The seed for the random number generator.",
        default_value=123456789,
        is_advanced=True,
    )


@knext.parameter_group(label="Cluster settings")
class ClusterSettings:
    class WeightModel(knext.EnumParameterOptions):
        QUEEN = (
            "Queen",
            "Queen contiguity weights.",
        )
        ROOK = ("Rook", "Rook contiguity weights.")

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to implement spatial clustering."
    )

    # expose seed

    bound_col = knext.ColumnParameter(
        "Bound column for minibound",
        "Select the bound column for clusters with minibound.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    # use column filter

    attribute_list = knext.ColumnFilterParameter(
        "Attribute columns for clustering",
        "Select columns for calculating attribute distance.",
        port_index=0,
        column_filter=knut.is_numeric,
    )
    minibound = knext.DoubleParameter(
        "Minimum total value for the bounding variable in each output cluster",
        "The sum of the bounding variable in each cluster must be greater than this minimum value.",
    )
    weight_mode = knext.EnumParameter(
        label="Spatial weight model",
        description="Input spatial weight mode.",
        default_value=WeightModel.QUEEN.name,
        enum=WeightModel,
    )

    def __weight(self, geoda_df):
        import pygeoda

        if self.weight_mode == self.WeightModel.QUEEN.name:
            w = pygeoda.queen_weights(geoda_df)
        else:
            w = pygeoda.rook_weights(geoda_df)
        return w

    def do_configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return input_schema.append(knext.Column(knext.int64(), name=_CLUSTER_ID))

    def do_clustering(self, exec_context: knext.ExecutionContext, input_1, fnc):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        import pygeoda

        geoda_df = pygeoda.open(gdf)
        attributelist = self.attribute_list.apply(input_1.schema).column_names
        data = geoda_df[attributelist]
        b_vals = geoda_df.GetRealCol(self.bound_col)
        w = self.__weight(geoda_df)
        m_bound = self.minibound

        knut.check_canceled(exec_context)
        cluster = fnc(w, data, b_vals, m_bound)
        gdf[_CLUSTER_ID] = cluster["Clusters"]
        gdf.reset_index(drop=True, inplace=True)
        return knut.to_table(gdf, exec_context)


############################################
# MeanCenter and StandardDistance
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
    description="Table with geometry column for center computation.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with mean center and standard distance.",
)
class MeanCenterNode:
    """
    Mean center and standard distance.

    This node uses the Pysal package [pointpats](http://pysal.org/pointpats/)
    to  measure the compactness of a spatial distribution of features around its mean center.
    Standard distance (or standard distance deviation) is usually represented as a circle where the radius of the circle is the standard distance.
    """

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute mean center and standard distance."
    )

    weight_col = knext.ColumnParameter(
        "Weight column",
        "Select the weighting column to calculate the mean center (Optional).",
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

        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
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
            {"X": [mc[0]], "Y": [mc[1]]},
            geometry=gp.GeoSeries.from_xy(x=[mc[0]], y=[mc[1]], crs=gdf.crs),
        )
        # Standard distance
        mcgdf["Distance"] = pointpats.std_distance(pp.points)
        mcgdf.reset_index(drop=True, inplace=True)
        return knut.to_table(mcgdf, exec_context)


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
    description="Table with geometry column for computing standard deviational ellipse.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with Standard deviational ellipse.",
)
class StandardEllipseNode:
    """
    Standard deviational ellipse.

    This node uses the Pysal package [pointpats](http://pysal.org/pointpats/)
    to  calculate parameters of standard deviational ellipse for a point pattern, which is a common way of measuring
    the trend for a set of points or areas. These measures define the axes of an ellipse (or ellipsoid) encompassing the distribution of features.
    """

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute standard deviational ellipse."
    )
    _COL_GEOMETRY = "geometry"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return knext.Schema.from_columns(
            [
                knext.Column(knut.TYPE_POLYGON, self._COL_GEOMETRY),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pointpats
        import numpy as np
        from matplotlib.patches import Ellipse
        from shapely.geometry import Polygon

        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
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
        return knut.to_table(gdf, exec_context)


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
    description="Geodata for spatial clustering implementation.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
@knut.geoda_node_description(
    short_description="Spatial C(K)luster Analysis by Tree Edge Removal(SKATER).",
    description="""The Spatial C(K)luster Analysis by Tree Edge Removal (SKATER). The Spatial C(K)luster Analysis by 
    Tree Edge Removal (SKATER) algorithm introduced by Assuncao et al. (2006)
    is based on the optimal pruning of a minimum spanning tree that reflects the contiguity structure among the observations.
    It provides an optimized algorithm to prune to tree into several clusters that their values of selected variables 
    are as similar as possible. 
    """,
    references={
        "Spatial Clustering": "https://geodacenter.github.io/pygeoda/spatial_clustering.html",
        "SKATER": "https://geodacenter.github.io/pygeoda/_modules/pygeoda/clustering/skater.html",
        "Spatially Constrained Clustering - Hierarchical Methods": "https://geodacenter.github.io/workbook/9c_spatial3/lab9c.html",
    },
)
class SKATERNode:
    cluster_settings = ClusterSettings()

    cluster_k = get_cluster_k()

    seed = get_seed()

    def configure(self, configure_context, input_schema):
        return self.cluster_settings.do_configure(configure_context, input_schema)

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pygeoda

        fnc = lambda w, data, b_vals, m_bound: pygeoda.skater(
            self.cluster_k,
            w,
            data,
            distance_method="euclidean",
            bound_variable=b_vals,
            min_bound=m_bound,
            random_seed=self.seed,
        )
        return self.cluster_settings.do_clustering(exec_context, input_1, fnc)


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
    description="Geodata for spatial clustering implementation.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
@knut.geoda_node_description(
    short_description="REDCAP (Regionalization with dynamically constrained agglomerative clustering and partitioning).",
    description="""It is developed by [D. Guo (2008).](https://doi.org/10.1080/13658810701674970) Like SKATER, REDCAP 
    starts from building a spanning tree with 4 different ways (single-linkage, average-linkage, 
    complete-linkage and wards-linkage). Then, REDCAP provides 2 different ways (first‐order and full-order constraining) 
    to prune the tree to find clusters. The first-order approach with a minimum spanning tree is exactly
    the same with SKATER. In GeoDa and pygeoda, the following methods are provided:
    
- First-order and Single-linkage: In this local approach, clusters are formed by considering only immediate neighbors, 
    and their distance is measured by the shortest distance between any pair of points from each cluster.
- Full-order and Complete-linkage: This method also considers all points in the dataset for clustering and calculates 
    the distance between clusters as the average distance between all pairs of points, one from each cluster.
- Full-order and Average-linkage: This method also considers all points in the dataset for clustering and calculates 
    the distance between clusters as the average distance between all pairs of points, one from each cluster.
- Full-order and Single-linkage: Using a global context, this approach considers all points in the dataset for clustering 
    and measures the distance between clusters as the shortest distance between any pair of points from each cluster.
- Full-order and Wards-linkage: All points in the dataset are considered for clustering, and the distance 
    between clusters is calculated in a way that minimizes the internal variance within each cluster.

    """,
    references={
        "Spatial Clustering": "https://geodacenter.github.io/pygeoda/spatial_clustering.html",
        "REDCAP": "https://geodacenter.github.io/pygeoda/_modules/pygeoda/clustering/redcap.html",
        "Spatially Constrained Clustering - Hierarchical Methods ": "https://geodacenter.github.io/workbook/9c_spatial3/lab9c.html",
    },
)
class REDCAPNode:
    class LinkageModes(knext.EnumParameterOptions):
        FIRSTORDER_SINGLELINKAGE = (
            "First-order and Single-linkage",
            """In this local approach, clusters are formed by considering only immediate neighbors, 
and their distance is measured by the shortest distance between any pair of points from each cluster.""",
        )
        FULLORDER_SINGLELINKAGE = (
            "Full-order and Single-linkage",
            """Using a global context, this approach considers all points in the dataset for clustering 
and measures the distance between clusters as the shortest distance between any pair of points from each cluster.""",
        )
        FULLORDER_COMPLETELINKAGE = (
            "Full-order and Complete-linkage",
            """This method also considers all points in the dataset for clustering and calculates 
the distance between clusters as the average distance between all pairs of points, one from each cluster.""",
        )
        FULLORDER_AVERAGELINKAGE = (
            "Full-order and Average-linkage",
            """This method also considers all points in the dataset for clustering and calculates 
the distance between clusters as the average distance between all pairs of points, one from each cluster.""",
        )
        FULLORDER_WARDLINKAGE = (
            "Full-order and Wards-linkage",
            """All points in the dataset are considered for clustering, and the distance 
between clusters is calculated in a way that minimizes the internal variance within each cluster.""",
        )

        @classmethod
        def get_default(cls):
            return cls.FULLORDER_COMPLETELINKAGE

    cluster_settings = ClusterSettings()

    cluster_k = get_cluster_k()

    seed = get_seed()

    link_mode = knext.EnumParameter(
        label="Linkage mode",
        description="Input linkage mode.",
        default_value=LinkageModes.get_default().name,
        enum=LinkageModes,
    )

    def configure(self, configure_context, input_schema):
        return self.cluster_settings.do_configure(configure_context, input_schema)

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        linkage = self.link_mode.lower().replace("_", "-")
        import pygeoda

        fnc = lambda w, data, b_vals, m_bound: pygeoda.redcap(
            self.cluster_k,
            w,
            data,
            linkage,
            bound_variable=b_vals,
            min_bound=m_bound,
            random_seed=self.seed,
        )
        return self.cluster_settings.do_clustering(exec_context, input_1, fnc)


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
    description="Geodata for spatial clustering implementation.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
@knut.geoda_node_description(
    short_description="SCHC (Spatially constrained hierarchical clustering).",
    description="""SCHC (Spatially constrained hierarchical clustering). It is a special form of constrained clustering, 
    where the constraint is based on contiguity (common borders). The method builds up the clusters using agglomerative 
    hierarchical clustering methods: single linkage, complete linkage, average linkage, and Ward’s method 
    (a special form of centroid linkage). Meanwhile, it also maintains the spatial contiguity when merging two clusters. 
    The method builds up the clusters using agglomerative hierarchical clustering methods:

- Single
- Complete
- Average
- Ward
    """,
    references={
        "Spatial Clustering": "https://geodacenter.github.io/pygeoda/spatial_clustering.html",
        "SCHC": "https://geodacenter.github.io/pygeoda/_modules/pygeoda/clustering/schc.html",
        "Spatially Constrained Clustering - Hierarchical Methods ": "https://geodacenter.github.io/workbook/9c_spatial3/lab9c.html",
    },
)
class SCHCNode:
    class LinkageModes(knext.EnumParameterOptions):
        SINGLE = (
            "Single linkage",
            """Forms clusters by linking geographical units based on the shortest distance between them, while maintaining spatial contiguity.""",
        )
        COMPLETE = (
            "Complete linkage",
            """Clusters geographical units based on the farthest distance within pairs, ensuring that each resulting cluster shares common borders.""",
        )
        AVERAGE = (
            "Average linkage",
            """Groups units into clusters by calculating the average distance among all possible pairs, while also preserving spatial contiguity.""",
        )
        WARD = (
            "Ward linkage",
            """Minimizes within-cluster variance by clustering units that have similar centroids, 
            all while maintaining common borders between units in each cluster.""",
        )

        @classmethod
        def get_default(cls):
            return cls.COMPLETE

    cluster_settings = ClusterSettings()

    cluster_k = get_cluster_k()

    link_mode = knext.EnumParameter(
        label="Linkage mode",
        description="Input linkage mode.",
        default_value=LinkageModes.get_default().name,
        enum=LinkageModes,
    )

    def configure(self, configure_context, input_schema):
        return self.cluster_settings.do_configure(configure_context, input_schema)

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        linkage = self.link_mode.lower()

        import pygeoda

        fnc = lambda w, data, b_vals, m_bound: pygeoda.schc(
            self.cluster_k,
            w,
            data,
            linkage,
            bound_variable=b_vals,
            min_bound=m_bound,
        )
        return self.cluster_settings.do_clustering(exec_context, input_1, fnc)


############################################
# MaxP Greedy
############################################
@knext.node(
    name="MaxP",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "maxp.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata for spatial clustering implementation.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
@knut.geoda_node_description(
    short_description="Max-P Greedy.",
    description="""Max-P Greedy. A greedy algorithm to solve the max-p-region problem.

The so-called max-p regions model (outlined in [Duque, Anselin, and Rey 2012](https://doi.org/10.1111/j.1467-9787.2011.00743.x))
uses a different approach and considers the regionalization problem as an application of
integer programming. In addition, the number of regions is determined endogenously.

The algorithm itself consists of a search process that starts with an initial feasible
solution and iteratively improves upon it while maintaining contiguity among the elements of each cluster.
    """,
    references={
        "Spatial Clustering": "https://geodacenter.github.io/pygeoda/spatial_clustering.html",
        "Max-P": "https://geodacenter.github.io/pygeoda/_modules/pygeoda/clustering/maxp.html",
        "Spatially Constrained Clustering - Partitioning Methods ": "https://geodacenter.github.io/workbook/9d_spatial4/lab9d.html",
    },
)
class MaxPgreedyNode:
    cluster_settings = ClusterSettings()

    seed = get_seed()

    def configure(self, configure_context, input_schema):
        return self.cluster_settings.do_configure(configure_context, input_schema)

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pygeoda

        fnc = lambda w, data, b_vals, m_bound: pygeoda.maxp_greedy(
            w, data, bound_variable=b_vals, min_bound=m_bound, random_seed=self.seed
        )
        return self.cluster_settings.do_clustering(exec_context, input_1, fnc)


############################################
# AZP Greedy
############################################
@knext.node(
    name="AZP",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "azp.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geodata table",
    description="Geodata for spatial clustering implementation.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
@knut.geoda_node_description(
    short_description="AZP Greedy.",
    description="""AZP Greedy. A greedy algorithm for automatic zoning procedure (AZP).

The automatic zoning procedure (AZP) was initially outlined in [Openshaw (1977)](https://doi.org/10.2307/622300) to 
address some of the consequences of the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic 
to find the best set of combinations of contiguous spatial units into p regions, minimizing the within sum of squares 
as a criterion of homogeneity. The number of regions needs to be specified beforehand, as in most other clustering 
methods considered so far.

    """,
    references={
        "Spatial Clustering": "https://geodacenter.github.io/pygeoda/spatial_clustering.html",
        "AZP": "https://geodacenter.github.io/pygeoda/_modules/pygeoda/clustering/azp.html",
        "Spatially Constrained Clustering - Partitioning Methods ": "https://geodacenter.github.io/workbook/9d_spatial4/lab9d.html",
    },
)
class AZPgreedyNode:
    cluster_settings = ClusterSettings()

    cluster_k = get_cluster_k()

    seed = get_seed()

    def configure(self, configure_context, input_schema):
        return self.cluster_settings.do_configure(configure_context, input_schema)

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pygeoda

        fnc = lambda w, data, b_vals, m_bound: pygeoda.azp_greedy(
            self.cluster_k,
            w,
            data,
            bound_variable=b_vals,
            min_bound=m_bound,
            random_seed=self.seed,
        )
        return self.cluster_settings.do_clustering(exec_context, input_1, fnc)


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
    description="Geodata for calculating Peano curve spatial order.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with Peano curve spatial order.",
)
class PeanoCurveNode:
    """
    Spatial Order by Peano Curve.

    Peano curves are space-filling curves first introduced by Italian mathematician Giuseppe Peano (1890).
    A variation and more complicated form is called Hilbert curves because Hilbert (1891) visualized the
    space filling idea described in Peano curves and later referred to it as “topological monsters” (Bartholdi and Platzman 1988).

    Peano and Hilbert curves have been used to find all-nearest-neighbors (Chen and Chang 2011) and spatial ordering of
    geographic data (Guo and Gahegan 2006). Conceptually, Peano curves use algorithms to assign spatial orders to
    points in 2D space and map the points onto one-dimensional (1D) space.

    Following the spatial order,  a point in 2D space can be mapped onto the 1D line underneath,
    and the connected line in 2D space (on the right) is the Peano curve. Following the spatial orders along the 1D line,
    spatial clustering can be achieved by classification with many methods.

    The value of the binary-digit scale (n) will be calculated as the grid scale value, 2^n.
    This value will then be used to generate a 2^n x 2^n grid.

    """

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to calculate spatial order."
    )

    grid_k = knext.IntParameter(
        "Binary-digit scale as powers of 2",
        "Grid scale of 2^n defined for spatial filling.",
        default_value=5,
        validator=lambda k: validate_k(
            k, 2, "Binary-digit scale must be greater than 2"
        ),
    )

    _PEANO_CURVE_ORDER = "Peano Order"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return input_schema.append(
            knext.Column(knext.double(), name=self._PEANO_CURVE_ORDER)
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import math
        import pandas as pd
        import numpy as np

        k = 2 ^ self.grid_k

        gdf0 = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        gdf = gp.GeoDataFrame(geometry=gdf0.geometry, crs=gdf0.crs)
        # Scale Coordinates
        gdf["ctrX"] = gdf.centroid.x
        gdf["ctrY"] = gdf.centroid.y
        Xmin, Ymin, Xmax, Ymax = gdf.total_bounds
        dx = Xmax - Xmin
        dy = Ymax - Ymin

        knut.check_canceled(exec_context)

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
        row_count = gdf.shape[0]
        for i in range(gdf.shape[0]):
            x = gdf.loc[(gdf.theid == i), "unitx"].item()
            y = gdf.loc[(gdf.theid == i), "unity"].item()
            gdf.loc[(gdf.theid == i), "peanoorder"] = Peano(x, y, self.grid_k)
            exec_context.set_progress(
                i / row_count, f"Processing row {i} of {row_count}"
            )
            knut.check_canceled(exec_context)

        gdf0[self._PEANO_CURVE_ORDER] = gdf.peanoorder
        gdf0.reset_index(drop=True, inplace=True)
        return knut.to_table(gdf0, exec_context)


@knext.parameter_group(label="Constraints")
class MSSCConstraints:
    columns = knext.ColumnFilterParameter(
        "Constraint columns for clustering",
        "The constraints columns used for the clustering.",
        port_index=0,
    )

    capacities = knext.StringParameter(
        "Minimum constraint values for clustering",
        "Input the capacity list with one value for each selected constrained column separated by semicolon.",
        "",
    )

    def get_capacities(self):
        if (
            not self.capacities or self.capacities.isspace()
        ):  # return empty array if capacities string is empty
            return []
        return [float(x) for x in self.capacities.split(";")]

    def get_columns(self, schema):
        return self.columns.apply(schema).column_names

    def validate_settings(self, schema):
        cols = len(self.get_columns(schema))
        caps = len(self.get_capacities())
        if cols != caps:
            raise knext.InvalidParametersError(
                f"Please provide a minimum constraint value for each selected constrained column ({caps} vs {cols})."
            )


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
    description="Geodata for implementing modified scale-space clustering.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with cluster tag.",
)
class MSSCNode:
    """
    Modified scale-space clustering(MSSC) Initialization.

    This node performs the initial MSSC clustering based on the input data.
    The MSSC model [(Mu and Wang 2008)](https://doi.org/10.1080/00045608.2014.968910) follows the values of spatial order along the Peano curve with breaking points
    that are defined by a threshold population size. Through iterations in programming, each cluster satisfies the
    criteria of ascending spatial order and aggregation volume minimum constraints. MSSC (Mu and Wang 2008) is
    developed based on scale-space theory, an earlier algorithm, and applications of the theory in remote sensing and GIS.
    Using analogies of solid melting and viewing images, scale-space theory treats “scale”—corresponding to temperature
    in solid melting or distance in viewing images—as a parameter in describing the processes and phenomena.
    With the increase of scale (as temperature in the melting algorithm), a piece of metal will melt into liquid
    but not evenly, showing a clustering pattern; with the increase of scale (as distance in the blurring algorithm),
    the same image can reveal different levels of generalizations and details, or different cluster centers.

    """

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to implement spatial clustering."
    )

    order_col = knext.ColumnParameter(
        "Weighted order column",
        "Select the order column for MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    constraints = MSSCConstraints()

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        self.constraints.validate_settings(input_schema)
        result = input_schema.append(
            [
                knext.Column(
                    knext.int64(),
                    knut.get_unique_column_name(_ISOLATED, input_schema),
                ),
                knext.Column(
                    knext.int64(),
                    knut.get_unique_column_name(_CLUSTER_ID, input_schema),
                ),
            ]
        )
        return result

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import libpysal

        constraintList = self.constraints.get_columns(input_1.schema)
        capacityList = self.constraints.get_capacities()

        newlist = constraintList + [self.order_col] + [self.geo_col]

        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)

        df = gdf[newlist].copy()
        df["OriginalID"] = list(range(df.shape[0]))
        df["theid"] = df[self.order_col].rank(method="first").astype(int) - 1
        df = df.set_index("theid", drop=False).sort_index().rename_axis(None)

        # Create spatial weight matrix
        wq = libpysal.weights.Rook.from_dataframe(df)
        w = wq.neighbors

        df["isolate"] = 0
        df["included"] = 0
        tmp_df = df.copy()

        class_value = 1
        count = 0
        classValueDict = {class_value: []}

        tot_rec = df.shape[0]
        roundCount = 0
        newRoundCount = tot_rec
        orphanIteration = 0

        knut.check_canceled(exec_context)

        acc_capacity = []
        for i in range(len(constraintList)):
            acc_capacity.append(0)

        # Major round of clustering starts

        while tmp_df.shape[0] > 0:
            tmp_df = tmp_df[tmp_df["included"] == 0]
            for indexs, row in tmp_df.iterrows():
                n_loop = tmp_df.shape[0]
                index = int(row["theid"])
                count += 1
                roundCount += 1
                exec_context.set_progress(
                    0.8 * roundCount / n_loop,
                    f"Batch {roundCount} of {n_loop} processed",
                )
                knut.check_canceled(exec_context)
                cur_capacity = []
                for i in range(len(constraintList)):
                    cur_capacity.append(row[constraintList[i]])

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
                    for i in range(len(constraintList)):
                        acc_capacity[i] += cur_capacity[i]
                    classValueDict[class_value].append(index)

                    satisfyAll = 1
                    for i in range(len(constraintList)):
                        if float(acc_capacity[i]) >= float(capacityList[i]):
                            satisfyAll = satisfyAll * 1
                        else:
                            satisfyAll = satisfyAll * 0

                    if satisfyAll == 1 and count < tot_rec:
                        for i in range(len(constraintList)):
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
                            for i in range(len(constraintList)):
                                acc_capacity[i] = 0
                            class_value += 1
                            classValueDict[class_value] = []

                        tmp_df = tmp_df[tmp_df["included"] == 0]
                        newRoundCount = tmp_df.shape[0]
                        roundCount = 0
                        break

        df["subclusid"] = 0
        for k in classValueDict.keys():
            df.loc[(df.theid.isin(classValueDict[k])), "subclusid"] = k
        df = df.set_index("OriginalID", drop=False).sort_index()

        gdf[_ISOLATED] = df.isolate.tolist()
        gdf[_CLUSTER_ID] = df.subclusid.tolist()
        gdf.reset_index(drop=True, inplace=True)
        return knut.to_table(gdf, exec_context)


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
    description="Geodata from the MSSC Initialization.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with modified cluster tag.",
)
class MSSCmodifierNode:
    """
    Modified scale-space clustering(MSSC) Refiner.

    This node takes the clusters generated by the MSSC Initialization node and applies additional
    steps to refine those clusters that do not do not meet certain criteria.
    The method combines spatial weight matrix with spatial orders to address the jumping problem— disconnected members
    in a cluster.

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
        "Isolate column from MSSC",
        "Select the column 'isolate' generated by MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    constraints = MSSCConstraints()

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        self.clusterid_col = knut.column_exists_or_preset(
            configure_context, _CLUSTER_ID, input_schema, knut.is_long
        )
        self.isolate_col = knut.column_exists_or_preset(
            configure_context, _ISOLATED, input_schema, knut.is_long
        )
        self.constraints.validate_settings(input_schema)
        return input_schema

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import libpysal
        import pandas as pd

        constraintList = self.constraints.get_columns(input_1.schema)
        capacityList = self.constraints.get_capacities()

        subclusid = self.clusterid_col
        theid = "theid"
        isolateid = self.isolate_col
        keycolum = constraintList + [subclusid] + [isolateid] + [self.geo_col]
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        df = gdf[keycolum].copy()

        df = df.rename(columns={isolateid: "isolate", subclusid: "subclusid"})
        df["theid"] = list(range(df.shape[0]))
        df = df.set_index("theid", drop=False).rename_axis(None)

        keycolum1 = ["subclusid"] + constraintList
        sum_tmp = df[keycolum1].groupby("subclusid").sum().reset_index()
        classValueDict = df.groupby("subclusid")["theid"].apply(list).to_dict()

        # Create spatial weight matrix
        wq = libpysal.weights.Rook.from_dataframe(df)
        w = wq.neighbors

        df_list = []
        # Loop through the constraint and capacity lists
        for j in range(len(constraintList)):
            df_filtered = sum_tmp[sum_tmp[constraintList[j]] < capacityList[j]]
            df_list.append(df_filtered)
        df2 = pd.concat(df_list).drop_duplicates()

        if df2.shape[0] > 0:
            clusAdjList = []
            isoClus = []
            rowCount = 0
            clusAdjDict = {}
            class_value = max(df2["subclusid"])
            for index, rowLocal in df2.iterrows():
                n_loop = df2.shape[0]
                rowCount += 1
                FoundIt = False
                foundCount = 0
                currentClus = rowLocal["subclusid"]
                theCluster = currentClus
                exec_context.set_progress(
                    0.8 * rowCount / n_loop,
                    f"Batch {rowCount} of {n_loop} processed",
                )
                knut.check_canceled(exec_context)

                capA = []
                for j in range(len(constraintList)):
                    capA.append(rowLocal[constraintList[j]])

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
                    foundRow = sum_tmp[sum_tmp.subclusid == i]
                    capB = []
                    for j in range(len(constraintList)):
                        capB.append(foundRow[constraintList[j]])
                    newcap = []
                    for j in range(len(constraintList)):
                        newcap.append(capA[j] + capB[j])

                    newSatisfyAll = 1
                    for j in range(len(constraintList)):
                        if float(newcap[j]) >= float(capacityList[j]):
                            newSatisfyAll = newSatisfyAll * 1
                        else:
                            newSatisfyAll = newSatisfyAll * 0

                    # if newcap1 >= capacity1 and newcap2 >= capacity2:
                    if newSatisfyAll == 1:
                        foundCount += 1
                        if foundCount == 1:
                            mincap = []
                            for k in range(len(constraintList)):
                                mincap.append(newcap[k])
                            FoundIt = True
                            theCluster = i
                            break
                        # elif foundCount > 1 and newcap1 < mincap1 and newcap2 < mincap2:
                        elif foundCount > 1:
                            minSatisfyAll = 1
                            for j in range(len(constraintList)):
                                if float(newcap[j]) < float(mincap[j]):
                                    minSatisfyAll = minSatisfyAll * 1
                                else:
                                    minSatisfyAll = minSatisfyAll * 0

                            if minSatisfyAll == 1:
                                for k in range(len(constraintList)):
                                    mincap[k] = newcap[k]
                                FoundIt = True
                                theCluster = i
                                break

                if FoundIt == True:
                    clusAdjDict[currentClus] = theCluster
                    df.loc[(df.subclusid == currentClus), "subclusid"] = theCluster
                    df.loc[(df.subclusid == currentClus), "isolate"] = 0
                else:
                    df.loc[(df.subclusid == currentClus), "isolate"] = 1
                    isoClus.append(currentClus)

                clusAdjList.append(FoundIt)

            clusAdjust = 0
            listInd = 0

            for indexs, rowLocal in df2.iterrows():
                if len(clusAdjList) > 0:
                    if clusAdjList[listInd] == True:
                        index = rowLocal["subclusid"]
                        df.loc[(df.subclusid > (index - clusAdjust)), "subclusid"] -= 1
                        clusAdjust += 1
                        listInd += 1

        # df = df.set_index("theid", drop=False).sort_index()
        gdf[subclusid] = df.subclusid.tolist()
        gdf[isolateid] = df.isolate.tolist()

        # gdf.reset_index(drop=True, inplace=True)
        return knut.to_table(gdf, exec_context)


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
    description="Geodata with tags of cluster and isoloate status.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with modified cluster tag.",
)
class MSSCisolationNode:
    """
    Modified scale-space clustering(MSSC).

    The node forces a cluster membership for unclaimed or loose-end units, especially useful for the enclave units in
    mixed level (MLR) regionalization model.

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
        "Isolate column from MSSC",
        "Select the column 'Isoloate' generated by MSSC clustering.",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    constraints = MSSCConstraints()

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        self.clusterid_col = knut.column_exists_or_preset(
            configure_context, _CLUSTER_ID, input_schema, knut.is_long
        )
        self.isolate_col = knut.column_exists_or_preset(
            configure_context, _ISOLATED, input_schema, knut.is_long
        )
        self.constraints.validate_settings(input_schema)
        return input_schema

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        # Copy input to output
        import libpysal
        import pandas as pd

        constraintList = self.constraints.get_columns(input_1.schema)
        capacityList = self.constraints.get_capacities()

        subclusid = self.clusterid_col
        isolateid = self.isolate_col
        keycolum = constraintList + [subclusid] + [isolateid] + [self.geo_col]
        gdf0 = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        gdf = gdf0[keycolum].copy()

        gdf = gdf.rename(
            columns={
                isolateid: "isolate",
                subclusid: "subclusid",
                self.geo_col: "geometry",
            }
        )
        gdf["FinalClus"] = gdf["subclusid"]

        knut.check_canceled(exec_context)

        # Create dictionary of Constraint list
        mixStatFlds = {}
        for i in constraintList:
            mixStatFlds[i] = "sum"
        mixStatFlds["subclusid"] = "first"
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
            lyr6 = tmpMixedClusFC[tmpMixedClusFC["isolate"] == 0]
            # get weight
            wq = libpysal.weights.Rook.from_dataframe(tmpMixedClusFC)
            w_mix = wq.neighbors

            lyr6select = lyr6[lyr6.index.isin(w_mix[theClus])]
            count = lyr6select.shape[0]
            if count > 0:
                # get the smaller group
                minClus = lyr6select[constraintList[0]].idxmin()
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
                    src_column="subclusid",
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
            process_counter = i + 1
            exec_context.set_progress(
                0.8 * process_counter / isoCount1,
                f"Batch {process_counter} of {isoCount1} processed",
            )
            knut.check_canceled(exec_context)

            theClusID, theNewClusID = MergeIsolated(tmpMixedClusFC, isoClusters[i])
            gdf.loc[(gdf["subclusid"] == theClusID), fldIso] = 0
            gdf.loc[(gdf["subclusid"] == theClusID), "subclusid"] = theNewClusID

        gdf0[subclusid] = gdf.subclusid.tolist()
        gdf0[isolateid] = gdf.isolate.tolist()

        gdf0.reset_index(drop=True, inplace=True)
        return knut.to_table(gdf0, exec_context)
