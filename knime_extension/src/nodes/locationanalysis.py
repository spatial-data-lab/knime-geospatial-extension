import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


# import pickle
# from io import StringIO
# import sys

__category = knext.category(
    path="/community/geo",
    level_id="LocationAnalysis",
    name="Location Analysis",
    description="Nodes that solve various location optimization problems.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/LocationAnalysisCategory.png",
    after="spatialmodels",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/LocationAnalysis/"

# Common used column names
_FACILITYID = "Facility ID"
_Chosen = "Chosen"


def validate_p(
    value: int,
    min_val: int = 1,
    msg: str = "Number of clusters must be no less than 1.",
) -> int:
    """
    Checks if the cluster k value is no less than 1.
    """
    if value < min_val:
        raise knext.InvalidParametersError(msg)
    return int(value)


def get_optimal_p():
    return knext.IntParameter(
        "Optimum p number",
        "The optimum number of facilities.",
        default_value=3,
        validator=validate_p,
    )


def get_id_col():
    return knext.ColumnParameter(
        "Demand ID column from demand table",
        "The column indicating the demand location IDs.",
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )


def get_required_id():
    return knext.ColumnParameter(
        "Distance column for nearest required facility",
        "The Column indicating distance to the closest required facility.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=True,
    )


def get_threshold():
    return knext.DoubleParameter(
        "Threshold",
        "Desired threshold distance or cost.",
    )


def get_candidates_dist():
    return knext.ColumnFilterParameter(
        "Distance columns for candidate facilities",
        "The columns representing the distance matrix between demand lcoation to candidate facilities.",
        column_filter=knut.is_numeric,
    )


def get_population_id():
    return knext.ColumnParameter(
        "Demand size column",
        "The column for demand size (e.g., population).",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )


# functions class for data preprocessing
class LocationData:
    # load and compare dataset
    def load(input_1, candidates_dist, required_id):
        # Get columns for distance matrix of candidate facilities
        col_list = candidates_dist.apply(input_1.schema).column_names
        if required_id in col_list:
            raise RuntimeError("Required column detected in the candidate columns ")
        if required_id is not None and required_id.lower() != "<none>":
            col_list = col_list + [required_id]
        df = input_1[col_list].to_pandas()
        return df

    # load and compare dataset
    def process(cost_matrix, required_id, p=None, nrequire=None):
        demand = cost_matrix.shape[0]
        candidate = cost_matrix.shape[1]
        nrequire = None
        # Whether exist required facilites
        if required_id is not None and required_id.lower() != "<none>":
            if p is not None:
                p = p + 1
            nrequire = candidate - 1
        # Return values based on the presence of p
        if p is None:
            return demand, candidate, nrequire
        else:
            return demand, candidate, p, nrequire

    # Define the decision variables
    def pulpxy(demand, candidate):
        import pulp

        x = pulp.LpVariable.dicts(
            "x", [(i, j) for i in range(demand) for j in range(candidate)], cat="Binary"
        )
        y = pulp.LpVariable.dicts("y", [j for j in range(candidate)], cat="Binary")
        return x, y


############################################
# Location-allocation Pmedian
############################################
@knext.node(
    name="P-median",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "pmedian.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input OD list with geometries ",
    description="Table with geometry information of demand and supply ",
)
@knext.output_table(
    name="Demand table with P-median result",
    description="Demand table with assigned supply point and link",
)
@knut.pulp_node_description(
    short_description="Solve P-median problem minimize total weighted spatial costs.",
    description="The p-median model, one of the most widely used location models of any kind,"
    + "locates p facilities and allocates demand nodes to them to minimize total weighted distance traveled. "
    + "The P-Median problem will be solved by PuLP package. ",
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class PmedianNode:
    DemandID = knext.ColumnParameter(
        "Serial id column for demand",
        "Integer id number starting with 0",
        # port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    SupplyID = knext.ColumnParameter(
        "Serial id column for supply",
        "Integer id number starting with 0",
        # port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    DemandPopu = knext.ColumnParameter(
        "Column for demand population",
        "The population of demand",
        # port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    DemandGeometry = knext.ColumnParameter(
        "Geometry column for demand points",
        "The Geometry column for demand points",
        # port_index=3,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    SupplyGeometry = knext.ColumnParameter(
        "Geometry column for supply points",
        "The Geometry column for supply points",
        # port_index=4,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    ODcost = knext.ColumnParameter(
        "Travel cost between supply and demand",
        "The travel cost between the points of supply and demand",
        # port_index=5,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    # CRSinfo = knext.StringParameter(
    #     "Input CRSinfo for geometries of demand and supply",
    #     "The CRS information for geometry columns",
    #     "epsg:3857",
    # )
    Pnumber = knext.IntParameter(
        "Input optimum p number of p-median", "The optimum p number of p-median", 5
    )

    def configure(self, configure_context, input_schema_1):
        # knut.geo_column_exists(self.geo_col, input_schema_1)
        # TODO Create combined schema
        return None
        # knext.Schema.from_columns(
        # [
        #     knext.Column(knext.(), self.DemandID),
        #     knext.Column(knext.double(), self.DemandPopu),
        #     knext.Column(knext.double(), "geometry"),
        #     knext.Column(knext.double(), "assignSID"),
        #     knext.Column(knext.double(), "SIDwkt"),
        #     knext.Column(knext.double(), "Linewkt"),
        # ])

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        df = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.DemandGeometry)
        # Sort with DID and SID
        df = df.sort_values(by=[self.DemandID, self.SupplyID]).reset_index(drop=True)

        # rebuild demand and supply data
        DemandPt = (
            df[[self.DemandID, self.DemandPopu, self.DemandGeometry]]
            .groupby([self.DemandID])
            .first()
            .reset_index()
            .rename(columns={"index": self.DemandID, self.DemandGeometry: "geometry"})
        )
        SupplyPt = (
            df[[self.SupplyID, self.SupplyGeometry]]
            .groupby([self.SupplyID])
            .first()
            .reset_index()
            .rename(columns={"index": self.SupplyID, self.SupplyGeometry: "geometry"})
        )
        DemandPt = gp.GeoDataFrame(DemandPt, geometry="geometry")
        SupplyPt = gp.GeoDataFrame(SupplyPt, geometry="geometry")
        DemandPt = DemandPt.set_crs(df.crs)
        SupplyPt = SupplyPt.set_crs(df.crs)

        # calculate parameter for matrix
        num_trt = DemandPt.shape[0]
        num_hosp = SupplyPt.shape[0]

        import pulp

        # define problem
        problem = pulp.LpProblem("pmedian", sense=pulp.LpMinimize)

        def getIndex(s):
            li = s.split()
            if len(li) == 2:
                return int(li[-1])
            else:
                return int(li[-2]), int(li[-1])

        X = [
            "X" + " " + str(i) + " " + str(j)
            for i in range(num_trt)
            for j in range(num_hosp)
        ]
        Y = ["Y" + " " + str(j) for j in range(num_hosp)]

        varX = pulp.LpVariable.dicts("x", X, cat="Binary")
        varY = pulp.LpVariable.dicts("y", Y, cat="Binary")

        object_list = []
        for it in X:
            i, j = getIndex(it)
            cost = df[(df[self.DemandID] == i) & (df[self.SupplyID] == j)][
                self.ODcost
            ].values[0]
            object_list.append(DemandPt[self.DemandPopu][i] * cost * varX[it])
        problem += pulp.lpSum(object_list)

        problem += pulp.lpSum([varY[it] for it in Y]) == self.Pnumber
        for i in range(num_trt):
            problem += (
                pulp.lpSum(
                    [varX["X" + " " + str(i) + " " + str(j)] for j in range(num_hosp)]
                )
                == 1
            )
        for i in range(num_trt):
            for j in range(num_hosp):
                problem += (
                    varX["X" + " " + str(i) + " " + str(j)] - varY["Y" + " " + str(j)]
                    <= 0
                )

        problem.solve()
        print(pulp.LpStatus[problem.status])

        # Extract result

        DemandPt["assignSID"] = None
        DemandPt["SIDcoord"] = None
        for it in X:
            v = varX[it]
            if v.varValue == 1:
                i, j = getIndex(it)
                DemandPt.at[i, "assignSID"] = SupplyPt.at[j, self.SupplyID]
                DemandPt.at[i, "SIDcoord"] = SupplyPt.at[j, "geometry"]

        from shapely.geometry import LineString

        DemandPt["linxy"] = [
            LineString(xy) for xy in zip(DemandPt["geometry"], DemandPt["SIDcoord"])
        ]
        SIDwkt = gp.GeoDataFrame(geometry=DemandPt.SIDcoord, crs=df.crs)
        Linewkt = gp.GeoDataFrame(geometry=DemandPt.linxy, crs=df.crs)
        DemandPt["SIDwkt"] = SIDwkt.geometry
        DemandPt["Linewkt"] = Linewkt.geometry
        DemandPt = DemandPt.drop(columns=["linxy", "SIDcoord"])
        DemandPt = DemandPt.reset_index(drop=True)
        return knext.Table.from_pandas(DemandPt)


############################################
# Location-allocation LSCP
############################################
@knext.node(
    name="LSCP",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "LSCP.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input OD list with geometries ",
    description="Table with geometry information of demand and supply ",
)
@knext.output_table(
    name="Demand table with LSCP result",
    description="Demand table with assigned supply point and link",
)
@knut.pulp_node_description(
    short_description="Solve LSCP problem to minimize the number of facilities.",
    description="The LSCP problem,location set covering problem (LSCP), aims to cover all demand points within a threshold distance,"
    + "with minimizing the number of facilities. "
    + "The LSCP problem will be solved by PuLP package. ",
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class LSCPNode:
    DemandID = knext.ColumnParameter(
        "Serial id column for demand",
        "Integer id number starting with 0",
        # port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    SupplyID = knext.ColumnParameter(
        "Serial id column for supply",
        "Integer id number starting with 0",
        # port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    DemandGeometry = knext.ColumnParameter(
        "Geometry column for demand points",
        "The Geometry column for demand points",
        # port_index=3,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    SupplyGeometry = knext.ColumnParameter(
        "Geometry column for supply points",
        "The Geometry column for supply points",
        # port_index=4,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    ODcost = knext.ColumnParameter(
        "Travel cost between supply and demand",
        "The travel cost between the points of supply and demand",
        # port_index=5,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    # CRSinfo = knext.StringParameter(
    #     "Input CRSinfo for geometries of demand and supply",
    #     "The CRS information for geometry columns",
    #     "epsg:3857",
    # )
    threshold = knext.DoubleParameter(
        "Input critical spatial cost", "Critical spatial cost to supply point", 13
    )

    def configure(self, configure_context, input_schema_1):
        # knut.geo_column_exists(self.geo_col, input_schema_1)
        # TODO Create combined schema
        return None
        # knext.Schema.from_columns(
        # [
        #     knext.Column(knext.(), self.DemandID),
        #     knext.Column(knext.double(), self.DemandPopu),
        #     knext.Column(knext.double(), "geometry"),
        #     knext.Column(knext.double(), "assignSID"),
        #     knext.Column(knext.double(), "SIDwkt"),
        #     knext.Column(knext.double(), "Linewkt"),
        # ])

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        df = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.DemandGeometry)
        # Sort with DID and SID
        df = df.sort_values(by=[self.DemandID, self.SupplyID]).reset_index(drop=True)

        # rebuild demand and supply data
        DemandPt = (
            df[[self.DemandID, self.DemandGeometry]]
            .groupby([self.DemandID])
            .first()
            .reset_index()
            .rename(columns={"index": self.DemandID, self.DemandGeometry: "geometry"})
        )
        SupplyPt = (
            df[[self.SupplyID, self.SupplyGeometry]]
            .groupby([self.SupplyID])
            .first()
            .reset_index()
            .rename(columns={"index": self.SupplyID, self.SupplyGeometry: "geometry"})
        )
        DemandPt = gp.GeoDataFrame(DemandPt, geometry="geometry")
        SupplyPt = gp.GeoDataFrame(SupplyPt, geometry="geometry")
        DemandPt = DemandPt.set_crs(df.crs)
        SupplyPt = SupplyPt.set_crs(df.crs)

        # calculate parameter for matrix
        num_trt = DemandPt.shape[0]
        num_hosp = SupplyPt.shape[0]

        df["within"] = 0
        df.loc[(df[self.ODcost] <= self.threshold), "within"] = 1

        import pulp

        # define problem
        problem = pulp.LpProblem("LSCPProblem", sense=pulp.LpMinimize)

        def getIndex(s):
            li = s.split()
            if len(li) == 2:
                return int(li[-1])
            else:
                return int(li[-2]), int(li[-1])

        X = ["X" + " " + str(j) for j in range(num_hosp)]
        varX = pulp.LpVariable.dicts("X", X, cat="Binary")

        object_list = []
        for it in X:
            object_list.append(varX[it])
        problem += pulp.lpSum(object_list)

        for i in range(num_trt):
            constrain = []
            for j in range(num_hosp):
                within = df[(df[self.DemandID] == i) & (df[self.SupplyID] == j)][
                    "within"
                ].values[0]
                constrain.append(within * varX["X" + " " + str(j)])
        problem += pulp.lpSum(constrain) >= 1

        problem.solve()
        # print(pulp.LpStatus[problem.status])

        DemandPt["assignSID"] = None
        DemandPt["SIDcoord"] = None

        import pandas as pd

        area = pd.DataFrame(columns=["OID"])
        for it in X:
            v = varX[it]
            if v.varValue == 1:
                print(1)
                j = getIndex(it)
                tempdict = {"OID": j}
                area = area.append(tempdict, ignore_index=True)

        for it in range(num_trt):
            index = pd.to_numeric(
                df[(df[self.DemandID] == it) & (df[self.SupplyID].isin(area["OID"]))][
                    self.ODcost
                ]
            ).idxmin()
            oid = int(df.at[index, self.SupplyID])
            DemandPt.at[it, "assignSID"] = SupplyPt.at[oid, self.SupplyID]
            DemandPt.at[it, "SIDcoord"] = SupplyPt.at[oid, "geometry"]

        from shapely.geometry import LineString

        DemandPt["linxy"] = [
            LineString(xy) for xy in zip(DemandPt["geometry"], DemandPt["SIDcoord"])
        ]
        SIDwkt = gp.GeoDataFrame(geometry=DemandPt.SIDcoord, crs=df.crs)
        Linewkt = gp.GeoDataFrame(geometry=DemandPt.linxy, crs=df.crs)
        DemandPt["SIDwkt"] = SIDwkt.geometry
        DemandPt["Linewkt"] = Linewkt.geometry
        DemandPt = DemandPt.drop(columns=["linxy", "SIDcoord"])
        DemandPt = DemandPt.reset_index(drop=True)
        return knext.Table.from_pandas(DemandPt)


############################################
# Location-allocation MCLP
############################################
@knext.node(
    name="MCLP",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "MCLP.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input OD list with geometries ",
    description="Table with geometry information of demand and supply.",
)
@knext.output_table(
    name="Demand table with MCLP result",
    description="Demand table with assigned supply point and link.",
)
@knut.pulp_node_description(
    short_description="Solve MCLP problem to Maximize Capacitated Coverage by setting an impedance cutoff.",
    description="""The MCLP model, maximum covering location problem (MCLP), aims to Locate p facilities,
and demand is covered if it is within a specified distance (time) of a facility. 
The MCLP problem will be solved by PuLP package. """,
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class MCLPNode:
    DemandID = knext.ColumnParameter(
        "Serial id column for demand",
        "Integer id number starting with 0",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    SupplyID = knext.ColumnParameter(
        "Serial id column for supply",
        "Integer id number starting with 0",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    DemandPopu = knext.ColumnParameter(
        "Column for demand population",
        "The population of demand",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    DemandGeometry = knext.ColumnParameter(
        "Geometry column for demand points",
        "The Geometry column for demand points",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    SupplyGeometry = knext.ColumnParameter(
        "Geometry column for supply points",
        "The Geometry column for supply points",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    ODcost = knext.ColumnParameter(
        "Travel cost between supply and demand",
        "The travel cost between the points of supply and demand",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Pnumber = knext.IntParameter(
        "Input optimum p number ", "The optimum p number of facilities", 3
    )

    threshold = knext.DoubleParameter(
        "Input critical spatial cost", "Critical spatial cost to supply point", 42000
    )

    def configure(self, configure_context, input_schema_1):
        # knut.geo_column_exists(self.geo_col, input_schema_1)
        # TODO Create combined schema
        return None
        # knext.Schema.from_columns(
        # [
        #     knext.Column(knext.(), self.DemandID),
        #     knext.Column(knext.double(), self.DemandPopu),
        #     knext.Column(knext.double(), "geometry"),
        #     knext.Column(knext.double(), "assignSID"),
        #     knext.Column(knext.double(), "SIDwkt"),
        #     knext.Column(knext.double(), "Linewkt"),
        # ])

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        df = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.DemandGeometry)
        # Sort with DID and SID
        df = df.sort_values(by=[self.DemandID, self.SupplyID]).reset_index(drop=True)

        # rebuild demand and supply data
        DemandPt = (
            df[[self.DemandID, self.DemandPopu, self.DemandGeometry]]
            .groupby([self.DemandID])
            .first()
            .reset_index()
            .rename(columns={"index": self.DemandID, self.DemandGeometry: "geometry"})
        )
        SupplyPt = (
            df[[self.SupplyID, self.SupplyGeometry]]
            .groupby([self.SupplyID])
            .first()
            .reset_index()
            .rename(columns={"index": self.SupplyID, self.SupplyGeometry: "geometry"})
        )
        DemandPt = gp.GeoDataFrame(DemandPt, geometry="geometry")
        SupplyPt = gp.GeoDataFrame(SupplyPt, geometry="geometry")
        DemandPt = DemandPt.set_crs(df.crs)
        SupplyPt = SupplyPt.set_crs(df.crs)

        # calculate parameter for matrix
        num_trt = DemandPt.shape[0]
        num_hosp = SupplyPt.shape[0]

        df["within"] = 0
        df.loc[(df[self.ODcost] <= self.threshold), "within"] = 1

        import pulp

        # define problem
        problem = pulp.LpProblem("MCLP", sense=pulp.LpMinimize)

        def getIndex(s):
            li = s.split()
            if len(li) == 2:
                return int(li[-1])
            else:
                return int(li[-2]), int(li[-1])

        X = ["X" + " " + str(j) for j in range(num_hosp)]
        Y = ["Y" + " " + str(i) for i in range(num_trt)]
        varX = pulp.LpVariable.dicts("assigned", X, cat="Binary")
        varY = pulp.LpVariable.dicts("covered", Y, cat="Binary")

        object_list = []
        for it in Y:
            i = getIndex(it)
            object_list.append(DemandPt[self.DemandPopu][i] * varY[it])
        problem += pulp.lpSum(object_list)

        problem += pulp.lpSum([varX[it] for it in X]) == self.Pnumber

        for i in range(num_trt):
            constrain = []
            for j in range(num_hosp):
                within = df[(df[self.DemandID] == i) & (df[self.SupplyID] == j)][
                    "within"
                ].values[0]
                constrain.append(within * varX["X" + " " + str(j)])
            constrain.append(varY["Y" + " " + str(i)])
            problem += pulp.lpSum(constrain) >= 1
        problem.solve()
        print(pulp.LpStatus[problem.status])

        DemandPt["assignSID"] = None
        DemandPt["SIDcoord"] = None

        import pandas as pd

        area = pd.DataFrame(columns=["OID"])
        for it in X:
            v = varX[it]
            if v.varValue == 1:
                j = getIndex(it)
                tempdict = {"OID": j}
                area = area.append(tempdict, ignore_index=True)

        for it in range(num_trt):
            index = pd.to_numeric(
                df[(df[self.DemandID] == it) & (df[self.SupplyID].isin(area["OID"]))][
                    self.ODcost
                ]
            ).idxmin()
            oid = int(df.at[index, self.SupplyID])
            DemandPt.at[it, "assignSID"] = SupplyPt.at[oid, self.SupplyID]
            DemandPt.at[it, "SIDcoord"] = SupplyPt.at[oid, "geometry"]

        from shapely.geometry import LineString

        DemandPt["linxy"] = [
            LineString(xy) for xy in zip(DemandPt["geometry"], DemandPt["SIDcoord"])
        ]
        SIDwkt = gp.GeoDataFrame(geometry=DemandPt.SIDcoord, crs=df.crs)
        Linewkt = gp.GeoDataFrame(geometry=DemandPt.linxy, crs=df.crs)
        DemandPt["SIDwkt"] = SIDwkt.geometry
        DemandPt["Linewkt"] = Linewkt.geometry
        DemandPt = DemandPt.drop(columns=["linxy", "SIDcoord"])
        DemandPt = DemandPt.reset_index(drop=True)
        return knext.Table.from_pandas(DemandPt)


############################################
# Location-allocation P-center Solver
############################################
@knext.node(
    name="P-center Solver",
    node_type=knext.NodeType.PREDICTOR,
    icon_path=__NODE_ICON_PATH + "pcenter.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input OD Matrix of demand locations",
    description="Distance matrix table (m rows X n columns) from demand locations(m) to (n-1) candidate facilities and the nearest required facilities.",
)
@knext.output_table(
    name="P-center result table",
    description="Candidate column names and chosen status (1/0).",
)
@knut.pulp_node_description(
    short_description="P-center (minimax) minimizes the maximum distance between demand points and their nearest facilities by locating p facilities.",
    description="""The [P-center model](https://en.wikipedia.org/wiki/Facility_location_problem) assists in selecting 
p facilities from a set of candidate facilities (n) for various locations (m), based on the OD matrix (m x n).
    
The input table should contain the columns for distance matrix to candidate facilities, and optionally a column for 
distance to closest required facilities. 

The P-center (Minimax) problem is solved using the PuLP package.
    """,
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class PcenterSolverNode:
    candidates_dist = get_candidates_dist()
    required_id = get_required_id()
    p_number = get_optimal_p()

    def configure(self, configure_context, input_schema_1):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), _FACILITYID),
                knext.Column(knext.int32(), _Chosen),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        cost_matrix = LocationData.load(input_1, self.candidates_dist, self.required_id)
        demand, candidate, p, nrequire = LocationData.process(
            cost_matrix, self.required_id, p=self.p_number
        )

        import pulp
        import pandas as pd

        # Define the optimization problem
        prob = pulp.LpProblem("p-center", pulp.LpMinimize)

        # Define the decision variables
        x, y = LocationData.pulpxy(demand, candidate)

        Z = pulp.LpVariable("Z", lowBound=0)

        # Define the objective function
        prob += Z

        # Define the constraints
        for i in range(demand):
            prob += pulp.lpSum([x[(i, j)] for j in range(candidate)]) == 1

        for j in range(candidate):
            for i in range(demand):
                prob += x[(i, j)] <= y[j]

        for i in range(demand):
            prob += (
                pulp.lpSum(
                    [x[(i, j)] * cost_matrix.iloc[i][j] for j in range(candidate)]
                )
                <= Z
            )

        prob += pulp.lpSum([y[j] for j in range(candidate)]) == p

        if self.required_id is not None and self.required_id.lower() != "<none>":
            prob += y[nrequire] == 1

        # Solve the problem
        prob.solve()

        status = pulp.LpStatus[prob.status]
        if status != "Optimal":
            raise RuntimeError(f"Solver failed with status: {status}")

        # Extract the value of y for each candidate facility
        dfy = pd.DataFrame({_FACILITYID: cost_matrix.columns})
        for j in range(candidate):
            dfy.loc[j, _Chosen] = pulp.value(y[j])
        dfy[_FACILITYID] = dfy[_FACILITYID].astype(str)
        dfy[_Chosen] = dfy[_Chosen].astype(int)
        return knext.Table.from_pandas(dfy)


############################################
# Location-allocation P-median Solver
############################################
@knext.node(
    name="P-median Solver",
    node_type=knext.NodeType.PREDICTOR,
    icon_path=__NODE_ICON_PATH + "pmedian.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input demand data",
    description="Table including a demand column and an optional distance column to the closest required facilities.",
)
@knext.output_table(
    name="P-median result table",
    description="Candidate columns name and chosen status (1/0).",
)
@knut.pulp_node_description(
    short_description="Solve the P-median problem to minimize total weighted spatial costs.",
    description="""The [P-median model,](https://en.wikipedia.org/wiki/Facility_location_problem) one of the most 
widely used location models, locates p facilities and allocates demand nodes to them in order to minimize the total 
weighted distance traveled. 

The input table should contain the three types of columns: demand size (e.g., population), 
distance matrix to candidate facilities, and optionally a column for distance to closest required facilities. 

The P-Median problem will be solved using the PuLP package.
    """,
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class PmedianSolverNode:
    candidates_dist = get_candidates_dist()
    population_id = get_population_id()
    required_id = get_required_id()
    p_number = get_optimal_p()

    def configure(self, configure_context, input_schema_1):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), _FACILITYID),
                knext.Column(knext.int32(), _Chosen),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        cost_matrix = LocationData.load(input_1, self.candidates_dist, self.required_id)
        demand, candidate, p, nrequire = LocationData.process(
            cost_matrix, self.required_id, p=self.p_number
        )

        popux = input_1[[self.population_id]].to_pandas()
        popu = popux[self.population_id].values.tolist()

        import pulp
        import pandas as pd

        # Define the optimization problem
        prob = pulp.LpProblem("p-media", pulp.LpMinimize)

        # Define the decision variables
        x, y = LocationData.pulpxy(demand, candidate)

        # Define the objective function
        prob += pulp.lpSum(
            [
                x[(i, j)] * cost_matrix.iloc[i][j] * popu[i]
                for i in range(demand)
                for j in range(candidate)
            ]
        )

        # Define the constraints
        for i in range(demand):
            prob += pulp.lpSum([x[(i, j)] for j in range(candidate)]) == 1

        for j in range(candidate):
            for i in range(demand):
                prob += x[(i, j)] <= y[j]

        prob += pulp.lpSum([y[j] for j in range(candidate)]) == p

        if self.required_id is not None and self.required_id.lower() != "<none>":
            prob += y[nrequire] == 1

        # Solve the problem
        prob.solve()

        status = pulp.LpStatus[prob.status]
        if status != "Optimal":
            raise RuntimeError(f"Solver failed with status: {status}")

        # Extract the value of y for each candidate facility
        dfy = pd.DataFrame({_FACILITYID: cost_matrix.columns})
        for j in range(candidate):
            dfy.loc[j, _Chosen] = pulp.value(y[j])
        dfy[_FACILITYID] = dfy[_FACILITYID].astype(str)
        dfy[_Chosen] = dfy[_Chosen].astype(int)
        return knext.Table.from_pandas(dfy)


############################################
# Location-allocation MCLP Solver
############################################
@knext.node(
    name="MCLP Solver",
    node_type=knext.NodeType.PREDICTOR,
    icon_path=__NODE_ICON_PATH + "MCLP.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input demand data",
    description="Table including a demand column and an optional distance column to the closest required facilities.",
)
@knext.output_table(
    name="MCLP result table",
    description="Candidate column names and chosen status (1/0).",
)
@knut.pulp_node_description(
    short_description="Solve MCLP problem to Maximize Capacitated Coverage by setting an impedance cutoff.",
    description="""The MCLP model, maximum covering location problem [(MCLP),](https://en.wikipedia.org/wiki/Maximum_coverage_problem)
aims to Locate p facilities, and demand is covered if it is within a specified distance (time) of a facility. 

The input table should contain the three types of columns: demand size (e.g., population), 
distance matrix to candidate facilities, and optionally a column for distance to closest required facilities. 

The MCLP problem will be solved by PuLP package. 
    """,
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class MCLPSolverNode:
    candidates_dist = get_candidates_dist()
    population_id = get_population_id()
    required_id = get_required_id()
    p_number = get_optimal_p()

    threshold = get_threshold()

    def configure(self, configure_context, input_schema_1):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), _FACILITYID),
                knext.Column(knext.int32(), _Chosen),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pulp
        import pandas as pd

        cost_matrix = LocationData.load(input_1, self.candidates_dist, self.required_id)
        demand, candidate, p, nrequire = LocationData.process(
            cost_matrix, self.required_id, p=self.p_number
        )

        popux = input_1[[self.population_id]].to_pandas()
        popu = popux[self.population_id].values.tolist()

        cutoff = self.threshold

        # Define the optimization problem
        prob = pulp.LpProblem("mclp", pulp.LpMinimize)

        # Define the decision variables
        x, y = LocationData.pulpxy(demand, candidate)

        # Define the objective function
        prob += pulp.lpSum(
            [x[(i, j)] * popu[i] for i in range(demand) for j in range(candidate)]
        )

        # Define the constraints
        for i in range(demand):
            prob += pulp.lpSum([x[(i, j)] for j in range(candidate)]) == 1

        for j in range(candidate):
            for i in range(demand):
                prob += x[(i, j)] <= y[j]

        # Define the threshold cutoff distance constraint
        for i in range(demand):
            prob += (
                pulp.lpSum(
                    [cost_matrix.iloc[i][j] * x[(i, j)] for j in range(candidate)]
                )
                <= cutoff
            )

        # Define the one facility required constraint
        prob += pulp.lpSum([y[j] for j in range(candidate)]) == p
        if self.required_id is not None and self.required_id.lower() != "<none>":
            prob += y[nrequire] == 1
        # Solve the problem
        prob.solve()

        status = pulp.LpStatus[prob.status]
        if status != "Optimal":
            raise RuntimeError(f"Solver failed with status: {status}")

        # Extract the value of y for each candidate facility
        dfy = pd.DataFrame({_FACILITYID: cost_matrix.columns})
        for j in range(candidate):
            dfy.loc[j, _Chosen] = pulp.value(y[j])
        dfy[_FACILITYID] = dfy[_FACILITYID].astype(str)
        dfy[_Chosen] = dfy[_Chosen].astype(int)
        return knext.Table.from_pandas(dfy)


############################################
# Location-allocation LSCP Solver
############################################
@knext.node(
    name="LSCP Solver",
    node_type=knext.NodeType.PREDICTOR,
    icon_path=__NODE_ICON_PATH + "LSCP.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input demand data",
    description="Table including OD Matrix for candidate facilities and an optional distance column to the closest required facilities.",
)
@knext.output_table(
    name="LSCP result table",
    description="candidate column names and chosen status (1/0).",
)
@knut.pulp_node_description(
    short_description="Solve LSCP problem to minimize the number of facilities within an impedance cutoff.",
    description="""The LSCP model, location set covering problem [(LSCP),](https://en.wikipedia.org/wiki/Facility_location_problem) 
aims to cover all demand points within a threshold distance. 

The input table should contain the two types of columns: distance matrix to candidate facilities, 
and optionally a column for distance to closest required facilities. 
    
The LSCP problem will be solved by PuLP package. 
    """,
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class LSCPSolverNode:
    candidates_dist = get_candidates_dist()
    required_id = get_required_id()
    threshold = get_threshold()

    def configure(self, configure_context, input_schema_1):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), _FACILITYID),
                knext.Column(knext.int32(), _Chosen),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        import pulp
        import pandas as pd

        cost_matrix = LocationData.load(input_1, self.candidates_dist, self.required_id)
        demand, candidate, nrequire = LocationData.process(
            cost_matrix, self.required_id
        )

        cutoff = self.threshold

        # transform costmatrix to binary matrix by cutoff
        cost_matrix = (cost_matrix <= cutoff).astype(int)

        # Define the optimization problem
        prob = pulp.LpProblem("lscp", pulp.LpMinimize)

        # Define the decision variables
        x, y = LocationData.pulpxy(demand, candidate)

        # Define the objective function to minimize the number of facilities
        prob += pulp.lpSum([y[j] for j in range(candidate)])

        # Define the constraints
        for i in range(demand):
            prob += (
                pulp.lpSum(
                    [x[(i, j)] * cost_matrix.iloc[i][j] for j in range(candidate)]
                )
                >= 1
            )

        for j in range(candidate):
            for i in range(demand):
                prob += x[(i, j)] <= y[j]

        if self.required_id is not None and self.required_id.lower() != "<none>":
            prob += y[nrequire] == 1
        # Solve the problem
        prob.solve()

        status = pulp.LpStatus[prob.status]
        if status != "Optimal":
            raise RuntimeError(f"Solver failed with status: {status}")

        # Extract the value of y for each candidate facility
        dfy = pd.DataFrame({_FACILITYID: cost_matrix.columns})
        for j in range(candidate):
            dfy.loc[j, _Chosen] = pulp.value(y[j])
        dfy[_FACILITYID] = dfy[_FACILITYID].astype(str)
        dfy[_Chosen] = dfy[_Chosen].astype(int)
        return knext.Table.from_pandas(dfy)


############################################
# Location-allocation MAEP Solver
############################################
@knext.node(
    name="MAEP Solver",
    node_type=knext.NodeType.PREDICTOR,
    icon_path=__NODE_ICON_PATH + "MAEP.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Demand table",
    description=""" Each row in this table represents a demand location (m). 
    It consists of three types of columns: demand size (e.g., population), 
    distance columns to the existing facilities (n) and new facilities (k).
    """,
)
@knext.input_table(
    name="Capacity table for existing facilities ",
    description="""This table provides information regarding the supply capacity of each existing facility.
    The values in the Facility ID column must exactly match the column names for facilities in the demand table.
    """,
)
@knext.output_table(
    name="MAEP result table",
    description="Facilities with assigned capacities.",
)
class MAEPSolverNode:
    """
    Maximal Accessibility Equality Problem (MAEP).

    The Maximal Accessibility Equality Problem (MAEP), originally introduced by
    [Jin et al.,](https://doi.org/10.1155/2017/2094654) addresses capacity adjustments for minimizing inequality in
    accessibility. It takes into account the match ratio between supply and demand, along with intricate spatial
    interactions.

    The optimization objective of MAEP aims to minimize inequality in facility accessibility, with a focus on
    reducing variance across geographic areas. This problem can be formulated as either a Nonlinear Programming (NLP)
    or a Quadratic Programming (QP) task.

    The top input table contains three types of columns: demand size (e.g., population), distance columns to the
    existing facilities (n) and new facilities (k).
    The bottom input table contains two columns: the ID and the capacity columns for existing facilities.

    The result table comprises three columns: Facility IDs and their assigned capacities under two scenarios, denoted
    by the 'All' and 'Fixed' columns:
    Global Optimization: New capacity is allocated to both existing and new facilities by adding the specified
    new capacity to the total capacity of existing facilities.
    Local Optimization: New capacity is exclusively assigned to new facilities based on the specified new capacity,
    with no impact on the capacities of existing facilities.

    To solve the MAEP problem, this node utilizes the [cvxopt package.](https://cvxopt.org/)

    """

    class DistDecayModes(knext.EnumParameterOptions):
        POWER = (
            "Power",
            """Apply the power function, power(distance, -n), to model distance decay effect.""",
        )
        EXPONENTIAL = (
            "Exponential",
            """Apply the exponential function, exp(-distance * n), to model distance decay effect.""",
        )
        BINARY = (
            "2SFCA binary",
            "Use a threshold to create a binary value for distance decay effect.",
        )

        @classmethod
        def get_default(cls):
            return cls.POWER

    # demand
    population_id = get_population_id()
    existing_dist = knext.ColumnFilterParameter(
        "Distance columns",
        "Distance columns representing the distances to the existing facilities.",
        column_filter=knut.is_numeric,
    )
    candidates_dist = get_candidates_dist()

    id_supply = knext.ColumnParameter(
        "Facility ID column ",
        "The column for facility IDs.",
        port_index=1,
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )
    supply = knext.ColumnParameter(
        "Capacity column of existing facilities",
        "The column representing the capacities of existing supply facilities.",
        port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    decay_model = knext.EnumParameter(
        label="Distance decay model",
        description="The model representing distance decay effect.",
        default_value=DistDecayModes.get_default().name,
        enum=DistDecayModes,
    )
    decay_parameter = knext.DoubleParameter(
        "Distance decay parameters",
        "It works for the parameters for the chosen corresponding models.",
    )
    new_capacity = knext.DoubleParameter(
        "Input new capacity",
        "Total capacity assigned to all new facilities.",
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), "Facility ID"),
                knext.Column(knext.double(), "All"),
                knext.Column(knext.double(), "Fixed"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2):
        # Extract data list
        existing_list = self.existing_dist.apply(input_1.schema).column_names
        candidate_list = self.candidates_dist.apply(input_1.schema).column_names
        df1 = input_1.to_pandas()
        # Get data for existing facility capacity
        df2 = input_2.to_pandas()
        # Align data for existing facility
        df_existing = df1[existing_list]
        df_existing = df_existing.reindex(sorted(df_existing.columns), axis=1)
        df2 = df2.sort_values(by=self.id_supply)
        # Check column names
        colname = df_existing.columns.tolist()
        rolname = df2[self.id_supply].tolist()
        if colname != rolname:
            raise RuntimeError(
                "Facility ID inconsistency between Distance Matrix and Facility Data. Ensure consistent IDs."
            )

        import numpy as np
        import math
        import cvxopt
        import pandas as pd

        D = df1[self.population_id]
        df_candidate = df1[candidate_list]

        dist = pd.concat([df_existing, df_candidate], axis=1)
        fixhosp = df2[self.supply]

        # Calcualte indicators
        fixh = fixhosp.shape[0]
        fixcapacity = fixhosp.sum()
        demand_popu = D.sum()
        Totalcapacity = fixcapacity + self.new_capacity
        ave_accessibility = Totalcapacity / demand_popu
        supplypt = dist.shape[1]
        demandpt = D.shape[0]
        newH = supplypt - fixh
        # Convert D dataframe to matrix
        dij = dist.values
        if self.decay_model == self.DistDecayModes.POWER.name:
            fij = np.power(dij, (-1 * self.decay_parameter))
        elif self.decay_model == self.DistDecayModes.EXPONENTIAL.name:
            math.exp(-1 * self.decay_parameter * dij)
            fij = np.power(dij, (-1 * self.decay_parameter))
        else:
            fij = np.where(fij <= self.decay_parameter, 1, 0)

        Dfij = fij * D.values.reshape(-1, 1)
        # Compute the row sums
        Dfij_sums = Dfij.sum(axis=0)
        # Divide each value by its corresponding row sum
        Fij = np.divide(fij, Dfij_sums)

        vec_D = D.values.squeeze()
        D = np.diag(vec_D)
        A = np.ones((demandpt, 1)) @ np.reshape(
            ave_accessibility, (1,)
        )  # compute average accessibility

        # Compute optimization problem parameters
        Dmat = Fij.T @ D @ Fij
        dvec = Fij.T @ D @ A
        dvec = dvec * -1

        # Ax = b, Gx ≤ h, here, B=E, x=S, d=supply_popu
        A1 = np.ones((1, supplypt))
        b1 = np.array([Totalcapacity])
        # AA = np.diag(np.ones(supplypt))
        G1 = np.identity(supplypt)
        np.fill_diagonal(G1, -1)
        h1 = np.zeros((supplypt, 1))

        # Define the quadratic objective function
        Q = cvxopt.matrix(Dmat)
        p = cvxopt.matrix(dvec)

        # Define the constraints
        G = cvxopt.matrix(G1)
        h = cvxopt.matrix(h1)
        A = cvxopt.matrix(A1)
        b = cvxopt.matrix(b1)

        # Solve the QP problem
        solution = cvxopt.solvers.qp(Q, p, G, h, A, b)

        # Extract the solution
        x = solution["x"]
        A2 = np.identity(supplypt)
        A2 = A2[:fixh, :]
        G2 = np.identity(supplypt)
        np.fill_diagonal(G2, -1)
        G2 = G2[fixh:, :]
        h2 = np.zeros((newH, 1))

        # print(G2.shape)

        # build constraints
        A3 = np.vstack((A1, A2))
        b2 = np.hstack((b1, fixhosp.values.flatten()))
        # Define the constraints
        Gn = cvxopt.matrix(G2)
        hn = cvxopt.matrix(h2)
        An = cvxopt.matrix(A3)
        bn = cvxopt.matrix(b2)

        # Solve the QP problem
        solution = cvxopt.solvers.qp(Q, p, Gn, hn, An, bn)
        # Extract the solution
        x1 = solution["x"]
        dff = pd.DataFrame({"Facility ID": dist.columns, "All": x, "Fixed": x1})

        return knext.Table.from_pandas(dff)
