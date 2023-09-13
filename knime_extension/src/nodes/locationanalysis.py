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
    description="Table with geometry information of demand and supply",
)
@knext.output_table(
    name="Demand table with MCLP result",
    description="Demand table with assigned supply point and link",
)
@knut.pulp_node_description(
    short_description="Solve MCLP problem to Maximize Capacitated Coverage by setting an impedance cutoff.",
    description="The MCLP model, maximum covering location problem (MCLP), aims to Locate p facilities,"
    + " and demand is covered if it is within a specified distance (time) of a facility. "
    + "The P-Median problem will be solved by PuLP package. ",
    references={
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",
    },
)
class MCLPNode:
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
    name="Input demand table",
    description=""" Each row in this table represents a demand location (m). It consists of three types of columns: 
    - demand ID
    - demand size (e.g., population)
    - distance columns, which represent the distances to the new facilities (k).
    """,
)
@knext.input_table(
    name="Input distance matrix to existing facilities ",
    description="""This table represents the Origin-Destination (OD) matrix, showing the distances between demand locations (m) and existing facilities (n).
    It MUST consist of m rows and (n+1) columns, where one column is dedicated to Demand ID, 
    and the remaining columns represent the distances between each demand location and an existing facility.
    """,
)
@knext.input_table(
    name="Input capacity table for existing facilities ",
    description="""This table provides information regarding the supply capacity of each existing facility.
    The values in the Facility ID column must exactly match the column names for facilityies in the distance matrix table.
    """,
)
@knext.output_table(
    name="MAEP result table",
    description="Facilities with assigned capacities",
)
class MAEPSolverNode:
    """
    Maximal Accessibility Equality Problem (MAEP).

    The Maximal Accessibility Equality Problem (MAEP), originally introduced by [Jin et al.](https://doi.org/10.1155/2017/2094654),
    addresses capacity adjustments for minimizing inequality in accessibility. It takes into account the match ratio between supply and demand, along with intricate spatial interactions.

    The optimization objective of MAEP aims to minimize inequality in facility accessibility,
    with a focus on reducing variance across geographic areas. This problem can be formulated as either a Nonlinear Programming (NLP) or a Quadratic Programming (QP) task.

    The result table comprises three columns: Facility IDs and their assigned capacities under two scenarios, denoted by the 'All' and 'Fixed' columns:
    Global Optimization: New capacity is allocated to both existing and new facilities by adding the specified new capacity to the total capacity of existing facilities.
    Local Optimization: New capacity is exclusively assigned to new facilities based on the specified new capacity, with no impact on the capacities of existing facilities.

    To solve the MAEP problem, this node utilizes the [cvxopt package](https://cvxopt.org/).

    """

    id_left = knext.ColumnParameter(
        "Demand ID column from demand table",
        "The column for demand IDs from the first table. ",
        port_index=0,
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    demand = knext.ColumnParameter(
        "Demand size column",
        "The column for demand size (e.g.,population). ",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    newdistance = knext.ColumnFilterParameter(
        "Distance columns",
        "Distance columns representing the distances to the new facilities ",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    id_right = knext.ColumnParameter(
        "Demand ID column from distance matrix table",
        "The column for demand IDs from the second input table. ",
        port_index=1,
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )
    id_supply = knext.ColumnParameter(
        "Facility ID column ",
        "The column for facility IDs from the third input table. ",
        port_index=2,
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )
    supply = knext.ColumnParameter(
        "Capacity column of existing facility",
        "The column representing the capacities of existing supply facilities.",
        port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    decaymodel = knext.StringParameter(
        label="Distance decay model",
        description="The model representing distance decay effect. ",
        default_value="Power",
        enum=[
            "Power",
            "Exponential",
            "2SFCA",
        ],
    )
    decaypara = knext.DoubleParameter(
        "Distance decay parameters",
        "It works for the parameters for the chosen corresponding models.",
    )
    newcapacity = knext.DoubleParameter(
        "Input new capacity",
        "Total capacility assigned to all new facilities.",
    )

    def configure(
        self, configure_context, input_schema_1, input_schema_2, input_schema_3
    ):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), "Facility ID"),
                knext.Column(knext.double(), "All"),
                knext.Column(knext.double(), "Fixed"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2, input_3):
        dfn1 = gp.GeoDataFrame(input_1.to_pandas()).reset_index(drop=True)
        df2 = gp.GeoDataFrame(input_2.to_pandas()).reset_index(drop=True)
        df3 = gp.GeoDataFrame(input_3.to_pandas()).reset_index(drop=True)
        import numpy as np
        import math
        import cvxopt
        import pandas as pd

        df1_list = self.newdistance.extend([self.id_left, self.demand])
        df1 = dfn1[[df1_list]]
        # TODO , use Class function from another pull request
        df1 = df1.sort_values(by=self.id_left).drop(columns=[self.id_left])
        df2 = df2.sort_values(by=self.id_right).drop(columns=[self.id_right])
        try:
            df2 = df2[df3[self.id_supply]]
        except:
            raise RuntimeError(
                "Facility ID inconsistency between Distance Matrix and Facility Data. Ensure consistent and aligned IDs."
            )

        D = df1[self.demand]

        dist = pd.concat([df2, df1.drop(columns=[self.demand])], axis=1)
        fixhosp = df3[self.supply]

        # Calcualte indicators
        fixh = fixhosp.shape[0]
        fixcapacity = fixhosp.sum()
        demand_popu = D.sum()
        Totalcapacity = fixcapacity + self.newcapacity
        ave_accessibility = Totalcapacity / demand_popu
        supplypt = dist.shape[1]
        demandpt = D.shape[0]
        newH = supplypt - fixh
        # Convert D dataframe to matrix
        dij = dist.values
        if self.decaymodel == "Power":
            fij = np.power(dij, (-1 * self.Decaypara))
        elif self.decaymodel == "Exponential":
            math.exp(-1 * self.decaypara * dij)
            fij = np.power(dij, (-1 * self.Decaypara))
        else:
            fij = np.where(fij <= self.Decaypara, 1, 0)

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
