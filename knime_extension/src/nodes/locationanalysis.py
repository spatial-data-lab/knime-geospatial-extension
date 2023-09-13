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
    name="Input tabale for demand and new facilities",
    description="""The table MUST contain m rows and (k+1) columns, 
    with the first column representing demand (e.g., population) 
    and the remaining k columns representing the distance 
    between each demand location and each new facility.
    """,
)
@knext.input_table(
    name="Input OD matrix for existing facilities ",
    description="""The table contains the origin-destination (OD) matrix for demand locations(m) to all existing facilities(n) .
    This table MUST have m rows and n columns, with each entry representing the distance between a demand location and an existing facility.
    """,
)
@knext.input_table(
    name="Input capacity table for existing facilities ",
    description="""The table provides information on the supply capacity of each existing facility.
    This table can be used in conjunction with the other input tables to determine the optimal location of 
    new facilities while taking into account the capacity constraints of existing facilities.
    Please ensure that the order of rows in the capacity column is consistent with the order of columns OD Matrix for Existing Facilities.
    """,
)
@knext.output_table(
    name="P-median result table",
    description="candidate column name and chosen status(1/0)",
)
@knut.pulp_node_description(
    short_description="maximal accessibility equality problem (MAEP)",
    description="""The optimization objective of MAEP is to minimize inequality in accessibility of facilities 
    and is currently formulated as minimal variance across geographic areas. Specifically, 
    it becomes a nonlinear programming (NLP) or quadratic programming (QP).
    The MAEP problem will be solved by cvxopt package. 
    """,
    references={
        "cvxopt": "https://cvxopt.org/",
    },
)
class MAEPSolverNode:

    Newcapacity = knext.DoubleParameter(
        "Input new capacity",
        "New resource or capacility for facilities .",
    )
    Demand = knext.ColumnParameter(
        "Demand column",
        "The column for demand(e.g.,population). ",
        port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Supply = knext.ColumnParameter(
        "Column for supply facility",
        "Travel cost column of  required facility or Travel cost column of nearest required facilities. ",
        port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    Decaymodel = knext.StringParameter(
        label="Distance decay model",
        description="The model representing distance decay effect. ",
        default_value="Power",
        enum=[
            "Power",
            "Exponential",
            "2SFCA",
        ],
    )
    Decaypara = knext.DoubleParameter(
        "Distance decay parameters",
        "It works for the parameters for the chosen corresponding models",
    )

    def configure(
        self, configure_context, input_schema_1, input_schema_2, input_schema_3
    ):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.string(), "FacilityID"),
                knext.Column(knext.double(), "All"),
                knext.Column(knext.double(), "Fixed"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1, input_2, input_3):
        df1 = gp.GeoDataFrame(input_1.to_pandas()).reset_index(drop=True)
        df2 = gp.GeoDataFrame(input_2.to_pandas()).reset_index(drop=True)
        df3 = gp.GeoDataFrame(input_3.to_pandas()).reset_index(drop=True)
        import numpy as np
        import math
        import cvxopt
        import pandas as pd

        D = df1[self.Demand]

        dist = pd.concat([df2, df1.drop(columns=[self.Demand])], axis=1)
        fixhosp = df3[self.Supply]

        # Calcualte indicators
        fixH = fixhosp.shape[0]
        fixcapacity = fixhosp.sum()
        demand_popu = D.sum()
        NewCapacity = 7200
        Totalcapacity = fixcapacity + self.Newcapacity
        ave_accessibility = Totalcapacity / demand_popu
        supplypt = dist.shape[1]
        demandpt = D.shape[0]
        newH = supplypt - fixH
        # Convert D dataframe to matrix
        dij = dist.values
        if self.Decaymodel == "Power":
            fij = np.power(dij, (-1 * self.Decaypara))
        elif self.Decaymodel == "Exponential":
            math.exp(-1 * self.Decaypara * dij)
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
        A2 = A2[:fixH, :]
        G2 = np.identity(supplypt)
        np.fill_diagonal(G2, -1)
        G2 = G2[fixH:, :]
        h2 = np.zeros((newH, 1))

        print(G2.shape)

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
        dff = pd.DataFrame({"FacilityID": dist.columns, "All": x, "Fixed": x1})

        return knext.Table.from_pandas(dff)