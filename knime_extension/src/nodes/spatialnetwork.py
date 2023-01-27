import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut

__category = knext.category(
    path="/community/geo",
    level_id="spatialnetwork",
    name="Spatial Network",
    description="Spatial network analysis.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialnetworkCategroy.png",
    after="opendataset",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Spatialnetwork/"


############################################
# Google Distance Matrix
############################################
@knext.node(
    name="Google Distance Matrix",
    node_type=knext.NodeType.MANIPULATOR,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GoogleDistMatrix.png",
)
@knext.input_table(
    name="Input Table as origin",
    description="Input origin table with geometry.",
)
@knext.input_table(
    name="Input Table as destination",
    description="Input destination table with geometry.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with travel Cost in Minutes and Meters.",
)
class GoogleDistanceMatrix:
    """
    Google Distance Matrix.
    The [Distance Matrix API](https://developers.google.com/maps/documentation/distance-matrix/overview)
    provides travel distance and time for a matrix of origins and destinations,
    and consists of rows containing duration and distance values for each pair.
    The API returns information based on the recommended route between start and end points.
    This node provides two travel modes: driving and transit. The units of distance and time are Meters and Minutes.
    If the input geometry is not point feature, the centroids will be used.
    The output includes numerical indices for the origin and destination data that will serve as a common key for merging the data.
    """

    O_geo_col = knext.ColumnParameter(
        "Origin Geometry Column",
        "Select the geometry column as origin.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    D_geo_col = knext.ColumnParameter(
        "Right geometry column",
        "Select the geometry column as destination.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    API_Key = knext.StringParameter(
        "Google API Key",
        "[Google API Key](https://developers.google.com/maps/documentation/distance-matrix/overview)  for distance matrix",
        "",
    )
    Travel_Mode = knext.StringParameter(
        "Travel Mode",
        "Set [Trave Mode](https://developers.google.com/maps/documentation/distance-matrix/distance-matrix) ",
        "Driving",
        enum=["Driving", "Transit"],
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.O_geo_col = knut.column_exists_or_preset(
            configure_context, self.O_geo_col, input_schema_1, knut.is_geo
        )
        self.D_geo_col = knut.column_exists_or_preset(
            configure_context, self.D_geo_col, input_schema_2, knut.is_geo
        )
        return knext.Schema.from_columns(
            [
                knext.Column(knext.int64(), "originid"),
                knext.Column(knext.int64(), "destinationid"),
                knext.Column(knext.double(), "duration"),
                knext.Column(knext.double(), "distance"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        # the output distance will have units of miles
        basestring = "https://maps.googleapis.com/maps/api/distancematrix/json?language=en&units=imperial&origins={0}&destinations={1}&key={2}&mode={3}"

        # define function to derive travel time and distance from Google Maps API
        import urllib.request as urllib2  # for Google Drive
        import json  # for OSRM

        def fetch_google_OD(download_link):
            nt_time = 0
            nt_dist = 0
            try:
                req = urllib2.urlopen(download_link)
                jsonout = json.loads(req.read())
                nt_time = jsonout["rows"][0]["elements"][0]["duration"][
                    "value"
                ]  # meters
                nt_dist = jsonout["rows"][0]["elements"][0]["distance"][
                    "value"
                ]  # seconds

                # transform seconds to minutes and meters to miles
                nt_time = round(nt_time / 60, 2)
                nt_dist = round(nt_dist, 2)
            except:
                print("Empty value generated")
            return [nt_time, nt_dist]

        O_gdf = left_input.to_pandas().rename(columns={self.O_geo_col: "geometry"})
        D_gdf = right_input.to_pandas().rename(columns={self.D_geo_col: "geometry"})
        O_gdf = gp.GeoDataFrame(O_gdf, geometry="geometry")
        D_gdf = gp.GeoDataFrame(D_gdf, geometry="geometry")

        # Set a lat\Lon CRS
        O_gdf = O_gdf.to_crs(4326)
        D_gdf = D_gdf.to_crs(4326)
        # Generate ID
        O_gdf["originid"] = range(1, (O_gdf.shape[0] + 1))
        D_gdf["destinationid"] = range(1, (D_gdf.shape[0] + 1))
        mergedf = O_gdf.merge(D_gdf, how="cross")
        mergedf_x = gp.GeoDataFrame(geometry=mergedf["geometry_x"])
        mergedf_y = gp.GeoDataFrame(geometry=mergedf["geometry_y"])

        distance_matrix = mergedf[["originid", "destinationid"]]
        distance_matrix["duration"] = 0
        distance_matrix["distance"] = 0
        # can process point or polygon features or both
        for i in range(0, len(distance_matrix)):
            Origins_latlon = "{},{}".format(
                mergedf_x.centroid.y[i], mergedf_x.centroid.x[i]
            )
            Destinations_latlon = "{},{}".format(
                mergedf_y.centroid.y[i], mergedf_y.centroid.x[i]
            )
            Google_Request_Link = basestring.format(
                Origins_latlon,
                Destinations_latlon,
                self.API_Key,
                self.Travel_Mode.lower(),
            )
            Google_Travel_Cost = fetch_google_OD(Google_Request_Link)
            distance_matrix.iloc[i, 2] = Google_Travel_Cost[0]
            distance_matrix.iloc[i, 3] = Google_Travel_Cost[1]
        distance_matrix["duration"] = distance_matrix.duration / 1.0
        distance_matrix["distance"] = distance_matrix.distance / 1.0
        return knext.Table.from_pandas(distance_matrix)


############################################
# OSRM
############################################
@knext.node(
    name="OSRM Matrix",
    node_type=knext.NodeType.MANIPULATOR,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "OSRMdistMatrix.png",
)
@knext.input_table(
    name="Input Table as origin",
    description="Input origin table with geometry.",
)
@knext.input_table(
    name="Input Table as destination",
    description="Input destination table with geometry.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with driving Cost in Minutes and Meters.",
)
class OSRMDriveMatrix:
    """
    This node can calculate the driving time, distance or route of the shortest paths between origin (left or upper input port)
    and destination (right or lower input port) points based on the Open Source Routing Machine or [OSRM](https://project-osrm.org/),
    which  is a C++ implementation of a high-performance routing engine for shortest paths in road networks. It combines sophisticated
    routing algorithms with the open and free road network data of the OpenStreetMap (OSM) project.
    If the input geometry is not point feature, the centroids will be used.
    The output includes numerical indices for the source and destination that will serve as a common key for merging the data.

    """

    # input parameters
    O_geo_col = knext.ColumnParameter(
        "Origin Geometry Column(Points)",
        "Select the geometry column as origin.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    D_geo_col = knext.ColumnParameter(
        "Right geometry column",
        "Select the geometry column as destination.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    osrmmodel = knext.StringParameter(
        label="Set Network Route Model",
        description="""Available Model are: 
        
        - **Travel Cost:** contains drive distance in meters and travel time in minutes.
        - **Route:** contains line feature for driving routes.
        - **Travel Cost and Route:** contains both travel cost and route.

""",
        default_value="Travel Cost",
        enum=[
            "Travel Cost",
            "Route",
            "Travel Cost and Route",
        ],
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        self.O_geo_col = knut.column_exists_or_preset(
            configure_context, self.O_geo_col, input_schema_1, knut.is_geo
        )
        self.D_geo_col = knut.column_exists_or_preset(
            configure_context, self.D_geo_col, input_schema_2, knut.is_geo
        )
        if self.osrmmodel == "Travel Cost":
            return knext.Schema.from_columns(
                [
                    knext.Column(knext.int64(), "originid"),
                    knext.Column(knext.int64(), "destinationid"),
                    knext.Column(knext.double(), "duration"),
                    knext.Column(knext.double(), "distance"),
                ]
            )
        elif self.osrmmodel == "Route":
            return knext.Schema.from_columns(
                [
                    knext.Column(knext.int64(), "originid"),
                    knext.Column(knext.int64(), "destinationid"),
                    knext.Column(knut.TYPE_LINE, "geometry"),
                ]
            )
        else:
            return knext.Schema.from_columns(
                [
                    knext.Column(knext.int64(), "originid"),
                    knext.Column(knext.int64(), "destinationid"),
                    knext.Column(knext.double(), "duration"),
                    knext.Column(knext.double(), "distance"),
                    knext.Column(knut.TYPE_LINE, "geometry"),
                ]
            )

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):

        import pandas as pd
        from shapely.geometry import LineString
        import requests  # for OSRM
        import json  # for OSRM
        import polyline

        # set digits for coordinates
        def RoundCoordList(coordlist, digits):
            coordlist = list(map(lambda x: [round(i, digits) for i in x], coordlist))
            coordlist = [tuple(i) for i in coordlist]
            return coordlist

        # extract route geometry
        def extractRoute(data):
            from shapely.geometry import LineString

            decodeline = polyline.decode(data["routes"][0]["geometry"])
            # Extract the location coordinates from the 'waypoints' field
            coordinates = [
                round(coord, 5)
                for waypoint in data["waypoints"]
                for coord in waypoint["location"]
            ]
            points = [
                (coordinates[i + 1], coordinates[i])
                for i in range(0, len(coordinates), 2)
            ]
            decodeline4 = RoundCoordList(decodeline, 4)
            points4 = RoundCoordList(points, 4)
            indexes = []
            tag = 0
            for i in points4:
                newline = decodeline4[tag:]
                for j, p in enumerate(newline):
                    if i == p:
                        tag = j + tag
                        break
                indexes.append(tag)
                tag = tag + 1
            redecodeline = [(y, x) for x, y in decodeline]
            getRoutes = [
                LineString(redecodeline[indexes[i] : (indexes[(i + 1)] + 1)])
                for i in range(0, len(indexes), 2)
            ]
            return getRoutes

        # update travel cost and route geometry
        def updatePart(df, ns, ne):
            dfs = df.copy().loc[ns:ne]
            dfs = dfs[["StartX", "StartY", "EndX", "EndY"]]
            dfs = dfs.astype(str)
            dfs["period"] = (
                dfs["StartX"]
                + ","
                + dfs["StartY"]
                + ";"
                + dfs["EndX"]
                + ","
                + dfs["EndY"]
            )
            Querylist = [";".join(dfs["period"])]
            address = osrm_route_service + Querylist[0]
            try:
                r = requests.get(
                    address, params={"continue_straight": "false"}, timeout=None
                )
                data = json.loads(r.text)
                if data["code"] == "Ok":
                    if self.osrmmodel != "Route":
                        dfr = pd.DataFrame(data["routes"][0]["legs"])[
                            ["duration", "distance"]
                        ].iloc[::2]
                        df.loc[ns:ne, "duration"] = dfr.duration.to_list()
                        df.loc[ns:ne, "distance"] = dfr.distance.to_list()
                    if self.osrmmodel != "Travel Cost":
                        # get route
                        temproute = extractRoute(data)
                        # get route
                        if len(temproute) == 1:
                            df.loc[ns:ne, "geometry"] = temproute[0]
                        else:
                            df.loc[ns:ne, "geometry"] = temproute
                else:
                    print("error from:{} to :{}".format(ns, ne))
            except:
                print("error from:{} to :{}".format(ns, ne))

        # Cross join Data
        O_gdf = left_input.to_pandas().rename(columns={self.O_geo_col: "geometry"})
        D_gdf = right_input.to_pandas().rename(columns={self.D_geo_col: "geometry"})
        O_gdf = gp.GeoDataFrame(O_gdf, geometry="geometry")
        D_gdf = gp.GeoDataFrame(D_gdf, geometry="geometry")
        # Set a lat\Lon CRS
        O_gdf = O_gdf.to_crs(4326)
        D_gdf = D_gdf.to_crs(4326)
        # Generate ID
        O_gdf["originid"] = range(1, (O_gdf.shape[0] + 1))
        D_gdf["destinationid"] = range(1, (D_gdf.shape[0] + 1))
        mergedf = O_gdf.merge(D_gdf, how="cross")
        mergedf_x = gp.GeoDataFrame(geometry=mergedf["geometry_x"])
        mergedf_y = gp.GeoDataFrame(geometry=mergedf["geometry_y"])

        df = mergedf[["originid", "destinationid"]]
        df["StartX"] = mergedf_x.centroid.x
        df["StartY"] = mergedf_x.centroid.y
        df["EndX"] = mergedf_y.centroid.x
        df["EndY"] = mergedf_y.centroid.y
        df = df.reset_index(drop=True)
        # set minimun set
        slnum = 50
        nlength = df.shape[0]
        nloop = nlength // slnum
        ntail = nlength % slnum
        if self.osrmmodel != "Route":
            df["duration"] = 0.0
            df["distance"] = 0.0
        if self.osrmmodel != "Travel Cost":
            df["geometry"] = LineString([(0, 0), (1, 1)])
        osrm_route_service = "http://router.project-osrm.org/route/v1/driving/"
        if nloop >= 1:
            for i in range(nloop):
                ns = slnum * i
                ne = ns + slnum - 1
                updatePart(df, ns, ne)
        if ntail > 0:
            ns = slnum * nloop
            ne = ns + ntail - 1
            updatePart(df, ns, ne)

        if self.osrmmodel == "Travel Cost":
            gdf = df[["originid", "destinationid", "duration", "distance"]]
            gdf["duration"] = gdf.duration / 60
        else:
            if self.osrmmodel == "Route":
                df1 = df[["originid", "destinationid"]]
            else:
                df1 = df[["originid", "destinationid", "duration", "distance"]]
                df1["duration"] = df1.duration / 60
            gdf = gp.GeoDataFrame(df1, geometry=df.geometry, crs=4326)
        return knext.Table.from_pandas(gdf)


############################################
# Road Network
############################################
@knext.node(
    name="Street Network Matrix",
    node_type=knext.NodeType.MANIPULATOR,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "OSMDistMatrix.png",
)
@knext.input_table(
    name="Input Table as origin",
    description="Input origin table with geometry.",
)
@knext.input_table(
    name="Input Table as destination",
    description="Input destination table with geometry.",
)
@knext.input_table(
    name="Input Table as road network",
    description="Input road network with LineString geometry.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with travel time and distance.",
)
class StreetNetworkMatrix:
    """
    This node can calculate the travel time, distance  between origin (left or upper input port)
    and destination (right or lower input port) points based on the graph derivred from the input road network and its speed column.
    It first snaps the origin and destination points to the road work, and uses the function
    [single_source_dijkstra_path_length](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.algorithms.shortest_paths.weighted.single_source_dijkstra_path_length.html)
    in [NetworkX](https://networkx.org/) to compute the shortest path length between source and all other reachable nodes for the weighted (time or distance) graph.

    The column for speed value in the road data is used to calculate travle time with the formula "length/speed", so the speed value should better be in meters/minute or meters/second.
    E.g., 1 mile/hour = 26.8224 meter/minute, 1 Kilometer/hour = 16.6667 meter/minute.

    The output table contains the distance between the orginal points and snap points along the road network, which can also be further calculated as additional travel time or distance with user-defined speed.
    The calcualtion depends on a projected coordinates systemn of road network datasets. If it doesn't have a projected CRS information, a default CRS"3857" will be applied.
    If the input geometry is not point feature, the centroids will be used.
    The output includes numerical indices for the origin and destination that will serve as a common key for merging the data.

    """

    # input parameters
    O_geo_col = knext.ColumnParameter(
        "Origin Geometry Column",
        "Select the geometry column as origin.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    D_geo_col = knext.ColumnParameter(
        "Destination geometry column",
        "Select the geometry column as destination.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    R_geo_col = knext.ColumnParameter(
        "Road Network geometry column",
        "Select the geometry column as road network.",
        # Allow only GeoValue compatible columns
        port_index=2,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    R_speed_col = knext.ColumnParameter(
        "Speed column from Road Network",
        "Select the speed column from road network.",
        # Allow only GeoValue compatible columns
        port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(
        self,
        configure_context,
        input_schema_1,
        input_schema_2,
        input_schema_3,
    ):
        self.O_geo_col = knut.column_exists_or_preset(
            configure_context, self.O_geo_col, input_schema_1, knut.is_geo
        )
        self.D_geo_col = knut.column_exists_or_preset(
            configure_context, self.D_geo_col, input_schema_2, knut.is_geo
        )
        self.R_geo_col = knut.column_exists_or_preset(
            configure_context, self.R_geo_col, input_schema_3, knut.is_geo
        )
        return knext.Schema.from_columns(
            [
                knext.Column(knext.int64(), "originid"),
                knext.Column(knext.int64(), "destinationid"),
                knext.Column(knext.double(), "duration"),
                knext.Column(knext.double(), "distance"),
                knext.Column(knext.double(), "snap_o_dist"),
                knext.Column(knext.double(), "snap_d_dist"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input1, input2, input3):

        import pandas as pd
        import networkx as nx
        from shapely.geometry import MultiPoint, LineString, Point
        import rtree
        import numpy as np
        from shapely.ops import snap, split
        import itertools
        from scipy.spatial import cKDTree
        from pyproj import CRS  # For CRS Units check
        import logging

        # Start ---------------Repeated code as those in Isochrone Map
        # Use momepy function: generate_primal,gdf_to_nx,primal_to_gdf,nx_to_gdf
        def generate_primal(G, gdf_network, fields, multigraph):
            G.graph["approach"] = "primal"
            key = 0
            for row in gdf_network.itertuples():

                first = row.geometry.coords[0]
                last = row.geometry.coords[-1]

                data = [r for r in row][1:]
                attributes = dict(zip(fields, data))
                if multigraph:
                    G.add_edge(first, last, key=key, **attributes)
                    key += 1
                else:
                    G.add_edge(first, last, **attributes)

        def gdf_to_nx(
            gdf_network,
            approach="primal",
            length="mm_len",
            multigraph=True,
            directed=False,
            angles=True,
            angle="angle",
        ):
            gdf_network = gdf_network.copy()
            if "key" in gdf_network.columns:
                gdf_network.rename(columns={"key": "__key"}, inplace=True)

            if multigraph and directed:
                net = nx.MultiDiGraph()
            elif multigraph and not directed:
                net = nx.MultiGraph()
            elif not multigraph and directed:
                net = nx.DiGraph()
            else:
                net = nx.Graph()

            net.graph["crs"] = gdf_network.crs
            gdf_network[length] = gdf_network.geometry.length
            fields = list(gdf_network.columns)
            generate_primal(net, gdf_network, fields, multigraph)

            return net

        def primal_to_gdf(net, points, lines, spatial_weights, nodeID):
            if points is True:
                gdf_nodes = points_to_gdf(net)

                if spatial_weights is True:
                    W = libpysal.weights.W.from_networkx(net)
                    W.transform = "b"

            if lines is True:
                gdf_edges = lines_to_gdf(net, points, nodeID)

            if points is True and lines is True:
                if spatial_weights is True:
                    return gdf_nodes, gdf_edges, W
                return gdf_nodes, gdf_edges
            if points is True and lines is False:
                if spatial_weights is True:
                    return gdf_nodes, W
                return gdf_nodes
            return gdf_edges

        def nx_to_gdf(
            net, points=True, lines=True, spatial_weights=False, nodeID="nodeID"
        ):
            # generate nodes and edges geodataframes from graph
            primal = True
            for nid, n in enumerate(net):
                net.nodes[n][nodeID] = nid
            return primal_to_gdf(
                net,
                points=points,
                lines=lines,
                spatial_weights=spatial_weights,
                nodeID=nodeID,
            )

        def points_to_gdf(net):
            node_xy, node_data = zip(*net.nodes(data=True))
            if isinstance(node_xy[0], int) and "x" in node_data[0].keys():
                geometry = [
                    Point(data["x"], data["y"]) for data in node_data
                ]  # osmnx graph
            else:
                geometry = [Point(*p) for p in node_xy]
            gdf_nodes = gp.GeoDataFrame(list(node_data), geometry=geometry)
            if "crs" in net.graph.keys():
                gdf_nodes.crs = net.graph["crs"]
            return gdf_nodes

        def lines_to_gdf(net, points, nodeID):
            starts, ends, edge_data = zip(*net.edges(data=True))
            gdf_edges = gp.GeoDataFrame(list(edge_data))

            if points is True:
                node_start = []
                node_end = []
                for s in starts:
                    node_start.append(net.nodes[s][nodeID])
                for e in ends:
                    node_end.append(net.nodes[e][nodeID])
                gdf_edges["node_start"] = node_start
                gdf_edges["node_end"] = node_end

            if "crs" in net.graph.keys():
                gdf_edges.crs = net.graph["crs"]

            return gdf_edges

        # End----------------Repeated code as those in Isochrone Map

        # Define locate nearest edge (kne) and projected point (pp)
        def find_kne(point, lines):
            dists = np.array(list(map(lambda l: l.distance(point), lines)))
            kne_pos = dists.argsort()[0]
            kne = lines.iloc[[kne_pos]]
            kne_idx = kne.index[0]
            return kne_idx, kne.values[0]

        def get_pp(point, line):
            """Get the projected point (pp) of 'point' on 'line'."""
            pp = line.interpolate(line.project(point))  # PP as a Point
            return pp

        # split line with projected points
        def split_line(line, pps):
            line = snap(line, pps, 1e-8)  # slow?
            try:
                new_lines = list(split(line, pps))  # split into segments
                return new_lines
            except TypeError as e:
                print("Error when splitting line: {}\n{}\n{}\n".format(e, line, pps))
                return []

        # for interpolation (split by pp): replicate old line
        def update_edges(edges, new_lines, replace):
            # for interpolation (split by pp): replicate old line
            if replace:
                # create a flattened gdf with all line segs and corresponding kne_idx
                kne_idxs = list(line_pps_dict.keys())
                lens = [len(item) for item in new_lines]
                new_lines_gdf = gp.GeoDataFrame(
                    {
                        "kne_idx": np.repeat(kne_idxs, lens),
                        "geometry": list(itertools.chain.from_iterable(new_lines)),
                    }
                )
                # merge to inherit the data of the replaced line
                cols = list(edges.columns)
                cols.remove("geometry")  # don't include the old geometry
                new_edges = new_lines_gdf.merge(
                    edges[cols], how="left", left_on="kne_idx", right_index=True
                )
                new_edges.drop("kne_idx", axis=1, inplace=True)
                new_lines = new_edges["geometry"]  # now a flatten list
                # for connection (to external poi): append new lines
            else:
                new_edges = gp.GeoDataFrame(
                    POI[[key_col]], geometry=new_lines, columns=[key_col, "geometry"]
                )

            # update features (a bit slow)
            new_edges["length"] = [l.length for l in new_lines]

            # remember to reindex to prevent duplication when concat
            start = edges.index[-1] + 1
            stop = start + len(new_edges)
            new_edges.index = range(start, stop)

            # for interpolation: remove existing edges
            if replace:
                edges = edges.drop(kne_idxs, axis=0)
            # for connection: filter invalid links
            else:
                valid_pos = np.where(new_edges["length"] <= threshold)[0]
                n = len(new_edges)
                n_fault = n - len(valid_pos)
                f_pct = n_fault / n * 100
                print(
                    "Remove faulty projections: {}/{} ({:.2f}%)".format(
                        n_fault, n, f_pct
                    )
                )
                new_edges = new_edges.iloc[valid_pos]  # use 'iloc' here

            # merge new edges
            dfs = [edges, new_edges]
            edges = gp.GeoDataFrame(
                pd.concat(dfs, ignore_index=False, sort=False), crs=dfs[0].crs
            )

            # all edges, newly added edges only
            return edges, new_edges

        # Difine function of nearest points
        def ckdnearest(gdA, gdB):

            nA = np.array(list(gdA.geometry.apply(lambda x: (x.x, x.y))))
            nB = np.array(list(gdB.geometry.apply(lambda x: (x.x, x.y))))
            btree = cKDTree(nB)
            dist, idx = btree.query(nA, k=1)
            gdB_nearest = gdB.iloc[idx].reset_index(drop=True)
            gdf = pd.concat(
                [
                    gdA.drop(columns="geometry").reset_index(drop=True),
                    gdB_nearest,
                    pd.Series(dist, name="dist"),
                ],
                axis=1,
            )
            return gdf

        # Cross join Data
        O_gdf = input1.to_pandas().rename(columns={self.O_geo_col: "geometry"})
        D_gdf = input2.to_pandas().rename(columns={self.D_geo_col: "geometry"})
        R_gdf = input3.to_pandas().rename(columns={self.R_geo_col: "geometry"})
        O_gdf = gp.GeoDataFrame(O_gdf, geometry="geometry")
        D_gdf = gp.GeoDataFrame(D_gdf, geometry="geometry")
        R_gdf = gp.GeoDataFrame(R_gdf, geometry="geometry")
        R_gdf = R_gdf.rename(columns={self.R_speed_col: "speed"})
        # Set a lat\Lon CRS

        crsinput = CRS.from_user_input(R_gdf.crs)
        if crsinput.is_geographic:
            logging.warning("Unit as Degree, Please use Projected CRS")
            R_gdf = R_gdf.to_crs(3857)

        O_gdf = O_gdf.to_crs(R_gdf.crs)
        D_gdf = D_gdf.to_crs(R_gdf.crs)
        O_gdf["geometry"] = O_gdf.geometry.centroid
        D_gdf["geometry"] = D_gdf.geometry.centroid

        O_gdf["key"] = list(range(1, (O_gdf.shape[0] + 1)))
        gdfmax = O_gdf.key.max()
        O_gdf["Category"] = "Origin"
        D_gdf["key"] = list(range((gdfmax + 1), (gdfmax + 1 + D_gdf.shape[0])))
        D_gdf["Category"] = "Destination"

        O_gdf = O_gdf[["geometry", "key", "Category"]]
        D_gdf = D_gdf[["geometry", "key", "Category"]]
        POI = pd.concat([O_gdf, D_gdf], ignore_index=True)
        # result = O_gdf.append(D_gdf)
        # result = result.reset_index(drop=True)
        # POI = result

        # Convert geoDataFrame to Graph with MomePy
        gdfroad = R_gdf
        graph = gdf_to_nx(gdfroad)
        nodes, edges = nx_to_gdf(graph, points=True, lines=True, spatial_weights=False)
        graph = gdf_to_nx(edges)

        # build rtree
        Rtree = rtree.index.Index()
        [Rtree.insert(fid, geom.bounds) for fid, geom in edges["geometry"].items()]

        knn = 5
        key_col = "key"
        threshold = 100000

        # locate nearest edge (kne) and projected point (pp)
        # Projecting POIs to the network...
        POI["near_idx"] = [
            list(Rtree.nearest(point.bounds, knn)) for point in POI["geometry"]
        ]  # slow
        POI["near_lines"] = [
            edges["geometry"][near_idx] for near_idx in POI["near_idx"]
        ]  # very slow
        POI["kne_idx"], knes = zip(
            *[
                find_kne(point, near_lines)
                for point, near_lines in zip(POI["geometry"], POI["near_lines"])
            ]
        )  # slow
        POI["pp"] = [get_pp(point, kne) for point, kne in zip(POI["geometry"], knes)]

        # 08-1: update internal edges (split line segments)
        line_pps_dict = {
            k: MultiPoint(list(v)) for k, v in POI.groupby(["kne_idx"])["pp"]
        }
        new_lines = [
            split_line(edges["geometry"][idx], pps)
            for idx, pps in line_pps_dict.items()
        ]  # bit slow
        edges, _ = update_edges(edges, new_lines, replace=True)

        edges["length"] = edges.length.astype(float)
        edges["time"] = edges.length / edges.speed

        # Convert geoDataFrame to Graph with MomePy
        graph = gdf_to_nx(edges)

        # Save geoDataFrame back to Points and Edges
        nodes1, edges1 = nx_to_gdf(
            graph, points=True, lines=True, spatial_weights=False
        )
        # Difine function of nearest points
        snappoint = ckdnearest(POI, nodes1)
        # set geommetry to Snapped nodes
        snappoint = snappoint.set_geometry("geometry")
        # edges = edges.reset_index(drop=True)
        snappoint = snappoint.reset_index(drop=True)
        x = snappoint[snappoint["Category"] == "Origin"]
        y = snappoint[snappoint["Category"] != "Origin"]
        dfx = x[["key", "Category", "nodeID", "dist"]]
        dfy = y[["key", "Category", "nodeID", "dist"]]
        dff = pd.DataFrame()
        for i in range(len(dfx)):
            a = dfx.iloc[i]
            b = dfy
            i1 = a["nodeID"]
            dist = []
            lent = []
            all_dist = nx.single_source_dijkstra_path_length(
                graph, list(graph.nodes)[i1], weight="time"
            )
            all_length = nx.single_source_dijkstra_path_length(
                graph, list(graph.nodes)[i1], weight="length"
            )
            for j in range(len(b)):
                i2 = b.iloc[j]["nodeID"]
                dist.append(all_dist[list(graph.nodes)[i2]])
                lent.append(all_length[list(graph.nodes)[i2]])
            data = {
                "originid": a["key"],
                "destinationid": b["key"],
                "duration": dist,
                "distance": lent,
                "snap_o_dist": a["dist"],
                "snap_d_dist": b["dist"],
            }
            df = pd.DataFrame(data)
            dff = pd.concat([df, dff], ignore_index=True)

        dff["destinationid"] = dff.destinationid - gdfmax
        dff.sort_values(by=["originid", "destinationid"], inplace=True)
        dff = dff.reset_index(drop=True)
        return knext.Table.from_pandas(dff)


############################################
# Isochrone map
############################################
@knext.node(
    name="Isochrone Map",
    node_type=knext.NodeType.MANIPULATOR,
    # node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "OSMisochrone.png",
)
@knext.input_table(
    name="Input Table as Center",
    description="Input table with geometry.",
)
@knext.input_table(
    name="Input Table as road network",
    description="Input road network with LineString geometry.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with isochrone geometry.",
)
class IsochroneMap:
    """
    This node can calculate the isochrone map for the input point  based on the input road network and its travel cost column.
    It first snaps input point to the road work, and uses the function
    [ego_graph](https://networkx.org/documentation/stable/reference/generated/networkx.generators.ego.ego_graph.html)
    in [NetworkX](https://networkx.org/) to isochrone for the weighted (time or distance) graph.

    The input value for interval list should be cautious, it should be in a reasonable boundary of the road network.

    The output table contains the two column, isochrone intervals and geometry.The calcualtion depends on a projected coordinates
    systemn of road network datasets. If it doesn't have a projected CRS information, a default CRS"3857" will be applied.
    If the input geometry is not point feature, the centroid will be used. If it contains multiple rows, the total centroid will be applied.
    """

    # input parameters
    c_geo_col = knext.ColumnParameter(
        "Point Geometry Column",
        "Select the geometry column as center.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    r_geo_col = knext.ColumnParameter(
        "Road Network geometry column",
        "Select the geometry column as road network.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    r_cost_col = knext.ColumnParameter(
        "Travel cost column from Road Network",
        "Select the travel cost column from road network.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    isolist = knext.StringParameter(
        "Set isochone interval list",
        "Input a interval list separated by comma ",
        "",
    )

    def configure(
        self,
        configure_context,
        input_schema_1,
        input_schema_2,
    ):
        self.c_geo_col = knut.column_exists_or_preset(
            configure_context, self.c_geo_col, input_schema_1, knut.is_geo
        )
        self.r_geo_col = knut.column_exists_or_preset(
            configure_context, self.r_geo_col, input_schema_2, knut.is_geo
        )

        return knext.Schema.from_columns(
            [
                knext.Column(knut.TYPE_MULTI_POLYGON, "geometry"),
                knext.Column(knext.int64(), "isochrone"),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input1, input2):

        import pandas as pd
        import networkx as nx
        from shapely.geometry import MultiPoint, LineString, Point, Polygon
        import numpy as np
        from pyproj import CRS  # For CRS Units check
        import logging
        import osmnx as ox

        # Start----------------Repeated code as those in Street Network Matrix
        # Use momepy function: generate_primal,gdf_to_nx,primal_to_gdf,nx_to_gdf
        def generate_primal(G, gdf_network, fields, multigraph):
            G.graph["approach"] = "primal"
            key = 0
            for row in gdf_network.itertuples():

                first = row.geometry.coords[0]
                last = row.geometry.coords[-1]

                data = [r for r in row][1:]
                attributes = dict(zip(fields, data))
                if multigraph:
                    G.add_edge(first, last, key=key, **attributes)
                    key += 1
                else:
                    G.add_edge(first, last, **attributes)

        def gdf_to_nx(
            gdf_network,
            approach="primal",
            length="mm_len",
            multigraph=True,
            directed=False,
            angles=True,
            angle="angle",
        ):
            gdf_network = gdf_network.copy()
            if "key" in gdf_network.columns:
                gdf_network.rename(columns={"key": "__key"}, inplace=True)
            import networkx as nx

            if multigraph and directed:
                net = nx.MultiDiGraph()
            elif multigraph and not directed:
                net = nx.MultiGraph()
            elif not multigraph and directed:
                net = nx.DiGraph()
            else:
                net = nx.Graph()

            net.graph["crs"] = gdf_network.crs
            gdf_network[length] = gdf_network.geometry.length
            fields = list(gdf_network.columns)
            generate_primal(net, gdf_network, fields, multigraph)

            return net

        def primal_to_gdf(net, points, lines, spatial_weights, nodeID):
            if points is True:
                gdf_nodes = points_to_gdf(net)

                if spatial_weights is True:
                    W = libpysal.weights.W.from_networkx(net)
                    W.transform = "b"

            if lines is True:
                gdf_edges = lines_to_gdf(net, points, nodeID)

            if points is True and lines is True:
                if spatial_weights is True:
                    return gdf_nodes, gdf_edges, W
                return gdf_nodes, gdf_edges
            if points is True and lines is False:
                if spatial_weights is True:
                    return gdf_nodes, W
                return gdf_nodes
            return gdf_edges

        def nx_to_gdf(
            net, points=True, lines=True, spatial_weights=False, nodeID="nodeID"
        ):
            # generate nodes and edges geodataframes from graph
            primal = True
            for nid, n in enumerate(net):
                net.nodes[n][nodeID] = nid
            return primal_to_gdf(
                net,
                points=points,
                lines=lines,
                spatial_weights=spatial_weights,
                nodeID=nodeID,
            )

        def points_to_gdf(net):
            node_xy, node_data = zip(*net.nodes(data=True))
            if isinstance(node_xy[0], int) and "x" in node_data[0].keys():
                geometry = [
                    Point(data["x"], data["y"]) for data in node_data
                ]  # osmnx graph
            else:
                geometry = [Point(*p) for p in node_xy]
            gdf_nodes = gp.GeoDataFrame(list(node_data), geometry=geometry)
            if "crs" in net.graph.keys():
                gdf_nodes.crs = net.graph["crs"]
            return gdf_nodes

        def lines_to_gdf(net, points, nodeID):
            starts, ends, edge_data = zip(*net.edges(data=True))
            gdf_edges = gp.GeoDataFrame(list(edge_data))

            if points is True:
                node_start = []
                node_end = []
                for s in starts:
                    node_start.append(net.nodes[s][nodeID])
                for e in ends:
                    node_end.append(net.nodes[e][nodeID])
                gdf_edges["node_start"] = node_start
                gdf_edges["node_end"] = node_end

            if "crs" in net.graph.keys():
                gdf_edges.crs = net.graph["crs"]

            return gdf_edges

        # End----------------Repeated code as those in Street Network Matrix

        def gdf_to_osmgraph(gdf):
            graph = gdf_to_nx(gdf)
            nodes, edges = nx_to_gdf(
                graph, points=True, lines=True, spatial_weights=False
            )
            nodes["x"] = nodes.geometry.x
            nodes["y"] = nodes.geometry.y
            nodes = nodes.rename(columns={"nodeID": "osmid"})
            edges = edges.rename(columns={"node_start": "u", "node_end": "v"})
            edges["k"] = list(range(1, (edges.shape[0] + 1)))
            edges = edges.set_index(["u", "v", "k"])
            G = ox.utils_graph.graph_from_gdfs(nodes, edges)
            return G

        # Cross join Data
        c_gdf = input1.to_pandas().rename(columns={self.c_geo_col: "geometry"})
        r_gdf = input2.to_pandas().rename(columns={self.r_geo_col: "geometry"})
        c_gdf = gp.GeoDataFrame(c_gdf, geometry="geometry")
        r_gdf = gp.GeoDataFrame(r_gdf, geometry="geometry")
        r_gdf = r_gdf.rename(columns={self.r_cost_col: "time"})
        # Set a lat\Lon CRS

        crsinput = CRS.from_user_input(r_gdf.crs)
        if crsinput.is_geographic:
            logging.warning("Unit as Degree, Please use Projected CRS")
            r_gdf = r_gdf.to_crs(3857)

        c_gdf = c_gdf.to_crs(r_gdf.crs)
        if c_gdf.shape[0] > 1:
            logging.warning("Only the total centroid was applied")
            c_gdf = gp.GeoDataFrame(geometry=gp.GeoSeries(c_gdf.unary_union.centroid))
        c_gdf["geometry"] = c_gdf.geometry.centroid

        # This example script simply outputs the node's input table.
        graph = gdf_to_osmgraph(r_gdf)
        nodes, edges = nx_to_gdf(graph, points=True, lines=True, spatial_weights=False)

        nearest_node = gp.sjoin_nearest(c_gdf, nodes)
        center_node = nearest_node["osmid"][0]

        trip_times = list(map(int, self.isolist.split(",")))
        trip_times = sorted(trip_times)

        isochrone_polys = []
        for trip_time in sorted(trip_times):
            subgraph = nx.ego_graph(
                graph, center_node, radius=trip_time, undirected=True, distance="time"
            )
            edgex = nx_to_gdf(subgraph, points=False, lines=True, spatial_weights=False)
            new_iso = Polygon(edgex.unary_union.buffer(50).exterior)
            isochrone_polys.append(new_iso)

        gdfx = gp.GeoDataFrame(geometry=gp.GeoSeries(isochrone_polys))
        gdfx["isochrone"] = trip_times
        if len(trip_times) > 1:
            for i in range(1, len(trip_times)):
                k = i - 1
                c0 = isochrone_polys[k]
                c1 = isochrone_polys[i]
                cd = c1.difference(c0)
                gdfx.at[i, "geometry"] = cd

        gdfx.crs = r_gdf.crs
        gdfx = gdfx.reset_index(drop=True)
        return knext.Table.from_pandas(gdfx)