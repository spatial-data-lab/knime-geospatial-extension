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

# short functions from momepy
class SimpleMomepy:
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
        import networkx as nx

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
        SimpleMomepy.generate_primal(net, gdf_network, fields, multigraph)

        return net

    def primal_to_gdf(net, points, lines, spatial_weights, nodeID):
        import libpysal

        if points is True:
            gdf_nodes = SimpleMomepy.points_to_gdf(net)

            if spatial_weights is True:
                W = libpysal.weights.W.from_networkx(net)
                W.transform = "b"

        if lines is True:
            gdf_edges = SimpleMomepy.lines_to_gdf(net, points, nodeID)

        if points is True and lines is True:
            if spatial_weights is True:
                return gdf_nodes, gdf_edges, W
            return gdf_nodes, gdf_edges
        if points is True and lines is False:
            if spatial_weights is True:
                return gdf_nodes, W
            return gdf_nodes
        return gdf_edges

    def nx_to_gdf(net, points=True, lines=True, spatial_weights=False, nodeID="nodeID"):
        # generate nodes and edges geodataframes from graph
        primal = True
        for nid, n in enumerate(net):
            net.nodes[n][nodeID] = nid
        return SimpleMomepy.primal_to_gdf(
            net,
            points=points,
            lines=lines,
            spatial_weights=spatial_weights,
            nodeID=nodeID,
        )

    def points_to_gdf(net):
        from shapely.geometry import Point

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


############################################
# Google Distance Matrix
############################################


class _GoogleTravelMode(knext.EnumParameterOptions):
    BICYCLING = (
        "Bicycling",
        """Requests bicycling directions or distance via bicycle paths & preferred streets (where available).
        Bicycling directions may sometimes not include clear bicycling paths.""",
    )
    DRIVING = (
        "Driving",
        "Indicates standard driving directions or distance using the road network.",
    )
    TRANSIT = (
        "Transit",
        """Requests directions or distance via public transit routes (where available). 
        The used departure time is now.""",
    )
    WALKING = (
        "Walking",
        """Requests walking directions or distance via pedestrian paths & sidewalks (where available).
        Walking directions may sometimes not include clear pedestrian paths.""",
    )

    @classmethod
    def get_default(cls):
        return cls.DRIVING


@knext.node(
    name="Google Distance Matrix",
    node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "GoogleDistMatrix.png",
)
@knext.input_table(
    name="Input table with origins",
    description="Input table with origin geometry and ID column.",
)
@knext.input_table(
    name="Input table with destinations",
    description="Input table with destination geometry and ID column.",
)
@knext.output_table(
    name="Output table",
    description="""Output table with the selected origin and destination ID columns and the corresponding travel costs
    in minutes and meters.""",
)
class GoogleDistanceMatrix:
    """
    This node uses the
    [Google Distance Matrix API](https://developers.google.com/maps/documentation/distance-matrix/overview)
    to create a distance matrix for the provided origins and destinations. The matrix is created by pairing each
    input origin with each input destination and will contain the travel distance and time for each pair.
    The distance unit is meter and the estimated travel time is returned in minutes.

    If the input geometry is not a point geometry, the centroids will be automatically computed and used.
    """

    o_geo_col = knext.ColumnParameter(
        "Origin geometry column",
        "Select the geometry column that describes the origins.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    o_id_col = knext.ColumnParameter(
        "Origin ID column",
        """Select the column which contains for each origin a unique ID. The selected column will be returned
        in the result table and can be used to link back to the original data.""",
        port_index=0,
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    d_geo_col = knext.ColumnParameter(
        "Destination geometry column",
        "Select the geometry column that describes the destinations.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    d_id_col = knext.ColumnParameter(
        "Destination ID column",
        """Select the column which contains for each destination a unique ID. The selected column will be returned
        in the result table and can be used to link back to the original data.""",
        port_index=1,
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )
    api_key = knext.StringParameter(
        "Google API key",
        """Click [here](https://developers.google.com/maps/documentation/distance-matrix/get-api-key) for details on
        how to obtain and use a Google API key for the Distance Matrix API.""",
        "",
    )
    travel_mode = knext.EnumParameter(
        "Travel mode",
        """The following 
        [travel modes](https://developers.google.com/maps/documentation/distance-matrix/distance-matrix#mode) 
        are supported: """,
        default_value=_GoogleTravelMode.get_default().name,
        enum=_GoogleTravelMode,
    )
    # Constant for distance matrix
    _COL_DURATION = "duration"
    _COL_DISTANCE = "distance"
    _BASE_URL = "https://maps.googleapis.com/maps/api/distancematrix/json?language=en&units=imperial&origins={0}&destinations={1}&key={2}&mode={3}"

    def configure(self, configure_context, o_schema, d_schema):
        self.o_geo_col = knut.column_exists_or_preset(
            configure_context, self.o_geo_col, o_schema, knut.is_geo
        )
        knut.column_exists(self.o_id_col, o_schema)
        o_id_type = o_schema[self.o_id_col].ktype

        self.d_geo_col = knut.column_exists_or_preset(
            configure_context, self.d_geo_col, d_schema, knut.is_geo
        )
        knut.column_exists(self.d_id_col, d_schema)
        d_id_type = d_schema[self.d_id_col].ktype

        return knext.Schema(
            [o_id_type, d_id_type, knext.double(), knext.int64()],
            [self.o_id_col, self.d_id_col, self._COL_DURATION, self._COL_DISTANCE],
        )

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        # define function to derive travel time and distance from Google Maps API
        import urllib.request as urllib2  # for Google Drive
        import json  # for OSRM

        knut.check_canceled(exec_context)

        def fetch_google_od(download_link):
            nt_duration = 0
            nt_distance = 0
            try:
                req = urllib2.urlopen(download_link)
                jsonout = json.loads(req.read())
                nt_duration = jsonout["rows"][0]["elements"][0]["duration"][
                    "value"
                ]  # meters
                nt_distance = jsonout["rows"][0]["elements"][0]["distance"][
                    "value"
                ]  # seconds

                # transform seconds to minutes
                nt_duration = round(nt_duration / 60, 2)
                nt_distance = round(nt_distance, 2)
            except Exception as err:
                knut.LOGGER.warning(f"Exception while calling Google Matrix API: {err}")
            return [nt_duration, nt_distance]

        o_gdf = knut.load_geo_data_frame(left_input, self.o_geo_col, exec_context)
        d_gdf = knut.load_geo_data_frame(right_input, self.d_geo_col, exec_context)
        o_gdf = knut.validify_id_column(o_gdf, self.o_id_col)
        d_gdf = knut.validify_id_column(d_gdf, self.d_id_col)

        # Set a lat\Lon CRS before renaming the geometry column
        o_gdf = o_gdf.to_crs(4326)
        d_gdf = d_gdf.to_crs(4326)

        # Filter all columns except the needed once and rename the geometry column to have a consistent result in merge
        o_gdf = o_gdf.filter(items=[self.o_geo_col, self.o_id_col]).rename(
            columns={self.o_geo_col: "geometry"}
        )
        d_gdf = d_gdf.filter(items=[self.d_geo_col, self.d_id_col]).rename(
            columns={self.d_geo_col: "geometry"}
        )

        # Generate origin destination matrix via cross join
        merge_df = o_gdf.merge(d_gdf, how="cross")
        merge_df_x = gp.GeoDataFrame(geometry=merge_df["geometry_x"])
        merge_df_y = gp.GeoDataFrame(geometry=merge_df["geometry_y"])

        # create the result matrix with the two id columns...
        distance_matrix = merge_df[[self.o_id_col, self.d_id_col]]
        # ... and default value 0 for the duration and distance column
        distance_matrix[self._COL_DURATION] = 0
        distance_matrix[self._COL_DISTANCE] = 0
        # can process any geometry since it computes the centroid first
        for i in range(0, len(distance_matrix)):
            origins_lat_lon = "{},{}".format(
                merge_df_x.centroid.y[i], merge_df_x.centroid.x[i]
            )
            destinations_lat_lon = "{},{}".format(
                merge_df_y.centroid.y[i], merge_df_y.centroid.x[i]
            )
            google_request_link = self._BASE_URL.format(
                origins_lat_lon,
                destinations_lat_lon,
                self.api_key,
                self.travel_mode.lower(),
            )
            google_travel_cost = fetch_google_od(google_request_link)
            # add duration result in minutes
            distance_matrix.iat[i, 2] = google_travel_cost[0]
            # add distance result in meters
            distance_matrix.iat[i, 3] = google_travel_cost[1]
        return knut.to_table(distance_matrix)


############################################
# OSRM
############################################


class _OSRMResultModel(knext.EnumParameterOptions):
    ROUTE = (
        "Route",
        "Returns only the travel route.",
    )
    TRAVEL = (
        "Travel cost",
        "Returns the drive distance in meters and travel time in minutes.",
    )
    TRAVEL_ROUTE = (
        "Travel cost and route",
        "Returns the drive distance in meters and travel time in minutes as well as the travel route.",
    )

    @classmethod
    def get_default(cls):
        return cls.TRAVEL

    def append_route(self) -> bool:
        return self is not _OSRMResultModel.TRAVEL

    def append_distance(self) -> bool:
        return self is not _OSRMResultModel.ROUTE


@knext.node(
    name="OSRM Distance Matrix",
    node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "OSRMdistMatrix.png",
)
@knext.input_table(
    name="Input table with origins",
    description="Input table with origin geometry and ID column.",
)
@knext.input_table(
    name="Input table with destinations",
    description="Input table with destination geometry and ID column.",
)
@knext.output_table(
    name="Output Table",
    description="""Output table with the selected origin and destination ID columns and the corresponding travel costs
    in minutes and meters as well as the travel route.""",
)
class OSRMDistanceMatrix:
    """
    This node uses the [Open Source Routing Machine (OSRM)](https://project-osrm.org/) to create a distance matrix
    for the provided origins and destinations. The matrix is created by pairing each input origin with each input
    destination and will contain the driving travel distance and time as well as the
    [route](http://project-osrm.org/docs/v5.5.1/api/?language=Python#route-service) for each pair.
    The travel distance unit is meter and the estimated drive time is returned in minutes.

    OSRM is a C++ implementation of a high-performance routing engine for shortest paths in road networks.
    It combines sophisticated routing algorithms with the open and free road network data of the
    [OpenStreetMap (OSM) project.](https://www.openstreetmap.org/about)

    If the input geometry is not a point geometry, the centroids will be automatically computed and used.

    #Usage Policy
    The general usage policy and limitations of the OSRM service that is used by this node can be found
    [here.](https://github.com/Project-OSRM/osrm-backend/wiki/Api-usage-policy) The current demo server is hosted
    by [FOSSGIS](https://www.fossgis.de/) with the following
    [terms of use (German).](https://www.fossgis.de/arbeitsgruppen/osm-server/nutzungsbedingungen/)

    #Copyright
    Data provided by [OpenStreetMap](https://www.openstreetmap.org/copyright)
    [(ODbl)](https://opendatacommons.org/licenses/odbl/index.html) under
    [CC-BY-SA](https://creativecommons.org/licenses/by-sa/2.0/).
    To report a problem and contribute to OpenStreetMap click [here.](https://www.openstreetmap.org/fixthemap)
    """

    # input parameters
    o_geo_col = knext.ColumnParameter(
        "Origin geometry column",
        "Select the geometry column that describes the origins.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    o_id_col = knext.ColumnParameter(
        "Origin ID column",
        """Select the column which contains for each origin a unique ID. The selected column will be returned
        in the result table and can be used to link back to the original data.""",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    d_geo_col = knext.ColumnParameter(
        "Destination geometry column",
        "Select the geometry column that describes the destinations.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    d_id_col = knext.ColumnParameter(
        "Destination ID column",
        """Select the column which contains for each destination a unique ID. The selected column will be returned
        in the result table and can be used to link back to the original data.""",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )
    result_model = knext.EnumParameter(
        label="Result model",
        description="Supports the following result models:",
        default_value=_OSRMResultModel.get_default().name,
        enum=_OSRMResultModel,
    )

    # Constant for distance matrix
    _COL_DURATION = "duration"
    _COL_DISTANCE = "distance"
    _COL_GEOMETRY = "route"

    # For details see: http://project-osrm.org/docs/v5.5.1/api/#route-service
    _BASE_URL = "https://router.project-osrm.org"
    # only supports car as profile: https://github.com/Project-OSRM/osrm-backend/issues/4034
    _PROFILE = "driving"
    # For details see: http://project-osrm.org/docs/v5.5.1/api/#route-service
    _REQUEST_PARAMETER = {"continue_straight": "false"}

    # number of pairs send per request
    _BATCH_SIZE = 50
    # Number of seconds to wait after each request
    _REQUEST_DELAY = 1
    # Request timeout
    # do not send more than 1 request per second https://github.com/Project-OSRM/osrm-backend/wiki/Demo-server
    _REQUEST_TIMEOUT = None

    def configure(self, configure_context, o_schema, d_schema):
        self.o_geo_col = knut.column_exists_or_preset(
            configure_context, self.o_geo_col, o_schema, knut.is_geo
        )
        knut.column_exists(self.o_id_col, o_schema)
        o_id_type = o_schema[self.o_id_col].ktype

        self.d_geo_col = knut.column_exists_or_preset(
            configure_context, self.d_geo_col, d_schema, knut.is_geo
        )
        knut.column_exists(self.d_id_col, d_schema)
        d_id_type = d_schema[self.d_id_col].ktype

        model = _OSRMResultModel[self.result_model]

        result_types = [o_id_type, d_id_type]
        result_names = [self.o_id_col, self.d_id_col]
        if model.append_distance():
            result_types += [knext.double(), knext.double()]
            result_names += [self._COL_DURATION, self._COL_DISTANCE]
        if model.append_route():
            result_types += [knut.TYPE_LINE]
            result_names += [self._COL_GEOMETRY]
        # check which model is selected
        return knext.Schema(
            result_types,
            result_names,
        )

    # set digits for coordinates
    def round_coord_list(self, coord_list, digits):
        coord_list = list(map(lambda x: [round(i, digits) for i in x], coord_list))
        coord_list = [tuple(i) for i in coord_list]
        return coord_list

    # extract route geometry
    def extract_route(self, data):

        import polyline
        from shapely.geometry import LineString

        decode_line = polyline.decode(data["routes"][0]["geometry"])
        # Extract the location coordinates from the 'waypoints' field
        coordinates = [
            round(coord, 5)
            for waypoint in data["waypoints"]
            for coord in waypoint["location"]
        ]
        points = [
            (coordinates[i + 1], coordinates[i]) for i in range(0, len(coordinates), 2)
        ]
        decode_line4 = self.round_coord_list(decode_line, 4)
        points4 = self.round_coord_list(points, 4)
        indexes = []
        tag = 0
        for i in points4:
            newline = decode_line4[tag:]
            for j, p in enumerate(newline):
                if i == p:
                    tag = j + tag
                    break
            indexes.append(tag)
            tag = tag + 1
        re_decode_line = [(y, x) for x, y in decode_line]
        routes = [
            LineString(re_decode_line[indexes[i] : (indexes[(i + 1)] + 1)])
            for i in range(0, len(indexes), 2)
        ]
        return routes

    # update travel cost and route geometry for the given batch of the given dataframe
    def update_part(self, model: _OSRMResultModel, df, ns, ne):
        import requests
        import time
        import json
        import pandas as pd

        df_batch = df.copy().loc[ns:ne]
        df_batch = df_batch[["StartX", "StartY", "EndX", "EndY"]]
        df_batch = df_batch.astype(str)
        df_batch["period"] = (
            df_batch["StartX"]
            + ","
            + df_batch["StartY"]
            + ";"
            + df_batch["EndX"]
            + ","
            + df_batch["EndY"]
        )
        # http://project-osrm.org/docs/v5.5.1/api/#route-service
        coordinate_query_list = [";".join(df_batch["period"])]
        request_url = (
            self._BASE_URL
            + "/route/v1/"
            + self._PROFILE
            + "/"
            + coordinate_query_list[0]
        )

        try:
            r = requests.get(
                request_url,
                params=self._REQUEST_PARAMETER,
                headers=knut.WEB_REQUEST_HEADER,
                timeout=self._REQUEST_TIMEOUT,
            )
            time.sleep(self._REQUEST_DELAY)
            data = json.loads(r.text)
            if data["code"] == "Ok":
                if model.append_distance():
                    dfr = pd.DataFrame(data["routes"][0]["legs"])[
                        [self._COL_DURATION, self._COL_DISTANCE]
                    ].iloc[::2]
                    # convert seconds to minutes
                    dfr.duration /= 60
                    df.loc[ns:ne, self._COL_DURATION] = dfr.duration.to_list()
                    df.loc[ns:ne, self._COL_DISTANCE] = dfr.distance.to_list()
                if model.append_route():
                    # get route
                    temp_route = self.extract_route(data)
                    # get route
                    if len(temp_route) == 1:
                        df.loc[ns:ne, self._COL_GEOMETRY] = temp_route[0]
                    else:
                        df.loc[ns:ne, self._COL_GEOMETRY] = temp_route
            else:
                knut.LOGGER.warning(f"No route found from:{ns} to :{ne}")
        except Exception as err:
            knut.LOGGER.warning(f"Error finding route from:{ns} to :{ne}. Error: {err}")

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):

        from shapely.geometry import LineString

        # Cross join Data
        o_gdf = knut.load_geo_data_frame(left_input, self.o_geo_col, exec_context)
        d_gdf = knut.load_geo_data_frame(right_input, self.d_geo_col, exec_context)
        o_gdf = knut.validify_id_column(o_gdf, self.o_id_col)
        d_gdf = knut.validify_id_column(d_gdf, self.d_id_col)

        # Set a lat\Lon CRS before renaming the geometry column
        o_gdf = o_gdf.to_crs(4326)
        d_gdf = d_gdf.to_crs(4326)

        # Filter all columns except the needed once and rename the geometry column to have a consistent result in merge
        o_gdf = o_gdf.filter(items=[self.o_geo_col, self.o_id_col]).rename(
            columns={self.o_geo_col: "geometry"}
        )
        d_gdf = d_gdf.filter(items=[self.d_geo_col, self.d_id_col]).rename(
            columns={self.d_geo_col: "geometry"}
        )

        # Generate origin destination matrix via cross join
        merge_df = o_gdf.merge(d_gdf, how="cross")
        merge_df_x = gp.GeoDataFrame(geometry=merge_df["geometry_x"])
        merge_df_y = gp.GeoDataFrame(geometry=merge_df["geometry_y"])
        df = merge_df[[self.o_id_col, self.d_id_col]]
        df["StartX"] = merge_df_x.centroid.x
        df["StartY"] = merge_df_x.centroid.y
        df["EndX"] = merge_df_y.centroid.x
        df["EndY"] = merge_df_y.centroid.y
        df = df.reset_index(drop=True)
        # compute the batches
        n_length = df.shape[0]
        n_loop = n_length // self._BATCH_SIZE
        n_tail = n_length % self._BATCH_SIZE

        model = _OSRMResultModel[self.result_model]

        if model.append_distance():
            df[self._COL_DURATION] = 0.0
            df[self._COL_DISTANCE] = 0.0
        if model.append_route():
            df[self._COL_GEOMETRY] = LineString([(0, 0), (1, 1)])

        # loop over the different batches
        if n_loop >= 1:
            for i in range(n_loop):
                ns = self._BATCH_SIZE * i
                ne = ns + self._BATCH_SIZE - 1
                self.update_part(model, df, ns, ne)
                # i starts with 0
                process_counter = i + 1
                exec_context.set_progress(
                    0.9 * process_counter / n_loop,
                    f"Batch {process_counter} of {n_loop} processed",
                )
                knut.check_canceled(exec_context)
        # process the remaining rows
        if n_tail > 0:
            exec_context.set_progress(0.95, "Processing left over batch")
            knut.check_canceled(exec_context)
            ns = self._BATCH_SIZE * n_loop
            ne = ns + n_tail - 1
            self.update_part(model, df, ns, ne)

        # remove the origin and destination columns
        rdf = df.loc[:, ~df.columns.isin(["StartX", "StartY", "EndX", "EndY"])]
        if model.append_route():
            gdf = gp.GeoDataFrame(rdf, geometry=self._COL_GEOMETRY, crs=4326)
        else:
            gdf = rdf

        return knut.to_table(gdf, exec_context)


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
    This node can calculate the travel time, and distance between the origin (left or upper input port)
    and destination (right or lower input port) points based on the graph derived from the input road network and its speed column.
    It first snaps the origin and destination points to the road work and uses the function
    [single_source_dijkstra_path_length](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.algorithms.shortest_paths.weighted.single_source_dijkstra_path_length.html)
    in [NetworkX](https://networkx.org/) to compute the shortest path length between source and all other reachable nodes for the weighted (time or distance) graph.

    The column for speed value in the road data is used to calculate travel time with the formula "length/speed", so the speed value should better be in meters/minute or meters/second.
    E.g., 1 mile/hour = 26.8224 meters/minute, 1 Kilometer/hour = 16.6667 meters/minute.

    The output table contains the distance between the original points and snap points along the road network, which can also be further calculated as additional travel time or distance with user-defined speed.
    The calculation depends on a projected coordinates system of road network datasets. If it is not in a projected CRS, a default CRS"3857" will be applied.
    If the input geometry is not a point feature, the centroids will be used.
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

    O_id_col = knext.ColumnParameter(
        "Unique origin ID Column",
        "Select the ID column for origin points.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_numeric_or_string,
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
    D_id_col = knext.ColumnParameter(
        "Unique destination ID column",
        "Select the ID  column for destination points.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_numeric_or_string,
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
        port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    speedunit = knext.StringParameter(
        label="Speed unit",
        description="The unit for speed column",
        default_value="Meters/second(m/s)",
        enum=[
            "Meter/second(m/s)",
            "Miles/hour(mph)",
            "kilometers/hour(km/h)",
            "None",
        ],
    )
    # Constant for distance matrix
    _OID = "originid"
    _DID = "destinationid"
    _TIMECOST = "duration"
    _DISTCOST = "distance"
    _OSNAP = "snap_o_dist"
    _DSNAP = "snap_d_dist"

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
        return None

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
        Or_gdf = knut.load_geo_data_frame(input1, self.O_geo_col, exec_context)
        De_gdf = knut.load_geo_data_frame(input2, self.D_geo_col, exec_context)
        R_gdf = knut.load_geo_data_frame(input3, self.R_geo_col, exec_context)
        # Valide Id column
        Or_gdf = knut.validify_id_column(Or_gdf, self.O_id_col)
        De_gdf = knut.validify_id_column(De_gdf, self.D_id_col)
        # create new data with key columns
        O_gdf = Or_gdf[[self.O_geo_col]].rename(columns={self.O_geo_col: "geometry"})
        D_gdf = De_gdf[[self.D_geo_col]].rename(columns={self.D_geo_col: "geometry"})
        R_gdf = R_gdf[[self.R_geo_col, self.R_speed_col]].rename(
            columns={self.R_geo_col: "geometry", self.R_speed_col: "speed"}
        )
        R_gdf = R_gdf.dropna(subset=["geometry"], how="any")
        if self.speedunit == "Meter/second(m/s)":
            R_gdf["speed"] = R_gdf.speed * 60
        elif self.speedunit == "Miles/hour(mph)":
            R_gdf["speed"] = R_gdf.speed * 26.8224
        elif self.speedunit == "kilometers/hour(km/h)":
            R_gdf["speed"] = R_gdf.speed * 16.6667
        else:
            R_gdf["speed"] = R_gdf.speed * 1
        # Set a lat\Lon CRS
        crsinput = CRS.from_user_input(R_gdf.crs)
        if crsinput.is_geographic:
            logging.warning(
                "Input CRS of road network has degree unit. Use default projected CRS epsg:3857"
            )
            R_gdf = R_gdf.to_crs(3857)

        O_gdf = O_gdf.to_crs(R_gdf.crs)
        D_gdf = D_gdf.to_crs(R_gdf.crs)
        # Generate ID
        O_gdf[self._OID] = Or_gdf[self.O_id_col].to_list()
        D_gdf[self._DID] = De_gdf[self.D_id_col].to_list()

        O_gdf["geometry"] = O_gdf.geometry.centroid
        D_gdf["geometry"] = D_gdf.geometry.centroid

        O_gdf["key"] = list(range(1, (O_gdf.shape[0] + 1)))
        gdfmax = O_gdf.key.max()
        O_gdf["Category"] = "Origin"
        D_gdf["key"] = list(range((gdfmax + 1), (gdfmax + 1 + D_gdf.shape[0])))
        D_gdf["Category"] = "Destination"

        O_gdf1 = O_gdf[["geometry", "key", "Category"]]
        D_gdf1 = D_gdf[["geometry", "key", "Category"]]
        POI = pd.concat([O_gdf1, D_gdf1], ignore_index=True)

        # Convert geoDataFrame to Graph with MomePy
        gdfroad = R_gdf
        graph = SimpleMomepy.gdf_to_nx(gdfroad)
        nodes, edges = SimpleMomepy.nx_to_gdf(
            graph, points=True, lines=True, spatial_weights=False
        )
        graph = SimpleMomepy.gdf_to_nx(edges)

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
        graph = SimpleMomepy.gdf_to_nx(edges)

        # Save geoDataFrame back to Points and Edges
        nodes1, edges1 = SimpleMomepy.nx_to_gdf(
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
                if list(graph.nodes)[i2] in all_dist:
                    dist.append(all_dist[list(graph.nodes)[i2]])
                    lent.append(all_length[list(graph.nodes)[i2]])
                else:
                    dist.append(999999)
                    lent.append(999999)
            data = {
                "originkey": a["key"],
                "destinationkey": b["key"],
                self._TIMECOST: dist,
                self._DISTCOST: lent,
                self._OSNAP: a["dist"],
                self._DSNAP: b["dist"],
            }
            df = pd.DataFrame(data)
            dff = pd.concat([df, dff], ignore_index=True)

        # dff["destinationid"] = dff.destinationid - gdfmax
        O_gdf2 = O_gdf[["key", self._OID]]
        D_gdf2 = D_gdf[["key", self._DID]]
        # # join the tables on the 'key' column
        result = pd.merge(dff, O_gdf2, left_on="originkey", right_on="key").drop(
            "key", axis=1
        )
        result = pd.merge(
            result, D_gdf2, left_on="destinationkey", right_on="key"
        ).drop("key", axis=1)
        result = result.copy()[
            [
                self._OID,
                self._DID,
                self._TIMECOST,
                self._DISTCOST,
                self._OSNAP,
                self._DSNAP,
            ]
        ]
        result.sort_values(by=[self._OID, self._DID], inplace=True)
        result.reset_index(drop=True, inplace=True)
        return knut.to_table(result)


############################################
# Isochrone map
############################################
@knext.node(
    name="Isochrone Map",
    node_type=knext.NodeType.MANIPULATOR,
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
    It first snaps  the input points to the road work, and uses the function
    [ego_graph](https://networkx.org/documentation/stable/reference/generated/networkx.generators.ego.ego_graph.html)
    in [NetworkX](https://networkx.org/) to isochrone for the weighted (time or distance) graph.

    The input value for the interval list should be selected carefully, it should be within a reasonable boundary of the road network.

    The output table contains the two columns, isochrone intervals, and geometry. The calculation depends on a projected coordinates
     system of road network datasets. If it is not in a projected CRS, a default CRS"3857" will be applied.
    If the input geometry is not a point feature, the centroid will be used. If it contains multiple rows, the total centroid will be applied.
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
        port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    isolist = knext.StringParameter(
        "Set isochone interval list",
        "Input a interval list separated by comma ",
        "5,10,15,20,25,30",
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

        return None

    def execute(self, exec_context: knext.ExecutionContext, input1, input2):

        import pandas as pd
        import networkx as nx
        from shapely.geometry import MultiPoint, LineString, Point, Polygon
        import numpy as np
        from pyproj import CRS  # For CRS Units check
        import logging
        import osmnx as ox

        def gdf_to_osmgraph(gdf):
            graph = SimpleMomepy.gdf_to_nx(gdf)
            nodes, edges = SimpleMomepy.nx_to_gdf(
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
        c_gdf = knut.load_geo_data_frame(input1, self.c_geo_col, exec_context)
        r_gdf = knut.load_geo_data_frame(input2, self.r_geo_col, exec_context)
        c_gdf = gp.GeoDataFrame(geometry=c_gdf[self.c_geo_col], crs=c_gdf.crs)
        r_gdf = r_gdf[[self.r_geo_col, self.r_cost_col]].rename(
            columns={self.r_geo_col: "geometry", self.r_cost_col: "time"}
        )
        # Set a lat\Lon CRS

        crsinput = CRS.from_user_input(r_gdf.crs)
        if crsinput.is_geographic:
            logging.warning(
                "Input CRS has degree unit. Use default projected CRS epsg:3857"
            )
            r_gdf = r_gdf.to_crs(3857)

        c_gdf = c_gdf.to_crs(r_gdf.crs)
        if c_gdf.shape[0] > 1:
            logging.warning("Only the total centroid was applied")
            c_gdf = gp.GeoDataFrame(geometry=gp.GeoSeries(c_gdf.unary_union.centroid))
        c_gdf["geometry"] = c_gdf.geometry.centroid

        # This example script simply outputs the node's input table.
        graph = gdf_to_osmgraph(r_gdf)
        nodes, edges = SimpleMomepy.nx_to_gdf(
            graph, points=True, lines=True, spatial_weights=False
        )

        nearest_node = gp.sjoin_nearest(c_gdf, nodes)
        center_node = nearest_node["osmid"][0]

        trip_times = list(map(int, self.isolist.split(",")))
        trip_times = sorted(trip_times)

        isochrone_polys = []
        for trip_time in sorted(trip_times):
            subgraph = nx.ego_graph(
                graph, center_node, radius=trip_time, undirected=True, distance="time"
            )
            edgex = SimpleMomepy.nx_to_gdf(
                subgraph, points=False, lines=True, spatial_weights=False
            )
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
        gdfx.reset_index(drop=True, inplace=True)
        return knut.to_table(gdfx)
