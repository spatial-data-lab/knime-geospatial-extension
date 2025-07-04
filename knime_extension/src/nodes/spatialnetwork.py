import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import datetime as dt

__category = knext.category(
    path="/community/geo",
    level_id="spatialnetwork",
    name="Spatial Network",
    description="Nodes that create distance matrices and isochrone maps using various routing engines.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/SpatialnetworkCategroy.png",
    after="opendataset",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Spatialnetwork/"

# Common used column names
_COL_O_ID = "Origin ID"
_COL_D_ID = "Destination ID"
_COL_DURATION = "Duration"
_COL_DURATION_TRAFFIC = "Duration in Traffic"
_COL_DISTANCE = "Distance"


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


class _GoogleTrafficModel(knext.EnumParameterOptions):
    BEST_GUESS = (
        "Best guess",
        """Indicates that the returned duration in traffic should be the best estimate of travel time given what is 
        known about both historical traffic conditions and live traffic. Live traffic becomes more important the closer
        the departure time is to now.""",
    )
    PESSIMISTIC = (
        "Pessimistic",
        """Indicates that the returned duration in traffic should be longer than the actual travel time on most days, 
        though occasional days with particularly bad traffic conditions may exceed this value.""",
    )
    OPTIMISTIC = (
        "Optimistic",
        """Indicates that the returned duration in traffic should be shorter than the actual travel time on most days, 
        though occasional days with particularly good traffic conditions may be faster than this value.""",
    )

    @classmethod
    def get_default(cls):
        return cls.BEST_GUESS


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
    This node uses the Google Distance Matrix API to create a distance matrix for the provided origins and destinations.

    This node uses the
    [Google Distance Matrix API](https://developers.google.com/maps/documentation/distance-matrix/overview)
    to create a distance matrix for the provided origins and destinations. The matrix is created by pairing each
    input origin with each input destination and will contain the travel distance and duration for each pair.
    The distance unit is meter and the duration is returned in minutes.

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
        label="Google API key",
        description="""Click [here](https://developers.google.com/maps/documentation/distance-matrix/get-api-key) for details on
        how to obtain and use a Google API key for the Distance Matrix API.""",
        default_value="your api key here",
        validator=knut.api_key_validator,
    )

    travel_mode = knext.EnumParameter(
        "Travel mode",
        """The following 
        [travel modes](https://developers.google.com/maps/documentation/distance-matrix/distance-matrix#mode) 
        are supported: """,
        default_value=_GoogleTravelMode.get_default().name,
        enum=_GoogleTravelMode,
    )

    consider_traffic = knext.BoolParameter(
        label="Consider traffic",
        description="""If checked, the travel time and distance will be computed considering the traffic conditions.
        If unchecked, the travel time and distance will be computed without considering the traffic conditions.""",
        default_value=False,
        since_version="1.3.0",
    )
    # add traffic models
    traffic_model = knext.EnumParameter(
        "Traffic model",
        """The [traffic model](https://developers.google.com/maps/documentation/distance-matrix/distance-matrix#traffic_model)
        specifies the assumptions to use when calculating time in traffic.
        """,
        default_value=_GoogleTrafficModel.get_default().name,
        since_version="1.3.0",
        enum=_GoogleTrafficModel,
    ).rule(knext.OneOf(consider_traffic, [True]), knext.Effect.SHOW)

    departure_time = knext.DateTimeParameter(
        label="Departure time",
        description="""The departure time may be specified in two cases:
                - For requests where the travel mode is transit: You can optionally specify departure_time or arrival_time. If None the departure_time defaults to now (that is, the departure time defaults to the current time).
                - For requests where the travel mode is driving: You can specify the departure_time to receive a route and trip duration (response field: duration_in_traffic) that take traffic conditions into account. The departure_time must be set to the current time or some time in the future. It cannot be in the past.
                
                Note: If departure time is not specified, choice of route and duration are based on road network and average time-independent traffic conditions. Results for a given request may vary over time due to changes in the road network, updated average traffic conditions, and the distributed nature of the service. Results may also vary between nearly-equivalent routes at any time or frequency.
                """,
        default_value=dt.datetime.now() + dt.timedelta(days=1),
        since_version="1.3.0",
        show_time=True,
        show_seconds=False,
    ).rule(knext.OneOf(consider_traffic, [True]), knext.Effect.SHOW)
    # Constant for distance matrix
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

        if self.api_key == "your api key here":
            configure_context.set_warning("Please provide a valid API key")

        if self.consider_traffic:
            return knext.Schema(
                [o_id_type, d_id_type, knext.double(), knext.double(), knext.int64()],
                [
                    _COL_O_ID,
                    _COL_D_ID,
                    _COL_DURATION_TRAFFIC,
                    _COL_DURATION,
                    _COL_DISTANCE,
                ],
            )
        else:
            return knext.Schema(
                [o_id_type, d_id_type, knext.double(), knext.int64()],
                [_COL_O_ID, _COL_D_ID, _COL_DURATION, _COL_DISTANCE],
            )

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        # define function to derive travel time and distance from Google Maps API

        import requests
        from pandas import json_normalize
        import numpy as np

        knut.check_canceled(exec_context)
        if self.api_key == "your api key here":
            raise ("Please provide a valid API key")

        def extract_coords(point):
            return point.centroid.y, point.centroid.x

        def update_distance_matrix(origins, destinations, start_index):
            origin_batch = "|".join([f"{lat},{lng}" for lat, lng in origins])
            destination_batch = "|".join([f"{lat},{lng}" for lat, lng in destinations])
            google_request_link = self._BASE_URL.format(
                origin_batch,
                destination_batch,
                self.api_key,
                self.travel_mode.lower(),
            )
            if self.consider_traffic:
                # handle departure time
                from datetime import datetime, timedelta

                departure_time_timestamp = int(self.departure_time.timestamp())
                departure_time = datetime.fromtimestamp(departure_time_timestamp)
                current_time = datetime.now()
                if departure_time < current_time + timedelta(minutes=10):
                    knut.LOGGER.warning(
                        "Departure time is in the past. Adjusting to the same time tomorrow."
                    )
                    departure_time = departure_time + timedelta(days=1)
                departure_time_timestamp = int(departure_time.timestamp())
                google_request_link = (
                    google_request_link
                    + "&traffic_model="
                    + "_".join(self.traffic_model.lower().split(" "))
                    + "&departure_time={}".format(departure_time_timestamp)
                )
            response = requests.get(google_request_link)
            data = response.json()
            if data["status"] == "OK":

                elements = json_normalize(data["rows"], record_path=["elements"])
                end_index = start_index + len(elements) - 1
                if self.consider_traffic:
                    duration_traffic_val = elements["duration_in_traffic.value"]
                    distance_matrix.loc[
                        start_index:end_index, _COL_DURATION_TRAFFIC
                    ] = np.array(duration_traffic_val / 60)
                duration_val = elements["duration.value"]
                distance_matrix.loc[start_index:end_index, _COL_DURATION] = np.array(
                    duration_val / 60
                )

                distance_matrix.loc[start_index:end_index, _COL_DISTANCE] = np.array(
                    elements["distance.value"]
                )
            else:
                knut.LOGGER.error(
                    f"Error fetching data: {data.get('error_message', 'No error message provided')}"
                )

        def calculate_od_batch(loop_coords, batch_coords, switch_od=False):
            for i in range(len(loop_coords)):
                process_counter = i + 1
                exec_context.set_progress(
                    0.9 * process_counter / len(loop_coords),
                    f"Batch {process_counter} of {len(loop_coords)} processed",
                )
                # Calculate the number of full batches
                jtime = len(batch_coords) // batchsize
                # Calculate the number of leftover destinations
                leftover = len(batch_coords) % batchsize
                for j in range(jtime):
                    batch_points = batch_coords[j * batchsize : (j + 1) * batchsize]
                    start_index = i * len(batch_coords) + j * batchsize
                    loop_points = loop_coords[i : i + 1]
                    if switch_od:
                        update_distance_matrix(loop_points, batch_points, start_index)
                    else:
                        update_distance_matrix(batch_points, loop_points, start_index)
                if leftover > 0:
                    left_start = len(batch_coords) - leftover
                    batch_points = batch_coords[left_start:]
                    start_index = (i + 1) * len(batch_coords) - leftover
                    loop_points = loop_coords[i : i + 1]
                    if switch_od:
                        update_distance_matrix(loop_points, batch_points, start_index)
                    else:
                        update_distance_matrix(batch_points, loop_points, start_index)

        o_gdf = knut.load_geo_data_frame(left_input, self.o_geo_col, exec_context)
        d_gdf = knut.load_geo_data_frame(right_input, self.d_geo_col, exec_context)

        # Set a lat\Lon CRS before renaming the geometry column
        o_gdf = o_gdf.to_crs(4326)
        d_gdf = d_gdf.to_crs(4326)

        # Filter all columns except the needed once and rename the geometry column to have a consistent result in merge
        o_gdf = o_gdf.filter(items=[self.o_geo_col, self.o_id_col]).rename(
            columns={self.o_geo_col: "geometry", self.o_id_col: _COL_O_ID}
        )
        d_gdf = d_gdf.filter(items=[self.d_geo_col, self.d_id_col]).rename(
            columns={self.d_geo_col: "geometry", self.d_id_col: _COL_D_ID}
        )

        # Create Coordinates
        o_gdf["origin_coords"] = o_gdf["geometry"].apply(extract_coords)
        d_gdf["destination_coords"] = d_gdf["geometry"].apply(extract_coords)
        # Flatten lists to send as OD pairs in batches
        origin_coords = list(o_gdf["origin_coords"])
        destination_coords = list(d_gdf["destination_coords"])
        # Generate origin destination matrix via cross join
        switch_od = (
            len(origin_coords) > 2 * len(destination_coords)
            and len(origin_coords) >= 25
        )
        if switch_od == False:
            merge_df = o_gdf.merge(d_gdf, how="cross")
        else:
            merge_df = d_gdf.merge(o_gdf, how="cross")

        # create the result matrix with the two id columns...
        distance_matrix = merge_df[[_COL_O_ID, _COL_D_ID]]
        # ... and default value 0 for the duration and distance column
        if self.consider_traffic:
            distance_matrix[_COL_DURATION_TRAFFIC] = 0
        distance_matrix[_COL_DURATION] = 0
        distance_matrix[_COL_DISTANCE] = 0

        batchsize = 25
        if len(distance_matrix) <= 25:
            update_distance_matrix(origin_coords, destination_coords, 0)
        else:
            if switch_od:
                calculate_od_batch(
                    destination_coords, origin_coords, switch_od=switch_od
                )
            else:
                calculate_od_batch(
                    origin_coords, destination_coords, switch_od=switch_od
                )

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
    This node uses the Open Source Routing Machine (OSRM) to create a distance matrix for the provided origins and destinations.

    This node uses the [Open Source Routing Machine (OSRM)](https://project-osrm.org/) to create a distance matrix
    for the provided origins and destinations. The matrix is created by pairing each input origin with each input
    destination and will contain the driving travel distance and time as well as the
    [route](http://project-osrm.org/docs/v5.5.1/api/?language=Python#route-service) for each pair.
    The travel distance unit is meter and the estimated drive time is returned in minutes.

    OSRM is a C++ implementation of a high-performance routing engine for shortest paths in road networks.
    It combines sophisticated routing algorithms with the open and free road network data of the
    [OpenStreetMap (OSM) project.](https://www.openstreetmap.org/about)

    If the input geometry is not a point geometry, the centroids will be automatically computed and used.

    ##Usage Policy
    Good practice and general limitations of the OSRM service that is used by this node can be found
    [here.](https://github.com/Project-OSRM/osrm-backend/wiki/Api-usage-policy)
    The current demo server is hosted by [FOSSGIS](https://www.fossgis.de/) which is subject to the usage policies and
    terms and conditions that can be found [here.](https://www.fossgis.de/arbeitsgruppen/osm-server/nutzungsbedingungen/)

    ##Note
    Data copyright by [OpenStreetMap](https://www.openstreetmap.org/copyright)
    [(ODbl)](https://opendatacommons.org/licenses/odbl/index.html) and provided under
    [CC-BY-SA.](https://creativecommons.org/licenses/by-sa/2.0/)
    To report a problem and contribute to OpenStreetMap click [here.](https://www.openstreetmap.org/fixthemap)
    Please note the licence/attribution guidelines as described [here.](https://wiki.osmfoundation.org/wiki/Licence/Attribution_Guidelines)
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
        label="Result mode",
        description="Supports the following result modes:",
        default_value=_OSRMResultModel.get_default().name,
        enum=_OSRMResultModel,
    )
    min_delay_seconds = knext.IntParameter(
        label="Minimum delay (seconds)",
        description="The minimum delay in seconds between two requests to the OSRM server.",
        default_value=1,
        since_version="1.2.0",
        is_advanced=True,
    )
    default_timeout = knext.IntParameter(
        label="Default timeout (seconds)",
        description="The default timeout in seconds for a request to the OSRM server.",
        default_value=10,
        since_version="1.2.0",
        is_advanced=True,
    )
    osrm_server = knext.StringParameter(
        label="OSRM server",
        description="""The URL of the OSRM server to use. The default value is the demo server hosted at https://map.project-osrm.org""",
        default_value="https://router.project-osrm.org",
        since_version="1.2.0",
        is_advanced=True,
    )
    # Constant for distance matrix
    _COL_GEOMETRY = "Route"

    # For details see: http://project-osrm.org/docs/v5.5.1/api/#route-service
    # _BASE_URL = "https://router.project-osrm.org"
    # only supports car as profile: https://github.com/Project-OSRM/osrm-backend/issues/4034
    _PROFILE = "driving"
    # For details see: http://project-osrm.org/docs/v5.5.1/api/#route-service
    _REQUEST_PARAMETER = {"continue_straight": "false"}

    # number of pairs send per request
    _BATCH_SIZE = 50
    # Number of seconds to wait after each request
    # _REQUEST_DELAY = 1
    # Request timeout
    # do not send more than 1 request per second https://github.com/Project-OSRM/osrm-backend/wiki/Demo-server
    # _REQUEST_TIMEOUT = None

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
        result_names = [_COL_O_ID, _COL_D_ID]
        if model.append_distance():
            result_types += [knext.double(), knext.double()]
            result_names += [_COL_DURATION, _COL_DISTANCE]
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
            self.osrm_server
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
                timeout=self.default_timeout,
            )
            time.sleep(self.min_delay_seconds)
            data = json.loads(r.text)
            if data["code"] == "Ok":
                if model.append_distance():
                    dfr = pd.DataFrame(data["routes"][0]["legs"])[
                        ["duration", "distance"]
                    ].iloc[::2]
                    # convert seconds to minutes
                    dfr.duration /= 60
                    df.loc[ns:ne, _COL_DURATION] = dfr.duration.to_list()
                    df.loc[ns:ne, _COL_DISTANCE] = dfr.distance.to_list()
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
            knut.LOGGER.warning(f"OSRM server not available. Error: {err}")

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        from shapely.geometry import LineString

        # Cross join Data
        o_gdf = knut.load_geo_data_frame(left_input, self.o_geo_col, exec_context)
        d_gdf = knut.load_geo_data_frame(right_input, self.d_geo_col, exec_context)

        # Set a lat\Lon CRS before renaming the geometry column
        o_gdf = o_gdf.to_crs(4326)
        d_gdf = d_gdf.to_crs(4326)

        # Filter all columns except the needed once and rename the geometry column to have a consistent result in merge
        o_gdf = o_gdf.filter(items=[self.o_geo_col, self.o_id_col]).rename(
            columns={self.o_geo_col: "geometry", self.o_id_col: _COL_O_ID}
        )
        d_gdf = d_gdf.filter(items=[self.d_geo_col, self.d_id_col]).rename(
            columns={self.d_geo_col: "geometry", self.d_id_col: _COL_D_ID}
        )

        # Generate origin destination matrix via cross join
        merge_df = o_gdf.merge(d_gdf, how="cross")
        merge_df_x = gp.GeoDataFrame(geometry=merge_df["geometry_x"], crs=4326)
        merge_df_y = gp.GeoDataFrame(geometry=merge_df["geometry_y"], crs=4326)
        df = merge_df[[_COL_O_ID, _COL_D_ID]]
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
            df[_COL_DURATION] = 0.0
            df[_COL_DISTANCE] = 0.0
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
    name="Road Network Distance Matrix",
    node_type=knext.NodeType.MANIPULATOR,
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
class RoadNetworkDistanceMatrix:
    """
    This node creates a distance matrix for the provided origins and destinations using the given road network.

    This node uses the [NetworkX library](https://networkx.org/) to create a distance matrix for the provided origins
    and destinations using the given road network. Prior computing the shortest path, origins and destinations
    are snapped to the closest point of the given road network. The matrix is then created by computing the
    [shortest path](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.algorithms.shortest_paths.weighted.single_source_dijkstra_path_length.html)
    between each snapped origin and all other reachable snapped destinations.

    The returned distance is in meters and the duration in minutes. In addition to the travel distance and duration,
    the output table contains the distance in meters between each origin and destination and its corresponding snap
    point along the road network, which can be incorporated into a total travel time and duration.

    The node projects the input coordinates to [EPSG:3857](https://epsg.io/3857) prior computing the length of each
    road and the snap distances of the origin and destinations.

    If the origin and destination geometries are not a point geometry, the centroids will be automatically computed
    and used.
    """

    class UnitModes(knext.EnumParameterOptions):
        METER_SECOND = (
            "Meters per second (m/s)",
            "The unit of the speed column is meters per second.",
        )
        MILE_HOUR = (
            "Miles per hour (mph)",
            "The unit of the speed column is miles per hour.",
        )
        KM_HOUR = (
            "Kilometers per hour (km/h)",
            "The unit of the speed column is kilometers per hour.",
        )
        DEFAULT_UNIT = (
            "Default unit",
            "If selected the speed column is used as is and assumed to be meters per minute.",
        )

        @classmethod
        def get_default(cls):
            return cls.KM_HOUR

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
        port_index=1,
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    r_geo_col = knext.ColumnParameter(
        "Road network geometry column",
        "Select the column which contains the road network data.",
        # Allow only GeoValue compatible columns
        port_index=2,
        column_filter=knut.is_geo_line,
        include_row_key=False,
        include_none_column=False,
    )

    r_speed_col = knext.ColumnParameter(
        "Road network speed column",
        "Select the column which contains the speed for the road network.",
        port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    speed_unit = knext.EnumParameter(
        label="Speed unit",
        description="The unit of the selected speed column.",
        default_value=UnitModes.get_default().name,
        enum=UnitModes,
    )
    # Constant for distance matrix
    _COL_O_SNAP = "Origin snap distance"
    _COL_D_SNAP = "Destination snap distance"

    _THRESHOLD = 100000

    # Define locate nearest edge (kne) and projected point (pp)
    def find_kne(self, point, lines):
        import numpy as np

        dists = np.array(list(map(lambda l: l.distance(point), lines)))
        kne_pos = dists.argsort()[0]
        kne = lines.iloc[[kne_pos]]
        kne_idx = kne.index[0]
        return kne_idx, kne.values[0]

    def get_pp(self, point, line):
        """Get the projected point (pp) of 'point' on 'line'."""
        pp = line.interpolate(line.project(point))  # PP as a Point
        return pp

    # split line with projected points
    def split_line(self, line, pps):
        from shapely.ops import snap, split

        line = snap(line, pps, 1e-8)  # slow?
        try:
            new_lines = list(split(line, pps).geoms)  # split into segments
            return new_lines
        except TypeError as e:
            print("Error when splitting line: {}\n{}\n{}\n".format(e, line, pps))
            return []

    # for interpolation (split by pp): replicate old line
    def update_edges(self, edges, new_lines, line_pps_dict, replace):
        import itertools
        import numpy as np
        import pandas as pd

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
            valid_pos = np.where(new_edges["length"] <= self._THRESHOLD)[0]
            n = len(new_edges)
            n_fault = n - len(valid_pos)
            f_pct = n_fault / n * 100
            print(
                "Remove faulty projections: {}/{} ({:.2f}%)".format(n_fault, n, f_pct)
            )
            new_edges = new_edges.iloc[valid_pos]  # use 'iloc' here

        # merge new edges
        dfs = [edges, new_edges]
        edges = gp.GeoDataFrame(
            pd.concat(dfs, ignore_index=False, sort=False), crs=dfs[0].crs
        )

        # all edges, newly added edges only
        return edges, new_edges

    # Define function of nearest points
    def ckd_nearest(self, gdA, gdB):
        import numpy as np
        import pandas as pd
        from scipy.spatial import cKDTree

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

    # check isolated node
    def connect_graph(self, G, threshold=0.1):
        import math
        import networkx as nx
        from shapely.geometry import LineString

        def distance(node1, node2):
            x1, y1 = node1
            x2, y2 = node2
            dx = x2 - x1
            dy = y2 - y1
            return math.sqrt(dx * dx + dy * dy)

        # check if the graph is connected
        if not nx.is_connected(G):
            print("The graph is not connected.")
            # find the connected components
            components = list(nx.connected_components(G))
            # add edges between the unconnected components
            for i in range(len(components)):
                for j in range(i + 1, len(components)):
                    component_i = components[i]
                    component_j = components[j]
                    edges = []
                    for node_i in component_i:
                        for node_j in component_j:
                            # check if the nodes are close to each other (optional)
                            if distance(node_i, node_j) < threshold:
                                line = LineString([node_i, node_j])
                                # add the edge with geometry and other attributes
                                G.add_edge(
                                    node_i,
                                    node_j,
                                    geometry=line,
                                    time=0,
                                    length=0,
                                    mm_len=0,
                                )
                                edges.append((node_i, node_j))
            # check if the graph is now connected
            if nx.is_connected(G):
                print("The graph is now connected.")
            else:
                print("The graph is still not connected.")
        return G

    def configure(
        self,
        configure_context,
        o_schema,
        d_schema,
        r_schema,
    ):
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

        self.r_geo_col = knut.column_exists_or_preset(
            configure_context, self.r_geo_col, r_schema, knut.is_geo
        )

        return knext.Schema(
            [
                o_id_type,
                d_id_type,
                knext.double(),
                knext.double(),
                knext.double(),
                knext.double(),
            ],
            [
                _COL_O_ID,
                _COL_D_ID,
                _COL_DURATION,
                _COL_DISTANCE,
                self._COL_O_SNAP,
                self._COL_D_SNAP,
            ],
        )

    def execute(self, exec_context: knext.ExecutionContext, input1, input2, input3):
        import networkx as nx
        from shapely.geometry import MultiPoint, LineString, Point
        import rtree
        from pyproj import CRS  # For CRS Units check
        import pandas as pd
        import logging

        # Cross join Data
        o_gdf = knut.load_geo_data_frame(input1, self.o_geo_col, exec_context)
        d_gdf = knut.load_geo_data_frame(input2, self.d_geo_col, exec_context)
        r_gdf = knut.load_geo_data_frame(input3, self.r_geo_col, exec_context)
        if not r_gdf.crs.is_projected:
            r_gdf = r_gdf.to_crs(3857)
            knut.LOGGER.warning("Road not projected. Using EPSG 3857 as default.")
        o_gdf = o_gdf.to_crs(r_gdf.crs)
        d_gdf = d_gdf.to_crs(r_gdf.crs)

        # Filter all columns except the needed once and rename geometry column since some methods expect it
        o_gdf = o_gdf.filter(items=[self.o_geo_col, self.o_id_col]).rename(
            columns={self.o_geo_col: "geometry", self.o_id_col: _COL_O_ID}
        )
        o_gdf.set_geometry("geometry", inplace=True)
        d_gdf = d_gdf.filter(items=[self.d_geo_col, self.d_id_col]).rename(
            columns={self.d_geo_col: "geometry", self.d_id_col: _COL_D_ID}
        )
        d_gdf.set_geometry("geometry", inplace=True)
        r_gdf = r_gdf.filter(items=[self.r_geo_col, self.r_speed_col]).rename(
            columns={self.r_geo_col: "geometry", self.r_speed_col: "speed"}
        )
        r_gdf.set_geometry("geometry", inplace=True)

        r_gdf = r_gdf.dropna(subset=["geometry"], how="any")

        # convert to meters per minute
        if self.speed_unit == self.UnitModes.METER_SECOND.name:
            r_gdf["speed"] = r_gdf.speed * 60
        elif self.speed_unit == self.UnitModes.MILE_HOUR.name:
            r_gdf["speed"] = r_gdf.speed * 26.8224
        elif self.speed_unit == self.UnitModes.KM_HOUR.name:
            r_gdf["speed"] = r_gdf.speed * 16.6667
        else:
            r_gdf["speed"] = r_gdf.speed * 1

        # ensure that origin and destination are points
        o_gdf["geometry"] = o_gdf.geometry.centroid
        d_gdf["geometry"] = d_gdf.geometry.centroid

        # generate unique key for all rows
        o_gdf["key"] = list(range(1, (o_gdf.shape[0] + 1)))
        gdf_max = o_gdf.key.max()
        o_gdf["Category"] = "Origin"
        d_gdf["key"] = list(range((gdf_max + 1), (gdf_max + 1 + d_gdf.shape[0])))
        d_gdf["Category"] = "Destination"

        # combine origin and destination to single table
        POI = pd.concat([o_gdf, d_gdf], ignore_index=True)

        # Convert geoDataFrame to Graph with MomePy
        graph = SimpleMomepy.gdf_to_nx(r_gdf)
        nodes, edges = SimpleMomepy.nx_to_gdf(
            graph, points=True, lines=True, spatial_weights=False
        )
        graph = SimpleMomepy.gdf_to_nx(edges)

        # build rtree
        r_tree = rtree.index.Index()
        [r_tree.insert(fid, geom.bounds) for fid, geom in edges["geometry"].items()]

        knn = 5

        # locate nearest edge (kne) and projected point (pp)
        # Projecting POIs to the network...
        POI["near_idx"] = [
            list(r_tree.nearest(point.bounds, knn)) for point in POI["geometry"]
        ]  # slow
        POI["near_lines"] = [
            edges["geometry"][near_idx] for near_idx in POI["near_idx"]
        ]  # very slow
        POI["kne_idx"], knes = zip(
            *[
                self.find_kne(point, near_lines)
                for point, near_lines in zip(POI["geometry"], POI["near_lines"])
            ]
        )  # slow
        POI["pp"] = [
            self.get_pp(point, kne) for point, kne in zip(POI["geometry"], knes)
        ]

        # 08-1: update internal edges (split line segments)
        line_pps_dict = {
            k: MultiPoint(list(v)) for k, v in POI.groupby(["kne_idx"])["pp"]
        }
        new_lines = [
            self.split_line(edges["geometry"][idx], pps)
            for idx, pps in line_pps_dict.items()
        ]  # bit slow
        edges, _ = self.update_edges(edges, new_lines, line_pps_dict, replace=True)

        edges["length"] = edges.length.astype(float)
        edges["time"] = edges.length / edges.speed

        # Convert geoDataFrame to Graph with MomePy
        graph = SimpleMomepy.gdf_to_nx(edges)
        graph = self.connect_graph(graph)
        # Save geoDataFrame back to Points and Edges
        nodes1, edges1 = SimpleMomepy.nx_to_gdf(
            graph, points=True, lines=True, spatial_weights=False
        )
        # Define function of nearest points
        snap_point = self.ckd_nearest(POI, nodes1)
        # set geometry to snapped nodes
        snap_point = snap_point.set_geometry("geometry")
        snap_point.reset_index(drop=True, inplace=True)
        x = snap_point[snap_point["Category"] == "Origin"]
        y = snap_point[snap_point["Category"] != "Origin"]

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
                _COL_DURATION: dist,
                _COL_DISTANCE: lent,
                self._COL_O_SNAP: a["dist"],
                self._COL_D_SNAP: b["dist"],
            }
            df = pd.DataFrame(data)
            dff = pd.concat([df, dff], ignore_index=True)

        o_gdf2 = o_gdf[["key", _COL_O_ID]]
        d_gdf2 = d_gdf[["key", _COL_D_ID]]
        # # join the tables on the 'key' column
        result = pd.merge(dff, o_gdf2, left_on="originkey", right_on="key").drop(
            "key", axis=1
        )
        result = pd.merge(
            result, d_gdf2, left_on="destinationkey", right_on="key"
        ).drop("key", axis=1)
        result = result.copy()[
            [
                _COL_O_ID,
                _COL_D_ID,
                _COL_DURATION,
                _COL_DISTANCE,
                self._COL_O_SNAP,
                self._COL_D_SNAP,
            ]
        ]
        result.sort_values(by=[_COL_O_ID, _COL_D_ID], inplace=True)
        result.reset_index(drop=True, inplace=True)
        return knut.to_table(result)


############################################
# Isochrone map
############################################
@knext.node(
    name="Road Network Isochrone Map",
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
class RoadNetworkIsochroneMap:
    """
    This node calculates the isochrone map for the input point based on the input road network and its travel cost column.

    This node calculates the [isochrone map](https://en.wikipedia.org/wiki/Isochrone_map) for the input point based on
    the input road network and its travel cost column.
    It first snaps the input points to the road network, and then uses the function
    [ego_graph](https://networkx.org/documentation/stable/reference/generated/networkx.generators.ego.ego_graph.html)
    in [NetworkX](https://networkx.org/) to isochrone for the weighted (time or distance) graph.

    The input value for the interval list should be selected carefully, it should be within a reasonable boundary of
    the road network.

    The output table contains the isochrone intervals and geometry. The calculation depends on a
    projected coordinates system of the input road network. If it is not in a projected CRS, it will be projected
    to [epsg:3857.](https://epsg.io/3857)
    If the input geometry is not a point feature, the centroid will be used.
    If it contains multiple rows, the total centroid will be applied.
    """

    # input parameters
    c_geo_col = knext.ColumnParameter(
        "Origin geometry column",
        "Select the geometry column that describes the origin.",
        # Allow only GeoValue compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    r_geo_col = knext.ColumnParameter(
        "Road network geometry column",
        "Select the column which contains the road network data.",
        # Allow only GeoValue compatible columns
        port_index=1,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    r_cost_col = knext.ColumnParameter(
        "Travel cost column from road network",
        "Select the column that contains the travel cost for the road network.",
        port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )
    iso_list = knext.StringParameter(
        "Isochrone intervals (comma separated)",
        "Input an interval list separated by comma e.g. 5,10,15,20,25,30",
    )

    _COL_GEOMETRY = "Geometry"
    _COL_ISOCHRONE = "Isochrone"

    def gdf_to_osmgraph(self, gdf):
        import osmnx as ox

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
        G = ox.convert.graph_from_gdfs(nodes, edges)
        return G

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

        return knext.Schema(
            [
                knut.TYPE_GEO,
                knext.int64(),
            ],
            [
                self._COL_GEOMETRY,
                self._COL_ISOCHRONE,
            ],
        )

    def execute(self, exec_context: knext.ExecutionContext, input1, input2):
        import pandas as pd
        import networkx as nx
        from shapely.geometry import MultiPoint, LineString, Point, Polygon
        import numpy as np
        from pyproj import CRS  # For CRS Units check
        import logging

        # Cross join Data
        c_gdf = knut.load_geo_data_frame(input1, self.c_geo_col, exec_context)
        r_gdf = knut.load_geo_data_frame(input2, self.r_geo_col, exec_context)
        c_gdf = gp.GeoDataFrame(geometry=c_gdf[self.c_geo_col], crs=c_gdf.crs)
        r_gdf = r_gdf[[self.r_geo_col, self.r_cost_col]].rename(
            columns={self.r_geo_col: "geometry", self.r_cost_col: "time"}
        )
        r_gdf.set_geometry("geometry", inplace=True)
        # Set a lat\Lon CRS

        crsinput = CRS.from_user_input(r_gdf.crs)
        if crsinput.is_geographic:
            r_gdf = r_gdf.to_crs(3857)
        c_gdf = c_gdf.to_crs(r_gdf.crs)
        if c_gdf.shape[0] > 1:
            # compute the global centroid
            c_gdf = gp.GeoDataFrame(geometry=gp.GeoSeries(c_gdf.unary_union.centroid))
        c_gdf[self.c_geo_col] = c_gdf.geometry.centroid

        # This example script simply outputs the node's input table.
        graph = self.gdf_to_osmgraph(r_gdf)
        nodes, edges = SimpleMomepy.nx_to_gdf(
            graph, points=True, lines=True, spatial_weights=False
        )

        nearest_node = gp.sjoin_nearest(c_gdf, nodes)
        center_node = nearest_node["osmid"][0]

        trip_times = list(map(int, self.iso_list.split(",")))
        trip_times = sorted(trip_times)

        len_trip_times = len(trip_times)
        isochrone_polys = []
        for trip_time in trip_times:
            subgraph = nx.ego_graph(
                graph, center_node, radius=trip_time, undirected=True, distance="time"
            )
            edge_x = SimpleMomepy.nx_to_gdf(
                subgraph, points=False, lines=True, spatial_weights=False
            )
            new_iso = Polygon(edge_x.unary_union.buffer(50).exterior)
            isochrone_polys.append(new_iso)
            exec_context.set_progress(
                0.7 * len(isochrone_polys) / float(len_trip_times),
                f"Isochrone {len(isochrone_polys)} of {len_trip_times} computed",
            )
            knut.check_canceled(exec_context)

        gd_fx = gp.GeoDataFrame(geometry=gp.GeoSeries(isochrone_polys), crs=r_gdf.crs)
        gd_fx[self._COL_ISOCHRONE] = trip_times

        if len(trip_times) > 1:
            # compute the difference between each isochrone and its predecessor
            for i in range(1, len(trip_times)):
                k = i - 1
                c0 = isochrone_polys[k]
                c1 = isochrone_polys[i]
                cd = c1.difference(c0)
                gd_fx.at[i, "geometry"] = cd
                exec_context.set_progress(
                    i / float(len_trip_times - 1),
                    f"Difference {i} of {len_trip_times - 1} computed",
                )
                knut.check_canceled(exec_context)

        gd_fx.rename(columns={"geometry": self._COL_GEOMETRY}, inplace=True)
        gd_fx.reset_index(drop=True, inplace=True)
        return knut.to_table(gd_fx)


############################################
# TomTom Isochrone Map Node
############################################
class _TomTomRouteType(knext.EnumParameterOptions):
    FASTEST = ("Fastest", "Focuses on reducing travel time.")
    SHORTEST = ("Shortest", "Prioritizes the shortest physical distance.")
    ECO = ("Eco", "Optimizes for fuel efficiency.")

    @classmethod
    def get_default(cls):
        return cls.FASTEST


class _TomTomTravelMode(knext.EnumParameterOptions):
    CAR = ("Car", "Car as vehicle type.")
    TRUCK = ("Truck", "Truck as vehicle type.")
    TAXI = ("Taxi", "Taxi as vehicle type.")
    BUS = ("Bus", "Bus as vehicle type.")
    VAN = ("Van", "Van as vehicle type.")
    MOTORCYCLE = ("Motorcycle", "Motorcycle as vehicle type.")

    @classmethod
    def get_default(cls):
        return cls.CAR


@knext.node(
    name="""TomTom Isochrone Map""",
    node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "TomTomIsochrone.png",
)
@knext.input_table(
    name="Input Table with the origin points.",
    description="Input table with point geometry representing the origin points.",
)
@knext.output_table(
    name="Output Table",
    description="Output table with isochrone geometry.",
)
class TomTomIsochroneMap:
    """This node calculates the isochrone map (reachable range) for a given geometric point using the
    Calculate Reachable Range service provided by TomTom.

    This node calculates the isochrone map (reachable range) for a list of given geometric points using the
    [Calculate Reachable Range service](https://developer.tomtom.com/routing-api/documentation/routing/calculate-reachable-range)
    of the [Routing service](https://www.tomtom.com/products/routing/)
    provided by [TomTom](https://www.tomtom.com/). It takes a geometry as origin and generates an isochrone map as output for each given geometry, illustrating
    areas reachable within a given time budget list. If the input geometry is not a point feature, the centroid will be used.

    Please note that this node requires a
    [TomTom API key](https://developer.tomtom.com/knowledgebase/platform/articles/how-to-get-an-tomtom-api-key/)
    that can be acquired for free by [registering here.](https://developer.tomtom.com/user/register)
    For more details about the number of free request and pricing go to the
    [TomTom pricing page.](https://developer.tomtom.com/store/maps-api)
    """

    # input parameters
    c_geo_col = knext.ColumnParameter(
        "Origin geometry column",
        """This parameter selects the geometry column from the input table that represents the origin point for 
        the isochrone calculation.""",
        # Allow only Geo compatible columns
        port_index=0,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    id_col = knext.ColumnParameter(
        "Origin ID column",
        """This parameter selects the column which contains for each origin a unique ID. The selected column will be
        returned in the result table and can be used to link back to the original data.""",
        port_index=0,
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    iso_time_budget_list = knext.StringParameter(
        "Isochrone time budget list",
        """Input an interval list in minutes separated by comma e.g. 5,10,15,20,25,30""",
        default_value="5,10,15,20,25,30",
    )

    depart_at = knext.DateTimeParameter(
        "Departure time",
        """The departure time for the isochrone calculation. This parameter can affect the isochrone map due to 
        varying traffic conditions at different times.""",
        default_value=None,
        show_time=True,
    )

    traffic = knext.BoolParameter(
        "Consider current traffic",
        """If selected the current traffic conditions is considered in the isochrone calculation. 
        Note that information on historic road speeds is always used.""",
        default_value=True,
    )

    route_type = knext.EnumParameter(
        "Route type",
        """Determines the type of route used for the isochrone calculation. 
        Different route types can result in different isochrones.""",
        default_value=_TomTomRouteType.get_default().name,
        enum=_TomTomRouteType,
        style=knext.EnumParameter.Style.DROPDOWN,
    )

    travel_mode = knext.EnumParameter(
        "Travel mode",
        """Specifies the mode of travel for the isochrone calculation, which can 
        significantly impact the shape and extent of the isochrone.""",
        default_value=_TomTomTravelMode.get_default().name,
        enum=_TomTomTravelMode,
    )

    tomtom_api_key = knext.StringParameter(
        "TomTom API Key",
        """The 
        [TomTom API key](https://developer.tomtom.com/knowledgebase/platform/articles/how-to-get-an-tomtom-api-key/)
        is required to authenticate requests to the 
        [Calculate Reachable Range service](https://developer.tomtom.com/routing-api/documentation/routing/calculate-reachable-range)
        which is part of the [Routing API](https://developer.tomtom.com/routing-api/documentation/routing/routing-service)
        provided by TomTom. 
        To get an API key, click
        [here.](https://developer.tomtom.com/knowledgebase/platform/articles/how-to-get-an-tomtom-api-key/)
        For details about the pricing go to the [TomTom pricing page.](https://developer.tomtom.com/store/maps-api)""",
        default_value="your api key here",
        validator=knut.api_key_validator,
    )

    timeout = knext.IntParameter(
        "Request timeout in seconds",
        "The maximum time in seconds to wait for the request to the TomTom API to succeed.",
        120,
        min_value=1,
        is_advanced=True,
    )

    _COL_GEOMETRY = knut.DEF_GEO_COL_NAME
    _COL_ISOCHRONE = "Time budget (Mins)"

    def configure(self, configure_context, input_schema_1):
        self.c_geo_col = knut.column_exists_or_preset(
            configure_context, self.c_geo_col, input_schema_1, knut.is_geo
        )
        self.id_col = knut.column_exists_or_preset(
            configure_context, self.id_col, input_schema_1, knut.is_numeric_or_string
        )

        return knext.Schema(
            [
                input_schema_1[self.id_col].ktype,
                knext.int64(),
                # input_schema_1[self.c_geo_col].ktype,
                knut.TYPE_POLYGON,
            ],
            [
                self.id_col,
                self._COL_ISOCHRONE,
                self._COL_GEOMETRY,
            ],
        )

    def execute(self, exec_context: knext.ExecutionContext, input1):

        tomtom_base_url = "https://api.tomtom.com/routing/1/calculateReachableRange/"
        import requests
        import json
        from shapely.geometry import Polygon

        c_gdf = knut.load_geo_data_frame(input1, self.c_geo_col, exec_context)
        iso_map_list = []
        loop_i = 1
        total_loops = len(c_gdf) * len(self.iso_time_budget_list.split(","))
        if self.tomtom_api_key == "your api key here" or self.tomtom_api_key == "":
            knut.LOGGER.error(
                "Please enter your TomTom API key. If you don't have one, you can get one [here](https://developer.tomtom.com/knowledgebase/platform/articles/how-to-get-an-tomtom-api-key/)."
            )
            raise ValueError(
                "Please enter your TomTom API key. If you don't have one, you can get one [here](https://developer.tomtom.com/knowledgebase/platform/articles/how-to-get-an-tomtom-api-key/)."
            )
        from datetime import datetime, timedelta

        current_time = datetime.now()
        depart_at_datetime = datetime.fromtimestamp(int(self.depart_at.timestamp()))
        if depart_at_datetime < current_time + timedelta(minutes=1):
            knut.LOGGER.warning(
                "Departure time is in the past. Adjusting to the same time tomorrow."
            )
            depart_at_datetime += timedelta(days=1)

        for k, row in c_gdf.iterrows():
            id_ = row[self.id_col]
            x = str(row[self.c_geo_col].centroid.x)
            y = str(row[self.c_geo_col].centroid.y)
            time_budgets = list(map(int, self.iso_time_budget_list.split(",")))

            for time_budget in time_budgets:
                URL = (
                    "%s%s,%s/json?timeBudgetInSec=%s&travelMode=%s&traffic=%s&key=%s&routeType=%s&departAt=%s"
                    % (
                        tomtom_base_url,
                        y,
                        x,
                        str(time_budget * 60),
                        self.travel_mode.lower(),
                        str(self.traffic).lower(),
                        self.tomtom_api_key,
                        self.route_type.lower(),
                        depart_at_datetime.isoformat()[0:19],
                    )
                )

                req = requests.get(URL, timeout=self.timeout)
                response_code = req.status_code
                if response_code != 200:
                    knut.LOGGER.error(f"Error! TomTom response code: {response_code} ")
                    raise ValueError(f"Error! TomTom response code: {response_code} ")
                data = json.loads(req.text)
                bounds = data["reachableRange"]["boundary"]
                bounds_polygon = Polygon(
                    [(x["longitude"], x["latitude"]) for x in bounds]
                )
                exec_context.set_progress(
                    0.9 * loop_i / float(total_loops),
                    f"Isochrone {loop_i} of {total_loops} computed",
                )
                knut.check_canceled(exec_context)
                loop_i += 1
                iso_map_list.append([id_, time_budget, bounds_polygon])
        gdf = gp.GeoDataFrame(
            iso_map_list, columns=[self.id_col, self._COL_ISOCHRONE, self._COL_GEOMETRY]
        )
        gdf.set_geometry(self._COL_GEOMETRY, inplace=True)
        gdf.crs = c_gdf.crs

        return knut.to_table(gdf)


############################################
# TomTom Distance Matrix Node
############################################


class _TomTomMatrixDepartTimeType(knext.EnumParameterOptions):
    ANY = ("Any", "Use the default departure time.")
    NOW = ("Now", "Use the current time as departure time.")
    DATETIME = ("Custom datetime", "Specify a custom departure time.")

    @classmethod
    def get_default(cls):
        return cls.ANY


class _TomTomMatrixRouteType(knext.EnumParameterOptions):
    FASTEST = ("Fastest", "Focuses on reducing travel time.")
    SHORTEST = ("Shortest", "Prioritizes the shortest physical distance.")

    @classmethod
    def get_default(cls):
        return cls.FASTEST


class _TomTomMatrixTravelMode(knext.EnumParameterOptions):
    CAR = ("Car", "Car as vehicle type.")
    TRUCK = ("Truck", "Truck as vehicle type.")
    PEDESTRIAN = ("Pedestrian", "Walking as travel mode.")

    @classmethod
    def get_default(cls):
        return cls.CAR


class _TomTomMatrixBatchSize(knext.EnumParameterOptions):
    N100 = ("100", "No limitations.")
    N200 = (
        "200",
        "All origins and destinations should be contained in an axis-aligned 400 km x 400 km bounding box. Otherwise, some matrix cells will be resolved as OUT_OF_REGION.",
    )
    N2500 = (
        "2500",
        "Route type will be FASTEST, traffic will be HISTORICAL, travel mode must be CAR or TRUCK, and departure time type will be ANY.",
    )

    @classmethod
    def get_default(cls):
        return cls.N100


@knext.node(
    name="""TomTom Distance Matrix""",
    node_type=knext.NodeType.MANIPULATOR,
    category=__category,
    icon_path=__NODE_ICON_PATH + "TomTomDistanceMatrix.png",
)
@knext.input_table(
    name="Origins Table",
    description="A table containing origin geometries and a unique ID column.",
)
@knext.input_table(
    name="Destinations Table",
    description="A table containing destination geometries and a unique ID column.",
)
@knext.output_table(
    name="Distance Matrix Table",
    description="""A table containing the origin and destination IDs along with the corresponding 
    travel distances (meters) and durations (minutes).""",
)
class TomTomDistanceMatrix:
    """Calculates travel distances and durations between origins and destinations using the [TomTom Routing API](https://www.tomtom.com/products/routing-apis/).

    This node uses the [TomTom Matrix Routing v2 API](https://developer.tomtom.com/matrix-routing-v2-api/documentation/product-information/introduction)
    to computes a distance matrix by pairing each origin with each destination, returning travel distances
    (in meters) and durations (in minutes) for each pair. If the input geometries are not point-based, their
    centroids are automatically used.

    **Note:** A [TomTom API key](https://developer.tomtom.com/knowledgebase/platform/articles/how-to-get-an-tomtom-api-key/)
    is required. You can obtain one for free by [registering](https://developer.tomtom.com/user/register) on the
    TomTom website. Refer to the [TomTom pricing page](https://developer.tomtom.com/pricing) for information about free request limits and costs.
    """

    # Origin parameters
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

    # Destination parameters
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

    # API key parameter
    tomtom_api_key = knext.StringParameter(
        "TomTom API key",
        """The 
        [TomTom API key](https://developer.tomtom.com/user/register)
        is required to authenticate requests to the 
        [Matrix Routing API V2](https://developer.tomtom.com/matrix-routing-v2-api/documentation)
        provided by TomTom. 
        To get an API key, click
        [here.](https://developer.tomtom.com/user/register)
        For details about the pricing go to the [TomTom pricing page.](https://developer.tomtom.com/matrix-routing-v2-api/documentation/pricing)""",
        validator=knut.api_key_validator,
    )

    # Routing parameters
    route_type = knext.EnumParameter(
        "Route type",
        """Determines the type of route used for the distance calculation. 
        Different route types can result in different route choices.""",
        default_value=_TomTomMatrixRouteType.get_default().name,
        enum=_TomTomMatrixRouteType,
        style=knext.EnumParameter.Style.DROPDOWN,
    )

    travel_mode = knext.EnumParameter(
        "Travel mode",
        """Specifies the mode of travel for the distance calculation, which can 
        significantly impact the route choice, distance and travel time.""",
        default_value=_TomTomMatrixTravelMode.get_default().name,
        enum=_TomTomMatrixTravelMode,
    )

    consider_traffic = knext.BoolParameter(
        "Consider current traffic",
        """If selected, the current traffic conditions are considered in the distance and duration calculations. 
        Note that information on historic road speeds is always used.
        This option is only available for CAR and TRUCK travel modes.""",
        default_value=True,
    ).rule(knext.OneOf(travel_mode, ["CAR", "TRUCK"]), knext.Effect.SHOW)

    depart_time_type = knext.EnumParameter(
        "Departure time type",
        """Specify how to set the departure time:
        - Any: Use the default departure time
        - Now: Use the current time as departure time
        - Custom datetime: Specify a custom departure time""",
        default_value=_TomTomMatrixDepartTimeType.get_default().name,
        enum=_TomTomMatrixDepartTimeType,
    ).rule(knext.OneOf(travel_mode, ["CAR", "TRUCK"]), knext.Effect.SHOW)

    depart_at = knext.DateTimeParameter(
        "Custom departure time",
        """Specify a custom departure time for the distance calculation.
        This time must be in the future.""",
        default_value=None,
        show_time=True,
        show_seconds=False,
    ).rule(knext.OneOf(depart_time_type, ["DATETIME"]), knext.Effect.SHOW)

    # Batch processing parameters
    batch_size = knext.EnumParameter(
        "Maximum matrix batch size",
        """The maximum number of cells (origin-destination pairs) to include in a single API request.
        Recommended values: 100 for unrestricted requests, up to 200 for geographical restrictions.
        For values as 2500, parameters will be automatically adjusted to comply with API restrictions.
        As billable requests are calculated as max(origins, destinations) × 5, the querying batch size 
        for origin and destination points uses the square root of the Maximum matrix size, e.g.,10, 14, 50.
        """,
        default_value=_TomTomMatrixBatchSize.get_default().name,
        enum=_TomTomMatrixBatchSize,
        is_advanced=True,
    )

    # Advanced parameters
    timeout = knext.IntParameter(
        "Request timeout in seconds",
        "The maximum time in seconds to wait for the request to the TomTom API to succeed.",
        120,
        min_value=1,
        is_advanced=True,
    )

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

        # Check for batch size restrictions
        if int(self.batch_size[1:]) > 200:
            configure_context.set_warning(
                "For batch size > 200, some parameters will be automatically adjusted: "
                + "route type = 'Fastest', traffic = 'historical', travel mode = 'Car/Truck', departure time = 'Any'"
            )

        # Standard schema with origin ID, destination ID, duration and distance
        return knext.Schema(
            [o_id_type, d_id_type, knext.double(), knext.int64()],
            [_COL_O_ID, _COL_D_ID, _COL_DURATION, _COL_DISTANCE],
        )

    def execute(self, exec_context: knext.ExecutionContext, left_input, right_input):
        import requests
        import json
        import pandas as pd
        import numpy as np
        import math
        import datetime as dt
        from datetime import datetime, timedelta

        # Check if API key is provided
        if self.tomtom_api_key is None or len(self.tomtom_api_key.strip()) == 0:
            knut.LOGGER.error(
                "Please enter your TomTom API key. If you don't have one, you can register [here](https://developer.tomtom.com/user/register)."
            )
            raise ValueError(
                "Please enter your TomTom API key. If you don't have one, you can register [here](https://developer.tomtom.com/user/register)."
            )

        knut.check_canceled(exec_context)

        # Load GeoDataFrames for origins and destinations
        o_gdf = knut.load_geo_data_frame(left_input, self.o_geo_col, exec_context)
        d_gdf = knut.load_geo_data_frame(right_input, self.d_geo_col, exec_context)

        # Set a lat/lon CRS before extracting coordinates
        o_gdf = o_gdf.to_crs(4326)
        d_gdf = d_gdf.to_crs(4326)

        # Filter to only needed columns and rename
        o_gdf = o_gdf.filter(items=[self.o_geo_col, self.o_id_col]).rename(
            columns={self.o_geo_col: "geometry", self.o_id_col: _COL_O_ID}
        )
        d_gdf = d_gdf.filter(items=[self.d_geo_col, self.d_id_col]).rename(
            columns={self.d_geo_col: "geometry", self.d_id_col: _COL_D_ID}
        )

        # Extract centroids and coordinates
        # Format coordinates correctly for TomTom API V2
        o_gdf["origin_points"] = o_gdf["geometry"].apply(
            lambda point: {
                "point": {
                    "latitude": float(point.centroid.y),
                    "longitude": float(point.centroid.x),
                }
            }
        )
        d_gdf["destination_points"] = d_gdf["geometry"].apply(
            lambda point: {
                "point": {
                    "latitude": float(point.centroid.y),
                    "longitude": float(point.centroid.x),
                }
            }
        )

        # Prepare origins and destinations lists for the API request
        origins = o_gdf["origin_points"].tolist()
        destinations = d_gdf["destination_points"].tolist()

        # Prepare the cross join result dataframe for the distance matrix
        merge_df = o_gdf[[_COL_O_ID]].merge(d_gdf[[_COL_D_ID]], how="cross")

        # Create the empty distance matrix with default values
        distance_matrix = merge_df.copy()
        distance_matrix[_COL_DURATION] = 0
        distance_matrix[_COL_DISTANCE] = 0

        # Prepare departure time if specified and traffic is considered
        departure_time = None
        if self.consider_traffic and self.travel_mode in ["CAR", "TRUCK"]:
            if self.depart_time_type == "DATETIME":
                depart_at_datetime = self.depart_at
                current_time = datetime.now()
                if depart_at_datetime < current_time + timedelta(minutes=1):
                    knut.LOGGER.warning(
                        "Departure time is in the past. Adjusting to the same time tomorrow."
                    )
                    depart_at_datetime += timedelta(days=1)
                departure_time = depart_at_datetime.isoformat()[0:19]
            elif self.depart_time_type == "NOW":
                departure_time = "now"
            # For ANY, we leave departure_time as None

        # Create a copy of selected options that might need adjustment
        effective_route_type = self.route_type
        effective_travel_mode = self.travel_mode
        effective_consider_traffic = self.consider_traffic
        effective_depart_time_type = self.depart_time_type

        # For large batch sizes, automatically adjust parameters if needed
        actual_batch_size = int(self.batch_size[1:])
        if actual_batch_size > 200:
            knut.LOGGER.info(
                "Using batch size > 200, automatically adjusting request parameters to comply with API restrictions"
            )
            # For large batch sizes, automatically adjust to required values
            effective_route_type = "FASTEST"
            # Only CAR or TRUCK allowed for large batches
            if effective_travel_mode not in ["CAR", "TRUCK"]:
                effective_travel_mode = "CAR"
                knut.LOGGER.info(
                    "Adjusted travel mode to 'CAR' for large batch processing"
                )
            effective_consider_traffic = False
            effective_depart_time_type = "ANY"

        # Always use batch processing to handle potential large matrices
        self.process_matrix_in_batches(
            exec_context,
            origins,
            destinations,
            o_gdf[_COL_O_ID].tolist(),
            d_gdf[_COL_D_ID].tolist(),
            distance_matrix,
            departure_time,
            effective_route_type,
            effective_travel_mode,
            effective_consider_traffic,
            effective_depart_time_type,
        )

        return knut.to_table(distance_matrix)

    def process_matrix_in_batches(
        self,
        exec_context,
        origins,
        destinations,
        origin_ids,
        dest_ids,
        distance_matrix,
        departure_time,
        effective_route_type=None,
        effective_travel_mode=None,
        effective_consider_traffic=None,
        effective_depart_time_type=None,
    ):
        """Process a large matrix by splitting it into smaller batches."""
        import math

        knut.LOGGER.info(
            f"Processing matrix with {len(origins)} origins and {len(destinations)} destinations in batches"
        )

        # Use effective parameters if provided, otherwise use the class parameters
        route_type = (
            effective_route_type
            if effective_route_type is not None
            else self.route_type
        )
        travel_mode = (
            effective_travel_mode
            if effective_travel_mode is not None
            else self.travel_mode
        )
        consider_traffic = (
            effective_consider_traffic
            if effective_consider_traffic is not None
            else self.consider_traffic
        )
        depart_time_type = (
            effective_depart_time_type
            if effective_depart_time_type is not None
            else self.depart_time_type
        )

        # Determine appropriate batch size based on parameters and API limits
        actual_batch_size = int(self.batch_size[1:])

        # Ensure neither dimension exceeds 1000 (API limit)
        origins_per_batch = int(math.sqrt(actual_batch_size))
        dests_per_batch = int(math.sqrt(actual_batch_size))

        # Calculate total number of batches
        num_o_batches = math.ceil(len(origins) / origins_per_batch)
        num_d_batches = math.ceil(len(destinations) / dests_per_batch)
        total_batches = num_o_batches * num_d_batches

        current_batch = 0

        # Process each batch
        for o_start in range(0, len(origins), origins_per_batch):
            o_end = min(o_start + origins_per_batch, len(origins))
            origin_batch = origins[o_start:o_end]

            for d_start in range(0, len(destinations), dests_per_batch):
                d_end = min(d_start + dests_per_batch, len(destinations))
                dest_batch = destinations[d_start:d_end]

                current_batch += 1
                exec_context.set_progress(
                    0.1 + 0.8 * current_batch / total_batches,
                    f"Processing batch {current_batch} of {total_batches}",
                )

                # Process this batch
                self.process_single_request(
                    exec_context,
                    origin_batch,
                    dest_batch,
                    distance_matrix,
                    departure_time,
                    o_start,
                    d_start,
                    o_end - o_start,
                    d_end - d_start,
                    route_type,
                    travel_mode,
                    consider_traffic,
                    depart_time_type,
                )

                knut.check_canceled(exec_context)

    def process_single_request(
        self,
        exec_context,
        origins,
        destinations,
        distance_matrix,
        departure_time,
        o_offset,
        d_offset,
        o_count,
        d_count,
        route_type=None,
        travel_mode=None,
        consider_traffic=None,
        depart_time_type=None,
    ):
        """Process a single API request for a matrix or sub-matrix."""
        import requests
        import json

        # Use provided parameters or fall back to class parameters
        route_type = route_type if route_type is not None else self.route_type
        travel_mode = travel_mode if travel_mode is not None else self.travel_mode
        consider_traffic = (
            consider_traffic if consider_traffic is not None else self.consider_traffic
        )
        depart_time_type = (
            depart_time_type if depart_time_type is not None else self.depart_time_type
        )

        # Prepare the API request
        tomtom_base_url = "https://api.tomtom.com/routing/matrix/2"

        # Prepare the request body according to the Matrix Routing API V2 documentation
        request_body = {
            "origins": origins,
            "destinations": destinations,
            "options": {
                "routeType": route_type.lower(),
                "travelMode": travel_mode.lower(),
            },
        }

        # Add traffic options if applicable
        if travel_mode in ["CAR", "TRUCK"]:
            if consider_traffic:
                request_body["options"]["traffic"] = "live"
                # When traffic is live, departAt must be specified
                if depart_time_type == "DATETIME" and departure_time:
                    request_body["options"]["departAt"] = departure_time
                else:
                    # Default to "now" if not specified or if ANY is selected
                    request_body["options"]["departAt"] = "now"
            else:
                request_body["options"]["traffic"] = "historical"
        else:
            # For PEDESTRIAN, traffic is not applicable
            if "traffic" in request_body["options"]:
                del request_body["options"]["traffic"]
            if "departAt" in request_body["options"]:
                del request_body["options"]["departAt"]

        # Log the request for debugging
        batch_info = f"{len(origins)}×{len(destinations)} matrix"
        if o_count < len(origins) or d_count < len(destinations):
            batch_info = (
                f"batch {o_offset}:{o_offset+o_count}, {d_offset}:{d_offset+d_count}"
            )

        knut.LOGGER.info(f"Making request to TomTom Matrix API for {batch_info}")
        knut.LOGGER.debug(f"Request body: {json.dumps(request_body)[:500]}...")

        # Make the API request
        url = f"{tomtom_base_url}?key={self.tomtom_api_key}"
        headers = {"Content-Type": "application/json"}

        try:
            response = requests.post(
                url,
                data=json.dumps(request_body),
                headers=headers,
                timeout=self.timeout,
            )

            # Check if the request was successful
            if response.status_code == 200:
                data = response.json()

                # Log the response for debugging
                knut.LOGGER.debug(f"TomTom API Response: {json.dumps(data)[:500]}...")

                # Process the response based on actual V2 API format
                # Extract data from the 'data' array which contains the matrix results
                matrix_data = data.get("data", [])

                # Process each origin-destination pair
                for route_data in matrix_data:
                    origin_index = route_data.get("originIndex", 0)
                    destination_index = route_data.get("destinationIndex", 0)

                    # Calculate the index in the full matrix
                    full_origin_index = o_offset + origin_index
                    full_dest_index = d_offset + destination_index

                    # Calculate the row index in our distance matrix
                    # This works based on how the cross join is created in pandas
                    global_dest_count = len(distance_matrix) // (
                        len(distance_matrix[_COL_O_ID].unique())
                    )
                    row_index = full_origin_index * global_dest_count + full_dest_index

                    # Get route summary
                    route_summary = route_data.get("routeSummary", {})

                    # Extract lengthInMeters and travelTimeInSeconds
                    distance_meters = route_summary.get("lengthInMeters", 0)
                    duration_seconds = route_summary.get("travelTimeInSeconds", 0)

                    # Store values in the distance matrix (convert seconds to minutes)
                    distance_matrix.loc[row_index, _COL_DISTANCE] = distance_meters
                    distance_matrix.loc[row_index, _COL_DURATION] = (
                        duration_seconds / 60
                    )
            else:
                error_message = f"TomTom API Error: Status code {response.status_code}"
                try:
                    error_json = response.json()
                    error_message += f" - Full response: {json.dumps(error_json)}"
                except:
                    error_message += f" - Response text: {response.text}"

                knut.LOGGER.error(error_message)
                raise ValueError(error_message)

        except requests.exceptions.Timeout:
            error_msg = f"Request to TomTom API timed out after {self.timeout} seconds"
            knut.LOGGER.error(error_msg)
            raise ValueError(error_msg)
        except Exception as e:
            error_msg = f"Error: {str(e)}"
            knut.LOGGER.error(error_msg)
            raise ValueError(error_msg)
