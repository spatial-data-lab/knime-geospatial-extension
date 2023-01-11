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
    Google Distance Matrix
    """

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
    API_Key = knext.StringParameter(
        "Google API Key",
        "Google API Key for distance matrix[FIPS](https://developers.google.com/maps/documentation/distance-matrix/overview) ",
        "",
    )
    Travel_Mode = knext.StringParameter(
        "Travel Mode",
        "Set Trave Mode[FIPS](https://developers.google.com/maps/documentation/distance-matrix/distance-matrix) ",
        "Driving",
        enum=["Driving", "Transit"],
    )

    def configure(self, configure_context, input_schema_1, input_schema_2):
        return knext.Schema.from_columns(
            [
                knext.Column(knext.int64(), "OriginID"),
                knext.Column(knext.int64(), "DestinationID"),
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

        O_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.O_geo_col)
        D_gdf = gp.GeoDataFrame(right_input.to_pandas(), geometry=self.D_geo_col)
        # Set a lat\Lon CRS
        O_gdf = O_gdf.to_crs(4326)
        D_gdf = D_gdf.to_crs(4326)
        # Generate ID
        O_gdf["OriginID"] = range(1, (O_gdf.shape[0] + 1))
        D_gdf["DestinationID"] = range(1, (D_gdf.shape[0] + 1))
        mergedf = O_gdf.merge(D_gdf, how="cross")
        mergedf_x = gp.GeoDataFrame(geometry=mergedf["geometry_x"])
        mergedf_y = gp.GeoDataFrame(geometry=mergedf["geometry_y"])

        distance_matrix = mergedf[["OriginID", "DestinationID"]]
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
    OSRM Distance Matrix
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
        # from shapely.geometry import LineString

        # if self.osrmmodel == "Travel Cost":
        #     return knext.Schema.from_columns(
        #         [
        #             knext.Column(knext.int64(), "OriginID"),
        #             knext.Column(knext.int64(), "DestinationID"),
        #             knext.Column(knext.double(), "duration"),
        #             knext.Column(knext.double(), "distance"),
        #         ]
        #     )
        # elif self.osrmmodel == "Route":
        #     return knext.Schema.from_columns(
        #         [
        #             knext.Column(knext.int64(), "OriginID"),
        #             knext.Column(knext.int64(), "DestinationID"),
        #             knext.Column(knext.logical(LineString), "geometry"),
        #         ]
        #     )
        # else:
        #     return knext.Schema.from_columns(
        #         [
        #             knext.Column(knext.int64(), "OriginID"),
        #             knext.Column(knext.int64(), "DestinationID"),
        #             knext.Column(knext.double(), "duration"),
        #             knext.Column(knext.double(), "distance"),
        #             knext.Column(knext.logical(LineString), "geometry"),
        #         ]
        #     )
        # TODO LineString can not be identified
        return None

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
                        df.loc[ns:ne, "geometry"] = extractRoute(data)
                else:
                    print("error from:{} to :{}".format(ns, ne))
            except:
                print("error from:{} to :{}".format(ns, ne))

        # Cross join Data
        O_gdf = gp.GeoDataFrame(left_input.to_pandas(), geometry=self.O_geo_col)
        D_gdf = gp.GeoDataFrame(right_input.to_pandas(), geometry=self.D_geo_col)
        # Set a lat\Lon CRS
        O_gdf = O_gdf.to_crs(4326)
        D_gdf = D_gdf.to_crs(4326)
        # Generate ID
        O_gdf["OriginID"] = range(1, (O_gdf.shape[0] + 1))
        D_gdf["DestinationID"] = range(1, (D_gdf.shape[0] + 1))
        mergedf = O_gdf.merge(D_gdf, how="cross")
        mergedf_x = gp.GeoDataFrame(geometry=mergedf["geometry_x"])
        mergedf_y = gp.GeoDataFrame(geometry=mergedf["geometry_y"])

        df = mergedf[["OriginID", "DestinationID"]]
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
            gdf = df[["OriginID", "DestinationID", "duration", "distance"]]
            gdf["duration"] = gdf.duration / 60
        else:
            if self.osrmmodel == "Route":
                df1 = df[["OriginID", "DestinationID"]]
            else:
                df1 = df[["OriginID", "DestinationID", "duration", "distance"]]
                df1["duration"] = df1.duration / 60
            gdf = gp.GeoDataFrame(df1, geometry=df.geometry, crs=4326)
        return knext.Table.from_pandas(gdf)
