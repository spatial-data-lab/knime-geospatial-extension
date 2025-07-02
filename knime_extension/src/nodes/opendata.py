import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


__category = knext.category(
    path="/community/geo",
    level_id="opendataset",
    name="Open Datasets",
    description="Nodes that provide access to various public datasets, such as OpenStreetMap and US census data.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/OpendatasetCategory.png",
    after="LocationAnalysis",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/OpenDataset/"


############################################
# US2020 TIGER/Line for states
############################################
@knext.node(
    name="US2020 TIGER Map",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "UStiger.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="TIGER Line Table ",
    description="Retrieved geodata from the [United States Census Bureau](www2.census.gov/geo/tiger/TIGER2020PL)",
)
@knut.census_node_description(
    short_description="Retrieve geospatial data from US Census TIGER/Line",
    description="""This node retrieves the specific geospatial boundaries for one specific state of the United States.
        The popular TIGER/Line levels are Block group, Roads, Blocks, Tracts.
        When the same State FIPS (2-digits) or * is used for County FIPS(3-digits), the geodata of all counties in the state will be retrieved,
        county10/20 and state10/20 can only be applicable in this case. 
        This node can help user to get the FIPS codes of target study area. Linking it to Geospatial View node will be more helpful. """,
    references={
        "TIGER/Line Shapefiles": "https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.2020.html",
        "FTP Archive by State": "https://www2.census.gov/geo/tiger/TIGER2020PL/STATE/",
        "FTP Archive by Layer": "https://www2.census.gov/geo/tiger/TIGER2020PL/LAYER/",
        "FIPS code list": "https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt",
    },
)
class US2020TIGERNode:
    StateFips = knext.StringParameter(
        label="State FIPS (2-digits)",
        description="The State to use [FIPS.](https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt)",
        default_value="",
    )

    County3Fips = knext.StringParameter(
        "County FIPS(3-digits)/ State FIPS or * ",
        "The County/State FIPS code to use [FIPS.](https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt)",
        "",
    )

    geofile = knext.StringParameter(
        label="TIGER/Line data type",
        description="""Available TIGER/Line are: 
        
        - **Block:** tabblock10, tabblock20, for US Census block (the minimum level) in 2010 and 2020.
        - **Block group:** bg10, bg20, for US Census block group in 2010 and 2020.
        - **Tract:** tract10, tract20, for US Census Tract in 2010 and 2020.
        - **County:** county10, county20, for US County in 2010 and 2020, only applicable for the setting when the input of County FIPs equals State FIPs.
        - **State:** state10, state20, for US State in 2010 and 2020,same as County.
""",
        default_value="",
        enum=[
            "tabblock10",
            "tabblock20",
            "bg10",
            "bg20",
            "tract10",
            "tract20",
            "county10",
            "county20",
            "state10",
            "state20",
            "roads",
        ],
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        import pandas as pd

        USdict = pd.DataFrame.from_dict(
            {
                "state": [
                    "01",
                    "02",
                    "04",
                    "05",
                    "06",
                    "08",
                    "09",
                    "10",
                    "11",
                    "12",
                    "13",
                    "15",
                    "16",
                    "17",
                    "18",
                    "19",
                    "20",
                    "21",
                    "22",
                    "23",
                    "24",
                    "25",
                    "26",
                    "27",
                    "28",
                    "29",
                    "30",
                    "31",
                    "32",
                    "33",
                    "34",
                    "35",
                    "36",
                    "37",
                    "38",
                    "39",
                    "40",
                    "41",
                    "42",
                    "44",
                    "45",
                    "46",
                    "47",
                    "48",
                    "49",
                    "50",
                    "51",
                    "53",
                    "54",
                    "55",
                    "56",
                    "72",
                ],
                "stateName": [
                    "ALABAMA",
                    "ALASKA",
                    "ARIZONA",
                    "ARKANSAS",
                    "CALIFORNIA",
                    "COLORADO",
                    "CONNECTICUT",
                    "DELAWARE",
                    "DISTRICT_OF_COLUMBIA",
                    "FLORIDA",
                    "GEORGIA",
                    "HAWAII",
                    "IDAHO",
                    "ILLINOIS",
                    "INDIANA",
                    "IOWA",
                    "KANSAS",
                    "KENTUCKY",
                    "LOUISIANA",
                    "MAINE",
                    "MARYLAND",
                    "MASSACHUSETTS",
                    "MICHIGAN",
                    "MINNESOTA",
                    "MISSISSIPPI",
                    "MISSOURI",
                    "MONTANA",
                    "NEBRASKA",
                    "NEVADA",
                    "NEW_HAMPSHIRE",
                    "NEW_JERSEY",
                    "NEW_MEXICO",
                    "NEW_YORK",
                    "NORTH_CAROLINA",
                    "NORTH_DAKOTA",
                    "OHIO",
                    "OKLAHOMA",
                    "OREGON",
                    "PENNSYLVANIA",
                    "RHODE_ISLAND",
                    "SOUTH_CAROLINA",
                    "SOUTH_DAKOTA",
                    "TENNESSEE",
                    "TEXAS",
                    "UTAH",
                    "VERMONT",
                    "VIRGINIA",
                    "WASHINGTON",
                    "WEST_VIRGINIA",
                    "WISCONSIN",
                    "WYOMING",
                    "PUERTO_RICO",
                ],
            }
        )

        USdict["path"] = USdict.state + "_" + USdict.stateName

        Statepath = USdict[USdict.state == self.StateFips].iloc[0, 2]

        County5Fips = self.StateFips + self.County3Fips

        base_url = "ftp://ftp2.census.gov/geo/tiger/TIGER2020PL/STATE/"

        if self.StateFips != self.County3Fips and self.County3Fips != "*":
            data_url = f"{base_url}{Statepath}/{County5Fips}/tl_2020_{County5Fips}_{self.geofile}.zip"
        else:
            County5Fips = self.StateFips
            if self.geofile == "roads":
                self.geofile = "prisecroads"
            data_url = f"{base_url}{Statepath}/{County5Fips}/tl_2020_{County5Fips}_{self.geofile}.zip"
        gdf = gp.read_file(data_url)
        gdf.reset_index(drop=True, inplace=True)
        return knext.Table.from_pandas(gdf)


############################################
# US2020 Census Data
############################################
@knext.node(
    name="US2020 Census Data",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "UScensus.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="Census data table",
    description="Retrieved geodata from 2020 Census Redistricting Data",
)
@knut.census_node_description(
    short_description="Retrieve US 2020 Census Redistricting Data for one specific state of United States.",
    description="""This node retrieves US 2020 Census Redistricting Data 
(Decennial Census P.L. 94-171 Redistricting Data). This node provides all variables such as population and 
household information from the US Census 2020 data base. The default variable names are GEO_ID (geography), 
P1_001N (Total Population), P1_003N (Population of one race:!!White alone), P1_004N (Black or African American alone), 
H1_001N (Total Housing Units), H1_002 (Total Occupied Housing Units). 

Only if county is chosen for geography, then * can be input in State FIPS (2-digits) to retrieve all the 
county level data of all states. Each query can include **at most 50 variables**.

You may use this node **without an API key**. However, anonymous access is **limited to 500 requests per day per IP**.  
To avoid throttling and improve performance (especially for batch queries), it is recommended to register and provide a Census API key.

Before using the node, user need to sign up and get a Census API key first by clicking 
[here.](https://api.census.gov/data/key_signup.html)
    """,
    references={
        "Census API Key Sign Up": "https://api.census.gov/data/key_signup.html",
        "Decennial Census P.L. 94-171 Redistricting Data": "https://www.census.gov/programs-surveys/decennial-census/about/rdo/summary-files.html",
        "Datasets and its descendants": "https://api.census.gov/data/2020/dec.html",
        "Geography": "https://api.census.gov/data/2020/dec/pl/geography.html",
        "Variables": "https://api.census.gov/data/2020/dec/pl/variables.html",
        "FIPS code list": "https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt",
        "Census API examples": "https://api.census.gov/data/2020/dec/pl/examples.html",
    },
)
class USCensus2020Node:
    StateFips = knext.StringParameter(
        "State FIPS (2-digits)",
        "The State [FIPS](https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt) to investigate, input * for all states (while choose county for geography).",
        "25",
    )

    County3Fips = knext.StringParameter(
        "County FIPS (3-digits)",
        "The County [FIPS](https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt) to investigate, input * for all counties.",
        "017",
    )

    Tract6Fips = knext.StringParameter(
        "Tract FIPS (6-digits)",
        "The Tract to investigate, input * for all tracts.",
        "*",
    )

    censusapikey = knext.StringParameter(
        label="US Census API key (optional)",
        description="Optional Census API key to increase query limits and reliability.",
        default_value="",
        is_advanced=True,
    )

    cols = knext.StringParameter(
        "US Census Variable Names",
        "The [US Census Variable](https://api.census.gov/data/2020/dec/pl/variables.html) list to use.",
        "GEO_ID,P1_001N,P1_003N,P1_004N,H1_001N,H1_002N",
    )

    geofile = knext.StringParameter(
        label="Geographic Level",
        description=" Available geographic level for this node are Block, Block group, Tract and County.",
        default_value="block group",
        enum=["block group", "block", "tract", "county"],
    )
    compatibleid = knext.BoolParameter(
        label="Make Census GEOIDs compatible to Tiger/Line GEOIDs",
        description="""If enabled, the first nine characters of the Census GEOID column are removed to be compatible 
with the GEOID of the Tiger/Line Shapefiles. This is useful if you want to join the result of this node with the 
result of the
[US2020 TIGER Map node.](https://hub.knime.com/spatialdatalab/extensions/sdl.harvard.features.geospatial/latest/org.knime.python3.nodes.extension.ExtensionNodeSetFactory$DynamicExtensionNodeFactory:6009c0c2/) 

For more details about the format and the conversion see 
[Understanding Geographic Identifiers.](https://www.census.gov/programs-surveys/geography/guidance/geo-identifiers.html)""",
        default_value=lambda v: False if v < knext.Version(1, 1, 0) else True,
        since_version="1.1.0",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        base_url = "https://api.census.gov/data/2020/dec/pl?get="

        if self.geofile == "county":
            data_url = f"{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&key={self.censusapikey}"
        elif self.geofile == "tract":
            data_url = f"{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&in=county:{self.County3Fips}&key={self.censusapikey}"
        else:
            data_url = f"{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&in=county:{self.County3Fips}&in=tract:{self.Tract6Fips}&key={self.censusapikey}"

        import requests

        response = requests.get(data_url)
        data = response.json()

        import pandas as pd

        gdf = pd.DataFrame(data[1:], columns=data[0])
        if "GEO_ID" in gdf.columns and self.compatibleid == True:
            gdf["GEO_ID"] = gdf["GEO_ID"].str.slice(start=9)
        return knext.Table.from_pandas(gdf)


############################################
# US Census ACS-5
############################################
@knext.node(
    name="US ACS 5-Year Estimates",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "CensusACS.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="US Census ACS 5-Year table",
    description="Retrieved data from Census ACS Datasets",
)
@knut.census_node_description(
    short_description="Retrieves American Community Survey 5-Year Data (2009-2020) of United States.",
    description="""This node retrieves American Community Survey 5-Year Data (2009-2020).
The American Community Survey (ACS) is an ongoing survey that provides data every year -- giving 
communities the current information they need to plan investments and services. The ACS covers a broad 
range of topics about social, economic, demographic, and housing characteristics of the U.S. population.
The 5-year estimates from the ACS are period estimates that represent data collected over a period. 
The primary advantage of using multiyear estimates is the increased statistical reliability of the data 
for less populated areas and small population subgroups. 
Only if county is chosen for geography, then * can be input in State FIPS (2-digits) to retrieve all the 
county level data of all states. Each query can include **at most 50 variables**.

You may use this node **without an API key**. However, anonymous access is **limited to 500 requests per day per IP**.  
To avoid throttling and improve performance (especially for batch queries), it is recommended to register and provide a Census API key.

The 5-year estimates are available for all geographies down to the block group level.
    """,
    references={
        "Census API Key Sign Up": "https://api.census.gov/data/key_signup.html",
        "American Community Survey 5-Year Data (2009-2020)": "https://www.census.gov/data/developers/data-sets/acs-5year.html",
        "Geography": "https://api.census.gov/data/2020/acs/acs5/geography.html",
        "Variables": "https://api.census.gov/data/2020/acs/acs5/variables.html",
        "Census API examples": "https://api.census.gov/data/2020/acs/acs5/subject/examples.html",
        "FIPS code list": "https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt",
    },
)
class UScensusACSNode:
    StateFips = knext.StringParameter(
        "State FIPS (2-digits)",
        "The State [FIPS](https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt) to investigate, input * for all states (while choose county for geography).",
        "25",
    )

    County3Fips = knext.StringParameter(
        "County FIPS (3-digits)",
        "The County [FIPS](https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt) to investigate, input * for all counties.",
        "017",
    )

    Tract6Fips = knext.StringParameter(
        "Tract FIPS (6-digits)",
        "The Tract to investigate, input * for all tracts.",
        "*",
    )

    censusapikey = knext.StringParameter(
        label="US Census API key (optional)",
        description="Optional Census API key to increase query limits and reliability.",
        default_value="",
        is_advanced=True,
    )

    cols = knext.StringParameter(
        "US Census ACS Variable Names",
        "The US Census Variable list to use.",
        "GEO_ID,B02001_001E,B02001_002E,B02001_003E",
    )

    geofile = knext.StringParameter(
        label="Geographic Level",
        description=" Available geographic Level for this node are Block, Block group, Tract and County.",
        default_value="block group",
        enum=[
            "block group",
            "tract",
            "county",
            "state",
        ],
    )

    year = knext.StringParameter(
        "US Census ACS5 Year Label", "The Year label of dataset.", "2020"
    )
    compatibleid = knext.BoolParameter(
        label="Make GEO_ID compatible to Tiger/Line GEOID ",
        description=" FIPS-based GEOID for Block, Block group, Tract and County.",
        default_value=lambda v: False if v <= knext.Version(1, 0, 0) else True,
        since_version="1.1.0",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        base_url = "https://api.census.gov/data/"
        Dataset = "acs/acs5"

        if self.geofile == "state":
            data_url = f"{base_url}{self.year}/{Dataset}?get={self.cols}&for=state:{self.StateFips}&key={self.censusapikey}"
        elif self.geofile == "county":
            data_url = f"{base_url}{self.year}/{Dataset}?get={self.cols}&for=county:{self.County3Fips}&in=state:{self.StateFips}&key={self.censusapikey}"
        elif self.geofile == "tract":
            data_url = f"{base_url}{self.year}/{Dataset}?get={self.cols}&for=tract:{self.Tract6Fips}&in=state:{self.StateFips}&in=county:{self.County3Fips}&key={self.censusapikey}"
        else:
            data_url = f"{base_url}{self.year}/{Dataset}?get={self.cols}&for=block%20group:*&in=state:{self.StateFips}&in=county:{self.County3Fips}&in=tract:{self.Tract6Fips}&key={self.censusapikey}"

        import requests

        response = requests.get(data_url)
        data = response.json()

        import pandas as pd

        gdf = pd.DataFrame(data[1:], columns=data[0])
        if "GEO_ID" in gdf.columns and self.compatibleid == True:
            gdf["GEO_ID"] = gdf["GEO_ID"].str.slice(start=9)
        return knext.Table.from_pandas(gdf)


############################################
# OSM nodes
############################################
# TODO:Use this global temp dir for the OSMNX cache file until we have a way to get the KNIME workspace location
osm_root_dir = None


def get_osmnx():
    """
    Initializes and returns the osmnx module
    """
    import os
    import osmnx as ox
    import tempfile

    # use global variable to reuse the cache in subsequent calls of the method
    global osm_root_dir
    if osm_root_dir is None:
        osm_root_dir = tempfile.gettempdir()
        cache_dir = os.path.join(osm_root_dir, "knime_osmnx_cache")
        ox.settings.use_cache = True
        ox.settings.cache_folder = cache_dir
    return ox


############################################
# OSM Data Tool
############################################
@knext.node(
    name="OSM POIs",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "OSMpoi.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Polygon Data",
    description="Table with geometry used for downloading POIs",
)
@knext.output_table(
    name="OSM POI data",
    description="POI Geodata from the Open Street Map",
)
@knut.osm_node_description(
    short_description="Get Points of Interests(POIs) from the Open Street Map.",
    description="""This node downloads geospatial entities’ geometries and attributes from 
[OpenStreetMap.](https://www.openstreetmap.org/about) Results returned are the union, not intersection of each 
individual tag. Each result matches at least one given tag. The place tags should be OSM tags, 
(e.g., building, landuse, highway, etc) and the value tags should be either *True* to retrieve all items with the given 
tag, or a single value to retrieve a single tag-value combination, or a comma separated list of values to get multiple 
values for the given tag. For example, *place tag=building, tag value=True* would return all building footprints in the 
area. *Place tag=landuse, tag value=retail, commercial* would return all retail and commercial landuses. For more 
details about a tag and how to find valid tags, please refer to the [OpenStreet TagFinder.](https://tagfinder.osm.ch/)

Please be aware that tags can change over time. For more details about changes in tags, please refer to the
[Changelog.](https://wiki.openstreetmap.org/wiki/Changelog)
""",
    references={
        "OpenStreet TagFinder": "https://tagfinder.osm.ch/",
        "OpenStreetMap Taginfo": "https://taginfo.openstreetmap.org/",
        "OSM Map Features": "https://wiki.openstreetmap.org/wiki/Map_features",
        "OSMnx": "https://github.com/gboeing/osmnx",
        "osmnx.geometries_from_place": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.geometries",
    },
)
class OSMdataNode:
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column as boundary to get POIs.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_polygon,
        include_row_key=False,
        include_none_column=False,
    )

    taginfo = knext.StringParameter(
        "Input place tags",
        "The value for tag to specify the type of facilities.",
        "amenity",
    )

    tagvalue = knext.StringParameter(
        "Input value tags",
        "The specific type of facilities, True for all.",
        "restaurant",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf.to_crs(4326, inplace=True)
        gdf_union = gdf.unary_union
        if self.tagvalue != "True":
            tags = {self.taginfo: self.tagvalue}
        else:
            tags = {self.taginfo: True}
        # tags = {self.taginfo: self.tagvalue}

        gdfpoi = get_osmnx().features.features_from_polygon(gdf_union, tags)
        gdfpoi = gdfpoi.reset_index(drop=True)
        return knext.Table.from_pandas(gdfpoi)


############################################
# OSM Network
############################################
@knext.node(
    name="OSM Road Network",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "OSMnetwork.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Polygon Data",
    description="Table with geometry used for downloading street network",
)
@knext.output_table(
    name="OSM Network data",
    description="Road Network Geodata from the Open Street Map",
)
@knut.osm_node_description(
    short_description="Get Road Network from the Open Street Map.",
    description="""This node downloads a geospatial network and its attributes from 
    [OpenStreetMap.](https://www.openstreetmap.org/about) 
    If the street network type 'drive' is selected, the node will append the speed information to the result table and
    the total travel time for each segment will be calculated. 
""",
    references={
        "OSMnx": "https://github.com/gboeing/osmnx",
        "osmnx.graph.graph_from_polygon": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.graph",
        "osmnx.speed.add_edge_speeds": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.speed ",
    },
)
class OSMnetworkNode:
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column as boundary to get POIs.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    networktype = knext.StringParameter(
        label="Street network type",
        description="Type of street network ",
        default_value="drive",
        enum=[
            "all",
            "all_private",
            "bike",
            "drive",
            "drive_service",
            "walk",
        ],
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf.to_crs(4326, inplace=True)
        gdf_union = gdf.unary_union
        # check geometry type
        geom_type = gdf_union.geom_type
        if geom_type not in ["Polygon", "MultiPolygon"]:
            raise RuntimeError("Input data must be Polygon or MultiPolygon")
        ox = get_osmnx()
        knut.check_canceled(exec_context)
        G = ox.graph.graph_from_polygon(gdf_union, network_type=self.networktype)
        edges = ox.convert.graph_to_gdfs(G, nodes=False)
        objcolumn = edges.select_dtypes(include=["object"]).columns.tolist()

        # Convert each element of the dataframe to a string but sort lists before doing so
        def convert_to_string(v):
            if type(v) == list:
                v.sort()
            return str(v)

        edges[objcolumn] = edges[objcolumn].applymap(convert_to_string)
        edges = edges.reset_index(drop=True)
        return knext.Table.from_pandas(edges)


############################################
# OSM Geocode Boundary
############################################
@knext.node(
    name="OSM Boundary Map",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "OSMboundary.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="OSM GeoBoundary data",
    description="Boundary of places from the Open Street Map",
)
@knut.osm_node_description(
    short_description="Get Boundary from OpenStreetMap with Geocoding.",
    description="""This node gets place boundary from [OpenStreetMap](https://www.openstreetmap.org/about) by the 
    geocoding place name. The queries you provide must be resolvable to places in the 
    [Nominatim database.](https://nominatim.org/) The resulting GeoDataFrame’s geometry column contains place 
    boundaries if they exist in OpenStreetMap.
""",
    references={
        "OSMnx": "https://github.com/gboeing/osmnx",
        "osmnx.geocoder.geocode_to_gdf": "https://osmnx.readthedocs.io/en/stable/osmnx.html?highlight=geocode_to_gdf#module-osmnx.geocoder",
    },
)
class OSMGeoBoundaryNode:
    placename = knext.StringParameter(
        label="Input place names",
        description="""Hierarchical place names from specific to more general delimited with commas, such as 
        Cambridge, MA, USA.  """,
        default_value="Cambridge, MA, USA",
    )
    ignore_error = knext.BoolParameter(
        "Ignore error",
        "If selected, it will return None for unknown places. "
        "If not selected, the node will fail when an unknown place is encountered.",
        default_value=False,
        since_version="1.4.0",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        import geopandas as gp

        try:
            gdf = get_osmnx().geocoder.geocode_to_gdf(self.placename)
        except Exception as e:
            if not self.ignore_error:
                raise RuntimeError(
                    f"Failed to process place name '{self.placename}': {str(e)}"
                )
            return knext.Table.from_pandas(
                gp.GeoDataFrame(geometry=[], crs="EPSG:4326")
            )

        if not gdf.empty:
            gdf = gdf.reset_index(drop=True)
            gdf = gp.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")
            return knext.Table.from_pandas(gdf)

        if not self.ignore_error:
            raise RuntimeError(f"No boundary found for place: {self.placename}")
        return knext.Table.from_pandas(gp.GeoDataFrame(geometry=[], crs="EPSG:4326"))


############################################
# OSM Geocode Boundary Table
############################################
@knext.node(
    name="OSM Boundary Map Table",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "OSMboundary.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input table",
    description="Table with place names to get boundaries for.",
)
@knext.output_table(
    name="OSM GeoBoundary data",
    description="Boundary of places from the Open Street Map",
)
@knut.osm_node_description(
    short_description="Get Boundary from OpenStreetMap with Geocoding.",
    description="""This node gets place boundaries from [OpenStreetMap](https://www.openstreetmap.org/about) by 
   geocoding place names from the input table. The queries must be resolvable to places in the
   [Nominatim database.](https://nominatim.org/) The resulting GeoDataFrame's geometry column contains place
   boundaries if they exist in OpenStreetMap.""",
    references={
        "OSMnx": "https://github.com/gboeing/osmnx",
        "osmnx.geocoder.geocode_to_gdf": "https://osmnx.readthedocs.io/en/stable/osmnx.html?highlight=geocode_to_gdf#module-osmnx.geocoder",
    },
)
class OSMGeoBoundaryTableNode:
    place_col = knext.ColumnParameter(
        "Place name column",
        "Select the column containing place names to geocode. "
        "Names should be hierarchical from specific to general, e.g. 'Cambridge, MA, USA'",
        column_filter=knut.is_string,
        include_row_key=False,
        include_none_column=False,
    )

    ignore_error = knext.BoolParameter(
        "Ignore error",
        "If selected, it will return a missing cell for unknown places. "
        "If not selected, the node will fail when an unknown place is encountered.",
        default_value=True,
    )

    def configure(self, configure_context, input_schema):
        self.place_col = knut.column_exists_or_preset(
            configure_context, self.place_col, input_schema, knut.is_string
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        import pandas as pd
        import geopandas as gp

        # Get original data
        input_df = input_table.to_pandas()
        ox = get_osmnx()

        # Create new DataFrame to store OSM results
        osm_results = []
        total = len(input_df)

        for idx, place_name in enumerate(input_df[self.place_col]):
            if pd.isna(place_name):
                osm_results.append(pd.Series())
                continue

            try:
                gdf = ox.geocoder.geocode_to_gdf(place_name)
                if not gdf.empty:
                    osm_results.append(gdf.iloc[0])
                else:
                    osm_results.append(pd.Series())
                    if not self.ignore_error:
                        raise RuntimeError(f"No boundary found for place: {place_name}")

            except Exception as e:
                osm_results.append(pd.Series())
                if not self.ignore_error:
                    raise RuntimeError(
                        f"Failed to process place name '{place_name}': {str(e)}"
                    )
                else:
                    knut.LOGGER.warning(
                        f"Failed to process place name '{place_name}': {str(e)}"
                    )

            exec_context.set_progress(
                0.9 * (idx + 1) / total, f"Processing {idx + 1} of {total}"
            )
            knut.check_canceled(exec_context)

        # Convert OSM results to DataFrame
        osm_df = pd.DataFrame(osm_results)

        # Reset index to ensure proper alignment
        input_df = input_df.reset_index(drop=True)
        osm_df = osm_df.reset_index(drop=True)

        # Combine original data with OSM results
        result_df = pd.concat([input_df, osm_df], axis=1)

        # Convert to GeoDataFrame
        result_gdf = gp.GeoDataFrame(result_df, geometry="geometry", crs="EPSG:4326")

        if len(result_gdf) == 0:
            raise RuntimeError(
                "No valid boundaries found for any of the input place names"
            )

        return knext.Table.from_pandas(result_gdf)


# moved from spatialdata extension


############################################
# GDELT nodes
############################################
@knext.node(
    name="GDELT Global Knowledge Graph",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "GDELT.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="GDELT GKG Data Table",
    description="Retrieved geodata from GDELT Global Knowledge Graph. For details on the result columns see the "
    + "[GDELT documentation.](https://blog.gdeltproject.org/announcing-our-first-api-gkg-geojson/)",
)
class GDELTGKGNode:
    """This node retrieves GDELT Global Knowledge Graph data.
    The GDELT Global Knowledge Graph (GKG) is a real-time knowledge graph of global human society for open research.
    The GKG is a massive archive of global news and translated into 65 languages, updated every 15 minutes.
    The GKG is a network diagram of the world's events and news coverage, containing more than 1.5 billion people,
    organizations, locations, themes, emotions, counts, quotes, images and events across the planet
    dating back to January 1, 1979, and updated every 15 minutes.

    The node generates queries in the following form:
    *https://api.gdeltproject.org/api/v1/gkg_geojson?QUERY=<QUERY>&TIMESPAN=<Last Hours \* 60>*.
    Please refer to the [GDELT documentation](https://blog.gdeltproject.org/announcing-our-first-api-gkg-geojson/)
    for more details about the supported input values and output field definitions.

    ##Note
    Data copyright by [GDELT Project](https://www.gdeltproject.org/) and provided under the following
    [Terms of Use](https://www.gdeltproject.org/about.html#termsofuse).
    """

    query = knext.StringParameter(
        label="Query",
        description="The query to search in GDELT GKG. This can be a single key word (click "
        + "[here](http://data.gdeltproject.org/documentation/GKG-MASTER-THEMELIST.TXT) for a complete list) or a "
        + "more complex query (click [here](https://blog.gdeltproject.org/announcing-our-first-api-gkg-geojson/) "
        + "for more info and examples).",
        default_value="FOOD_SECURITY",
    )

    last_hours = knext.IntParameter(
        label="Last Hours",
        description="The last hours to search in GDELT GKG.",
        default_value=1,
        min_value=1,
        max_value=24,
    )

    timeout = knext.IntParameter(
        label="Request timeout in seconds",
        description="The timeout in seconds for the request for GDELT GKG.",
        default_value=120,
        min_value=1,
        is_advanced=True,
    )

    def configure(self, configure_context):
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        import geopandas as gp
        import requests

        url = "https://api.gdeltproject.org/api/v1/gkg_geojson?QUERY=%s&TIMESPAN=%d" % (
            self.query,
            self.last_hours * 60,
        )
        response = requests.get(url, timeout=self.timeout)
        data = response.json()
        gdf = gp.GeoDataFrame.from_features(data, crs="EPSG:4326")
        return knext.Table.from_pandas(gdf)


############################################
# Open Sky Network Data Node
############################################
@knext.node(
    name="Open Sky Network Data",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "OpenSkyNetwork.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="Open Sky Network Data Table",
    description="Retrieved [state vectors](https://openskynetwork.github.io/opensky-api/index.html#state-vectors) "
    + "with geodata from Open Sky Network Data. For details on the result columns see the "
    + "[OpenSky documentation.](https://openskynetwork.github.io/opensky-api/rest.html#response)",
)
class OpenSkyNetworkDataNode:
    """This node retrieves Open Sky Network Data.
    It returns live airspace information in the form of
    [state vectors](https://openskynetwork.github.io/opensky-api/index.html#state-vectors)
    for __research and non-commercial__ purposes. It does not provide commercial flight data such as airport schedules,
    delays or similar information that cannot be derived from ADS-B data contents!

    The OpenSky Network is a non-profit association based in Switzerland that operates a crowd sourced global
    database of air traffic control data. The network consists of thousands of sensors connected to the Internet
    by volunteers, whose main purpose is to measure the radio signals emitted by aircraft to track their position.
    Please refer to [Open Sky Network homepage](https://opensky-network.org/) for more details.

    ##Note
    Data copyright by [The OpenSky Network](https://opensky-network.org/) and provided under the following
    [Terms of Use](https://opensky-network.org/index.php/about/terms-of-use).
    The following section is copied from the terms of use:

    *OpenSky Network’s authorization to access the data grants You a limited, non-exclusive, non-transferable,
    non-assignable, and terminable license to copy, modify, and use the data in accordance with this AGREEMENT
    __solely for the purpose of non-profit research, non-profit education, or for government purposes__.
    No license is granted for any other purpose and there are no implied licenses in this AGREEMENT.
    In particular, any use by a for-profit entity requires written permission by the OpenSky Network.*
    """

    user = knext.StringParameter(
        label="User (optional)",
        description="The optional user name to access Open Sky Network Data. "
        + "If not provided [limitations](https://openskynetwork.github.io/opensky-api/rest.html#limitations) apply.",
        default_value="",
    )

    password = knext.StringParameter(
        label="Password (optional)",
        description="The optional password to access Open Sky Network Data."
        + "If not provided [limitations](https://openskynetwork.github.io/opensky-api/rest.html#limitations) apply.",
        default_value="",
    )

    timeout = knext.IntParameter(
        label="Request timeout in seconds",
        description="The timeout in seconds for the request for GDELT GKG.",
        default_value=120,
        min_value=1,
        is_advanced=True,
    )

    def configure(self, configure_context):
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        import geopandas as gp
        import pandas as pd
        import requests

        # API documentation https://openskynetwork.github.io/opensky-api/rest.html

        url = "https://opensky-network.org/api/states/all"
        kws = {"url": url, "timeout": self.timeout}
        if len(self.user) != 0 and len(self.password) != 0:
            kws["auth"] = (self.user, self.password)

        response = requests.get(**kws)
        json_data = response.json()
        states = pd.DataFrame(
            json_data["states"],
            columns=[
                "icao24",
                "callsign",
                "origin_country",
                "time_position",
                "last_contact",
                "longitude",
                "latitude",
                "baro_altitude",
                "on_ground",
                "velocity",
                "true_track",
                "vertical_rate",
                "sensors",
                "geo_altitude",
                "squawk",
                "spi",
                "position_source",
            ],
        )
        gdf = gp.GeoDataFrame(
            states,
            geometry=gp.points_from_xy(states.longitude, states.latitude),
            crs="EPSG:4326",
        )
        return knext.Table.from_pandas(gdf)


############################################
# NATURAL EARTH Node Dictionary
############################################


# Data dictionaries for file information
ne_dic = {
    "all_scale": {
        "countries": {
            "file": "admin_0_countries",
            "display": "Country boundaries",
            "description": "Countries as polygons with attributes",
        },
        "countries_lakes": {
            "file": "admin_0_countries_lakes",
            "display": "Country boundaries with lakes",
            "description": "Countries with major lakes as polygons",
        },
        "sovereignty": {
            "file": "admin_0_sovereignty",
            "display": "Sovereignty boundaries",
            "description": "Sovereign states as polygons",
        },
        "map_units": {
            "file": "admin_0_map_units",
            "display": "Map units boundaries",
            "description": "Map units as polygons",
        },
        "scale_rank": {
            "file": "admin_0_scale_rank",
            "display": "Scale rank boundaries",
            "description": "Scale ranked boundaries as polygons",
        },
        "boundary_lines_land": {
            "file": "admin_0_boundary_lines_land",
            "display": "Boundary lines land",
            "description": "Country boundary lines on land",
        },
        "states_provinces": {
            "file": "admin_1_states_provinces",
            "display": "States and provinces",
            "description": "States and provinces as polygons",
        },
        "states_provinces_lakes": {
            "file": "admin_1_states_provinces_lakes",
            "display": "States and provinces with lakes",
            "description": "States and provinces with lakes as polygons",
        },
        "states_provinces_lines": {
            "file": "admin_1_states_provinces_lines",
            "display": "States and provinces lines",
            "description": "State and province boundary lines",
        },
        "populated_places": {
            "file": "populated_places",
            "display": "Populated places",
            "description": "Major populated places as points",
        },
        "populated_places_simple": {
            "file": "populated_places_simple",
            "display": "Populated places (simple)",
            "description": "Simplified populated places as points",
        },
        "coastline": {
            "file": "coastline",
            "display": "Coastlines",
            "description": "Coastlines as lines",
        },
        "land": {
            "file": "land",
            "display": "Land areas",
            "description": "Land masses as polygons",
        },
        "ocean": {
            "file": "ocean",
            "display": "Ocean areas",
            "description": "Ocean areas as polygons",
        },
        "rivers_lake_centerlines": {
            "file": "rivers_lake_centerlines",
            "display": "Rivers and lake centerlines",
            "description": "Major rivers and lake centerlines as lines",
        },
        "lakes": {
            "file": "lakes",
            "display": "Lakes",
            "description": "Major lakes as polygons",
        },
    },
    "pov_10m": {
        "argentina": {
            "file": "admin_0_countries_arg",
            "display": "Argentina",
            "description": "Argentina country boundary and attributes",
        },
        "bangladesh": {
            "file": "admin_0_countries_bdg",
            "display": "Bangladesh",
            "description": "Bangladesh country boundary and attributes",
        },
        "brazil": {
            "file": "admin_0_countries_bra",
            "display": "Brazil",
            "description": "Brazil country boundary and attributes",
        },
        "china": {
            "file": "admin_0_countries_chn",
            "display": "China",
            "description": "China country boundary and attributes",
        },
        "egypt": {
            "file": "admin_0_countries_egy",
            "display": "Egypt",
            "description": "Egypt country boundary and attributes",
        },
        "france": {
            "file": "admin_0_countries_fra",
            "display": "France",
            "description": "France country boundary and attributes",
        },
        "germany": {
            "file": "admin_0_countries_deu",
            "display": "Germany",
            "description": "Germany country boundary and attributes",
        },
        "greece": {
            "file": "admin_0_countries_grc",
            "display": "Greece",
            "description": "Greece country boundary and attributes",
        },
        "indonesia": {
            "file": "admin_0_countries_idn",
            "display": "Indonesias",
            "description": "Indonesia country boundary and attributes",
        },
        "india": {
            "file": "admin_0_countries_ind",
            "display": "India",
            "description": "India country boundary and attributes",
        },
        "iso": {
            "file": "admin_0_countries_iso",
            "display": "ISO standard",
            "description": "ISO standard country boundaries and attributes",
        },
        "israel": {
            "file": "admin_0_countries_isr",
            "display": "Israel",
            "description": "Israel country boundary and attributes",
        },
        "italy": {
            "file": "admin_0_countries_ita",
            "display": "Italy",
            "description": "Italy country boundary and attributes",
        },
        "japan": {
            "file": "admin_0_countries_jpn",
            "display": "Japan",
            "description": "Japan country boundary and attributes",
        },
        "korea": {
            "file": "admin_0_countries_kor",
            "display": "Korea",
            "description": "Korea country boundary and attributes",
        },
        "morocco": {
            "file": "admin_0_countries_mar",
            "display": "Morocco",
            "description": "Morocco country boundary and attributes",
        },
        "nepal": {
            "file": "admin_0_countries_nep",
            "display": "Nepal",
            "description": "Nepal country boundary and attributes",
        },
        "netherlands": {
            "file": "admin_0_countries_nld",
            "display": "Netherlands",
            "description": "Netherlands country boundary and attributes",
        },
        "pakistan": {
            "file": "admin_0_countries_pak",
            "display": "Pakistan",
            "description": "Pakistan country boundary and attributes",
        },
        "poland": {
            "file": "admin_0_countries_pol",
            "display": "Poland",
            "description": "Poland country boundary and attributes",
        },
        "portugal": {
            "file": "admin_0_countries_prt",
            "display": "Portugal",
            "description": "Portugal country boundary and attributes",
        },
        "palestine": {
            "file": "admin_0_countries_pse",
            "display": "Palestine",
            "description": "Palestine country boundary and attributes",
        },
        "russia": {
            "file": "admin_0_countries_rus",
            "display": "Russia",
            "description": "Russia country boundary and attributes",
        },
        "saudi_arabia": {
            "file": "admin_0_countries_sau",
            "display": "Saudi Arabia",
            "description": "Saudi Arabia country boundary and attributes",
        },
        "spain": {
            "file": "admin_0_countries_esp",
            "display": "Spain",
            "description": "Spain country boundary and attributes",
        },
        "sweden": {
            "file": "admin_0_countries_swe",
            "display": "Sweden",
            "description": "Sweden country boundary and attributes",
        },
        "turkey": {
            "file": "admin_0_countries_tur",
            "display": "Turkey boundaries",
            "description": "Turkey country boundary and attributes",
        },
        "united_kingdom": {
            "file": "admin_0_countries_gbr",
            "display": "United Kingdom",
            "description": "United Kingdom country boundary and attributes",
        },
        "usa": {
            "file": "admin_0_countries_usa",
            "display": "United States",
            "description": "United States country boundary and attributes",
        },
        "ukraine": {
            "file": "admin_0_countries_ukr",
            "display": "Ukraine",
            "description": "Ukraine country boundary and attributes",
        },
        "vietnam": {
            "file": "admin_0_countries_vnm",
            "display": "Vietnam boundaries",
            "description": "Vietnam country boundary and attributes",
        },
    },
    "special_10m": {
        "roads": {
            "file": "roads",
            "display": "Global roads",
            "description": "Major roads and highways as linestrings globally",
        },
        "roads_north_america": {
            "file": "roads_north_america",
            "display": "North American roads",
            "description": "Detailed roads as linestrings in North America",
        },
        "railroads": {
            "file": "railroads",
            "display": "Global railroads",
            "description": "Major rail networks as linestrings globally",
        },
        "railroads_north_america": {
            "file": "railroads_north_america",
            "display": "North American railroads",
            "description": "Detailed rail networks as linestrings in North America",
        },
        "airports": {
            "file": "airports",
            "display": "Airports",
            "description": "Major airports as points with attributes",
        },
        "ports": {
            "file": "ports",
            "display": "Ports",
            "description": "Major ports as points with attributes",
        },
        "urban_areas": {
            "file": "urban_areas",
            "display": "Urban areas",
            "description": "Urban areas as polygons with population data",
        },
        "parks_and_protected_lands": {
            "file": "parks_and_protected_lands",
            "display": "Parks and protected lands",
            "description": "Protected areas and parks as polygons",
        },
    },
}

############################################
# NATURAL EARTH Node
############################################

# Constants
BASE_URL = "https://naciscdn.org/naturalearth/"


def is_physical_data(file_key):
    """Determine if the file is physical data based on its key."""
    physical_files = ["coastline", "land", "ocean", "rivers_lake_centerlines", "lakes"]
    return file_key in physical_files


class _DataCategory(knext.EnumParameterOptions):
    """Enumeration for Natural Earth data categories."""

    ALL_SCALE = (
        "All scale",
        "General scale data available in 1:10m, 1:50m, and 1:110m resolutions",
    )
    POV_10M = ("10M-POV", "Point of view (country-specific) data at 1:10m resolution")
    SPECIAL_10M = ("10M-Special", "Special features data at 1:10m resolution")

    @classmethod
    def get_default(cls):
        return cls.ALL_SCALE


class _ScaleOptions(knext.EnumParameterOptions):
    """Enumeration for scale options in All Scale category."""

    SCALE_10M = ("1:10m", "1:10 million scale, most detailed")
    SCALE_50M = ("1:50m", "1:50 million scale, medium detail")
    SCALE_110M = ("1:110m", "1:110 million scale, least detailed")

    @classmethod
    def get_default(cls):
        return cls.SCALE_10M


class _AllScaleFiles(knext.EnumParameterOptions):
    """Enumeration of all scale file options."""

    exec(
        "\n".join(
            [
                f"{key.upper()} = (ne_dic['all_scale']['{key}']['display'], ne_dic['all_scale']['{key}']['description'])"
                for key in ne_dic["all_scale"].keys()
            ]
        )
    )

    @classmethod
    def get_default(cls):
        return "COUNTRIES"


class _POVFiles(knext.EnumParameterOptions):
    """Enumeration of POV file options."""

    exec(
        "\n".join(
            [
                f"{key.upper()} = (ne_dic['pov_10m']['{key}']['display'], ne_dic['pov_10m']['{key}']['description'])"
                for key in ne_dic["pov_10m"].keys()
            ]
        )
    )

    @classmethod
    def get_default(cls):
        return "ARGENTINA"


class _SpecialFiles(knext.EnumParameterOptions):
    """Enumeration of special file options."""

    exec(
        "\n".join(
            [
                f"{key.upper()} = (ne_dic['special_10m']['{key}']['display'], ne_dic['special_10m']['{key}']['description'])"
                for key in ne_dic["special_10m"].keys()
            ]
        )
    )

    @classmethod
    def get_default(cls):
        return "ROADS"


@knext.node(
    name="Natural Earth Global Data",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "NaturalEarth.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="Natural Earth Data",
    description="Contains geospatial vector data retrieved from the Natural Earth dataset.",
)
class NaturalEarthNode:
    """Node for downloading and processing Natural Earth geographic data.

    This node provides direct access to [Natural Earth data,](https://www.naturalearthdata.com/)
    a public domain dataset offering detailed global geographic information at multiple map scales.
    It enables users to download and import vector data—such as countries, states, provinces, and other geographic
    features—directly into their KNIME workflows.

    This is particularly useful for creating choropleth maps, for example, to visualize population data
    across countries or states.

    All users of Natural Earth are highly encouraged to read about data sources and manipulation in the
    [Data Creation](https://www.naturalearthdata.com/about/data-creation/) section of the Natural Earth homepage.

    Natural Earth data is hosted by the
    [North American Cartographic Information Society (NACIS)](https://nacis.org/initiatives/natural-earth/),
    and is free for use in any type of project (see
    [Natural Earth Terms of Use](https://www.naturalearthdata.com/about/terms-of-use/) for more details).

    ##Note
    Made with Natural Earth. Free vector and raster map data @ naturalearthdata.com.
    """

    category = knext.EnumParameter(
        "Data category",
        "Select the category of Natural Earth data to download",
        default_value=_DataCategory.ALL_SCALE.name,
        enum=_DataCategory,
    )

    scale = knext.EnumParameter(
        "Scale",
        "Select the scale of the data (only applicable for All Scale category)",
        default_value=_ScaleOptions.SCALE_10M.name,
        enum=_ScaleOptions,
    ).rule(knext.OneOf(category, [_DataCategory.ALL_SCALE.name]), knext.Effect.SHOW)

    all_scale_file = knext.EnumParameter(
        "All scale file type",
        "Select the file to download for all scale category",
        default_value=next(iter(_AllScaleFiles)).name,
        enum=_AllScaleFiles,
    ).rule(knext.OneOf(category, [_DataCategory.ALL_SCALE.name]), knext.Effect.SHOW)

    pov_file = knext.EnumParameter(
        "10M-POV file type",
        "Select the country-specific file to download",
        default_value=next(iter(_POVFiles)).name,
        enum=_POVFiles,
    ).rule(knext.OneOf(category, [_DataCategory.POV_10M.name]), knext.Effect.SHOW)

    special_file = knext.EnumParameter(
        "10M-Special file type",
        "Select the special feature file to download",
        default_value=next(iter(_SpecialFiles)).name,
        enum=_SpecialFiles,
    ).rule(knext.OneOf(category, [_DataCategory.SPECIAL_10M.name]), knext.Effect.SHOW)

    def _build_download_url(self):
        """Build the download URL for Natural Earth data."""
        category_map = {
            _DataCategory.ALL_SCALE.name: ("all_scale", self.all_scale_file),
            _DataCategory.POV_10M.name: ("pov_10m", self.pov_file),
            _DataCategory.SPECIAL_10M.name: ("special_10m", self.special_file),
        }

        category_key, file_param = category_map[self.category]
        file_key = file_param.lower()
        file_info = ne_dic[category_key][file_key]

        if self.category == _DataCategory.ALL_SCALE.name:
            scale_prefix = self.scale.replace("SCALE_", "").lower()
        else:
            scale_prefix = "10m"

        subfolder = (
            "physical/"
            if (
                self.category == _DataCategory.ALL_SCALE.name
                and is_physical_data(file_key)
            )
            else "cultural/"
        )

        # Construct filename and URL
        filename = f"ne_{scale_prefix}_{file_info['file']}"
        url = f"{BASE_URL}{scale_prefix}/{subfolder}{filename}.zip"

        return url

    def configure(self, configure_context):
        """Configure the output table structure."""
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        """Execute the node to download and process Natural Earth data."""
        try:
            url = self._build_download_url()
            gdf = gp.read_file(url, encoding="utf-8")
            exec_context.set_progress(
                0.5,
                f"Downloading Natural Earth data from {url}. "
                + "Progress won't update until file is downloaded. Please wait...",
            )
            gdf.reset_index(drop=True, inplace=True)
            return knext.Table.from_pandas(gdf)
        except Exception as e:
            raise RuntimeError(
                f"Failed to download or process Natural Earth data: {str(e)}{url}"
            )


############################################
# Dataverse File Downloader
############################################


def validate_path(path: str) -> None:
    # no path check
    pass


class ExistingFile(knext.EnumParameterOptions):
    FAIL = (
        "Fail",
        "Will issue an error during the node's execution (to prevent unintentional overwrite).",
    )
    OVERWRITE = (
        "Overwrite",
        "Will replace any existing file.",
    )


@knext.node(
    name="Dataverse File Downloader",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "dataverse.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="Downloader File Path",
    description="Retrieved  data from Dataverse",
)
class DataverseFileDownloaderNode:
    """Downloads a file from a Dataverse repository.

    This node downloads a file from a Dataverse repository based on the provided File ID.
    The user must specify: Dataverse server URL (e.g., "https://dataverse.harvard.edu");
    File ID (a numerical identifier from Dataverse); Save location (local file path to save the downloaded file)
    The node retrieves the file and stores it at the user-defined folder.
    """

    server_url = knext.StringParameter(
        label="Dataverse Server URL",
        description="Base URL of the Dataverse server (e.g., https://dataverse.harvard.edu).",
        default_value="https://dataverse.harvard.edu",
        is_advanced=True,
    )

    file_id = knext.StringParameter(
        label="File ID",
        description="The unique file identifier in Dataverse.",
        default_value="",
    )

    save_path = knext.LocalPathParameter(
        label="Save Path",
        description="Select the directory to save the downloaded file.",
        placeholder_text="Select output directory...",
        validator=validate_path,
    )

    timeout = knext.IntParameter(
        label="Request Timeout (seconds)",
        description="Maximum time to wait for the server response.",
        default_value=120,
        min_value=1,
        is_advanced=True,
    )

    existing_file = knext.EnumParameter(
        "If exists:",
        "Specify the behavior of the node in case the output file already exists.",
        lambda v: (
            ExistingFile.OVERWRITE.name
            if v < knext.Version(1, 2, 0)
            else ExistingFile.FAIL.name
        ),
        enum=ExistingFile,
    )

    def configure(self, configure_context):
        return knext.Schema.from_columns([knext.Column(knext.string(), "File Path")])

    def execute(self, exec_context: knext.ExecutionContext):
        import requests
        import os
        import pandas as pd

        base_url = self.server_url.rstrip("/")
        download_url = f"{base_url}/api/access/datafile/{self.file_id}"
        self.__check_overwrite(self.save_path)
        try:
            save_dir = os.path.dirname(self.save_path)
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)

            response = requests.get(download_url, timeout=self.timeout)
            response.raise_for_status()

            with open(self.save_path, "wb") as file:
                file.write(response.content)

            output_table = pd.DataFrame({"File Path": [self.save_path]})

            return knext.Table.from_pandas(output_table)

        except Exception as e:
            raise ValueError(f"Download Error: {str(e)}")

    def __check_overwrite(self, fileurl):
        if self.existing_file == ExistingFile.FAIL.name:
            import os.path

            if os.path.exists(fileurl):
                raise knext.InvalidParametersError(
                    "File already exists and should not be overwritten."
                )


@knext.node(
    name="Dataverse Search",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "dataverse.png",
    category=__category,
    after="",
)
@knext.output_table(
    name="Search Results",
    description="Retrieved data from Dataverse",
)
class DataverseSearchNode:
    """Search for datasets and files in Dataverse repositories.

    This node allows you to search Dataverse using various parameters.

    Query Syntax Examples:
    1. Simple keyword search: "climate change"
       - Searches for items containing these terms anywhere

    2. Field-specific search: "title:climate"
       - Searches only in the title field
       - Other fields: author, description, keywords

    3. Boolean operators: "climate AND temperature"
       - AND: Both terms must be present
       - OR: Either term must be present
       - NOT: Exclude items with the term

    4. Combining operators: "climate AND (temperature OR rainfall)"
       - Use parentheses for complex queries

    Search results are returned as a table with all available metadata from the API.
    """

    server_url = knext.StringParameter(
        label="Dataverse Server URL",
        description="Base URL of the Dataverse server.",
        default_value="https://dataverse.harvard.edu",
        is_advanced=True,
    )

    query = knext.StringParameter(
        label="Search Query",
        description="Search keywords. Examples: climate, title:climate, climate AND temperature",
        default_value="",
    )

    search_type = knext.StringParameter(
        label="Search Type",
        description="Limit the type of objects to search for.",
        default_value="file",
        enum=["dataverse", "dataset", "file", "all"],
        is_advanced=True,
    )

    max_results = knext.IntParameter(
        label="Maximum Results",
        description="Maximum number of results to return",
        default_value=100,
        min_value=1,
        max_value=1000,
        is_advanced=True,
    )

    timeout = knext.IntParameter(
        label="Timeout (seconds)",
        description="Request timeout in seconds.",
        default_value=30,
        min_value=1,
        is_advanced=True,
    )

    def configure(self, configure_context):
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        import requests
        import pandas as pd

        if not self.query:
            raise ValueError("Search query must be provided")

        base_url = self.server_url.rstrip("/")
        search_url = f"{base_url}/api/search"

        params = {"q": self.query, "per_page": 20}
        if self.search_type != "all":
            params["type"] = self.search_type

        all_results = []
        start = 0

        try:
            while len(all_results) < self.max_results:
                params["start"] = start

                query_params = "&".join([f"{k}={v}" for k, v in params.items()])
                full_url = f"{search_url}?{query_params}"

                response = requests.get(full_url, timeout=self.timeout)

                if response.status_code != 200:
                    break

                data = response.json()
                items = data.get("data", {}).get("items", [])

                if not items:
                    break

                all_results.extend(items)

                start += params["per_page"]

                if len(all_results) >= self.max_results:
                    break

            all_results = all_results[: self.max_results]

            results_df = pd.DataFrame(all_results)

            if results_df.empty:
                results_df = pd.DataFrame({"no_results_found": []})

            return knext.Table.from_pandas(results_df)

        except requests.exceptions.RequestException as e:
            raise ValueError(f"Search error: {str(e)}")
        except Exception as e:
            raise ValueError(f"Processing error: {str(e)}")
