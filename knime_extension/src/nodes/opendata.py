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

        base_url = "https://www2.census.gov/geo/tiger/TIGER2020PL/STATE/"

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
county level data of all states. 

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
        "US Census APIkey", "The Census API key to use.", "Input Your Census APIkey "
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
        The 5-year estimates are available for all geographies down to the block group level.""",
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
        "US Census APIkey",
        "The APIkey to use.",
        "Input Your Census APIkey ",
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
    description="""This node downloads geospatial entities’ geometries and attributes from [OpenStreetMap.](https://www.openstreetmap.org/about)
Results returned are the union, not intersection of each individual tag. Each result matches at least one given tag. 
The place tags should be OSM tags, (e.g., building, landuse, highway, etc) and the value tags should be either True 
to retrieve all items with the given tag, or a single value to retrieve a single tag-value combination, or a 
comma separated list of values to get multiple values for the given tag. For example, 
place tag=building, tag value=True would return all building footprints in the area.
Place tag=landuse, tag value=retail, commercial would return all retail and commercial landuses.
""",
    references={
        "OpenStreet TagFinder": "https://tagfinder.osm.ch/",
        "OpenStreetMap Taginfo": "https://taginfo.openstreetmap.org/",
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
        G = ox.graph.graph_from_polygon(gdf_union, self.networktype)
        edges = ox.utils_graph.graph_to_gdfs(G, nodes=False)
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

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        gdf = get_osmnx().geocoder.geocode_to_gdf(self.placename)
        gdf = gdf.reset_index(drop=True)
        return knext.Table.from_pandas(gdf)
