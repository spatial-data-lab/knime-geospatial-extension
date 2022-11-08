from typing import Callable
from wsgiref.util import shift_path_info
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import requests
import osmnx as ox

__category = knext.category(
    path="/geo",
    level_id="opendataset",
    name="Open Datasets",
    description="Nodes for providing open geospatial dataset.",
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
)
@knext.output_table(
    name="TIGER Line Table ",
    description="Retrieved geodata from www2.census.gov/geo/tiger/TIGER2020PL",
)
@knut.census_node_description(
    short_description="Retrieve geospatial data from US Census TIGER/Line",
    description="This node Retrieve the specific geospatial boundaries for one specific state of United States."
    + "The popular TIGER/Line are Block group, Roads, Blocks,Tracts."
    + "While input the same State FIPS (2-digits) for County FIPS(3-digits), the geodata of all counties in the state will be retrieved,"
    + "county10/20 and state10/20 can only be applicable in this case. ",
    references={
        "TIGER/Line Shapefiles": "https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.2020.html",
        "FTP Archive by State": "https://www2.census.gov/geo/tiger/TIGER2020PL/STATE/",
        "FTP Archive by Layer": "https://www2.census.gov/geo/tiger/TIGER2020PL/LAYER/",
    },
)
class US2020TIGERNode:

    StateFips = knext.StringParameter(
        label="State FIPS (2-digits)",
        description="The State to use",
        default_value="25",
    )

    County3Fips = knext.StringParameter(
        "County FIPS(3-digits)/ or State FIPS ",
        "The County/State FIPS code to use",
        "017",
    )

    geofile = knext.StringParameter(
        label="TIGER/Line data type",
        description="""<p>Available TIGER/Line are: 
<ul>
<li><b>Block:</b> tabblock10,tabblock20, for US Census block (the minimun level) in 2010 and 2020.</li>
<li><b>Block group:</b> bg10, bg20, for US Census block group in 2010 and 2020.</li>
<li><b>Tract:</b> tract10,tract20, for US Census Tract in 2010 and 2020.</li>
<li><b>County:</b> county10,county20, for US County in 2010 and 2020, only applicable for the setting when the input of County FIPs equals State FIPs.</li>
<li><b>State:</b> state10,state20, for US State in 2010 and 2020,same as County.</li>
</ul></p>
""",
        default_value="bg20",
        enum=[
            "bg10",
            "bg20",
            "roads",
            "tabblock10",
            "tabblock20",
            "tract10",
            "tract20",
            "county10",
            "county20",
            "state10",
            "state20",
        ],
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
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
                    "COLUMBIA",
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
                    "HAMPSHIRE",
                    "JERSEY",
                    "MEXICO",
                    "YORK",
                    "CAROLINA",
                    "DAKOTA",
                    "OHIO",
                    "OKLAHOMA",
                    "OREGON",
                    "PENNSYLVANIA",
                    "ISLAND",
                    "CAROLINA",
                    "DAKOTA",
                    "TENNESSEE",
                    "TEXAS",
                    "UTAH",
                    "VERMONT",
                    "VIRGINIA",
                    "WASHINGTON",
                    "VIRGINIA",
                    "WISCONSIN",
                    "WYOMING",
                    "RICO",
                ],
            }
        )

        USdict["path"] = USdict.state + "_" + USdict.stateName

        Statepath = USdict[USdict.state == self.StateFips].iloc[0, 2]

        County5Fips = self.StateFips + self.County3Fips

        base_url = "https://www2.census.gov/geo/tiger/TIGER2020PL/STATE/"

        if self.StateFips != self.County3Fips:
            data_url = f"{base_url}{Statepath}/{County5Fips}/tl_2020_{County5Fips}_{self.geofile}.zip"
        else:
            County5Fips = self.StateFips
            if self.geofile == "roads":
                self.geofile = "prisecroads"
            data_url = f"{base_url}{Statepath}/{County5Fips}/tl_2020_{County5Fips}_{self.geofile}.zip"
        gdf = gp.read_file(data_url)
        gdf1 = gdf.reset_index(drop=True)
        # gdf1=gdf1[["GEOID20","geometry"]]
        return knext.Table.from_pandas(gdf1)


############################################
# US2020 Census Data
############################################
@knext.node(
    name="US2020 Census Data",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "UScensus.png",
    category=__category,
)
@knext.output_table(
    name="Census data table",
    description="Retrieved geodata from 2020 Census Redistricting Data",
)
@knut.census_node_description(
    short_description="Retrieve US 2020 Census Redistricting Data for one specific state of United States.",
    description="This node retrieve US 2020 Census Redistricting Data (Decennial Census P.L. 94-171 Redistricting Data)."
    + "This node privde all variables such as population and household inforamtion in US Census 2020 data base."
    + "The default variable names are GEO_ID (geography), ,P1_001N(Total Population),P1_003N(Population of one race:!!White alone),"
    + "P1_004N(Black or African American alone),H1_001N(Total Housing Units),H1_002 (Total Occupied Housing Units)."
    + "Only if county is choosed for geography, then * can be input in  State FIPS (2-digits) to retrieve all the county level data of all states."
    + "Before using the node, user need to sign up and get a census api key first by clicking the following hyperlink",
    references={
        "Census API Key Sign Up": "https://api.census.gov/data/key_signup.html",
        "Decennial Census P.L. 94-171 Redistricting Data": "https://www.census.gov/programs-surveys/decennial-census/about/rdo/summary-files.html",
        "Datasets and its descendants": "https://api.census.gov/data/2020/dec.html",
        "Geography": "https://api.census.gov/data/2020/dec/pl/geography.html",
        "Variables": "https://api.census.gov/data/2020/dec/pl/variables.html",
        "Census API examples": "https://api.census.gov/data/2020/dec/pl/examples.html",
    },
)
class USCensus2020Node:

    StateFips = knext.StringParameter(
        "State FIPS (2-digits)",
        "The State to investigate, input * for all states (while choose county for geography)",
        "25",
    )

    County3Fips = knext.StringParameter(
        "County FIPS (3-digits)",
        "The County to investigate, input * for all counties",
        "017",
    )

    Tract6Fips = knext.StringParameter(
        "Tract FIPS (6-digits)", "The Tract to investigate, input * for all tracts", "*"
    )

    censusapikey = knext.StringParameter(
        "US Census APIkey", "The APIkey  to use", "Input Your Census APIkey "
    )

    cols = knext.StringParameter(
        "US Census Variable Names",
        "The US Census Variable list to use",
        "GEO_ID,P1_001N,P1_003N,P1_004N,H1_001N,H1_002N",
    )

    geofile = knext.StringParameter(
        label="Geographic Level",
        description=" Available geographic Level  for this node are Block,Block group,Tract and Couty",
        default_value="block group",
        enum=["block group", "block", "tract", "county"],
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

        response = requests.get(data_url)
        data = response.json()
        gdf = pd.DataFrame(data[1:], columns=data[0])
        return knext.Table.from_pandas(gdf)


############################################
# US Census ACS-5
############################################
@knext.node(
    name="US ACS 5-Year Estimates",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "CensusACS.png",
    category=__category,
)
@knext.output_table(
    name="US Census ACS 5-Year table",
    description="Retrieved data from Census ACS Datasets",
)
@knut.census_node_description(
    short_description="Retrieve American Community Survey 5-Year Data (2009-2020) of United States.",
    description="This node retrieve American Community Survey 5-Year Data (2009-2020)."
    + "The American Community Survey (ACS) is an ongoing survey that provides data every year -- giving "
    + "communities the current information they need to plan investments and services. The ACS covers a broad "
    + "range of topics about social, economic, demographic, and housing characteristics of the U.S. population."
    + "The 5-year estimates from the ACS are period estimates that represent data collected over a period of "
    + "time. The primary advantage of using multiyear estimates is the increased statistical reliability of the data "
    + "for less populated areas and small population subgroups. "
    + " The 5-year estimates are available for all geographies down to the block group level.",
    references={
        "Census API Key Sign Up": "https://api.census.gov/data/key_signup.html",
        "American Community Survey 5-Year Data (2009-2020)": "https://www.census.gov/data/developers/data-sets/acs-5year.html",
        "Geography": "https://api.census.gov/data/2020/acs/acs5/geography.html",
        "Variables": "https://api.census.gov/data/2020/acs/acs5/variables.html",
        "Census API examples": "https://api.census.gov/data/2020/acs/acs5/subject/examples.html",
    },
)
class UScensusACSNode:

    StateFips = knext.StringParameter(
        "State FIPS (2-digits)",
        "The State to investigate, input * for all states (while choose county for geography)",
        "25",
    )

    County3Fips = knext.StringParameter(
        "County FIPS (3-digits)",
        "The County to investigate, input * for all counties",
        "017",
    )

    Tract6Fips = knext.StringParameter(
        "Tract FIPS (6-digits)", "The Tract to investigate, input * for all tracts", "*"
    )

    censusapikey = knext.StringParameter(
        "US Census APIkey", "The APIkey  to use", "Input Your Census APIkey "
    )

    cols = knext.StringParameter(
        "US Census ACS Variable Names",
        "The US Census Variable list to use",
        "GEO_ID,B02001_001E,B02001_002E,B02001_003E",
    )

    geofile = knext.StringParameter(
        label="Geographic Level",
        description=" Available geographic Level  for this node are Block,Block group,Tract and Couty",
        default_value="block group",
        enum=["block group", "tract", "county", "state"],
    )

    year = knext.StringParameter(
        "US Census ACS5 Year Label", "The Year label of dataset", "2020"
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

        response = requests.get(data_url)
        data = response.json()
        gdf = pd.DataFrame(data[1:], columns=data[0])
        return knext.Table.from_pandas(gdf)


############################################
# OSM Data Tool
############################################
@knext.node(
    name="OSM POIs",
    node_type=knext.NodeType.SOURCE,
    icon_path=__NODE_ICON_PATH + "OSMpoi.png",
    category=__category,
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
    description="This node Download geospatial entities’ geometries and attributes from OpenStreetMap."
    + "Results returned are the union, not intersection of each individual tag. Each result matches at least one given tag. "
    + "The dict keys should be OSM tags, (e.g., building, landuse, highway, etc) and the dict values should be either True "
    + "to retrieve all items with the given tag, or a string to get a single tag-value combination, or a list of strings to get"
    + " multiple values for the given tag. For example, tags = {‘building’: True} would return all building footprints in the area."
    + " tags = {‘amenity’:True, ‘landuse’:[‘retail’,’commercial’], ‘highway’:’bus_stop’} would return all amenities, landuse=retail, "
    + "landuse=commercial, and highway=bus_stop.",
    references={
        "osmnx.geometries_from_place": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.geometries",
        "Open Street Map Taginfo": "https://taginfo.openstreetmap.org/ ",
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
        "Input places tags",
        "The value for tag to specify the type of facilities",
        "amenity",
    )

    tagvalue = knext.StringParameter(
        "Input tags value",
        "The specific type of facilities, True for all",
        "restaurant",
    )

    def configure(self, configure_context, input_schema_1):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf_union = gdf.to_crs(4326).unary_union
        if self.tagvalue != "True":
            tags = {self.taginfo: self.tagvalue}
        else:
            tags = {self.taginfo: True}
        # tags = {self.taginfo: self.tagvalue}
        gdfpoi = ox.geometries.geometries_from_polygon(gdf_union, tags)
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
    short_description="Get Road Netwrok from the Open Street Map.",
    description="This node Download geospatial network and its attributes from OpenStreetMap."
    + "If the network type of drive is used for query, it will use the module add_edge_sppeds to append the speed information directly. "
    + "The total travel time for the segment will be calculated as well with the module add_edge_travel_times. ",
    references={
        "osmnx.graph.graph_from_polygon": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.graph",
        "osmnx.speed.add_edge_speeds": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.speed ",
    },
)
class OSMnetworkNode:

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column as boundary to get POIs.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_polygon,
        include_row_key=False,
        include_none_column=False,
    )

    networktype = knext.StringParameter(
        label="Street network type",
        description="Type of street network ",
        default_value="drive",
        enum=["all_private", "all", "bike", "drive", "drive_service", "walk"],
    )

    def configure(self, configure_context, input_schema_1):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf_union = gdf.to_crs(4326).unary_union
        G = ox.graph.graph_from_polygon(gdf_union, self.networktype)
        if self.networktype == "drive":
            # impute speed on all edges missing data
            G = ox.add_edge_speeds(G)
            # calculate travel time (seconds) for all edges
            G = ox.add_edge_travel_times(G)
        edges = ox.utils_graph.graph_to_gdfs(G, nodes=False)
        objcolumn = edges.select_dtypes(include=["object"]).columns.tolist()
        edges[objcolumn] = edges[objcolumn].astype("string")
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
)
@knext.output_table(
    name="OSM GeoBoundary data",
    description="BOundary of places from the Open Street Map",
)
@knut.osm_node_description(
    short_description="Get Boundary from the Open Street Map with Geocoding.",
    description="This node get place boudnary  from OpenStreetMap by the geocoding place name."
    + "If query argument is a list, then which_result should be either a single value or a list with the same length as query. . "
    + "The queries you provide must be resolvable to places in the Nominatim database. "
    + " The resulting GeoDataFrame’s geometry column contains place boundaries if they exist in OpenStreetMap.",
    references={
        "osmnx.geocoder.geocode_to_gdf": "https://osmnx.readthedocs.io/en/stable/osmnx.html?highlight=geocode_to_gdf#module-osmnx.geocoder",
    },
)
class OSMGeoBoundaryNode:

    placename = knext.StringParameter(
        label="Input Place Names",
        description="Hierachial place Names deliminated with commas, such as Cambridge, MA, USA ",
        default_value="Cambridge, MA, USA",
    )

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        gdf = ox.geocode_to_gdf(self.placename)
        gdf = gdf.reset_index(drop=True)
        return knext.Table.from_pandas(gdf)
