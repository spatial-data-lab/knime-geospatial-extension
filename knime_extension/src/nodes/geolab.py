from typing import Callable
from wsgiref.util import shift_path_info
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import requests
import osmnx as ox
import pulp
from shapely.geometry import  LineString
from keplergl import KeplerGl

__category = knext.category(
    path="/geo",
    level_id="geolab",
    name="Geospatial Lab",
    description="Nodes that for testing and future exploration.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/GeolabCategroy.png",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Geolab/"

############################################
# US2020 TIGER/Line for states
############################################
@knext.node(
    name="US2020 TIGER",
    node_type=knext.NodeType.MANIPULATOR,
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
        "017"
    )
    
    geofile =knext.StringParameter(
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
        enum=["bg10","bg20","roads","tabblock10","tabblock20","tract10","tract20","county10","county20","state10","state20"],
    )        
 

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        USdict = pd.DataFrame.from_dict({
        'state': ['01', '02', '04', '05', '06', '08', '09', '10', '11', '12', '13', '15', '16', '17', '18', '19', '20', 
        '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39',
         '40', '41', '42', '44', '45', '46', '47', '48', '49', '50', '51', '53', '54', '55', '56', '72'],
        'stateName': ['ALABAMA', 'ALASKA', 'ARIZONA', 'ARKANSAS', 'CALIFORNIA', 'COLORADO', 'CONNECTICUT', 'DELAWARE',
         'COLUMBIA', 'FLORIDA', 'GEORGIA', 'HAWAII', 'IDAHO', 'ILLINOIS', 'INDIANA', 'IOWA', 'KANSAS', 'KENTUCKY', 'LOUISIANA', 
         'MAINE', 'MARYLAND', 'MASSACHUSETTS', 'MICHIGAN', 'MINNESOTA', 'MISSISSIPPI', 'MISSOURI', 'MONTANA', 'NEBRASKA', 'NEVADA',
         'HAMPSHIRE', 'JERSEY', 'MEXICO', 'YORK', 'CAROLINA', 'DAKOTA', 'OHIO', 'OKLAHOMA', 'OREGON', 'PENNSYLVANIA', 'ISLAND',
         'CAROLINA', 'DAKOTA', 'TENNESSEE', 'TEXAS', 'UTAH', 'VERMONT', 'VIRGINIA', 'WASHINGTON', 'VIRGINIA', 'WISCONSIN', 'WYOMING', 'RICO']})
 
        USdict['path']=USdict.state+"_"+USdict.stateName

        Statepath=USdict[USdict.state==self.StateFips].iloc[0,2] 

        County5Fips=self.StateFips+self.County3Fips

        base_url ='https://www2.census.gov/geo/tiger/TIGER2020PL/STATE/'
        
        if self.StateFips!=self.County3Fips: 
            data_url = f'{base_url}{Statepath}/{County5Fips}/tl_2020_{County5Fips}_{self.geofile}.zip'
        else:
            County5Fips=self.StateFips
            if self.geofile == "roads": self.geofile = "prisecroads"
            data_url = f'{base_url}{Statepath}/{County5Fips}/tl_2020_{County5Fips}_{self.geofile}.zip'
        gdf=gp.read_file(data_url)
        gdf1=gdf.reset_index(drop=True)
        #gdf1=gdf1[["GEOID20","geometry"]]
        return knext.Table.from_pandas(gdf1)


############################################
# US2020 Census Data
############################################
@knext.node(
    name="US2020 Census",
    node_type=knext.NodeType.MANIPULATOR,
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
    
    StateFips = knext.StringParameter("State FIPS (2-digits)", "The State to investigate, input * for all states (while choose county for geography)", "25")

    County3Fips = knext.StringParameter("County FIPS (3-digits)", "The County to investigate, input * for all counties", "017")

    Tract6Fips = knext.StringParameter("Tract FIPS (6-digits)", "The County to investigate, input * for all tracts", "*")

    censusapikey = knext.StringParameter("US Census APIkey", "The APIkey  to use", "Input Your Census APIkey ")

    cols=knext.StringParameter("US Census Variable Names", "The US Census Variable list to use", "GEO_ID,P1_001N,P1_003N,P1_004N,H1_001N,H1_002N")
    
    geofile =knext.StringParameter(
        label="Geographic Level",
        description=" Available geographic Level  for this node are Block,Block group,Tract and Couty",
        default_value="block group",
        enum=["block group","block","tract","county"],
    )  
 
    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):
        
        base_url = "https://api.census.gov/data/2020/dec/pl?get="

        if self.geofile =="county":
            data_url= f'{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&key={self.censusapikey}'
        elif self.geofile =="tract":
            data_url= f'{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&in=county:{self.County3Fips}&key={self.censusapikey}'
        else:
            data_url= f'{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&in=county:{self.County3Fips}&in=tract:{self.Tract6Fips}&key={self.censusapikey}'    

        response=requests.get(data_url)
        data=response.json()
        gdf=pd.DataFrame(data[1:], columns=data[0])
        return knext.Table.from_pandas(gdf)





############################################
# OSM Data Tool
############################################
@knext.node(
    name="OSM POIs",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "OSMdata.png",
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
    +"Results returned are the union, not intersection of each individual tag. Each result matches at least one given tag. "
    +"The dict keys should be OSM tags, (e.g., building, landuse, highway, etc) and the dict values should be either True "
    +"to retrieve all items with the given tag, or a string to get a single tag-value combination, or a list of strings to get"
    +" multiple values for the given tag. For example, tags = {‘building’: True} would return all building footprints in the area."
    +" tags = {‘amenity’:True, ‘landuse’:[‘retail’,’commercial’], ‘highway’:’bus_stop’} would return all amenities, landuse=retail, "
    +"landuse=commercial, and highway=bus_stop.",
    references={        
        "osmnx.geometries_from_place": "https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.geometries",  
        "Open Street Map Taginfo": "https://taginfo.openstreetmap.org/ ",    
    },
)
class OSMdataNode:
    
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_polygon,
        include_row_key=False,
        include_none_column=False,
    )

    taginfo = knext.StringParameter("Input places tags", "The value for tag to specify the type of facilities", "amenity")

    tagvalue = knext.StringParameter("Input tags value", "The specific type of facilities, True for all", "restaurant")

    def configure(self, configure_context, input_schema_1):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        gdf_union = gdf.to_crs(4326).unary_union
        if self.tagvalue != "True" :
            tags = {self.taginfo: self.tagvalue}        
        else: 
            tags = {self.taginfo: True} 
        #tags = {self.taginfo: self.tagvalue}  
        gdfpoi=ox.geometries.geometries_from_polygon(gdf_union,tags)
        gdfpoi=gdfpoi.reset_index(drop=True)
        return knext.Table.from_pandas(gdfpoi )

############################################
# Parquet Reader
############################################
@knext.node(
    name="Parquet Reader",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "parquet.png",
    category=__category,
)

@knext.output_table(
    name="Dataframe table ",
    description="Dataframefrom the input Parquet file",
)

@knut.pd_node_description(
    short_description="Read Parquet data",
    description="This node Load a parquet object from the file path, returning a DataFrame.",
    references={
        "pandas.read_parquet": "https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html",
    },
)
class ParquetReaderNode:
    
    data_url = knext.StringParameter("Input File Path", "The file path of Parquet file", "Input Parquet file")

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):        
        gdf=pd.read_parquet(self.data_url)
        gdf=gdf.reset_index(drop=True)
        return knext.Table.from_pandas(gdf)

############################################
# GeoFile Reader
############################################
@knext.node(
    name="GeoFile Reader",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoFileReader.png",
    category=__category,
)

@knext.output_table(
    name="Geodata table",
    description="Geodata from the input file path",
)
@knut.geo_node_description(
    short_description="Read single layer GeoFile.",
    description="This node read the Shapefile, zipped  Shapefile or Geojson with geopandas.read_file().",
    references={
        "Reading Spatial Data": "https://geopandas.org/en/stable/docs/user_guide/io.html",
    },
)

class GeoFileReaderNode:    
    data_url = knext.StringParameter(
        "Input File Path", 
        "The file path for reading data", 
        "https://raw.githubusercontent.com/UrbanGISer/Test/main/JsonMap/countries.geojson")

    def configure(self, configure_context):
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext):        
        gdf=gp.read_file(self.data_url)
        gdf=gdf.reset_index(drop=True)
        return knext.Table.from_pandas(gdf)


############################################
# Location-allocation Pmedian
############################################
@knext.node(
    name="P-media",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "pmedian.png",
    category=__category,
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
    short_description="Solve P-median problem for location-allocation.",
    description="The p-median model, one of the most widely used location models of any kind,"
    +"locates p facilities and allocates demand nodes to them to minimize total weighted distance traveled. "
    +"The P-Median problem will be solved by PuLP package. ",
    references={        
        "Pulp.Solver": "https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html",  
    },
)
class PmediaNode: 

    DemandID  = knext.ColumnParameter(
        "Serial id column for demand",
        "Integer id number starting with 0",
        #port_index=0,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )   

    SupplyID  = knext.ColumnParameter(
        "Serial id column for supply",
        "Integer id number starting with 0",
        #port_index=1,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )  

    DemandPopu  = knext.ColumnParameter(
        "Column for demand population",
        "The populaiton of demand",
        #port_index=2,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    ) 

    DemandGeometry  = knext.ColumnParameter(
        "Geometry column for demand points",
        "The Geometry column for demand points",
        #port_index=3,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    ) 

    SupplyGeometry= knext.ColumnParameter(
        "Geometry column for supply points",
        "The Geometry column for supply points",
        #port_index=4,
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )  
    ODcost= knext.ColumnParameter(
        "Travel cost between supply and demand",
        "The travel cost between the points of supply and demand",
        #port_index=5,
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    ) 
    CRSinfo = knext.StringParameter("Input CRSinfo for geometries of demand and supply", "The CRS information for geometry columns", "epsg:3857")
    Pnumber = knext.IntParameter("Input optimum p number of p-median", "The optimum p number of p-median", 5)

    def configure(self, configure_context, input_schema_1):
        #knut.geo_column_exists(self.geo_col, input_schema_1)
        # TODO Create combined schema
        return None
            #knext.Schema.from_columns(
            # [
            #     knext.Column(knext.(), self.DemandID),
            #     knext.Column(knext.double(), self.DemandPopu),
            #     knext.Column(knext.double(), "geometry"),
            #     knext.Column(knext.double(), "assignSID"),
            #     knext.Column(knext.double(), "SIDwkt"),
            #     knext.Column(knext.double(), "Linewkt"),
            # ])


    def execute(self, exec_context: knext.ExecutionContext,input_1):        
        df = gp.GeoDataFrame(input_1.to_pandas())
        # Sort with DID and SID
        df=df.sort_values(by=[self.DemandID, self.SupplyID]).reset_index(drop=True)

        # rebuild demand and supply data
        DemandPt=df[[self.DemandID,self.DemandPopu,self.DemandGeometry]].groupby([self.DemandID]).first().reset_index().rename(columns={"index": self.DemandID,self.DemandGeometry:'geometry'})
        SupplyPt=df[[self.SupplyID,self.SupplyGeometry]].groupby([self.SupplyID]).first().reset_index().rename(columns={"index": self.SupplyID,self.SupplyGeometry:'geometry'})
        DemandPt=gp.GeoDataFrame(DemandPt,geometry='geometry')
        SupplyPt=gp.GeoDataFrame(SupplyPt,geometry='geometry')
        DemandPt = DemandPt.set_crs(self.CRSinfo)
        SupplyPt=SupplyPt.set_crs(self.CRSinfo)

        # calculate parameter for matrix
        num_trt=DemandPt.shape[0] 
        num_hosp=SupplyPt.shape[0] 

        # define problem
        problem=pulp.LpProblem("MyProblem",sense=pulp.LpMinimize)

        def getIndex(s):
            li=s.split()
            if len(li)==2 :
                return int(li[-1])
            else :
                return int(li[-2]),int(li[-1])

        X=['X'+' '+str(i)+' '+str(j) for i in range(num_trt) for j in range(num_hosp)]
        Y=['Y'+' '+str(j) for j in range(num_hosp) ]

        varX=pulp.LpVariable.dicts("x",X,cat="Binary")
        varY=pulp.LpVariable.dicts("y",Y,cat="Binary")

        object_list=[]
        for it in X:
            i,j=getIndex(it)
            cost = df[(df[self.DemandID]== i) & (df[self.SupplyID]== j)][self.ODcost].values[0]
            object_list.append(DemandPt[self.DemandPopu][i]*cost*varX[it])
        problem+=pulp.lpSum(object_list)

        problem+=pulp.lpSum([varY[it] for it in Y])==self.Pnumber
        for i in range(num_trt):
            problem+=pulp.lpSum([varX['X'+' '+str(i)+' '+str(j)] for j in range(num_hosp) ])==1
        for i in range(num_trt):
            for j in range(num_hosp):
                problem+=(varX['X'+' '+str(i)+' '+str(j)]-varY['Y'+' '+str(j)]<=0)

        problem.solve()
        print(pulp.LpStatus[problem.status])

        # Extract result             

        DemandPt['assignSID']=None
        DemandPt['SIDcoord']=None
        for it in X:
            v = varX[it]
            if v.varValue == 1:
                i,j=getIndex(it)
                DemandPt.at[i,'assignSID']=SupplyPt.at[j,self.SupplyID]
                DemandPt.at[i,'SIDcoord']=SupplyPt.at[j,'geometry']
        DemandPt['linxy']=[LineString(xy) for xy in zip(DemandPt['geometry'],DemandPt['SIDcoord'])]
        DemandPt['SIDwkt']=DemandPt.set_geometry('SIDcoord').geometry.to_wkt()
        DemandPt['Linewkt']=DemandPt.set_geometry('linxy').geometry.to_wkt()
        DemandPt=DemandPt.drop(columns=['linxy','SIDcoord'])
        DemandPt=DemandPt.reset_index(drop=True)
        #DemandPt=DemandPt.rename(columns={'geometry':'DIDgeometry'})
        gdf = gp.GeoDataFrame(DemandPt, geometry="geometry", crs=self.CRSinfo)
        return knext.Table.from_pandas(gdf)




############################################
# Kepler.gl 
############################################
@knext.node(
    name="Kepler.gl Geoview ",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "Kepler.gl.png",
    category=__category,
)
@knext.input_table(
    name="Geospatial table to visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial view", description="Showing a map with the geospatial data"
)
class ViewNodeKepler:
    """
    This node will visualize the given geometric elements on a map.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,  
        include_row_key=False,
        include_none_column=False,
    )

    save_config = knext.BoolParameter(
        "Save config",
        "Save the config for the map",
        default_value=False,
    )

    load_config = knext.BoolParameter(
        "Load config",
        "Load the config for the map",
        default_value=False,
    )


    def configure(self, configure_context, input_schema):
        knut.columns_exist([ self.geo_col], input_schema)
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)
        
        map_1 = KeplerGl(show_docs=False)
        map_1.add_data(data=gdf.copy(), name="state")
        config = {}
        if self.save_config:
            # Save map_1 config to a file
            # config_str = json.dumps(map_1.config)
            # if type(config) == str:
            #     config = config.encode("utf-8")
            with open('kepler_config.json', 'w') as f:
                f.write(json.dumps(map_1.config))
        
        
        if self.load_config:
            with open('kepler_config.json', 'r') as f:
                config = json.loads(f.read())
        map_1.config = config

        # map_1.add_data(data=data.copy(), name="haha")
        html = map_1._repr_html_()
        html = html.decode("utf-8")
        # knext.view_html(html)
        # knut.check_canceled(exec_context)
        # cx.add_basemap(map, crs=gdf.crs.to_string(), source=cx.providers.flatten()[self.base_map])

        return knext.view_html(html)



# @knext.node(
#     name="Echarts Geoview ",
#     node_type=knext.NodeType.VISUALIZER,
#     icon_path="icons/icon/Visulization/StaticMap.png",
#     category=category,
# )
# @knext.input_table(
#     name="Geospatial table to visualize",
#     description="Table with geospatial data to visualize",
# )
# @knext.output_view(
#     name="Geospatial view", description="Showing a map with the geospatial data"
# )
# class ViewNodeEcharts:
#     """
#     This node will visualize the given geometric elements on a map.
#     """

#     geo_col = knext.ColumnParameter(
#         "Geometry column",
#         "Select the geometry column to visualize.",
#         column_filter=knut.is_geo,  
#         include_row_key=False,
#         include_none_column=False,
#     )

#     save_config = knext.BoolParameter(
#         "Save config",
#         "Save the config for the map",
#         default_value=False,
#     )

#     load_config = knext.BoolParameter(
#         "Load config",
#         "Load the config for the map",
#         default_value=False,
#     )


#     def configure(self, configure_context, input_schema):
#         knut.columns_exist([ self.geo_col], input_schema)
#         # if self.name_cols is None:
#         #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
#         return None

#     def execute(self, exec_context: knext.ExecutionContext, input_table):
#         gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)
        
#         bar = Bar()
#         bar.add_xaxis(["衬衫", "羊毛衫", "雪纺衫", "裤子", "高跟鞋", "袜子"])
#         bar.add_yaxis("商家A", [5, 20, 36, 10, 75, 90])
#         # render 会生成本地 HTML 文件，默认会在当前目录生成 render.html 文件
#         # 也可以传入路径参数，如 bar.render("mycharts.html")
#         # html = bar.render_notebook()._repr_html_()
#         html = bar.render_embed()

#         return knext.view_html(html)

