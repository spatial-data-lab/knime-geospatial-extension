from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import requests

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
    name="US2020TIGER",
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
<li><b>Block group":</b> bg10, bg20, for US Census block group in 2010 and 2020.</li>
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

        return knext.Table.from_pandas(gdf)


############################################
# US2020 Census Data
############################################
@knext.node(
    name="US2020Census",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "UScensus.png",
    category=__category,
)

@knext.input_table(
    name="Census data table",
    description="Retrieved geodata from census.gov/data/2020",
)

@knext.output_table(
    name="Table with coordinates",
    description="Table with the coordinates for all point geometries",
)

@knut.geo_node_description(
    short_description="Extracts the XYZ coordinates.",
    description="This node extracts the XYZ coordinates for all given point objects."
    + "The coordinates are appended to the input table as xyz coordinate columns.",
    references={
        "X coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.x.html",
        "Y coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.y.html",
        "Z coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.z.html",
    },
)
class USCensus2020Node:
    
    input_column = knext.ColumnParameter(
        label="Any Column",
        description="[Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column to convert",
        include_row_key=False,
        include_none_column=False,
    )
    
    StateFips = knext.StringParameter("State FIPS (2-digits)", "The State to use", "48")

    County3Fips = knext.StringParameter("US County 3-digits FIPS", "The County to use", "015")

    Tract6Fips = knext.StringParameter("US Tract 6-digits FIPS or *", "The County  to use", "*")

    censusapikey = knext.StringParameter("US census apikey", "The Apikey  to use", "761438c28421cc4b9d179d23f7531489a4cf6718")

    cols=knext.StringParameter("US Census Variable Name *", "The US Census Variable Name to use", "P1_001N,H1_001N")
    
    geofile =knext.StringParameter(
        "Selection",
        "The geographic data to choose from.",
        "block",
        enum=["block","block group"],
    )  
 
    def configure(self, configure_context, input_schema):
        knut.column_exists(self.input_column, input_schema)
        # TODO Create combined schema
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        
        base_url = "https://api.census.gov/data/2020/dec/pl?get="
        data_url= f'{base_url}{self.cols}&for={self.geofile}:*&in=state:{self.StateFips}&in=county:{self.County3Fips}&in=tract:{self.Tract6Fips}&key={self.censusapikey}'
        response=requests.get(data_url)
        data=response.json()
        gdf=pd.DataFrame(data[1:], columns=data[0])
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