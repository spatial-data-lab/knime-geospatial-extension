from doctest import debug_script
import logging
import plotly.express as px
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
# import contextily as cx
import folium
from io import StringIO
import matplotlib.pyplot as plt
import hvplot.pandas
import hvplot
import plotly.io as io
import altair as alt
from vega_datasets import data
from keplergl import KeplerGl
from pyecharts.charts import Bar

import pydeck as pdk
import pandas as pd

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/geo",
    level_id="viz",
    name="Views",
    description="Spatial view nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/VisulizationCategory.png",
)


@knext.node(
    name="Geospatial View",
    node_type=knext.NodeType.VISUALIZER,
    icon_path="icons/icon/Visulization/InteractiveMap.png",
    category=category,
)
@knext.input_table(
    name="Geospatial table to visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial view", description="Showing a interactive map with the geospatial data"
)
class ViewNode:
    """
    This node cteates a interative map view based on the selected geometric elements of the input table.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        # "geometry",
        column_filter=knut.is_geo,  # Allows all geo columns
        include_row_key=False,
        include_none_column=False, # must contains a geometry column
    )

    color_col = knext.ColumnParameter(
        "Marker color column",
        "Select marker color column to be plotted.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    color_map = knext.StringParameter(
        "Color map",
        "Select the color map to use for the color column. `xxx_r` mean reverse of the `xxx` colormap. See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html)",
        default_value="viridis",
        enum=["viridis", "plasma", "inferno", "magma", "cividis",
                'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper',

                "viridis_r", "plasma_r", "inferno_r", "magma_r", "cividis_r",
                'Greys_r', 'Purples_r', 'Blues_r', 'Greens_r', 'Oranges_r', 'Reds_r',
                'YlOrBr_r', 'YlOrRd_r', 'OrRd_r', 'PuRd_r', 'RdPu_r', 'BuPu_r',
                'GnBu_r', 'PuBu_r', 'YlGnBu_r', 'PuBuGn_r', 'BuGn_r', 'YlGn_r',
                'binary_r', 'gist_yarg_r', 'gist_gray_r', 'gray_r', 'bone_r',
                'pink_r', 'spring_r', 'summer_r', 'autumn_r', 'winter_r', 'cool_r',
                'Wistia_r', 'hot_r', 'afmhot_r', 'gist_heat_r', 'copper_r'],
        
    )

    base_map = knext.StringParameter(
        "Base map",
        "Select the base map to use for the visualization. See [Folium base maps](https://python-visualization.github.io/folium/quickstart.html#Tiles).",
        default_value="OpenStreetMap",
        enum=["OpenStreetMap", "Stamen Terrain", "Stamen Toner", "Stamen Watercolor" "CartoDB positron", "CartoDB dark_matter"]
    )

    use_classify = knext.BoolParameter(
        "Use classification",
        "If checked, the color column will be classified using the selected classification method. The `Number of classes` will be used to determine the number of classes.",
        default_value=True,
    )

    classification_method = knext.StringParameter(
        "Classification method",
        "Select the classification method to use for the color column.",
        default_value="EqualInterval",
        enum=['BoxPlot', 'EqualInterval', 'FisherJenks', 'FisherJenksSampled', 'HeadTailBreaks', 'JenksCaspall', 'JenksCaspallForced', 'JenksCaspallSampled', 'MaxP', 'MaximumBreaks', 'NaturalBreaks', 'Quantiles', 'Percentiles', 'StdMean']
        
    )

    classification_bins = knext.IntParameter(
        "Number of classes",
        "Select the number of classes of the classification method.",
        default_value=5,
        min_value=1,
        max_value=50,
    )

    size_col = knext.ColumnParameter(
        "Marker size column",
        "Select marker size column. The size is fixed by default. If a size column is selected, the size will be scaled by the values of the column. For point features, the size is the radius of the circle. For line features, the size is the width of the line. For polygon features, the size is the radius of the centroid of the ploygon.",
        column_filter=knut.is_numeric,
        include_none_column=True,
    )


    name_cols = knext.MultiColumnParameter(
        "Tooltip columns",
        "Select columns which should be shown in the marker tooltip.",
        column_filter=knut.is_string,
    )

    popup_cols = knext.MultiColumnParameter(
        "Popup columns",
        "Select columns which should be shown in the marker popup.",
        column_filter=knut.is_string,
    )

    plot_legend = knext.BoolParameter(
        "Show legend",
        "If checked, a legend will be shown in the plot.",
        default_value=True,

    )

    legend_caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the legend. By default, the caption is the name of the selected color column.",
        default_value="",
    )


    def configure(self, configure_context, input_schema):
        knut.columns_exist([ self.geo_col], input_schema)
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        if (self.legend_caption is None) or (self.legend_caption == ""):
            self.legend_caption = self.color_col
        if "none" not in str(self.size_col).lower():

            
            max_pop_est = gdf[self.size_col].max()
            min_pop_est = gdf[self.size_col].min()

            # check whether is line 
            geo_types = gdf["geometry"].geom_type.unique()
            if  ("LineString" in geo_types) or ("MultiLineString" in geo_types):
                size = "weight"
                max_size = 8
                m = None
            # else:
            #     size = "radius"
            #     max_size = 30
            #     # m = 1
            elif ("Polygon" in geo_types) or ("MultiPolygon" in geo_types):
                size = "radius"
                m = gdf.explore()
                gdf["geometry"] = gdf.centroid
                max_size = 30
            else:
                size = "radius"
                max_size = 30
                m = None

            if self.use_classify:
                map = gdf.explore(
                    column=self.color_col, 
                    cmap=self.color_map,
                    tooltip=self.name_cols,
                    tiles=self.base_map,
                    popup=self.popup_cols,
                    scheme=self.classification_method,
                    k=self.classification_bins,
                    legend=self.plot_legend,
                    m=m,
                    # marker_kwds={"radius":self.size_col},
                    # mape POP_EST to 1-60
                    style_kwds={ 
                        "style_function": lambda x: {
                            size: (x["properties"][self.size_col] - min_pop_est) / (max_pop_est - min_pop_est) * max_size
                            }
                        },

                    legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":20,"colorbar":False}
                )
            else:
                map = gdf.explore(
                    column=self.color_col, 
                    cmap=self.color_map,
                    tooltip=self.name_cols,
                    tiles=self.base_map,
                    popup=self.popup_cols,
                    legend=self.plot_legend,
                    m=m,
                    style_kwds={ 
                        "style_function": lambda x: {
                            size: (x["properties"][self.size_col] - min_pop_est) / (max_pop_est - min_pop_est) * 30
                            }
                        },
                    legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":3,"colorbar":True}
                )
        else:

            if self.use_classify:
                map = gdf.explore(
                    column=self.color_col, 
                    cmap=self.color_map,
                    tooltip=self.name_cols,
                    tiles=self.base_map,
                    popup=self.popup_cols,
                    scheme=self.classification_method,
                k=self.classification_bins,
                legend=self.plot_legend,
                legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":20,"colorbar":False}
                )
            else:
                map = gdf.explore(
                    column=self.color_col, 
                    cmap=self.color_map,
                    tooltip=self.name_cols,
                    tiles=self.base_map,
                    popup=self.popup_cols,
                    legend=self.plot_legend,
                    # marker_kwds={"radius":self.size_col},
                    legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":3,"colorbar":True}
                )                                                           
    
        # knut.check_canceled(exec_context)
        return knext.view(map)

# geo view static
# TODO:
# - add legend location
# - add legend caption
@knext.node(
    name="Geospatial View Static",
    node_type=knext.NodeType.VISUALIZER,
    icon_path="icons/icon/Visulization/StaticMap.png",
    category=category,
)
@knext.input_table(
    name="Geospatial table to visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial view", description="Showing a map with the geospatial data"
)
class ViewNodeStatic:
    """
    This node will visualize the given geometric elements on a static map.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,  
        include_row_key=False,
        include_none_column=False,
    )

    color_col = knext.ColumnParameter(
        "Marker color column",
        "Select marker color column. The column must contain the color name e.g. red, green, blue, etc.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    color_map = knext.StringParameter(
        "Color map",
        "Select the color map to use for the color column. See https://matplotlib.org/stable/tutorials/colors/colormaps.html",
        default_value="viridis",
        enum=["viridis", "plasma", "inferno", "magma", "cividis"],
        
    )

    # base_map = knext.StringParameter(
    #     "Base map",
    #     "Select the base map to use for the visualization. See https://contextily.readthedocs.io/en/latest/providers_deepdive.html",
    #     default_value="OpenStreetMap",
    #     enum=['OpenStreetMap.Mapnik',
    #          'OpenTopoMap',
    #          'Stamen.Toner',
    #          'Stamen.TonerLite',
    #          'Stamen.Terrain',
    #          'Stamen.TerrainBackground',
    #          'Stamen.Watercolor',
    #          'NASAGIBS.ViirsEarthAtNight2012',
    #          'CartoDB.Positron',
    #          'CartoDB.Voyager'
    #         ]
    # )

    use_classify = knext.BoolParameter(
        "Use classification",
        "If checked, the color column will be classified using the selected classification method.",
        default_value=True,
    )

    classification_method = knext.StringParameter(
        "Classification method",
        "Select the classification method to use for the color column.",
        default_value="EqualInterval",
        enum=['BoxPlot', 'EqualInterval', 'FisherJenks', 'FisherJenksSampled', 'HeadTailBreaks', 'JenksCaspall', 'JenksCaspallForced', 'JenksCaspallSampled', 'MaxP', 'MaximumBreaks', 'NaturalBreaks', 'Quantiles', 'Percentiles', 'StdMean']
    )

    classification_bins = knext.IntParameter(
        "Number of classes",
        "Select the number of classes to use for the color column.",
        default_value=5,
        min_value=1,
        max_value=10,
    )

    figure_title = knext.StringParameter(
        "Figure title",
        "Set the title of the figure.",
        default_value="",
    )

    figure_title_size = knext.IntParameter(
        "Figure title size",
        "Set the size of the figure title.",
        default_value=10,
        min_value=1,
        max_value=100,
    )

    plot_legend = knext.BoolParameter(
        "Show legend",
        "If checked, a legend will be shown in the plot.",
        default_value=True,

    )

    # size_col = knext.ColumnParameter(
    #     "Marker size column",
    #     "Select marker size column. The column must contain the size value.",
    #     column_filter=knut.is_numeric,
    # )

    legend_caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the legend.",
        default_value="",
        # default_value=color_col,
    )
    
    legend_caption_fontsize = knext.IntParameter(
        "Legend caption font size",
        "Set the font size for the legend caption.",
        default_value=10,
        min_value=1,
        max_value=100,
    )

    legend_expand = knext.BoolParameter(
        "Expand legend",
        "If checked, the legend will be horizontally expanded to fill the axes area",
        default_value=False,
    )

    legend_location = knext.StringParameter(
        "Legend location",
        "Select the location for the legend.",
        default_value="best",
        enum=['best', 
            'upper right', 
            'upper left', 
            'lower left', 
            'lower right', 
            'right', 
            'center left', 
            'center right', 
            'lower center', 
            'upper center', 
            'center',
            'outside_top',
            'outside_bottom',
            ]
    )

    legend_columns = knext.IntParameter(
        "Legend columns",
        "Select the number of columns for the legend.",
        default_value=1,
        min_value=1,
        max_value=30,
    )

    legend_size = knext.IntParameter(
        "Legend size",
        "Select the size for the legend.",
        default_value=8,
        min_value=1,
        max_value=30,
    )

    legend_fontsize = knext.IntParameter(
        "Legend font size",
        "Select the font size for the legend.",
        default_value=10,
        min_value=1,
        max_value=30,
    )

    legend_labelcolor = knext.StringParameter(
        "Legend label color",
        "Select the label color for the legend.",
        default_value="black",
        enum=['black', 'red', 'green', 'blue', 'yellow', 'purple', 'orange', 'white']
    )

    legend_frame = knext.BoolParameter(
        "Show legend frame",
        "If checked, a frame will be shown in the legend.",
        default_value=True,

    )

    legend_framealpha = knext.DoubleParameter(
        "Legend frame alpha",
        "Select the alpha value for the legend frame.",
        default_value=1.0,
        min_value=0.0,
        max_value=1.0,
    )

    legend_borderpad = knext.DoubleParameter(
        "Legend border pad",
        "Select the border pad for the legend.",
        default_value=0.5,
        min_value=0.0,
        max_value=3.0,
    )

    legend_labelspacing = knext.DoubleParameter(
        "Legend label spacing",
        "Select the label spacing for the legend.",
        default_value=0.5,
        min_value=0.0,
        max_value=1.0,
    )

    legend_orientation = knext.StringParameter(
        "Legend orientation",
        "Select the orientation for the legend.",
        default_value="vertical",
        enum=['vertical', 'horizontal']
    )



    def configure(self, configure_context, input_schema):
        knut.columns_exist([ self.geo_col], input_schema)
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        if (self.legend_caption is None) or (self.legend_caption == ""):
            self.legend_caption = self.color_col
        legend_bbox_to_anchor = None
        if self.legend_location == "outside_top":
            self.legend_location = "lower right"
            legend_bbox_to_anchor = (0.0, 1.02, 1.0, 0.102)
        if self.legend_location == "outside_bottom":
            self.legend_location = "upper right"
            legend_bbox_to_anchor =  (0.0, -0.2, 1.0, 0.102)
        if self.legend_expand:
            legend_expand = "expand"
        else:
            legend_expand = None

        # if self.size_col is not None:
        #     gdf["geometry"] = gdf.centroid

        if self.use_classify:
            map = gdf.plot(
                column=self.color_col, 
                cmap=self.color_map,
                # tiles=self.base_map,
                alpha=1,
                scheme=self.classification_method,
                k=self.classification_bins,
                legend=self.plot_legend,
                legend_kwds={
                    # 'shrink': 0.56,
                    'fmt':"{:.0f}",
                    'loc': self.legend_location,
                    # "bbox_to_anchor":(0.47, -0.1),
                    "title": self.legend_caption,
                    'ncols': self.legend_columns,
                    'prop': {'size': self.legend_size},
                    'fontsize': self.legend_fontsize,
                    'bbox_to_anchor': legend_bbox_to_anchor, 
                    'labelcolor': self.legend_labelcolor,
                    'frameon': self.legend_frame,
                    'framealpha': self.legend_framealpha,
                    'fancybox': True,
                    'mode': legend_expand,
                    'alignment': "left",
                    'title': "Population",
                    'title_fontsize': self.legend_caption_fontsize,
                    'labelspacing': self.legend_labelspacing,
                    # 'handletextpad': 0.5,
                    # 'handlelength': 1.5,
                    # 'handleheight': 1,
                    'borderaxespad':self.legend_borderpad,
                    # 'columnspacing' : 4,
                    },
                )
            map.set_title(self.figure_title, fontsize=self.figure_title_size)
        else:
            map = gdf.plot(
                column=self.color_col, 
                cmap=self.color_map,
                # tiles=self.base_map,
                alpha=0.7,
                legend=self.plot_legend,
                # marker_kwds={"radius":self.size_col},
                legend_kwds={
                    'shrink': 0.56,
                    'orientation': self.legend_orientation,
                    }
                # legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":3,"colorbar":True}
            )
    
        # knut.check_canceled(exec_context)
        # return knext.
        # cx.add_basemap(map, crs=gdf.crs.to_string(), source=cx.providers.flatten()[self.base_map])

        return knext.view_matplotlib(map.get_figure())


# geo view html
# geo view static
# TODO:
# - add legend location
# - add legend caption
@knext.node(
    name="Geospatial View html",
    node_type=knext.NodeType.VISUALIZER,
    icon_path="icons/icon/Visulization/StaticMap.png",
    category=category,
)
@knext.input_table(
    name="Geospatial table to visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial view", description="Showing a map with the geospatial data"
)
class ViewNodeStaticHTML:
    """
    This node will visualize the given geometric elements on a static map.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,  
        include_row_key=False,
        include_none_column=False,
    )

    color_col = knext.ColumnParameter(
        "Marker color column",
        "Select marker color column. The column must contain the color name e.g. red, green, blue, etc.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    color_map = knext.StringParameter(
        "Color map",
        "Select the color map to use for the color column. See https://matplotlib.org/stable/tutorials/colors/colormaps.html",
        default_value="viridis",
        enum=["viridis", "plasma", "inferno", "magma", "cividis"],
        
    )

    # base_map = knext.StringParameter(
    #     "Base map",
    #     "Select the base map to use for the visualization. See https://contextily.readthedocs.io/en/latest/providers_deepdive.html",
    #     default_value="OpenStreetMap",
    #     enum=['OpenStreetMap.Mapnik',
    #          'OpenTopoMap',
    #          'Stamen.Toner',
    #          'Stamen.TonerLite',
    #          'Stamen.Terrain',
    #          'Stamen.TerrainBackground',
    #          'Stamen.Watercolor',
    #          'NASAGIBS.ViirsEarthAtNight2012',
    #          'CartoDB.Positron',
    #          'CartoDB.Voyager'
    #         ]
    # )

    use_classify = knext.BoolParameter(
        "Use classification",
        "If checked, the color column will be classified using the selected classification method.",
        default_value=True,
    )

    classification_method = knext.StringParameter(
        "Classification method",
        "Select the classification method to use for the color column.",
        default_value="EqualInterval",
        enum=['BoxPlot', 'EqualInterval', 'FisherJenks', 'FisherJenksSampled', 'HeadTailBreaks', 'JenksCaspall', 'JenksCaspallForced', 'JenksCaspallSampled', 'MaxP', 'MaximumBreaks', 'NaturalBreaks', 'Quantiles', 'Percentiles', 'StdMean']
    )

    classification_bins = knext.IntParameter(
        "Number of classes",
        "Select the number of classes to use for the color column.",
        default_value=5,
        min_value=1,
        max_value=10,
    )

    plot_legend = knext.BoolParameter(
        "Show legend",
        "If checked, a legend will be shown in the plot.",
        default_value=True,

    )

    # size_col = knext.ColumnParameter(
    #     "Marker size column",
    #     "Select marker size column. The column must contain the size value.",
    #     column_filter=knut.is_numeric,
    # )

    legend_caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the legend.",
        default_value="",
        # default_value=color_col,
    )
    
    legend_location = knext.StringParameter(
        "Legend location",
        "Select the location for the legend.",
        default_value="best",
        enum=['best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center']
    )

    legend_orientation = knext.StringParameter(
        "Legend orientation",
        "Select the orientation for the legend.",
        default_value="vertical",
        enum=['vertical', 'horizontal']
    )


    def configure(self, configure_context, input_schema):
        knut.columns_exist([ self.geo_col], input_schema)
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        if (self.legend_caption is None) or (self.legend_caption == ""):
            self.legend_caption = self.color_col


# hvplot
        # m= gdf.hvplot(geo=True,  tiles='OSM', cmap='viridis', c=self.color_col, edgecolor='black')
        
        # html =StringIO()
        # hvplot.save(m,html)
        # html.seek(0)
        # html = html.read()


        # plotly
        # df = px.data.gapminder()
        # fig = px.scatter_geo(df, locations="iso_alpha", color="continent", hover_name="country", size="pop",
        #             animation_frame="year", projection="natural earth")
        # # fig.show()
        # html = io.to_html(fig, full_html=False)


# pydeck
        # CPU_GRID_LAYER_DATA = (
        #     "https://raw.githubusercontent.com/uber-common/" "deck.gl-data/master/website/sf-bike-parking.json"
        # )
        # df = pd.read_json(CPU_GRID_LAYER_DATA)

        # # Define a layer to display on a map

        # layer = pdk.Layer(
        #     "GridLayer",
        #     df,
        #     pickable=True,
        #     extruded=True,
        #     cell_size=200,
        #     elevation_scale=4,
        #     get_position="COORDINATES",
        # )

        # view_state = pdk.ViewState(latitude=37.7749295, longitude=-122.4194155, zoom=11, bearing=0, pitch=45)

        # # Render
        # r = pdk.Deck(
        #     layers=[layer],
        #     initial_view_state=view_state,
        #     tooltip={"text": "{position}\nCount: {count}"},
        # )

        # html = r._repr_html_()

# keplrgl

        # map_1 = KeplerGl()
        # html = map_1._repr_html_()
        # html = html.decode("utf-8")

# pyecharts


        bar = Bar()
        bar.add_xaxis(["衬衫", "羊毛衫", "雪纺衫", "裤子", "高跟鞋", "袜子"])
        bar.add_yaxis("商家A", [5, 20, 36, 10, 75, 90])
        # render 会生成本地 HTML 文件，默认会在当前目录生成 render.html 文件
        # 也可以传入路径参数，如 bar.render("mycharts.html")
        # html = bar.render_notebook()._repr_html_()
        html = bar.render_embed()

        return knext.view_html(html)


# TODO:
# make point view node interactive and static map
# set size of the point, set color of the point
# make polygon view node interactive and static map
# only set the color of the polygon
# multi-layer map
# get two layers of data, and show them on the map
# make dynamic map
# not support yet
# density map
# line view
