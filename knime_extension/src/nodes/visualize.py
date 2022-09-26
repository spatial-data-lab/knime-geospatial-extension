from doctest import debug_script
import logging

import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
# import contextily as cx
import folium
import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/geo",
    level_id="viz",
    name="Views",
    description="Spatial view nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/visualize.png",
)


@knext.node(
    name="Geospatial View",
    node_type=knext.NodeType.VISUALIZER,
    icon_path="icons/visualize.png",
    category=category,
)
@knext.input_table(
    name="Geospatial table to visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial view", description="Showing a map with the geospatial data"
)
class ViewNode:
    """
    This node will visualize the given geometric elements on a map.
    """

    # geo_col = knext.ColumnParameter(
    #     "Geometry column",
    #     "Select the geometry column to visualize.",
    #     column_filter=knut.is_geo_point,  # Allows only GeoPoints
    #     include_row_key=False,
    #     include_none_column=False,
    # )

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

    base_map = knext.StringParameter(
        "Base map",
        "Select the base map to use for the visualization. See https://python-visualization.github.io/folium/quickstart.html#Tiles",
        default_value="OpenStreetMap",
        enum=["OpenStreetMap", "Stamen Terrain", "Stamen Toner", "Stamen Watercolor" "CartoDB positron", "CartoDB dark_matter"]
    )

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

    size_col = knext.ColumnParameter(
        "Marker size column",
        "Select marker size column. The column must contain the size value.",
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
        "Set the caption for the legend.",
        default_value="Legend",
    )


    def configure(self, configure_context, input_schema):
        # knut.columns_exist([ self.color_col,self.name_cols], input_schema)
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry="geometry")

        if (self.legend_caption is None) or (self.legend_caption == ""):
            self.legend_caption = self.color_col
        if "none" not in self.size_col:

            
            max_pop_est = gdf[self.size_col].max()
            min_pop_est = gdf[self.size_col].min()

            # check whether is line 
            geo_types = gdf["geometry"].geom_type.unique()
            if  ("LineString" in geo_types) or ("MultiLineString" in geo_types):
                size = "weight"
                max_size = 8
            else:
                size = "radius"
                gdf["geometry"] = gdf.centroid
                max_size = 30

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
@knext.node(
    name="Geospatial View Static",
    node_type=knext.NodeType.VISUALIZER,
    icon_path="icons/visualize.png",
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
    This node will visualize the given geometric elements on a map.
    """

    # geo_col = knext.ColumnParameter(
    #     "Geometry column",
    #     "Select the geometry column to visualize.",
    #     column_filter=knut.is_geo_point,  # Allows only GeoPoints
    #     include_row_key=False,
    #     include_none_column=False,
    # )

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

    base_map = knext.StringParameter(
        "Base map",
        "Select the base map to use for the visualization. See https://contextily.readthedocs.io/en/latest/providers_deepdive.html",
        default_value="OpenStreetMap",
        enum=['OpenStreetMap.Mapnik',
             'OpenTopoMap',
             'Stamen.Toner',
             'Stamen.TonerLite',
             'Stamen.Terrain',
             'Stamen.TerrainBackground',
             'Stamen.Watercolor',
             'NASAGIBS.ViirsEarthAtNight2012',
             'CartoDB.Positron',
             'CartoDB.Voyager'
            ]
    )

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
        # default_value=color_col,
    )

    def configure(self, configure_context, input_schema):
        # knut.columns_exist([ self.color_col,self.name_cols], input_schema)
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry="geometry")

        # if self.size_col is not None:
        #     gdf["geometry"] = gdf.centroid

        if self.use_classify:
            map = gdf.plot(
                column=self.color_col, 
                cmap=self.color_map,
                # tiles=self.base_map,
                alpha=0.7,
                scheme=self.classification_method,
                k=self.classification_bins,
                legend=self.plot_legend,
                # marker_kwds={"radius":self.size_col},
                # legend_kwds={'shrink': 0.6}
                # legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":20,"colorbar":False}
                )
        else:
            map = gdf.plot(
                column=self.color_col, 
                cmap=self.color_map,
                # tiles=self.base_map,
                alpha=0.7,
                legend=self.plot_legend,
                # marker_kwds={"radius":self.size_col},
                # legend_kwds={'shrink': 0.6}
                # legend_kwds={"caption": self.legend_caption,"scale":False,"max_labels":3,"colorbar":True}
            )
    
        # knut.check_canceled(exec_context)
        # return knext.
        # cx.add_basemap(map, crs=gdf.crs.to_string(), source=cx.providers.flatten()[self.base_map])

        return knext.view_matplotlib(map.get_figure())



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
