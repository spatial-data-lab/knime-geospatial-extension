import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
from keplergl import KeplerGl
import json
from folium import plugins
import folium


category = knext.category(
    path="/geo",
    level_id="viz",
    name="Spatial Visualization",  # Spatial Visualization
    description="Spatial view nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/VisulizationCategory.png",
    after="conversion",
)
# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Visulization/"


@knext.node(
    name="Geospatial View",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "InteractiveMap.png",
    category=category,
)
@knext.input_table(
    name="Geospatial table to visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial view",
    description="Showing a interactive map with the geospatial data",
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
        include_none_column=False,  # must contains a geometry column
    )

    color_col = knext.ColumnParameter(
        "Marker color column",
        "Select marker color column to be plotted.",
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    color_map = knext.StringParameter(
        "Color map",
        "Select the color map to use for the color column. `xxx_r` mean reverse of the `xxx` colormap. See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html)",
        default_value="viridis",
        enum=[
            "viridis",
            "plasma",
            "inferno",
            "magma",
            "cividis",
            "Greys",
            "Purples",
            "Blues",
            "Greens",
            "Oranges",
            "Reds",
            "YlOrBr",
            "YlOrRd",
            "OrRd",
            "PuRd",
            "RdPu",
            "BuPu",
            "GnBu",
            "PuBu",
            "YlGnBu",
            "PuBuGn",
            "BuGn",
            "YlGn",
            "binary",
            "gist_yarg",
            "gist_gray",
            "gray",
            "bone",
            "pink",
            "spring",
            "summer",
            "autumn",
            "winter",
            "cool",
            "Wistia",
            "hot",
            "afmhot",
            "gist_heat",
            "copper",
            "PiYG",
            "PRGn",
            "BrBG",
            "PuOr",
            "RdGy",
            "RdBu",
            "RdYlBu",
            "RdYlGn",
            "Spectral",
            "coolwarm",
            "bwr",
            "seismic",
            "viridis_r",
            "plasma_r",
            "inferno_r",
            "magma_r",
            "cividis_r",
            "Greys_r",
            "Purples_r",
            "Blues_r",
            "Greens_r",
            "Oranges_r",
            "Reds_r",
            "YlOrBr_r",
            "YlOrRd_r",
            "OrRd_r",
            "PuRd_r",
            "RdPu_r",
            "BuPu_r",
            "GnBu_r",
            "PuBu_r",
            "YlGnBu_r",
            "PuBuGn_r",
            "BuGn_r",
            "YlGn_r",
            "binary_r",
            "gist_yarg_r",
            "gist_gray_r",
            "gray_r",
            "bone_r",
            "pink_r",
            "spring_r",
            "summer_r",
            "autumn_r",
            "winter_r",
            "cool_r",
            "Wistia_r",
            "hot_r",
            "afmhot_r",
            "gist_heat_r",
            "copper_r",
            "PiYG_r",
            "PRGn_r",
            "BrBG_r",
            "PuOr_r",
            "RdGy_r",
            "RdBu_r",
            "RdYlBu_r",
            "RdYlGn_r",
            "Spectral_r",
            "coolwarm_r",
            "bwr_r",
            "seismic_r",
        ],
    )

    stroke = knext.BoolParameter(
        "Stroke",
        "Whether to draw strokes along the polygon boundary.",
        default_value=True,
    )

    base_map = knext.StringParameter(
        "Base map",
        "Select the base map to use for the visualization. See [Folium base maps](https://python-visualization.github.io/folium/quickstart.html#Tiles).",
        default_value="OpenStreetMap",
        enum=[
            "OpenStreetMap",
            "Stamen Watercolor",
            "Stamen Toner",
            "Stamen TonerBackground",
            "Stamen TonerHybrid",
            "Stamen TonerLines",
            "Stamen TonerLabels",
            "Stamen TonerLite",
            "Stamen Watercolor",
            "Stamen Terrain",
            "Stamen TerrainBackground",
            "Stamen TerrainLabels",
            "Stamen TopOSMRelief",
            "Stamen TopOSMFeatures",
            "CartoDB Positron",
            "CartoDB PositronNoLabels",
            "CartoDB PositronOnlyLabels",
            "CartoDB DarkMatter",
            "CartoDB DarkMatterNoLabels",
            "CartoDB DarkMatterOnlyLabels",
            "CartoDB Voyager",
            "CartoDB VoyagerNoLabels",
            "CartoDB VoyagerOnlyLabels",
            "CartoDB VoyagerLabelsUnder",
            "NASAGIBS ModisTerraTrueColorCR",
            "NASAGIBS ModisTerraBands367CR",
            "NASAGIBS ViirsEarthAtNight2012",
            "NASAGIBS ModisTerraLSTDay",
            "NASAGIBS ModisTerraSnowCover",
            "NASAGIBS ModisTerraAOD",
            "NASAGIBS ModisTerraChlorophyll",
            "NASAGIBS ModisTerraBands721CR",
            "NASAGIBS ModisAquaTrueColorCR",
            "NASAGIBS ModisAquaBands721CR",
            "NASAGIBS ViirsTrueColorCR",
            "NASAGIBS BlueMarble3413",
            "NASAGIBS BlueMarble3031",
            "NASAGIBS BlueMarble",
            "NASAGIBS ASTER_GDEM_Greyscale_Shaded_Relief",
            "Esri WorldStreetMap",
            "Esri DeLorme",
            "Esri WorldTopoMap",
            "Esri WorldImagery",
            "Esri WorldTerrain",
            "Esri WorldShadedRelief",
            "Esri WorldPhysical",
            "Esri OceanBasemap",
            "Esri NatGeoWorldMap",
            "Esri WorldGrayCanvas",
            "Gaode Normal",
            "Gaode Satellite",
            "OpenRailwayMap",
            "Strava All",
            "Strava Ride",
            "Strava Run",
            "Strava Water",
            "Strava Winter",
        ],
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
        enum=[
            "BoxPlot",
            "EqualInterval",
            "FisherJenks",
            "FisherJenksSampled",
            "HeadTailBreaks",
            "JenksCaspall",
            "JenksCaspallForced",
            "JenksCaspallSampled",
            "MaxP",
            "MaximumBreaks",
            "NaturalBreaks",
            "Quantiles",
            "Percentiles",
            "StdMean",
        ],
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

    size_scale = knext.IntParameter(
        "Size scale",
        "Select the size scale of the markers.",
        default_value=1,
        min_value=1,
        max_value=100,
    )

    name_cols = knext.MultiColumnParameter(
        "Tooltip columns",
        "Select columns which should be shown in the marker tooltip.",
        column_filter=knut.is_numeric_or_string,
    )

    popup_cols = knext.MultiColumnParameter(
        "Popup columns",
        "Select columns which should be shown in the marker popup.",
        column_filter=knut.is_numeric_or_string,
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
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        kws = {
            # "column":self.color_col,
            # "cmap":self.color_map,
            "tooltip": self.name_cols,
            "tiles": self.base_map,
            "popup": self.popup_cols,
            "legend": self.plot_legend,
            "m": None,
            "legend_kwds": {
                "caption": self.legend_caption,
                "scale": False,
                "max_labels": 3,
                "colorbar": True,
            },
        }

        if "none" not in str(self.color_col).lower():
            kws["column"] = self.color_col
            kws["cmap"] = self.color_map

            if type(gdf[self.color_col].values[0]) == str:
                kws["categorical"] = True
                kws["legend_kwds"]["colorbar"] = False

            if (self.legend_caption is None) or (self.legend_caption == ""):
                self.legend_caption = self.color_col
        else:
            self.plot_legend = False

        if self.use_classify:
            kws["scheme"] = self.classification_method
            kws["k"] = self.classification_bins
            kws["legend_kwds"]["colorbar"] = False
            kws["legend_kwds"]["max_labels"] = 20

        if "none" not in str(self.size_col).lower():

            max_pop_est = gdf[self.size_col].max()
            min_pop_est = gdf[self.size_col].min()

            # check whether is line
            # FIXME: change it to use utlis.is_line
            geo_types = gdf[self.geo_col].geom_type.unique()
            if ("LineString" in geo_types) or ("MultiLineString" in geo_types):
                max_size = 8
                kws["style_kwds"] = {
                    "style_function": lambda x: {
                        "weight": (x["properties"][self.size_col] - min_pop_est)
                        / (max_pop_est - min_pop_est)
                        * max_size
                    }
                }
            elif ("Polygon" in geo_types) or ("MultiPolygon" in geo_types):
                max_size = 30
                kws["style_kwds"] = {
                    "style_function": lambda x: {
                        "radius": (x["properties"][self.size_col] - min_pop_est)
                        / (max_pop_est - min_pop_est)
                        * max_size
                    }
                }
                kws["m"] = gdf.explore(tiles=self.base_map)
                gdf[self.geo_col] = gdf.centroid

            else:
                max_size = 30
                kws["style_kwds"] = {
                    "style_function": lambda x: {
                        "radius": (x["properties"][self.size_col] - min_pop_est)
                        / (max_pop_est - min_pop_est)
                        * max_size
                    }
                }
        if "none" not in str(self.size_scale).lower():
            geo_types = gdf[self.geo_col].geom_type.unique()
            if ("LineString" in geo_types) or ("MultiLineString" in geo_types):
                kws["style_kwds"] = {"weight": self.size_scale}
            elif ("Polygon" in geo_types) or ("MultiPolygon" in geo_types):
                pass
            else:
                kws["style_kwds"] = {"radius": self.size_scale}
        if not self.stroke:
            if "style_kwds" in kws:
                kws["style_kwds"]["stroke"] = False
            else:
                kws["style_kwds"] = {"stroke": False}
        map = gdf.explore(**kws)
        # knut.check_canceled(exec_context)
        return knext.view(map)


# geo view static


@knext.node(
    name="Geospatial View Static",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "StaticMap.png",
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
        "Select one column to map to the marker color. If you select none, it will not map any column to color.",
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=True,
    )

    color = knext.StringParameter(
        "Marker color",
        "Select marker color. Can select none if you don't want to set a unified marker color ",
        default_value="none",
        enum=[
            "none",
            "red",
            "blue",
            "green",
            "orange",
            "purple",
            "darkred",
            "lightred",
            "beige",
            "darkblue",
            "darkgreen",
            "cadetblue",
            "darkpurple",
            "white",
            "pink",
            "lightblue",
            "lightgreen",
            "gray",
            "black",
            "lightgray",
        ],
    )

    color_map = knext.StringParameter(
        "Color map",
        "Select the color map to use for the color column. See https://matplotlib.org/stable/tutorials/colors/colormaps.html",
        default_value="viridis",
        enum=["viridis", "plasma", "inferno", "magma", "cividis"],
    )

    edge_color = knext.StringParameter(
        "Edge color",
        "Set the edge color. See https://matplotlib.org/stable/tutorials/colors/colormaps.html",
        default_value="none",
        enum=[
            "none",
            "black",
            "red",
            "blue",
            "green",
            "orange",
            "purple",
            "darkred",
            "lightred",
            "beige",
            "darkblue",
            "darkgreen",
            "cadetblue",
            "darkpurple",
            "white",
            "pink",
            "lightblue",
            "lightgreen",
            "gray",
            "black",
            "lightgray",
        ],
    )

    size_col = knext.ColumnParameter(
        "Marker size column",
        "Select marker size column. The size is fixed by default. If a size column is selected, the size will be scaled by the values of the column. For point features, the size is the radius of the circle. For line features, the size is the width of the line. For polygon features, the size is the radius of the centroid of the ploygon.",
        column_filter=knut.is_numeric,
        include_none_column=True,
    )

    line_width_col = knext.ColumnParameter(
        "Line width column",
        "Select line width column. The width is fixed by default. If a width column is selected, the width will be scaled by the values of the column.",
        column_filter=knut.is_numeric,
        include_none_column=True,
    )

    line_width = knext.IntParameter(
        "Line width",
        "Select a unified line width, can be set to none",
        default_value=1,
        min_value=1,
        max_value=10,
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
        enum=[
            "BoxPlot",
            "EqualInterval",
            "FisherJenks",
            "FisherJenksSampled",
            "HeadTailBreaks",
            "JenksCaspall",
            "JenksCaspallForced",
            "JenksCaspallSampled",
            "MaxP",
            "MaximumBreaks",
            "NaturalBreaks",
            "Quantiles",
            "Percentiles",
            "StdMean",
        ],
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
        default_value="lower right",
        enum=[
            "best",
            "upper right",
            "upper left",
            "lower left",
            "lower right",
            "right",
            "center left",
            "center right",
            "lower center",
            "upper center",
            "center",
            "outside_top",
            "outside_bottom",
        ],
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
        enum=["black", "red", "green", "blue", "yellow", "purple", "orange", "white"],
    )

    legend_frame = knext.BoolParameter(
        "Show legend frame",
        "If checked, a frame will be shown in the legend.",
        default_value=True,
    )

    legend_framealpha = knext.DoubleParameter(
        "Legend frame alpha",
        "Select the transparent value for the legend frame.",
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

    legend_colorbar_shrink = knext.DoubleParameter(
        "Colorbar legend shrink",
        "Select the shrink value for the colorbar legend. Only work for colorbar",
        default_value=1.0,
        min_value=0.0,
        max_value=1.0,
    )

    legend_colorbar_pad = knext.DoubleParameter(
        "Colorbar legend pad",
        "Select the pad value for the colorbar legend. Only work for colorbar",
        default_value=0.1,
        min_value=0.0,
        max_value=0.99,
    )

    set_axis_off = knext.BoolParameter(
        "Set axis off",
        "If checked, the axis will be set off.",
        default_value=False,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        # check legend caption
        if (self.legend_caption is None) or (self.legend_caption == ""):
            if "none" not in str(self.color_col).lower():
                self.legend_caption = self.color_col

        #  set legend location
        if self.legend_location == "outside_top":
            colorbar_legend_location = "top"
        elif self.legend_location == "outside_bottom":
            colorbar_legend_location = "bottom"
        else:
            colorbar_legend_location = "right"

        legend_bbox_to_anchor = None
        if self.legend_location == "outside_top":
            self.legend_location = "lower right"
            legend_bbox_to_anchor = (0.0, 1.02, 1.0, 0.102)
        if self.legend_location == "outside_bottom":
            self.legend_location = "upper right"
            legend_bbox_to_anchor = (0.0, -0.2, 1.0, 0.102)

        if self.legend_expand:
            legend_expand = "expand"
        else:
            legend_expand = None

        kws = {"alpha": 1, "legend": self.plot_legend, "aspect": 1}

        if "none" not in str(self.edge_color):
            kws["edgecolor"] = self.edge_color
        if "none" not in str(self.color_col).lower():
            kws["column"] = self.color_col
            kws["cmap"] = self.color_map
        if "none" not in str(self.color).lower():
            kws["color"] = self.color

        if ("none" not in str(self.edge_color)) or (
            "none" not in str(self.color_col).lower()
        ):
            if self.use_classify:
                kws["legend_kwds"] = {
                    "fmt": "{:.0f}",
                    "loc": self.legend_location,
                    "title": self.legend_caption,
                    "ncols": self.legend_columns,
                    "prop": {"size": self.legend_size},
                    "fontsize": self.legend_fontsize,
                    "bbox_to_anchor": legend_bbox_to_anchor,
                    "labelcolor": self.legend_labelcolor,
                    "frameon": self.legend_frame,
                    "framealpha": self.legend_framealpha,
                    "fancybox": True,
                    "mode": legend_expand,
                    "alignment": "left",
                    # "title": "Population",
                    "title_fontsize": self.legend_caption_fontsize,
                    "labelspacing": self.legend_labelspacing,
                    "borderaxespad": self.legend_borderpad,
                }
                kws["scheme"] = self.classification_method
                kws["k"] = self.classification_bins
            else:
                kws["legend_kwds"] = {
                    "shrink": self.legend_colorbar_shrink,
                    "fmt": "{:.0f}",
                    "location": colorbar_legend_location,
                    "pad": self.legend_colorbar_pad,
                }

        geo_types = gdf[self.geo_col].geom_type.unique()
        if not (("LineString" in geo_types) or ("MultiLineString" in geo_types)):
            if "none" not in str(self.size_col).lower():
                max_point_size = 2000
                max_val = gdf[self.size_col].max()
                min_val = gdf[self.size_col].min()
                normal_base = (gdf[self.size_col] - min_val) / max_val
                kws["makersize"] = normal_base * max_point_size
        if "none" not in str(self.line_width_col).lower():
            max_line_width = 5
            max_val = gdf[self.line_width_col].max()
            min_val = gdf[self.line_width_col].min()
            normal_base = (gdf[self.line_width_col] - min_val) / max_val
            kws["linewidth"] = normal_base * max_line_width

        map = gdf.plot(**kws)
        map.set_title(self.figure_title, fontsize=self.figure_title_size)
        if self.set_axis_off:
            map.set_axis_off()
        # knut.check_canceled(exec_context)
        # cx.add_basemap(map, crs=gdf.crs.to_string(), source=cx.providers.flatten()[self.base_map])

        return knext.view_matplotlib(map.get_figure())


############################################
# Kepler.gl
############################################
@knext.node(
    name="Kepler.gl Geoview ",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "Kepler.gl.png",
    category=category,
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
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        df = input_table.to_pandas()
        df.rename(columns={self.geo_col: "geometry"}, inplace=True)
        gdf = gp.GeoDataFrame(df, geometry="geometry")
        
        map_1 = KeplerGl(show_docs=False)
        map_1.add_data(data=gdf.copy(), name="state")
        config = {}
        if self.save_config:
            # Save map_1 config to a file
            # config_str = json.dumps(map_1.config)
            # if type(config) == str:
            #     config = config.encode("utf-8")
            with open("kepler_config.json", "w") as f:
                f.write(json.dumps(map_1.config))

        if self.load_config:
            with open("kepler_config.json", "r") as f:
                config = json.loads(f.read())
        map_1.config = config

        # map_1.add_data(data=data.copy(), name="haha")
        html = map_1._repr_html_()
        html = html.decode("utf-8")
        # knext.view_html(html)
        # knut.check_canceled(exec_context)
        # cx.add_basemap(map, crs=gdf.crs.to_string(), source=cx.providers.flatten()[self.base_map])

        return knext.view_html(html)


############################################
# heatmap node
############################################


@knext.node(
    name="Heatmap",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "Heatmap.png",
    category=category,
)
@knext.input_table(
    name="Table to visualize",
    description="Table with data to visualize",
)
@knext.output_view(name="Heatmap view", description="Showing a heatmap with the data")
class ViewNodeHeatmap:
    """
    This node will visualize the given data on a heatmap. More details about heatmap [here](https://www.gislounge.com/heat-maps-in-gis/).
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    weight_col = knext.ColumnParameter(
        "Weight column",
        "Select the weight column to visualize.",
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=True,
    )

    min_opacity = knext.DoubleParameter(
        "Minimum opacity",
        "The minimum opacity the heat will start at.",
        default_value=0.5,
    )

    max_zoom = knext.IntParameter(
        "Maximum zoom",
        "Zoom level where the points reach maximum intensity (as intensity scales with zoom), equals maxZoom of the map by default",
        default_value=18,
    )

    radius = knext.IntParameter(
        "Radius",
        "Radius of each “point” of the heatmap",
        default_value=25,
    )

    blur = knext.IntParameter(
        "Blur",
        "Amount of blur",
        default_value=15,
    )

    def configure(self, configure_context, input_schema):

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):

        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        map = gdf.explore()
        gdf["lon"] = gdf.geometry.centroid.x
        gdf["lat"] = gdf.geometry.centroid.y

        if "none" not in str(self.weight_col).lower():
            heat_data = gdf[["lat", "lon", self.weight_col]]
        else:
            heat_data = gdf[["lat", "lon"]]

        plugins.HeatMap(
            heat_data,
            name="heatmap",
            min_opacity=self.min_opacity,
            max_zoom=self.max_zoom,
            radius=self.radius,
            blur=self.blur,
        ).add_to(map)

        folium.LayerControl().add_to(map)

        return knext.view(map)
