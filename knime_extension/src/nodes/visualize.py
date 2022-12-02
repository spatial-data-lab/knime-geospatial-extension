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


# helper class for geoview
# TODO: add this class
# TODO: add this class
# class __GeoView:
#     def __init__(self, geo_df, color_col, color_map, geo_col):
#         self.geo_df = geo_df
#         self.color_col = color_col
#         self.color_map = color_map
#         self.geo_col = geo_col

#     def __repr__(self):
#         return self.geo_df.to_json()

#     def __str__(self):
#         return self.geo_df.to_json()

#     def _repr_html_(self):
#         return self.geo_df._repr_html_()

#     def _repr_mimebundle_(self, include=None, exclude=None):
#         return self.geo_df._repr_mimebundle_(include, exclude)


def replace_external_js_css_paths(
    replacement: str,
    html: str,
    regex_patter: str = '(<script src="|<link( rel="stylesheet")? href=")https?[^"]*\/([^"]*)"( rel="stylesheet")?',
) -> str:
    """
    Uses a regular expression to find all script and stylesheet tags in a given html page.
    The first matching group is either the script ar stylesheet part up until the opening " of the
    URL. The second matching group is the file name.
    The method will also add the closing ".

    For example if the html code is <script src="https://cdn.jsdelivr.net/npm/leaflet@1.6.0/dist/leaflet.js">
    the first group is <script src=" and the second group is leaflet.js so using the following replacement
    r'\1./libs/kepler/2.5.5/\3\"\4' will lead to this URL: <script src="./libs/leaflet/1.6.0/leaflet.js">.
    """
    import re

    result = re.sub(
        regex_patter,
        replacement,
        html,
    )
    return result


@knext.parameter_group(label="Coloring Settings")
class ColorSettings:
    """
    Group of settings that define coloring of the geometric objects. The color column can be either nominal or numerical.
    If a numerical column is selected you might want to enable the classification of the numeric values to group
    them into bins prior assigning a color to each bin using the color map information.
    """

    color_col = knext.ColumnParameter(
        "Marker color column",
        """Select marker color column to be plotted. If numerical you might want to adapt the classification 
        settings accordingly.""",
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=True,
    )

    color_map = knext.StringParameter(
        "Color map",
        """Select the color map to use for the color column. `xxx_r` mean reverse of the `xxx` colormap. 
        See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html)""",
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

    use_classify = knext.BoolParameter(
        "Classify numerical marker color columns",
        """If checked, a numerical marker color column will be classified using the selected classification method. 
        The 'Number of classes' will be used to determine the number of bins.""",
        default_value=True,
    )

    classification_method = knext.StringParameter(
        "Classification method",
        "Select the classification method to use for the selected numerical marker color column.",
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
        "Select the number of classes used by the classification method.",
        default_value=5,
        min_value=1,
        max_value=50,
    )


@knext.parameter_group(label="Legend Settings")
class LegendSettings:
    """
    Group of settings that define if a legend is shown on the map and if so how it should be formatted.
    """

    plot = knext.BoolParameter(
        "Show legend",
        "If checked, a legend will be shown in the plot.",
        default_value=True,
    )

    caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the legend. By default, the caption is the name of the selected color column.",
        default_value="",
    )


@knext.node(
    name="Geospatial View",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "InteractiveMap.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geospatial Table to Visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial View",
    description="Showing a interactive map with the geospatial data",
    static_resources="libs/leaflet/1.6.0",
)
class ViewNode:
    """Creates an interactive map view based on the selected geometric elements of the input table.
    This node creates an interactive map view based on the selected geometric elements of the input table.
    It provides various dialog options to modify the appearance of teh view e.g. the base map, shape color and
    size.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        # "geometry",
        column_filter=knut.is_geo,  # Allows all geo columns
        include_row_key=False,
        include_none_column=False,  # must contains a geometry column
    )

    base_map = knext.StringParameter(
        "Base map",
        """Select the base map to use for the visualization. 
        See [Folium base maps](https://python-visualization.github.io/folium/quickstart.html#Tiles).""",
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

    stroke = knext.BoolParameter(
        "Stroke",
        "Whether to draw strokes along the polygon boundary.",
        default_value=True,
    )

    size_col = knext.ColumnParameter(
        "Marker size column",
        """Select marker size column. The size is fixed by default. If a size column is selected, the size will be 
        scaled by the values of the column. For point features, the size is the radius of the circle. 
        For line features, the size is the width of the line. For polygon features, the size is the radius of 
        the centroid of the polygon.""",
        column_filter=knut.is_numeric,
        include_none_column=True,
    )

    size_scale = knext.IntParameter(
        "Marker size scale",
        "Select the size scale of the markers.",
        default_value=1,
        min_value=1,
        max_value=100,
    )

    name_cols = knext.MultiColumnParameter(
        "Marker tooltip columns",
        "Select columns which should be shown in the marker tooltip.",
        column_filter=knut.is_numeric_or_string,
    )

    popup_cols = knext.MultiColumnParameter(
        "Marker popup columns",
        "Select columns which should be shown in the marker popup.",
        column_filter=knut.is_numeric_or_string,
    )

    color_settings = ColorSettings()

    legend_settings = LegendSettings()

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
            "legend": self.legend_settings.plot,
            "m": None,
            "legend_kwds": {
                "caption": self.legend_settings.caption,
                "scale": False,
                "max_labels": 3,
                "colorbar": True,
            },
        }

        if "none" not in str(self.color_settings.color_col).lower():
            kws["column"] = self.color_settings.color_col
            kws["cmap"] = self.color_settings.color_map

            if type(gdf[self.color_settings.color_col].values[0]) == str:
                kws["categorical"] = True
                kws["legend_kwds"]["colorbar"] = False

            if (self.legend_settings.caption is None) or (
                self.legend_settings.caption == ""
            ):
                self.legend_settings.caption = self.color_settings.color_col
        else:
            self.legend_settings.plot = False

        if self.color_settings.use_classify:
            kws["scheme"] = self.color_settings.classification_method
            kws["k"] = self.color_settings.classification_bins
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
        elif "none" not in str(self.size_scale).lower():
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

        # replace css and JavaScript paths
        html = map.get_root().render()
        html = replace_external_js_css_paths(
            r"\1./libs/leaflet/1.6.0/\3\"\4",
            html,
        )

        return knext.view(html)


# geo view static


@knext.parameter_group(label="Coloring Settings")
class StaticColorSettings:
    """
    Group of settings that define coloring of the geometric objects. The color column can be either nominal or numerical.
    If a numerical column is selected you might want to enable the classification of the numeric values to group
    them into bins prior assigning a color to each bin using the color map information.
    """

    color_col = knext.ColumnParameter(
        "Marker color column",
        """Select marker color column to be plotted. If numerical you might want to adapt the classification 
        settings accordingly.""",
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=True,
    )

    color = knext.StringParameter(
        "Marker color",
        "Select marker color. Select none if you don't want to set a unified marker color.",
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
        """Select the color map to use for the color column. 
        See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html).""",
        default_value="viridis",
        enum=["viridis", "plasma", "inferno", "magma", "cividis"],
    )

    use_classify = knext.BoolParameter(
        "Classify numerical marker color columns",
        """If checked, the numerical marker color column will be classified using the selected classification method. 
        # The 'Number of classes' will be used to determine the number of bins.""",
        default_value=True,
    )

    classification_method = knext.StringParameter(
        "Classification method",
        "Select the classification method to use for the selected numerical marker color column.",
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
        "Select the number of classes used by the classification method.",
        default_value=5,
        min_value=1,
        max_value=50,
    )

    edge_color = knext.StringParameter(
        "Edge color",
        "Set the edge color. See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html).",
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


@knext.parameter_group(label="Legend Settings")
class StaticLegendSettings:
    """
    Group of settings that define if a legend is shown on the map and if so how it should be formatted.
    """

    plot = knext.BoolParameter(
        "Show legend",
        "If checked, a legend will be shown in the plot.",
        default_value=True,
    )

    caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the legend.",
        default_value="",
        # default_value=color_col,
    )

    caption_fontsize = knext.IntParameter(
        "Legend caption font size",
        "Set the font size for the legend caption.",
        default_value=10,
        min_value=1,
        max_value=100,
    )

    expand = knext.BoolParameter(
        "Expand legend",
        "If checked, the legend will be horizontally expanded to fill the axes area",
        default_value=False,
    )

    location = knext.StringParameter(
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

    columns = knext.IntParameter(
        "Legend columns",
        "Select the number of columns for the legend.",
        default_value=1,
        min_value=1,
        max_value=30,
    )

    size = knext.IntParameter(
        "Legend size",
        "Select the size for the legend.",
        default_value=8,
        min_value=1,
        max_value=30,
    )

    fontsize = knext.IntParameter(
        "Legend font size",
        "Select the font size for the legend.",
        default_value=10,
        min_value=1,
        max_value=30,
    )

    labelcolor = knext.StringParameter(
        "Legend label color",
        "Select the label color for the legend.",
        default_value="black",
        enum=["black", "red", "green", "blue", "yellow", "purple", "orange", "white"],
    )

    frame = knext.BoolParameter(
        "Show legend frame",
        "If checked, a frame will be shown in the legend.",
        default_value=True,
    )

    framealpha = knext.DoubleParameter(
        "Legend frame alpha",
        "Select the transparent value for the legend frame.",
        default_value=1.0,
        min_value=0.0,
        max_value=1.0,
    )

    borderpad = knext.DoubleParameter(
        "Legend border pad",
        "Select the border pad for the legend.",
        default_value=0.5,
        min_value=0.0,
        max_value=3.0,
    )

    labelspacing = knext.DoubleParameter(
        "Legend label spacing",
        "Select the label spacing for the legend.",
        default_value=0.5,
        min_value=0.0,
        max_value=1.0,
    )

    colorbar_shrink = knext.DoubleParameter(
        "Colorbar legend shrink",
        "Select the shrink value for the colorbar legend. Only work for colorbar",
        default_value=1.0,
        min_value=0.0,
        max_value=1.0,
    )

    colorbar_pad = knext.DoubleParameter(
        "Colorbar legend pad",
        "Select the pad value for the colorbar legend. Only work for colorbar",
        default_value=0.1,
        min_value=0.0,
        max_value=0.99,
    )


@knext.node(
    name="Geospatial View Static",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "StaticMap.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geospatial Table to Visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial View", description="Showing a map with the geospatial data"
)
class ViewNodeStatic:
    """Creates a static visualization of the geometric elements.
    This node will visualize the given geometric elements using the [matplotlib}(https://matplotlib.org/stable/).
    It can be used to create [Choropleth Maps](https://en.wikipedia.org/wiki/Choropleth_map) by assigning
    a marker color column. The node further supports various settings to adapt the legend to your needs.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    size_col = knext.ColumnParameter(
        "Marker size column",
        """Select marker size column. The size is fixed by default. If a size column is selected, the size will be 
        scaled by the values of the column. For point features, the size is the radius of the circle. 
        For line features, the size is the width of the line. For polygon features, the size is the radius of 
        the centroid of the polygon.""",
        column_filter=knut.is_numeric,
        include_none_column=True,
    )

    line_width_col = knext.ColumnParameter(
        "Line width column",
        """Select line width column. The width is fixed by default. If a width column is selected, the width will 
        be scaled by the values of the column.""",
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

    # size_col = knext.ColumnParameter(
    #     "Marker size column",
    #     "Select marker size column. The column must contain the size value.",
    #     column_filter=knut.is_numeric,
    # )

    set_axis_off = knext.BoolParameter(
        "Set axis off",
        "If checked, the axis will be set off.",
        default_value=False,
    )

    color_settings = StaticColorSettings()

    legend_settings = StaticLegendSettings()

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
        if (self.legend_settings.caption is None) or (
            self.legend_settings.caption == ""
        ):
            if "none" not in str(self.color_settings.color_col).lower():
                self.legend_settings.caption = self.color_settings.color_col

        #  set legend location
        if self.legend_settings.location == "outside_top":
            colorbar_legend_location = "top"
        elif self.legend_settings.location == "outside_bottom":
            colorbar_legend_location = "bottom"
        else:
            colorbar_legend_location = "right"

        legend_bbox_to_anchor = None
        if self.legend_settings.location == "outside_top":
            self.legend_settings.location = "lower right"
            legend_bbox_to_anchor = (0.0, 1.02, 1.0, 0.102)
        if self.legend_settings.location == "outside_bottom":
            self.legend_settings.location = "upper right"
            legend_bbox_to_anchor = (0.0, -0.2, 1.0, 0.102)

        if self.legend_settings.expand:
            legend_expand = "expand"
        else:
            legend_expand = None

        kws = {"alpha": 1, "legend": self.legend_settings.plot, "aspect": 1}

        if "none" not in str(self.color_settings.edge_color):
            kws["edgecolor"] = self.color_settings.edge_color
        if "none" not in str(self.color_settings.color_col).lower():
            kws["column"] = self.color_settings.color_col
            kws["cmap"] = self.color_settings.color_map
        if "none" not in str(self.color_settings.color).lower():
            kws["color"] = self.color_settings.color

        if ("none" not in str(self.color_settings.edge_color)) or (
            "none" not in str(self.color_settings.color_col).lower()
        ):
            if self.color_settings.use_classify:
                kws["legend_kwds"] = {
                    "fmt": "{:.0f}",
                    "loc": self.legend_settings.location,
                    "title": self.legend_settings.caption,
                    "ncols": self.legend_settings.columns,
                    "prop": {"size": self.legend_settings.size},
                    "fontsize": self.legend_settings.fontsize,
                    "bbox_to_anchor": legend_bbox_to_anchor,
                    "labelcolor": self.legend_settings.labelcolor,
                    "frameon": self.legend_settings.frame,
                    "framealpha": self.legend_settings.framealpha,
                    "fancybox": True,
                    "mode": legend_expand,
                    "alignment": "left",
                    # "title": "Population",
                    "title_fontsize": self.legend_settings.caption_fontsize,
                    "labelspacing": self.legend_settings.labelspacing,
                    "borderaxespad": self.legend_settings.borderpad,
                }
                kws["scheme"] = self.color_settings.classification_method
                kws["k"] = self.color_settings.classification_bins
            else:
                kws["legend_kwds"] = {
                    "shrink": self.legend_settings.colorbar_shrink,
                    "fmt": "{:.0f}",
                    "location": colorbar_legend_location,
                    "pad": self.legend_settings.colorbar_pad,
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
    after="",
)
@knext.input_table(
    name="Geospatial Table to Visualize",
    description="Table with geospatial data to visualize",
)
@knext.output_view(
    name="Geospatial View",
    description="Showing a map with the geospatial data",
    static_resources="libs/kepler/2.5.5",
)
class ViewNodeKepler:
    """Visualizes given geometric elements on a map.

    This node will visualize the given geometric elements on a map using the [kepler.gl](https://kepler.gl/)
    visualization framework. This view is highly interactive and allows you to change various aspects of the view
    within the visualization itself e.g. adding [layers](https://docs.kepler.gl/docs/user-guides/c-types-of-layers)
    and [filters](https://docs.kepler.gl/docs/user-guides/e-filters). It also allows you to create a filter that
    creates and animation for a given timeseries column. For more information about the supported interactions
    see the [kepler.gl user guides](https://docs.kepler.gl/docs/user-guides).

    This node uses the [Mapbox GL JS API](https://www.mapbox.com/pricing#map-loads-for-web) which for commercial
    usage might require an [access token](https://docs.mapbox.com/help/glossary/access-token/).
    If you want to use a different base map, you can configure it inside the interactive
    view with Kepler.gl's UI. You can also configure the
    [Mapbox style](https://docs.kepler.gl/docs/user-guides/f-map-styles#custom-map-styles) you want to use and
    the access token.


    By default, it takes all column information that is included inside the input table.
    If you want to limit the amount of information send to th e node view you can use the one of the
    [column filter](https://kni.me/n/DOkyMaii62U05xZ1) nodes to filter the input table.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    # save_config = knext.BoolParameter(
    #     "Save config",
    #     "Save the config for the map",
    #     default_value=False,
    # )

    # load_config = knext.BoolParameter(
    #     "Load config",
    #     "Load the config for the map",
    #     default_value=False,
    # )

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
        # if self.save_config:
        #     # Save map_1 config to a file
        #     # config_str = json.dumps(map_1.config)
        #     # if type(config) == str:
        #     #     config = config.encode("utf-8")
        #     with open("kepler_config.json", "w") as f:
        #         f.write(json.dumps(map_1.config))

        # if self.load_config:
        #     with open("kepler_config.json", "r") as f:
        #         config = json.loads(f.read())
        map_1.config = config

        html = map_1._repr_html_()
        html = html.decode("utf-8")

        # replace css and JavaScript paths
        html = replace_external_js_css_paths(
            r"\1./libs/kepler/2.5.5/\3\"\4",
            html,
        )
        # replace any stylesheet links that are dynamically created
        html = replace_external_js_css_paths(
            r"\1./libs/kepler/2.5.5/\2\3",
            html,
            """(createElement\("link",\{rel:"stylesheet",href:")[^"']*\/([^"']*)("\}\))""",
        )

        # remove all Google Analytics scripts to prevent sending any information
        html = replace_external_js_css_paths(
            "", html, "<script>[^<]*www\.google-analytics\.com[^<]*<\/script>"
        )
        html = replace_external_js_css_paths(
            "",
            html,
            """s\.a\.createElement\("script"[^)]*googletagmanager\.com[^\)]*\),""",
        )
        html = replace_external_js_css_paths(
            "",
            html,
            """s\.a\.createElement\("script",null,"[^"]*gtag\([^"]*"\)""",
        )

        return knext.view_html(html)


############################################
# Spatial heatmap node
############################################


@knext.node(
    name="Spatial Heatmap",
    node_type=knext.NodeType.VISUALIZER,
    icon_path=__NODE_ICON_PATH + "SpatialHeatmap.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Table to Visualize",
    description="Table with data to visualize",
)
@knext.output_view(
    name="Heatmap View",
    description="Showing a heatmap with the data",
    static_resources="libs/leaflet/1.6.0",
)
class ViewNodeHeatmap:
    """This node will visualize the given data on a heatmap.

    This node will visualize the given data as interactive heatmap. The selected weight column defines
    the "heat" of each data point which is visualized on a world map.

    Please find more information about the heatmap [here](https://www.gislounge.com/heat-maps-in-gis/).
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
        "The minimum opacity the lowest value in the heatmap will have.",
        default_value=0.5,
    )

    max_zoom = knext.IntParameter(
        "Maximum zoom",
        """Zoom level where the points reach maximum intensity (as intensity scales with zoom), 
        equals maxZoom of the map by default""",
        default_value=18,
    )

    radius = knext.IntParameter(
        "Radius",
        "Radius of each datapoint of the heatmap",
        default_value=25,
    )

    blur = knext.IntParameter(
        "Blur",
        """The blur factor that will be applied to all datapoints. 
        The higher the blur factor is, the smoother the gradients will be.""",
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

        # replace css and JavaScript paths
        html = map.get_root().render()
        html = replace_external_js_css_paths(
            r"\1./libs/leaflet/1.6.0/\3\"\4",
            html,
        )

        return knext.view(html)
