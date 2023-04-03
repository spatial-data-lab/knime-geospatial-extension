import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


category = knext.category(
    path="/community/geo",
    level_id="viz",
    name="Spatial Visualization",  # Spatial Visualization
    description="Spatial view nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/VisualizationCategory.png",
    after="conversion",
)
# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/Visualization/"


def replace_external_js_css_paths(
    replacement: str,
    html: str,
    regex_patter: str = '(<script src="|<link( rel="stylesheet")? href=")https?[^"]*\/([^"]*)"( rel="stylesheet")?',
) -> str:
    """
    Uses a regular expression to find all script and stylesheet tags in a given HTML page.
    The first matching group is either the script or stylesheet part up until the opening " of the
    URL. The second matching group is the file name.
    The method will also add the closing ".

    For example if the HTML code is <script src="https://cdn.jsdelivr.net/npm/leaflet@1.9.3/dist/leaflet.js">
    the first group is <script src=" and the second group is leaflet.js so using the following replacement
    r'\1./libs/leaflet/1.9.3/\3"\4' will lead to this URL: <script src="./libs/leaflet/1.9.3/leaflet.js">.
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
    Group of settings that define the coloring of the geometric objects. The color column can be either nominal or numerical.
    If a numerical column is selected you might want to enable the classification of the numeric values to group
    them into bins prior to assigning a color to each bin using the color map information.
    If a nominal column is selected the color map will be used to assign a color to each unique value in the column.
    Noticed that if a nominal column is selected the classification settings will be ignored.
    """

    color_col = knext.ColumnParameter(
        "Marker color column",
        """Select the marker color column to be plotted. If numerical you might want to adapt the classification 
        settings accordingly.""",
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=True,
    )

    color_map = knext.StringParameter(
        "Color map",
        """Select the color map to use for the color column. `xxx_r` mean the reverse of the `xxx` color map. 
        See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html)""",
        default_value="viridis",
        enum=[
            "afmhot",
            "afmhot_r",
            "autumn",
            "autumn_r",
            "binary",
            "binary_r",
            "Blues",
            "Blues_r",
            "bone",
            "bone_r",
            "BrBG",
            "BrBG_r",
            "brg",
            "BuGn",
            "BuGn_r",
            "BuPu",
            "BuPu_r",
            "bwr",
            "bwr_r",
            "cividis",
            "cividis_r",
            "CMRmap",
            "cool",
            "cool_r",
            "coolwarm",
            "coolwarm_r",
            "copper",
            "copper_r",
            "cubehelix",
            "flag",
            "gist_earth",
            "gist_gray",
            "gist_gray_r",
            "gist_heat",
            "gist_heat_r",
            "gist_ncar",
            "gist_rainbow",
            "gist_stern",
            "gist_yarg",
            "gist_yarg_r",
            "GnBu",
            "GnBu_r",
            "gnuplot",
            "gnuplot2",
            "gray",
            "gray_r",
            "Greens",
            "Greens_r",
            "Greys",
            "Greys_r",
            "hot",
            "hot_r",
            "hsv",
            "inferno",
            "inferno_r",
            "jet",
            "magma",
            "magma_r",
            "nipy_spectral",
            "ocean",
            "Oranges",
            "Oranges_r",
            "OrRd",
            "OrRd_r",
            "pink",
            "pink_r",
            "PiYG",
            "PiYG_r",
            "plasma",
            "plasma_r",
            "PRGn",
            "PRGn_r",
            "prism",
            "PuBu",
            "PuBu_r",
            "PuBuGn",
            "PuBuGn_r",
            "PuOr",
            "PuOr_r",
            "PuRd",
            "PuRd_r",
            "Purples",
            "Purples_r",
            "rainbow",
            "RdBu",
            "RdBu_r",
            "RdGy",
            "RdGy_r",
            "RdPu",
            "RdPu_r",
            "RdYlBu",
            "RdYlBu_r",
            "RdYlGn",
            "RdYlGn_r",
            "Reds",
            "Reds_r",
            "seismic",
            "seismic_r",
            "Spectral",
            "Spectral_r",
            "spring",
            "spring_r",
            "summer",
            "summer_r",
            "terrain",
            "turbo",
            "twilight",
            "twilight_shifted",
            "viridis",
            "viridis_r",
            "winter",
            "winter_r",
            "Wistia",
            "Wistia_r",
            "YlGn",
            "YlGn_r",
            "YlGnBu",
            "YlGnBu_r",
            "YlOrBr",
            "YlOrBr_r",
            "YlOrRd",
            "YlOrRd_r",
        ],
    )

    use_classify = knext.BoolParameter(
        "Classify numerical marker color columns",
        """If checked, a numerical marker color column will be classified using the selected classification method. 
        The 'Number of classes' will be used to determine the number of bins.""",
        default_value=False,
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


@knext.parameter_group(label="Color Legend Settings")
class LegendSettings:
    """
    Group of settings that define if a color legend is shown on the map and if so how it should be formatted.
    The color legend is only shown if you have selected a color column.
    """

    plot = knext.BoolParameter(
        "Show color legend",
        "If checked, the color legend will be shown in the plot.",
        default_value=False,
    )

    caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the color legend. By default, the caption is the name of the selected color column or "
        + "empty for heat map.",
        default_value="",
    )


@knext.parameter_group(label="Size Settings")
class SizeSettings:
    """
    Group of settings that define the size of the geometric objects. The size column should be numerical.
    The size is fixed by default. If the `Marker size column` is selected, the `Marker size scale` option will
    be ignored, and size will be scaled by the values of the column. For point features, the size is the radius
    of the circle. For line features, the size is the width of the line. For polygon features, the size is the
    radius of the centroid of the polygon.
    """

    size_col = knext.ColumnParameter(
        "Marker size column",
        """Select marker size column. 
        """,
        column_filter=knut.is_numeric,
        include_none_column=True,
    )

    size_scale = knext.DoubleParameter(
        "Marker size scale",
        """Select the size scale of the markers. 
        If the Marker size column is selected, this option will be ignored.
        Noticed that the size scale only works for point features.
        """,
        default_value=float(1.0),
        min_value=None,
        max_value=None,
    )


@knext.parameter_group(label="Base Map Setting")
class BaseMapSettings:
    """Base map setting for the visualization."""

    base_map = knext.StringParameter(
        "Base map",
        """Select the base map to use for the visualization. If choose `Don't show base map`, the base map will be hidden.
        The default base map is `OpenStreetMap`.
        See [Folium base maps](https://python-visualization.github.io/folium/quickstart.html#Tiles).""",
        default_value="OpenStreetMap",
        enum=[
            "CartoDB DarkMatter",
            "CartoDB DarkMatterNoLabels",
            "CartoDB DarkMatterOnlyLabels",
            "CartoDB Positron",
            "CartoDB PositronNoLabels",
            "CartoDB PositronOnlyLabels",
            "CartoDB Voyager",
            "CartoDB VoyagerLabelsUnder",
            "CartoDB VoyagerNoLabels",
            "CartoDB VoyagerOnlyLabels",
            "Esri DeLorme",
            "Esri NatGeoWorldMap",
            "Esri OceanBasemap",
            "Esri WorldGrayCanvas",
            "Esri WorldImagery",
            "Esri WorldPhysical",
            "Esri WorldShadedRelief",
            "Esri WorldStreetMap",
            "Esri WorldTerrain",
            "Esri WorldTopoMap",
            "Gaode Normal",
            "Gaode Satellite",
            "NASAGIBS ASTER_GDEM_Greyscale_Shaded_Relief",
            "NASAGIBS BlueMarble",
            "NASAGIBS BlueMarble3031",
            "NASAGIBS BlueMarble3413",
            "NASAGIBS ModisAquaBands721CR",
            "NASAGIBS ModisAquaTrueColorCR",
            "NASAGIBS ModisTerraAOD",
            "NASAGIBS ModisTerraBands367CR",
            "NASAGIBS ModisTerraBands721CR",
            "NASAGIBS ModisTerraChlorophyll",
            "NASAGIBS ModisTerraLSTDay",
            "NASAGIBS ModisTerraSnowCover",
            "NASAGIBS ModisTerraTrueColorCR",
            "NASAGIBS ViirsEarthAtNight2012",
            "NASAGIBS ViirsTrueColorCR",
            "OpenRailwayMap",
            "OpenStreetMap",
            "Stamen Terrain",
            "Stamen TerrainBackground",
            "Stamen TerrainLabels",
            "Stamen Toner",
            "Stamen TonerBackground",
            "Stamen TonerHybrid",
            "Stamen TonerLabels",
            "Stamen TonerLines",
            "Stamen TonerLite",
            "Stamen TopOSMFeatures",
            "Stamen TopOSMRelief",
            "Stamen Watercolor",
            "Stamen Watercolor",
            "Strava All",
            "Strava Ride",
            "Strava Run",
            "Strava Water",
            "Strava Winter",
            "Don't show base map",
        ],
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
    description="Showing an interactive map with the geospatial data",
    static_resources="libs/leaflet/1.9.3",
)
class ViewNode:
    """Creates an interactive map view based on the selected geometric elements of the input table.
    This node creates an interactive map view based on the selected geometric elements of the input table.
    It provides various dialog options to modify the appearance of the view e.g. the base map, shape color and
    size.
    The geometric elements are drawn in the order they appear in the input table. For example, if you want to show
    points within a polygon you want to have the points drawn last on top of the polygon. To do so sort the input
    table to have polygons as first rows followed by the points.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        # "geometry",
        column_filter=knut.is_geo,  # Allows all geo columns
        include_row_key=False,
        include_none_column=False,  # must contains a geometry column
    )

    stroke = knext.BoolParameter(
        "Stroke",
        "Whether to draw strokes along the polygon boundary.",
        default_value=True,
    )

    name_cols = knext.MultiColumnParameter(
        "Marker tooltip columns",
        "Select columns which should be shown in the marker tooltips.",
        column_filter=knut.negate(knut.is_geo),  # Filter out all geo columns
    )

    popup_cols = knext.MultiColumnParameter(
        "Marker popup columns",
        "Select columns which should be shown in the marker popups.",
        column_filter=knut.negate(knut.is_geo),  # Filter out all geo columns
    )

    size_settings = SizeSettings()

    color_settings = ColorSettings()

    basemap_setting = BaseMapSettings()

    legend_settings = LegendSettings()

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        # if self.name_cols is None:
        #     self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # keep only the selected columns
        selected_col_names = {self.geo_col}
        if self.name_cols is not None:
            selected_col_names.update(self.name_cols)
        if self.popup_cols is not None:
            selected_col_names.update(self.popup_cols)
        if "none" not in str(self.color_settings.color_col).lower():
            selected_col_names.add(self.color_settings.color_col)
        if "none" not in str(self.size_settings.size_col).lower():
            selected_col_names.add(self.size_settings.size_col)
        filtered_table = input_table[list(selected_col_names)]

        gdf = gp.GeoDataFrame(filtered_table.to_pandas(), geometry=self.geo_col)
        # convert all remaining none string_numeric columns to string
        schema = filtered_table.schema
        for c in schema:
            if not knut.is_numeric_or_string(c) and not knut.is_geo(c):
                gdf[c.name] = gdf[c.name].apply(str)
        if self.basemap_setting.base_map == "Don't show base map":
            base_map = None
        else:
            base_map = self.basemap_setting.base_map
        kws = {
            # "column":self.color_col,
            # "cmap":self.color_map,
            "tooltip": self.name_cols,
            "tiles": base_map,
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

        if "none" not in str(self.size_settings.size_col).lower():
            max_pop_est = gdf[self.size_settings.size_col].max()
            min_pop_est = gdf[self.size_settings.size_col].min()

            # check whether is line
            geo_types = gdf[self.geo_col].geom_type.unique()
            if ("LineString" in geo_types) or ("MultiLineString" in geo_types):
                max_size = 8
                kws["style_kwds"] = {
                    "style_function": lambda x: {
                        "weight": (
                            x["properties"][self.size_settings.size_col] - min_pop_est
                        )
                        / (max_pop_est - min_pop_est)
                        * max_size
                    }
                }
            elif ("Polygon" in geo_types) or ("MultiPolygon" in geo_types):
                max_size = 30
                kws["style_kwds"] = {
                    "style_function": lambda x: {
                        "radius": (
                            x["properties"][self.size_settings.size_col] - min_pop_est
                        )
                        / (max_pop_est - min_pop_est)
                        * max_size
                    }
                }
                kws["m"] = gdf.explore(tiles=base_map)
                gdf[self.geo_col] = gdf.centroid

            else:
                max_size = 30
                kws["style_kwds"] = {
                    "style_function": lambda x: {
                        "radius": (
                            x["properties"][self.size_settings.size_col] - min_pop_est
                        )
                        / (max_pop_est - min_pop_est)
                        * max_size
                    }
                }
        elif "none" not in str(self.size_settings.size_scale).lower():
            geo_types = gdf[self.geo_col].geom_type.unique()
            if ("LineString" in geo_types) or ("MultiLineString" in geo_types):
                kws["style_kwds"] = {"weight": self.size_settings.size_scale}
            elif ("Polygon" in geo_types) or ("MultiPolygon" in geo_types):
                kws["style_kwds"] = {"weight": self.size_settings.size_scale}
            else:
                kws["style_kwds"] = {"radius": self.size_settings.size_scale}
        if not self.stroke:
            if "style_kwds" in kws:
                kws["style_kwds"]["stroke"] = False
            else:
                kws["style_kwds"] = {"stroke": False}
        map = gdf.explore(**kws)

        # replace css and JavaScript paths
        html = map.get_root().render()
        html = replace_external_js_css_paths(
            r'\1./libs/leaflet/1.9.3/\3"\4',
            html,
        )

        return knext.view(html)


# geo view static


@knext.parameter_group(label="Coloring Settings")
class StaticColorSettings:
    """
    Group of settings that define the coloring of the geometric objects. The color column can be either nominal or numerical.
    If a numerical column is selected you might want to enable the classification of the numeric values to group
    them into bins prior to assigning a color to each bin using the color map information.
    If a nominal column is selected the color map will be used to assign a color to each unique value in the column.
    Noticed that if a nominal column is selected the classification settings will be ignored.
    """

    color_col = knext.ColumnParameter(
        "Marker color column",
        """Select marker color column to be plotted. If numerical you might want to adapt the classification 
        settings accordingly. If nominal the color map will be used to assign a color to each unique value in the column.
        Noticed that if a nominal column is selected the classification settings will be ignored.
        Select 'None' if you want a unified marker color.""",
        column_filter=knut.is_numeric_or_string,
        include_row_key=False,
        include_none_column=True,
    )

    color = knext.StringParameter(
        "Marker color",
        """Select marker color. It will assign a unified color for all features. If the a `Marker color column` 
        is selected and not `None` option, this option will be ignored.
        Select none if you don't want to set a unified marker color.""",
        default_value="none",
        enum=[
            "beige",
            "black",
            "blue",
            "cadetblue",
            "darkblue",
            "darkgreen",
            "darkpurple",
            "darkred",
            "gray",
            "green",
            "lightblue",
            "lightgray",
            "lightgreen",
            "lightred",
            "none",
            "orange",
            "pink",
            "purple",
            "red",
            "white",
        ],
    )

    color_map = knext.StringParameter(
        "Color map",
        """Select the color map to use for the color column. 
        See [Colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html).""",
        default_value="viridis",
        enum=[
            "cividis",
            "inferno",
            "magma",
            "plasma",
            "viridis",
        ],
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
            "beige",
            "black",
            "black",
            "blue",
            "cadetblue",
            "darkblue",
            "darkgreen",
            "darkpurple",
            "darkred",
            "gray",
            "green",
            "lightblue",
            "lightgray",
            "lightgreen",
            "lightred",
            "none",
            "orange",
            "pink",
            "purple",
            "red",
            "white",
        ],
    )


@knext.parameter_group(label="Color Legend Settings")
class StaticLegendSettings:
    """
    Group of settings that define if a color legend is shown on the map and if so how it should be formatted.
    The color legend is only shown if you have selected a color column.

    More details about the legend settings can be found [matplotlib.pyplot.legend](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html#matplotlib.pyplot.legend)
    and [matplotlib.pyplot.colorbar](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html#matplotlib.pyplot.colorbar).
    """

    plot = knext.BoolParameter(
        "Show color legend",
        "If checked, the color legend will be shown in the plot.",
        default_value=True,
    )

    caption = knext.StringParameter(
        "Legend caption",
        "Set the caption for the legend. By default, the caption is the name of the selected color column or "
        + "empty for heat map.",
        default_value="",
        # default_value=color_col,
    )

    caption_fontsize = knext.IntParameter(
        "Caption font size",
        "Set the font size for the legend caption.",
        default_value=10,
        min_value=1,
        max_value=100,
    )

    expand = knext.BoolParameter(
        "Expand legend",
        "If checked, the legend will be horizontally expanded to fill the axes area.",
        default_value=False,
    )

    location = knext.StringParameter(
        "Legend location",
        "Select the location for the legend.",
        default_value="lower right",
        enum=[
            "best",
            "center left",
            "center right",
            "center",
            "lower center",
            "lower left",
            "lower right",
            "outside_bottom",
            "outside_top",
            "right",
            "upper center",
            "upper left",
            "upper right",
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
        enum=["black", "blue", "green", "orange", "purple", "red", "white", "yellow"],
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
        "Color bar legend shrink",
        "Select the shrinking value for the color bar legend. Only works for color bar.",
        default_value=1.0,
        min_value=0.0,
        max_value=1.0,
    )

    colorbar_pad = knext.DoubleParameter(
        "Color bar legend pad",
        "Select the pad value for the color bar legend. Only works for color bar",
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
    This node will visualize the given geometric elements using the [matplotlib](https://matplotlib.org/stable/).
    It can be used to create [Choropleth Maps](https://en.wikipedia.org/wiki/Choropleth_map) by assigning
    a marker color column. The node further supports various settings to adapt the legend to your needs.
    For more information on the settings, please refer to the [geopandas.GeoDataFrame.plot](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.plot.html).
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    line_width_col = knext.ColumnParameter(
        "Line width column",
        """Select line width column. The width is fixed by default. If a line width column is selected, the width will 
        be scaled by the values of the column.""",
        column_filter=knut.is_numeric,
        include_none_column=True,
    )

    line_width = knext.IntParameter(
        "Line width",
        "Select a unified line width, can be set to none.",
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

    set_axis_off = knext.BoolParameter(
        "Set axis off",
        "If checked, the axis will be set off.",
        default_value=False,
    )

    size_settings = SizeSettings()

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

        # set color
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
            if type(gdf[self.color_settings.color_col].values[0]) == str:
                kws["legend_kwds"] = {
                    # "fmt": "{:.0f}",
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
            elif self.color_settings.use_classify:
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

        # set size
        geo_types = gdf[self.geo_col].geom_type.unique()
        if not (("LineString" in geo_types) or ("MultiLineString" in geo_types)):
            if "none" not in str(self.size_settings.size_col).lower():
                max_point_size = 2000
                max_val = gdf[self.size_settings.size_col].max()
                min_val = gdf[self.size_settings.size_col].min()
                normal_base = (gdf[self.size_settings.size_col] - min_val) / max_val
                kws["markersize"] = normal_base * max_point_size
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
    creates an animation for a given time series column. For more information about the supported interactions
    see the [kepler.gl user guides](https://docs.kepler.gl/docs/user-guides).

    This node uses the [Mapbox GL JS API](https://www.mapbox.com/pricing#map-loads-for-web) which for commercial
    usage might require an [access token](https://docs.mapbox.com/help/glossary/access-token/).
    If you want to use a different base map, you can configure it inside the interactive
    view with Kepler.gl's UI. You can also configure the
    [Mapbox style](https://docs.kepler.gl/docs/user-guides/f-map-styles#custom-map-styles) you want to use and
    the access token.


    By default, it takes all column information that is included inside the input table.
    If you want to limit the amount of information sent to the node view you can use one of the
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

        # convert all none string_numeric and geometry columns to string
        schema = input_table.schema
        for c in schema:
            if not knut.is_numeric_or_string(c) and not knut.is_geo(c):
                gdf[c.name] = gdf[c.name].apply(str)

        from keplergl import KeplerGl

        map_1 = KeplerGl(show_docs=False)
        map_1.add_data(data=gdf.copy(), name="state")
        # config = {}

        # import json
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
        # map_1.config = config

        html = map_1._repr_html_(center_map=True)
        html = html.decode("utf-8")

        # replace css and JavaScript paths
        html = replace_external_js_css_paths(
            r'\1./libs/kepler/2.5.5/\3"\4',
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
    static_resources="libs/leaflet/1.9.3",
)
class ViewNodeHeatmap:
    """This node will visualize the given data on a heatmap.

    This node will visualize the given data as an interactive heatmap. The selected weight column defines
    the "heat" of each data point which is visualized on a world map.

    The geometric elements are drawn in the order they appear in the input table. For example, if you want to show
    points within a polygon you want to have the points drawn last on top of the polygon. To do so sort the input
    table to have polygons as first rows followed by the points.

    Please make sure the input table does not contain any rows with missing values.

    Find more information about the supported options in [folium.plugins.HeatMap](https://python-visualization.github.io/folium/plugins.html#folium.plugins.HeatMap).

    Find more information about the heatmap algorithm [here](https://www.gislounge.com/heat-maps-in-gis/).
    """

    import branca.colormap

    _color_bars = {
        "Accent_03": branca.colormap.linear.Accent_03,
        "Accent_04": branca.colormap.linear.Accent_04,
        "Accent_05": branca.colormap.linear.Accent_05,
        "Accent_06": branca.colormap.linear.Accent_06,
        "Accent_07": branca.colormap.linear.Accent_07,
        "Accent_08": branca.colormap.linear.Accent_08,
        "Blues_03": branca.colormap.linear.Blues_03,
        "Blues_04": branca.colormap.linear.Blues_04,
        "Blues_05": branca.colormap.linear.Blues_05,
        "Blues_06": branca.colormap.linear.Blues_06,
        "Blues_07": branca.colormap.linear.Blues_07,
        "Blues_08": branca.colormap.linear.Blues_08,
        "Blues_09": branca.colormap.linear.Blues_09,
        "BrBG_03": branca.colormap.linear.BrBG_03,
        "BrBG_04": branca.colormap.linear.BrBG_04,
        "BrBG_05": branca.colormap.linear.BrBG_05,
        "BrBG_06": branca.colormap.linear.BrBG_06,
        "BrBG_07": branca.colormap.linear.BrBG_07,
        "BrBG_08": branca.colormap.linear.BrBG_08,
        "BrBG_09": branca.colormap.linear.BrBG_09,
        "BrBG_10": branca.colormap.linear.BrBG_10,
        "BrBG_11": branca.colormap.linear.BrBG_11,
        "BuGn_03": branca.colormap.linear.BuGn_03,
        "BuGn_04": branca.colormap.linear.BuGn_04,
        "BuGn_05": branca.colormap.linear.BuGn_05,
        "BuGn_06": branca.colormap.linear.BuGn_06,
        "BuGn_07": branca.colormap.linear.BuGn_07,
        "BuGn_08": branca.colormap.linear.BuGn_08,
        "BuGn_09": branca.colormap.linear.BuGn_09,
        "BuPu_03": branca.colormap.linear.BuPu_03,
        "BuPu_04": branca.colormap.linear.BuPu_04,
        "BuPu_05": branca.colormap.linear.BuPu_05,
        "BuPu_06": branca.colormap.linear.BuPu_06,
        "BuPu_07": branca.colormap.linear.BuPu_07,
        "BuPu_08": branca.colormap.linear.BuPu_08,
        "BuPu_09": branca.colormap.linear.BuPu_09,
        "Dark2_03": branca.colormap.linear.Dark2_03,
        "Dark2_04": branca.colormap.linear.Dark2_04,
        "Dark2_05": branca.colormap.linear.Dark2_05,
        "Dark2_06": branca.colormap.linear.Dark2_06,
        "Dark2_07": branca.colormap.linear.Dark2_07,
        "Dark2_08": branca.colormap.linear.Dark2_08,
        "GnBu_03": branca.colormap.linear.GnBu_03,
        "GnBu_04": branca.colormap.linear.GnBu_04,
        "GnBu_05": branca.colormap.linear.GnBu_05,
        "GnBu_06": branca.colormap.linear.GnBu_06,
        "GnBu_07": branca.colormap.linear.GnBu_07,
        "GnBu_08": branca.colormap.linear.GnBu_08,
        "GnBu_09": branca.colormap.linear.GnBu_09,
        "Greens_03": branca.colormap.linear.Greens_03,
        "Greens_04": branca.colormap.linear.Greens_04,
        "Greens_05": branca.colormap.linear.Greens_05,
        "Greens_06": branca.colormap.linear.Greens_06,
        "Greens_07": branca.colormap.linear.Greens_07,
        "Greens_08": branca.colormap.linear.Greens_08,
        "Greens_09": branca.colormap.linear.Greens_09,
        "Greys_03": branca.colormap.linear.Greys_03,
        "Greys_04": branca.colormap.linear.Greys_04,
        "Greys_05": branca.colormap.linear.Greys_05,
        "Greys_06": branca.colormap.linear.Greys_06,
        "Greys_07": branca.colormap.linear.Greys_07,
        "Greys_08": branca.colormap.linear.Greys_08,
        "Greys_09": branca.colormap.linear.Greys_09,
        "Oranges_03": branca.colormap.linear.Oranges_03,
        "Oranges_04": branca.colormap.linear.Oranges_04,
        "Oranges_05": branca.colormap.linear.Oranges_05,
        "Oranges_06": branca.colormap.linear.Oranges_06,
        "Oranges_07": branca.colormap.linear.Oranges_07,
        "Oranges_08": branca.colormap.linear.Oranges_08,
        "Oranges_09": branca.colormap.linear.Oranges_09,
        "OrRd_03": branca.colormap.linear.OrRd_03,
        "OrRd_04": branca.colormap.linear.OrRd_04,
        "OrRd_05": branca.colormap.linear.OrRd_05,
        "OrRd_06": branca.colormap.linear.OrRd_06,
        "OrRd_07": branca.colormap.linear.OrRd_07,
        "OrRd_08": branca.colormap.linear.OrRd_08,
        "OrRd_09": branca.colormap.linear.OrRd_09,
        "Paired_03": branca.colormap.linear.Paired_03,
        "Paired_04": branca.colormap.linear.Paired_04,
        "Paired_05": branca.colormap.linear.Paired_05,
        "Paired_06": branca.colormap.linear.Paired_06,
        "Paired_07": branca.colormap.linear.Paired_07,
        "Paired_08": branca.colormap.linear.Paired_08,
        "Paired_09": branca.colormap.linear.Paired_09,
        "Paired_10": branca.colormap.linear.Paired_10,
        "Paired_11": branca.colormap.linear.Paired_11,
        "Paired_12": branca.colormap.linear.Paired_12,
        "Pastel1_03": branca.colormap.linear.Pastel1_03,
        "Pastel1_04": branca.colormap.linear.Pastel1_04,
        "Pastel1_05": branca.colormap.linear.Pastel1_05,
        "Pastel1_06": branca.colormap.linear.Pastel1_06,
        "Pastel1_07": branca.colormap.linear.Pastel1_07,
        "Pastel1_08": branca.colormap.linear.Pastel1_08,
        "Pastel1_09": branca.colormap.linear.Pastel1_09,
        "Pastel2_03": branca.colormap.linear.Pastel2_03,
        "Pastel2_04": branca.colormap.linear.Pastel2_04,
        "Pastel2_05": branca.colormap.linear.Pastel2_05,
        "Pastel2_06": branca.colormap.linear.Pastel2_06,
        "Pastel2_07": branca.colormap.linear.Pastel2_07,
        "Pastel2_08": branca.colormap.linear.Pastel2_08,
        "PiYG_03": branca.colormap.linear.PiYG_03,
        "PiYG_04": branca.colormap.linear.PiYG_04,
        "PiYG_05": branca.colormap.linear.PiYG_05,
        "PiYG_06": branca.colormap.linear.PiYG_06,
        "PiYG_07": branca.colormap.linear.PiYG_07,
        "PiYG_08": branca.colormap.linear.PiYG_08,
        "PiYG_09": branca.colormap.linear.PiYG_09,
        "PiYG_10": branca.colormap.linear.PiYG_10,
        "PiYG_11": branca.colormap.linear.PiYG_11,
        "PRGn_03": branca.colormap.linear.PRGn_03,
        "PRGn_04": branca.colormap.linear.PRGn_04,
        "PRGn_05": branca.colormap.linear.PRGn_05,
        "PRGn_06": branca.colormap.linear.PRGn_06,
        "PRGn_07": branca.colormap.linear.PRGn_07,
        "PRGn_08": branca.colormap.linear.PRGn_08,
        "PRGn_09": branca.colormap.linear.PRGn_09,
        "PRGn_10": branca.colormap.linear.PRGn_10,
        "PRGn_11": branca.colormap.linear.PRGn_11,
        "PuBu_03": branca.colormap.linear.PuBu_03,
        "PuBu_04": branca.colormap.linear.PuBu_04,
        "PuBu_05": branca.colormap.linear.PuBu_05,
        "PuBu_06": branca.colormap.linear.PuBu_06,
        "PuBu_07": branca.colormap.linear.PuBu_07,
        "PuBu_08": branca.colormap.linear.PuBu_08,
        "PuBu_09": branca.colormap.linear.PuBu_09,
        "PuBuGn_03": branca.colormap.linear.PuBuGn_03,
        "PuBuGn_04": branca.colormap.linear.PuBuGn_04,
        "PuBuGn_05": branca.colormap.linear.PuBuGn_05,
        "PuBuGn_06": branca.colormap.linear.PuBuGn_06,
        "PuBuGn_07": branca.colormap.linear.PuBuGn_07,
        "PuBuGn_08": branca.colormap.linear.PuBuGn_08,
        "PuBuGn_09": branca.colormap.linear.PuBuGn_09,
        "PuOr_03": branca.colormap.linear.PuOr_03,
        "PuOr_04": branca.colormap.linear.PuOr_04,
        "PuOr_05": branca.colormap.linear.PuOr_05,
        "PuOr_06": branca.colormap.linear.PuOr_06,
        "PuOr_07": branca.colormap.linear.PuOr_07,
        "PuOr_08": branca.colormap.linear.PuOr_08,
        "PuOr_09": branca.colormap.linear.PuOr_09,
        "PuOr_10": branca.colormap.linear.PuOr_10,
        "PuOr_11": branca.colormap.linear.PuOr_11,
        "PuRd_03": branca.colormap.linear.PuRd_03,
        "PuRd_04": branca.colormap.linear.PuRd_04,
        "PuRd_05": branca.colormap.linear.PuRd_05,
        "PuRd_06": branca.colormap.linear.PuRd_06,
        "PuRd_07": branca.colormap.linear.PuRd_07,
        "PuRd_08": branca.colormap.linear.PuRd_08,
        "PuRd_09": branca.colormap.linear.PuRd_09,
        "Purples_03": branca.colormap.linear.Purples_03,
        "Purples_04": branca.colormap.linear.Purples_04,
        "Purples_05": branca.colormap.linear.Purples_05,
        "Purples_06": branca.colormap.linear.Purples_06,
        "Purples_07": branca.colormap.linear.Purples_07,
        "Purples_08": branca.colormap.linear.Purples_08,
        "Purples_09": branca.colormap.linear.Purples_09,
        "RdBu_03": branca.colormap.linear.RdBu_03,
        "RdBu_04": branca.colormap.linear.RdBu_04,
        "RdBu_05": branca.colormap.linear.RdBu_05,
        "RdBu_06": branca.colormap.linear.RdBu_06,
        "RdBu_07": branca.colormap.linear.RdBu_07,
        "RdBu_08": branca.colormap.linear.RdBu_08,
        "RdBu_09": branca.colormap.linear.RdBu_09,
        "RdBu_10": branca.colormap.linear.RdBu_10,
        "RdBu_11": branca.colormap.linear.RdBu_11,
        "RdGy_03": branca.colormap.linear.RdGy_03,
        "RdGy_04": branca.colormap.linear.RdGy_04,
        "RdGy_05": branca.colormap.linear.RdGy_05,
        "RdGy_06": branca.colormap.linear.RdGy_06,
        "RdGy_07": branca.colormap.linear.RdGy_07,
        "RdGy_08": branca.colormap.linear.RdGy_08,
        "RdGy_09": branca.colormap.linear.RdGy_09,
        "RdGy_10": branca.colormap.linear.RdGy_10,
        "RdGy_11": branca.colormap.linear.RdGy_11,
        "RdPu_03": branca.colormap.linear.RdPu_03,
        "RdPu_04": branca.colormap.linear.RdPu_04,
        "RdPu_05": branca.colormap.linear.RdPu_05,
        "RdPu_06": branca.colormap.linear.RdPu_06,
        "RdPu_07": branca.colormap.linear.RdPu_07,
        "RdPu_08": branca.colormap.linear.RdPu_08,
        "RdPu_09": branca.colormap.linear.RdPu_09,
        "RdYlBu_03": branca.colormap.linear.RdYlBu_03,
        "RdYlBu_04": branca.colormap.linear.RdYlBu_04,
        "RdYlBu_05": branca.colormap.linear.RdYlBu_05,
        "RdYlBu_06": branca.colormap.linear.RdYlBu_06,
        "RdYlBu_07": branca.colormap.linear.RdYlBu_07,
        "RdYlBu_08": branca.colormap.linear.RdYlBu_08,
        "RdYlBu_09": branca.colormap.linear.RdYlBu_09,
        "RdYlBu_10": branca.colormap.linear.RdYlBu_10,
        "RdYlBu_11": branca.colormap.linear.RdYlBu_11,
        "RdYlGn_03": branca.colormap.linear.RdYlGn_03,
        "RdYlGn_04": branca.colormap.linear.RdYlGn_04,
        "RdYlGn_05": branca.colormap.linear.RdYlGn_05,
        "RdYlGn_06": branca.colormap.linear.RdYlGn_06,
        "RdYlGn_07": branca.colormap.linear.RdYlGn_07,
        "RdYlGn_08": branca.colormap.linear.RdYlGn_08,
        "RdYlGn_09": branca.colormap.linear.RdYlGn_09,
        "RdYlGn_10": branca.colormap.linear.RdYlGn_10,
        "RdYlGn_11": branca.colormap.linear.RdYlGn_11,
        "Reds_03": branca.colormap.linear.Reds_03,
        "Reds_04": branca.colormap.linear.Reds_04,
        "Reds_05": branca.colormap.linear.Reds_05,
        "Reds_06": branca.colormap.linear.Reds_06,
        "Reds_07": branca.colormap.linear.Reds_07,
        "Reds_08": branca.colormap.linear.Reds_08,
        "Reds_09": branca.colormap.linear.Reds_09,
        "Set1_03": branca.colormap.linear.Set1_03,
        "Set1_04": branca.colormap.linear.Set1_04,
        "Set1_05": branca.colormap.linear.Set1_05,
        "Set1_06": branca.colormap.linear.Set1_06,
        "Set1_07": branca.colormap.linear.Set1_07,
        "Set1_08": branca.colormap.linear.Set1_08,
        "Set1_09": branca.colormap.linear.Set1_09,
        "Set2_03": branca.colormap.linear.Set2_03,
        "Set2_04": branca.colormap.linear.Set2_04,
        "Set2_05": branca.colormap.linear.Set2_05,
        "Set2_06": branca.colormap.linear.Set2_06,
        "Set2_07": branca.colormap.linear.Set2_07,
        "Set2_08": branca.colormap.linear.Set2_08,
        "Set3_03": branca.colormap.linear.Set3_03,
        "Set3_04": branca.colormap.linear.Set3_04,
        "Set3_05": branca.colormap.linear.Set3_05,
        "Set3_06": branca.colormap.linear.Set3_06,
        "Set3_07": branca.colormap.linear.Set3_07,
        "Set3_08": branca.colormap.linear.Set3_08,
        "Set3_09": branca.colormap.linear.Set3_09,
        "Set3_10": branca.colormap.linear.Set3_10,
        "Set3_11": branca.colormap.linear.Set3_11,
        "Set3_12": branca.colormap.linear.Set3_12,
        "Spectral_03": branca.colormap.linear.Spectral_03,
        "Spectral_04": branca.colormap.linear.Spectral_04,
        "Spectral_05": branca.colormap.linear.Spectral_05,
        "Spectral_06": branca.colormap.linear.Spectral_06,
        "Spectral_07": branca.colormap.linear.Spectral_07,
        "Spectral_08": branca.colormap.linear.Spectral_08,
        "Spectral_09": branca.colormap.linear.Spectral_09,
        "Spectral_10": branca.colormap.linear.Spectral_10,
        "Spectral_11": branca.colormap.linear.Spectral_11,
        "viridis": branca.colormap.linear.viridis,
        "YlGn_03": branca.colormap.linear.YlGn_03,
        "YlGn_04": branca.colormap.linear.YlGn_04,
        "YlGn_05": branca.colormap.linear.YlGn_05,
        "YlGn_06": branca.colormap.linear.YlGn_06,
        "YlGn_07": branca.colormap.linear.YlGn_07,
        "YlGn_08": branca.colormap.linear.YlGn_08,
        "YlGn_09": branca.colormap.linear.YlGn_09,
        "YlGnBu_03": branca.colormap.linear.YlGnBu_03,
        "YlGnBu_04": branca.colormap.linear.YlGnBu_04,
        "YlGnBu_05": branca.colormap.linear.YlGnBu_05,
        "YlGnBu_06": branca.colormap.linear.YlGnBu_06,
        "YlGnBu_07": branca.colormap.linear.YlGnBu_07,
        "YlGnBu_08": branca.colormap.linear.YlGnBu_08,
        "YlGnBu_09": branca.colormap.linear.YlGnBu_09,
        "YlOrBr_03": branca.colormap.linear.YlOrBr_03,
        "YlOrBr_04": branca.colormap.linear.YlOrBr_04,
        "YlOrBr_05": branca.colormap.linear.YlOrBr_05,
        "YlOrBr_06": branca.colormap.linear.YlOrBr_06,
        "YlOrBr_07": branca.colormap.linear.YlOrBr_07,
        "YlOrBr_08": branca.colormap.linear.YlOrBr_08,
        "YlOrBr_09": branca.colormap.linear.YlOrBr_09,
        "YlOrRd_03": branca.colormap.linear.YlOrRd_03,
        "YlOrRd_04": branca.colormap.linear.YlOrRd_04,
        "YlOrRd_05": branca.colormap.linear.YlOrRd_05,
        "YlOrRd_06": branca.colormap.linear.YlOrRd_06,
        "YlOrRd_07": branca.colormap.linear.YlOrRd_07,
        "YlOrRd_08": branca.colormap.linear.YlOrRd_08,
        "YlOrRd_09": branca.colormap.linear.YlOrRd_09,
    }
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )

    color_map = knext.StringParameter(
        "Color map",
        "Select the color map to use for the heatmap. See [branca](https://python-visualization.github.io/branca/colormap.html) for more information.",
        default_value="YlOrRd_09",
        enum=list(_color_bars.keys()),
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
        equals maxZoom of the map by default.""",
        default_value=18,
    )

    radius = knext.IntParameter(
        "Radius",
        "Radius of each datapoint of the heatmap.",
        default_value=25,
    )

    blur = knext.IntParameter(
        "Blur",
        """The blur factor that will be applied to all data points. 
        The higher the blur factor is, the smoother the gradients will be.""",
        default_value=15,
    )

    basemap_settings = BaseMapSettings()

    legend_settings = LegendSettings()

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)
        # convert to WGS84
        gdf = gdf.to_crs(4326)

        if self.basemap_settings.base_map == "Don't show base map":
            base_map = None
        else:
            base_map = self.basemap_settings.base_map

        map = gdf.explore(tiles=base_map)
        gdf["lon"] = gdf.geometry.centroid.x
        gdf["lat"] = gdf.geometry.centroid.y

        if "none" not in str(self.weight_col).lower():
            heat_data = gdf[["lat", "lon", self.weight_col]]
        else:
            heat_data = gdf[["lat", "lon"]]

        # add color bar on the top of the map
        steps = 20
        linear_colormap = self._color_bars[self.color_map]
        colormap = linear_colormap.scale(0, 1).to_step(steps)
        colormap.caption = self.legend_settings.caption

        from collections import defaultdict

        gradient_map = defaultdict(dict)
        for i in range(steps):
            gradient_map[1 / steps * i] = colormap.rgb_hex_str(1 / steps * i)
        if self.legend_settings.plot:
            colormap.add_to(map)  # add color bar at the top of the map

        from folium import plugins
        import folium

        plugins.HeatMap(
            heat_data,
            name="heatmap",
            min_opacity=self.min_opacity,
            max_zoom=self.max_zoom,
            radius=self.radius,
            blur=self.blur,
            gradient=gradient_map,
        ).add_to(map)

        folium.LayerControl().add_to(map)

        # replace css and JavaScript paths
        html = map.get_root().render()
        html = replace_external_js_css_paths(
            r'\1./libs/leaflet/1.9.3/\3"\4',
            html,
        )

        return knext.view(html)
