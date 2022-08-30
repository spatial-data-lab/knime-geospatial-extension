from doctest import debug_script
import logging

import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut

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


@knext.node(name="Geospatial View", node_type=knext.NodeType.VISUALIZER, icon_path="icons/icon.png", category=category)
@knext.input_table(name="Geospatial table to visualize", description="Table with geospatial data to visualize")
@knext.output_view(name="Geospatial view", description="Showing a map woth the selected geo points.")
class ViewNode:
    """
    This node will visualize the given geometric elements on a map.
    """
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to visualize.",
        column_filter=knut.is_geo_point, # Allows only GeoPoints
        include_row_key=False,
        include_none_column=False
    )

    name_cols: knext.MultiColumnParameter = knext.MultiColumnParameter(
        "Tooltip columns",
        "Select columns which should be shown in the marker tooltip.",
        column_filter=knut.is_string
    )

    color_col = knext.ColumnParameter(
        "Marker color column",
        "Select marker color column. The column must contain the color name e.g. red, green, blue, etc.",
        column_filter=knut.is_string,
        include_row_key=False,
        include_none_column=False
    )

    def configure(self, configure_context, input_schema):
        knut.columns_exist([self.geo_col, self.color_col], input_schema)
        if self.name_cols is None:
            self.name_cols = [c.name for c in input_schema if knut.is_string(c)]
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        gdf = gp.GeoDataFrame(input_table.to_pandas(), geometry=self.geo_col)

        exec_context.set_progress(0.1, "Checking CRS")
        knut.check_canceled(exec_context)
        
        #project the gdf into the default CRS used by folium if necessary
        #TODO: DOes not seem to work like this
        # if gdf.crs != "EPSG:3857":
        #     gdf = gdf.to_crs("EPSG:3857")

        exec_context.set_progress(0.3, "Computing map center")
        knut.check_canceled(exec_context)
        
        #compute the mean coordinate for the tile center
        gdf['Lat'] = gdf[self.geo_col].y
        gdf['Long'] = gdf[self.geo_col].x
        mean_lat_long = gdf[['Lat', 'Long']].mean().values.tolist()
       
        map =  folium.Map(location = mean_lat_long, tiles = "OpenStreetMap", zoom_start = 3)

        exec_context.set_progress(0.5, "Create markers")
        knut.check_canceled(exec_context)
        # Create a geometry list from the GeoDataFrame
        geo_df_list = [[point.xy[1][0], point.xy[0][0]] for point in gdf.geometry ]
        # Iterate through list and add a marker for each volcano, color-coded by its type.
        i = 0
        for coordinates in geo_df_list:
            # Place the markers with the popup labels and data
            tooltip = ""
            for col_name in self.name_cols:
                tooltip += col_name + ": " + str(gdf[col_name][i]) + '<br>'
                
            map.add_child(folium.Marker(location = coordinates,
                                    popup = tooltip + "Coordinates: " + str(geo_df_list[i]),
                                    icon = folium.Icon(color = "%s" % gdf[self.color_col][i])))
            i = i + 1
            
        exec_context.set_progress(0.8, "Fit map view to markers")
        knut.check_canceled(exec_context)
        #calculate the corner coordinates and fit the view port
        sw = gdf[['Lat', 'Long']].min().values.tolist()
        ne = gdf[['Lat', 'Long']].max().values.tolist()
        map.fit_bounds([sw, ne]) 

        exec_context.set_progress(0.95, "Create view")
        knut.check_canceled(exec_context)
        return knext.view(map)
