from re import S
from typing import Callable
import pandas as pd
import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut
import pyproj as pyp
import geopy
from geopy.extra.rate_limiter import RateLimiter



__category = knext.category(
    path="/community/geo",
    level_id="conversion",
    name="Spatial Conversion",
    description="Nodes that perform conversions from/to geometric objects.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/ConversionCategory.png",
    after="transform",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/GeometryConversion/"

########################################################################################
# Conversion helper methods
########################################################################################


def crs_input_parameter(
    label: str = "CRS",
    description: str = knut.DEF_CRS_DESCRIPTION,
) -> knext.StringParameter:
    """Returns a CRS (coordinate reference system) string input parameter."""
    return knext.StringParameter(
        label=label,
        description=description,
        default_value=knut.DEFAULT_CRS,
    )


def validate_crs(crs: str) -> None:
    """Validates the input crs and throws an appropriate exception if it is invalid."""
    try:
        parse_crs(crs)
    except Exception as error:
        raise knext.InvalidParametersError(str(error))


def parse_crs(crs: str) -> pyp.CRS:
    """Parses the input crs into a CRS object and throws an exception if it is invalid."""
    return pyp.CRS.from_user_input(crs)


########################################################################################
# Converter that converts a single input column to a geometry column
########################################################################################
class _ToGeoConverter:
    """
    Helper class for conversion that creates a new geometry column from a single input column.
    """

    DEF_GEO_COL_NAME = "geometry"

    def __init__(
        self,
        converter: Callable[[pd.DataFrame, knext.Column], gp.GeoSeries],
        column_type_check: Callable[[knext.Column], bool] = None,
    ):
        self.converter = converter
        self.column_type_check = column_type_check

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Data table",
                description=f"Input table to append the generated geometry column to.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description="Input table with additional geometry column in the units of the provided CRS.",
            )
        ]

    def configure(self, configure_context, input_schema):
        # input_column and crs need to be defined in the child class
        self.input_column = knut.column_exists_or_preset(
            configure_context, self.input_column, input_schema, self.column_type_check
        )
        validate_crs(self.crs)

        import knime.types.geospatial as gt

        return input_schema.append(
            knext.Column(
                ktype=knext.logical(gt.GeoValue),
                name=knut.get_unique_column_name(self.DEF_GEO_COL_NAME, input_schema),
            )
        )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        df = input_table.to_pandas()
        try:
            # input_column and crs need to be defined in the child class
            result_col_name = knut.get_unique_column_name(
                self.DEF_GEO_COL_NAME, input_table.schema
            )
            df[result_col_name] = self.converter(df, self.input_column)
            gdf = gp.GeoDataFrame(df, geometry=result_col_name, crs=self.crs)
        except Exception as error:
            raise RuntimeError(
                f"Error converting '{self.input_column}' column to geometry: {str(error)}"
            )
        # execute the function and append the result
        return knut.to_table(gdf, exec_context)


############################################
# WKT to Geometry
############################################
@knext.node(
    name="WKT to Geometry",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "WKTtoGeo.png",
    category=__category,
    after="",
)
@knut.geo_node_description(
    short_description="Converts the input Well-known-text (WKT) column to a geometry column.",
    description="""This node converts the selected 
    [Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) input column to 
    a geometry column in the units of the provided CRS.
    """,
    references={
        "From WKT": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.from_wkt.html",
        "CRS parser from pyproj": "https://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input",
    },
)
class WKTtoGeoNode(_ToGeoConverter):

    input_column = knext.ColumnParameter(
        label="WKT column",
        description="[Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column to convert",
        column_filter=knut.is_string,
        include_row_key=False,
        include_none_column=False,
    )

    crs = crs_input_parameter()

    def __init__(self):
        super().__init__(
            lambda df, col: gp.GeoSeries.from_wkt(df[col].astype("str")), knut.is_string
        )


############################################
# WKB to Geometry
############################################
# binary data cells not supported yet in the KNIME Python Extension
# @knext.node(
#     name="WKB to Geometry",
#     node_type=knext.NodeType.MANIPULATOR,
#     icon_path=__NODE_ICON_PATH + "WKBtoGeo.png",
#     category=__category,
#     after=""
# )
# @knut.geo_node_description(
#     short_description="Converts the input Well-known-binary (WKB) column to a geometry column.",
#     description="""This node converts the selected
#     [Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) input column to
#     a geometry column in the units of the provided CRS.
#     """,
#     references={
#         "From WKB": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.from_wkb.html",
#         "CRS parser from pyproj": "https://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input",
#     },
# )
# class WKBtoGeoNode(_ToGeoConverter):

#     input_column = knext.ColumnParameter(
#         label="WKB column",
#         description="[Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column to convert",
#         column_filter=knut.is_binary,
#         include_row_key=False,
#         include_none_column=False,
#     )

#     crs = crs_input_parameter()

#     def __init__(self):
#         super().__init__(lambda df, col: gp.GeoSeries.from_wkb(df[col]), knut.is_binary)


############################################
# WKT to Geometry
############################################
@knext.node(
    name="GeoJSON to Geometry",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoJSONtoGeo.png",
    category=__category,
    after="",
)
@knut.geo_node_description(
    short_description="Converts the input GeoJSON column to a geometry column.",
    description="""This node converts the selected [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON) input column to 
    a geometry column in the units of the provided CRS.
    """,
    references={
        "GeoJSON cells are read using the read_file function": "https://geopandas.org/en/stable/docs/reference/api/geopandas.read_file.html#geopandas.read_file",
        "CRS parser from pyproj": "https://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input",
    },
)
class GeoJSONtoGeoNode(_ToGeoConverter):

    input_column = knext.ColumnParameter(
        label="GeoJSON formatted string column",
        description="String column with a [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON) string to convert",
        column_filter=knut.is_string,
        include_row_key=False,
        include_none_column=False,
    )

    crs = crs_input_parameter()

    def __init__(self):
        super().__init__(
            lambda df, col: df,
            knut.is_string,
        )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        df = input_table.to_pandas()

        result_col_name = knut.get_unique_column_name(
            self.DEF_GEO_COL_NAME, input_table.schema
        )
        try:
            # read each GeoJson cell and append the result to the table
            df[result_col_name] = df[self.input_column].apply(
                lambda row: gp.read_file(row, driver="GeoJSON").geometry
            )
            gdf = gp.GeoDataFrame(df, geometry=result_col_name, crs=self.crs)
        except Exception as error:
            raise RuntimeError(
                f"Error converting '{self.input_column}' column to geometry: {str(error)}"
            )
        # execute the function and append the result
        return knut.to_table(gdf, exec_context)


############################################
# Lat/Lon to Geometry
############################################
@knext.node(
    name="Lat/Lon to Geometry",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "LatLongToGeo.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Table with latitude and longitude",
    description="Table with columns for latitude and longitude to convert to a geometric point.",
)
@knext.output_table(
    name="Geo table",
    description="Input table with additional geometry column in the units of the provided CRS.",
)
@knut.geo_node_description(
    short_description="Converts the given latitude and longitude column to a geometric point column.",
    description="This node converts the given latitude and longitude into a geometric point."
    + "The geometric point is appended to the input table as geometry column.",
    references={
        "Points from XY": "https://geopandas.org/en/stable/docs/reference/api/geopandas.points_from_xy.html",
    },
)
class LatLongToGeoNode:

    lat_col = knext.ColumnParameter(
        label="Latitude column",
        description="Please select the latitude column",
        column_filter=knut.is_numeric,
    )

    lon_col = knext.ColumnParameter(
        label="Longitude column",
        description="Please select the longitude column",
        column_filter=knut.is_numeric,
    )

    crs = crs_input_parameter()

    def configure(self, configure_context, input_schema):
        knut.column_exists(self.lat_col, input_schema, knut.is_numeric)
        knut.column_exists(self.lon_col, input_schema, knut.is_numeric)
        knut.fail_if_column_exists(
            _ToGeoConverter.DEF_GEO_COL_NAME,
            input_schema,
            f"Output column '{_ToGeoConverter.DEF_GEO_COL_NAME}'  exists in input table",
        )
        validate_crs(self.crs)

        from shapely.geometry import Point

        return input_schema.append(
            knext.Column(
                ktype=knext.logical(Point),
                name=knut.get_unique_column_name(
                    _ToGeoConverter.DEF_GEO_COL_NAME, input_schema
                ),
            )
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # https://geopandas.org/en/stable/docs/user_guide/projections.html#the-axis-order-of-a-crs
        # In GeoPandas the coordinates are always stored as (x, y), and thus as (lon, lat) order, regardless of the CRS
        df = input_table.to_pandas()
        result_col_name = knut.get_unique_column_name(
            _ToGeoConverter.DEF_GEO_COL_NAME, input_table.schema
        )
        try:
            # input_column and crs need to be defined in the child class
            df[result_col_name] = gp.points_from_xy(df[self.lon_col], df[self.lat_col])
            gdf = gp.GeoDataFrame(df, geometry=result_col_name, crs=self.crs)
        except Exception as error:
            raise RuntimeError(
                f"Error converting '{self.lat_col}' and '{self.lon_col}' column to geometry: {str(error)}"
            )
        # execute the function and append the result
        return knut.to_table(gdf, exec_context)


########################################################################################
# Converter that converts a geometry column to another column
########################################################################################
class _FromGeoConverter:
    """
    Helper class for conversion that creates a new column from a geometry column.
    """

    def __init__(
        self,
        column_name: str,
        column_type: knext.KnimeType,
        converter: Callable[[gp.GeoSeries], gp.GeoSeries],
    ):
        self.column_name = column_name
        self.column_type = column_type
        self.converter = converter

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Data table",
                description=f"Input table with geometry column to extract the {self.column_name} column from.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description=f"Input table with additional {self.column_name} column extracted from the specified geometry column.",
            )
        ]

    def configure(self, configure_context, input_schema):
        # geo_column needs to be defined in the child class
        self.geo_column = knut.column_exists_or_preset(
            configure_context, self.geo_column, input_schema, knut.is_geo
        )
        return input_schema.append(
            knext.Column(
                self.column_type,
                knut.get_unique_column_name(self.column_name, input_schema),
            )
        )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        gdf = knut.load_geo_data_frame(input_table, self.geo_column, exec_context)
        # execute the function and append the result
        result_col_name = knut.get_unique_column_name(
            self.column_name, input_table.schema
        )
        try:
            result = self.converter(gdf[self.geo_column])
            gdf[result_col_name] = result
        except Exception as error:
            raise knext.InvalidParametersError(
                f"Error converting to {result_col_name}: {str(error)}"
            )
        return knut.to_table(gdf, exec_context)


############################################
# Geometry to WKT
############################################
@knext.node(
    name="Geometry to WKT",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoToWKT.png",
    category=__category,
    after="",
)
@knut.geo_node_description(
    short_description="Converts the input geometry column to a Well-known-text (WKT) column.",
    description="""This node converts the selected geometry column into a
    [Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry).
    """,
    references={
        "To WKT": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_wkt.html"
    },
)
class GeoToWKTNode(_FromGeoConverter):

    geo_column = knut.geo_col_parameter(
        description="Geometry column to convert to [Well-known-text (WKT)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry)"
    )

    def __init__(self):
        super().__init__(
            "WKT",
            knext.string(),
            lambda gs: gs.to_wkt(),
        )


############################################
# Geometry to WKB
############################################
# binary data cells not supported yet in the KNIME Python Extension
# @knext.node(
#     name="Geometry to WKB",
#     node_type=knext.NodeType.MANIPULATOR,
#     icon_path=__NODE_ICON_PATH + "GeoToWKB.png",
#     category=__category,
#     after=""
# )
# @knut.geo_node_description(
#     short_description="Converts the input geometry column to a Well-known-binary (WKB) column.",
#     description="""This node converts the selected geometry column into a
#     [Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) column.
#     """,
#     references={
#         "To WKB": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_wkb.html"
#     },
# )
# class GeoToWKBNode(_FromGeoConverter):

#     geo_column = knut.geo_col_parameter(
#         description="Geometry column to convert to [Well-known-binary (WKB)](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry)"
#     )

#     def __init__(self):
#         super().__init__(
#             "WKB",
#             knext.blob(),
#             lambda gs: gs.to_wkb(),
#         )


############################################
# Geometry to GeoJSON
############################################
@knext.node(
    name="Geometry to GeoJSON",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoToGeoJSON.png",
    category=__category,
    after="",
)
@knut.geo_node_description(
    short_description="Converts the input geometry column to a Well-known-binary (WKB) column.",
    description="""This node converts the selected geometry column into a 
    [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON) column.
    """,
    references={
        "To JSON": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.to_json.html"
    },
)
class GeoToGeoJSONNode(_FromGeoConverter):

    geo_column = knut.geo_col_parameter(
        description="Geometry column to convert to [GeoJSON](https://en.wikipedia.org/wiki/GeoJSON)"
    )

    def __init__(self):
        super().__init__(
            "GeoJSON",
            # TODO: JSON
            knext.string(),
            lambda gs: gs.to_json(),
        )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column
        gdf = knut.load_geo_data_frame(input_table, self.geo_column, exec_context)
        # execute the function and append the result
        result_col_name = knut.get_unique_column_name(
            self.column_name, input_table.schema
        )
        try:
            gdf[result_col_name] = gdf[self.geo_column].apply(
                lambda e: gp.GeoSeries(e).to_json()
            )
        except Exception as error:
            raise knext.InvalidParametersError(
                f"Error converting to {result_col_name}: {str(error)}"
            )
        return knut.to_table(gdf, exec_context)


############################################
# Geometry to Lat/Lon
############################################
@knext.node(
    name="Geometry to Lat/Long",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoToLatLon.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with a point geometry column to extract the latitude and longitude from.",
)
@knext.output_table(
    name="Table with latitude and longitude",
    description="Table with the latitude and longitude for all point geometries",
)
@knut.geo_node_description(
    short_description="Extracts the latitude and longitude.",
    description="This node extracts the latitude and longitude for all given point objects."
    + "The longitude and latitude columns are appended to the input table as lat and long column.",
    references={
        "X coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.x.html",
        "Y coordinate": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.y.html",
    },
)
class GeoToLatLongNode:

    geo_col = knut.geo_point_col_parameter(
        description="Select the point geometry column to extract the latitude and longitude from."
    )

    col_lat = "lat"
    col_lon = "lon"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo_point
        )
        return input_schema.append(
            [
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self.col_lat, input_schema),
                ),
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self.col_lon, input_schema),
                ),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = knut.load_geo_data_frame(input, self.geo_col, exec_context)
        gs = gdf.loc[:, self.geo_col]
        # https://geopandas.org/en/stable/docs/user_guide/projections.html#the-axis-order-of-a-crs
        # In GeoPandas the coordinates are always stored as (x, y), and thus as (lon, lat) order, regardless of the CRS
        # However the standard order for EPSG:4326 is latitude, longitude which is why we use this order in the output
        # Also the ISO 6709 standard (https://en.wikipedia.org/wiki/ISO_6709) uses this order
        gdf[knut.get_unique_column_name(self.col_lat, input.schema)] = gs.y
        gdf[knut.get_unique_column_name(self.col_lon, input.schema)] = gs.x
        return knut.to_table(gdf, exec_context)

############################################
# geocoding (address to geometry)
############################################

@knext.node(
    name="Geocoding",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoGeocoding.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input table",
    description="Table with an address column to geocode.",
)
@knext.output_table(
    name="Table with geometry",
    description="Table with the geometry for all given addresses.",
)
@knut.geo_node_description(
    short_description="Geocodes the given addresses.",
    description="This node geocodes the given addresses and appends the geometry to the input table.",
    references={
        "Geocoding": "https://en.wikipedia.org/wiki/Geocoding",
        "Geopy": "https://geopy.readthedocs.io/en/stable/",
    },
)
class GeoGeocodingNode:
    
        address_col = knext.ColumnParameter(
            "Address column",
            "Select the address column to geocode." 
            + "The column must contain a string with the address to geocode.",
             column_filter=knut.is_string,
            include_row_key=False,
        include_none_column=False,
            )

        service_provider = knext.StringParameter(
            "Service provider",
            "Select the service provider to use for geocoding.",
            default_value="Nominatim",
            enum=list(geopy.geocoders.SERVICE_TO_GEOCODER.keys())
        )

        api_key = knext.StringParameter(
            "API key",
            "Enter the API key for the service provider.",
            default_value=""
        )

        min_delay_seconds = knext.IntParameter(
            "Minimum delay (seconds)",
            "Enter the minimum delay in seconds between two geocoding requests.",
            default_value=1
        )

        default_timeout = knext.IntParameter(
            "Default timeout (seconds)",
            "Enter the default timeout in seconds for geocoding requests.",
            default_value=10
        )

    
        def configure(self, configure_context, input_schema):
           
            return None
    
        def execute(self, exec_context: knext.ExecutionContext, input_table):
            df = input_table.to_pandas()
            # geopy.geocoders.options.default_user_agent = 'my_app/1'
            geopy.geocoders.options.default_timeout = self.default_timeout

            service_provider = geopy.geocoders.SERVICE_TO_GEOCODER[self.service_provider] 
            geolocator = service_provider(api_key=self.api_key)

            geocode = RateLimiter(geolocator.geocode, min_delay_seconds=self.min_delay_seconds)
            df['latitude'] = df[self.address_col].apply(lambda x: geocode(x).latitude)
            df['longitude'] = df[self.address_col].apply(lambda x: geocode(x).longitude)
            gdf = gp.GeoDataFrame(df, geometry=gp.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")

            return knut.to_table(df)
        

############################################
# reverse geocoding (geometry to address)
############################################

@knext.node(
    name="Reverse Geocoding",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "GeoReverseGeocoding.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Input table",
    description="Table with a geometry column to reverse geocode.",
)
@knext.output_table(
    name="Table with address",
    description="Table with the address for all given geometries.",
)
@knut.geo_node_description(
    short_description="Reverse geocodes the given geometries.",
    description="This node reverse geocodes the given geometries and appends the address to the input table.",
    references={
        "Reverse Geocoding": "https://en.wikipedia.org/wiki/Reverse_geocoding",
        "Geopy": "https://geopy.readthedocs.io/en/stable/",
    },
)
class GeoReverseGeocodingNode:
        
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to reverse geocode." 
        + "The column must contain a geometry.",
        column_filter=knut.is_geo,
        include_row_key=False,
    include_none_column=False,
        )

    service_provider = knext.StringParameter(
        "Service provider",
        "Select the service provider to use for reverse geocoding.",
        default_value="Nominatim",
        enum=list(geopy.geocoders.SERVICE_TO_GEOCODER.keys())
    )

    api_key = knext.StringParameter(
        "API key",
        "Enter the API key for the service provider.",
        default_value=""
    )

    min_delay_seconds = knext.IntParameter(
        "Minimum delay (seconds)",
        "Enter the minimum delay in seconds between two reverse geocoding requests.",
        default_value=1
    )

    default_timeout = knext.IntParameter(
        "Default timeout (seconds)",
        "Enter the default timeout in seconds for reverse geocoding requests.",
        default_value=10
    )


    def configure(self, configure_context, input_schema):
        
        return None
    
    def execute(self, exec_context: knext.ExecutionContext, input_table):
        df = input_table.to_pandas()
        df.rename(columns={self.geo_col: "geometry"}, inplace=True)
        gdf = gp.GeoDataFrame(df, geometry="geometry")
        # geopy.geocoders.options.default_user_agent = 'my_app/1'
        geopy.geocoders.options.default_timeout = self.default_timeout

        service_provider = geopy.geocoders.SERVICE_TO_GEOCODER[self.service_provider] 
        geolocator = service_provider(api_key=self)

        reverse = RateLimiter(geolocator.reverse, min_delay_seconds=self.min_delay_seconds)
        gdf['address'] = gdf['geometry'].apply(lambda x: reverse((x.y,x.x)).address)
        

        return knut.to_table(gdf)