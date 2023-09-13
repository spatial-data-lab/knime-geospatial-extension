import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


__category = knext.category(
    path="/community/geo",
    level_id="calculation",
    name="Spatial Calculation",
    description="Nodes that facilitate the calculation and representation of spatial properties, enveloping geometries, and complex geometric unions.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/CalculationCategory.png",
    after="io",
)

# Root path for all node icons in this file
__NODE_ICON_PATH = "icons/icon/GeometryCalculation/"


############################################
# Simple transformation helper class
############################################
class _SingleCalculator:
    """
    Helper class for transformations that append a single result column to an input table with geo column.
    """

    from typing import Callable

    def __init__(
        self,
        func: Callable[[gp.GeoDataFrame], gp.GeoDataFrame],
        col_name: str,
        col_type: knext.KnimeType = knext.double(),
    ):
        self.func = func
        self.col_name = col_name
        self.coltype = col_type

        # standard input port description can be overwritten in the child classes
        self.input_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name="Geo table",
                description=f"Input table with a geospatial column to compute the '{self.col_name}' for.",
            )
        ]
        # standard output port description can be overwritten in the child classes
        self.output_ports = [
            knext.Port(
                type=knext.PortType.TABLE,
                name=f"Geo table with {self.col_name}",
                description=f"""Geo input table with additional '{self.col_name}' column containing the
                {self.col_name} of each input geometry in the units of the CRS.
                """,
            )
        ]

    def configure(self, configure_context, input_schema):
        # geo_col and result_settings needs to be defined in the child class

        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        if self.coltype is not None:
            return self.result_settings.get_result_schema(
                configure_context,
                input_schema,
                self.geo_col,
                self.coltype,
            )

    def execute(self, exec_context, input_table):
        # create GeoDataFrame with the selected column as active geometry column.
        # The geo_col and result_settings needs to be defined in the child class

        return self.result_settings.get_computed_result_table(
            exec_context, input_table, self.geo_col, self.func
        )


############################################
# Area node
############################################
@knext.node(
    name="Area",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Area.png",
    category=__category,
    after="",
)
@knut.geo_node_description(
    short_description="Calculates the area of geometric objects.",
    description="""This node will add a new column named 'area' to the input data with the calculated area of
    each geometry in the units of the CRS.  
    Area may be invalid for a geographic CRS using degrees as units. Use the Projection node to project 
    geometries to a planar CRS before using this node.
    This operation is planar, i.e. the potential third dimension is not taken into account.
    """,
    references={
        "Area": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.area.html"
    },
)
class AreaNode(_SingleCalculator):
    __COL_NAME = "area"

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the area."
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name=__COL_NAME,
        since_version="1.2.0",
    )

    def __init__(self):
        super().__init__(lambda gdf: gdf.area, self.__COL_NAME)


############################################
# Length node
############################################
@knext.node(
    name="Length",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Length.png",
    category=__category,
    after="",
)
@knut.geo_node_description(
    short_description="Calculates the length of geometric objects.",
    description="""This node will add a new column named 'length' to the input data
    with the calculated length of each geometry in the units of the CRS.  
    Length may be invalid for a geographic CRS using degrees as units. Use the Projection node to project 
    geometries to a planar CRS before using this node.
    This operation is planar, i.e. the potential third dimension is not taken into account.
    """,
    references={
        "Length": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.length.html"
    },
)
class LengthNode(_SingleCalculator):
    __COL_NAME = "length"

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the length."
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name=__COL_NAME,
        since_version="1.2.0",
    )

    def __init__(self):
        super().__init__(lambda gdf: gdf.length, self.__COL_NAME)


############################################
# Bounding Box
############################################
@knext.node(
    name="Bounding Box",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "BoundingBox.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
@knut.geo_node_description(
    short_description="This node generate rectangles representing the envelope of each geometry.",
    description="""This node generate rectangles representing the envelope of each geometry. 
    That is, the point or smallest rectangular polygon (with sides parallel to the coordinate axes) 
    that contains each of the geometries.""",
    references={
        "Envelop": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.envelope.html",
    },
)
class BoundingBoxNode(_SingleCalculator):
    __COL_NAME = "box"

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the bounding box."
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name=__COL_NAME,
        since_version="1.2.0",
    )

    def __init__(self):
        super().__init__(lambda gdf: gdf.envelope, self.__COL_NAME, knut.TYPE_GEO)


############################################
# Convex Hull
############################################
@knext.node(
    name="Convex Hull",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "ConvexHull.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed Geo table",
    description="Transformed Geo input table",
)
@knut.geo_node_description(
    short_description="This node generate the smallest convex Polygon containing all the points in each geometry.",
    description="""This node generate the smallest convex Polygon containing all the points in each geometry.
    The convex hull of a geometry is the smallest convex Polygon containing all the points in each geometry, 
    unless the number of points in the geometric object is less than three. For two points, the convex hull 
    collapses to a LineString; for 1, a Point.""",
    references={
        "Convex hull": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.convex_hull.html",
    },
)
class ConvexHullNode(_SingleCalculator):
    __COL_NAME = "hull"

    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the convex hull."
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name=__COL_NAME,
        since_version="1.2.0",
    )

    def __init__(self):
        super().__init__(lambda gdf: gdf.convex_hull, self.__COL_NAME, knut.TYPE_GEO)


############################################
# Coordinates node
############################################
@knext.node(
    name="Coordinates XYZ",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "XYZcoordinates.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with a point geometry column to extract the coordinates from.",
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
class CoordinatesNode:
    geo_col = knut.geo_point_col_parameter(
        description="Select the point geometry column to extract the coordinates from."
    )

    add_z = knext.BoolParameter(
        "Extract Z coordinate",
        "Select this option to extract the Z coordinate. If not selected only the X and Y coordinate are extracted",
        False,
    )

    _COL_X = "x"
    _COL_Y = "y"
    _COL_Z = "z"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo_point
        )
        result = input_schema_1.append(
            [
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_X, input_schema_1),
                ),
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_Y, input_schema_1),
                ),
            ]
        )
        if self.add_z:
            result = result.append(
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_Z, input_schema_1),
                )
            )
        return result

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        gs = gdf.loc[:, self.geo_col]
        gdf[knut.get_unique_column_name(self._COL_X, input_1.schema)] = gs.x
        gdf[knut.get_unique_column_name(self._COL_Y, input_1.schema)] = gs.y
        if self.add_z:
            gdf[knut.get_unique_column_name(self._COL_Z, input_1.schema)] = gs.z
        return knut.to_table(gdf, exec_context)


############################################
# Bounds node
############################################
@knext.node(
    name="Bounds",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "Bounds.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to compute the bounds per geometric object",
)
@knext.output_table(
    name="Geo table with bounds",
    description="Input geo table with the bounds per geometric object",
)
@knut.geo_node_description(
    short_description="Computes the bound for geometric objects.",
    description="""This node computes the bounds for each geometry object. The bounds are appended to the input 
    table as x/y coordinate columns.""",
    references={
        "Bounds": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.bounds.html"
    },
)
class BoundsNode:
    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the bounds for."
    )
    _COL_MIN_X = "minx"
    _COL_MIN_Y = "miny"
    _COL_MAX_X = "maxx"
    _COL_MAX_Y = "maxy"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return input_schema_1.append(
            [
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_MIN_X, input_schema_1),
                ),
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_MIN_Y, input_schema_1),
                ),
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_MAX_X, input_schema_1),
                ),
                knext.Column(
                    knext.double(),
                    knut.get_unique_column_name(self._COL_MAX_Y, input_schema_1),
                ),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        bounds = gdf.bounds
        bounds.rename(
            columns={
                self._COL_MIN_X: knut.get_unique_column_name(
                    self._COL_MIN_X, input_1.schema
                ),
                self._COL_MIN_Y: knut.get_unique_column_name(
                    self._COL_MIN_Y, input_1.schema
                ),
                self._COL_MAX_X: knut.get_unique_column_name(
                    self._COL_MAX_X, input_1.schema
                ),
                self._COL_MAX_Y: knut.get_unique_column_name(
                    self._COL_MAX_Y, input_1.schema
                ),
            },
            inplace=True,
        )

        import pandas as pd

        gdf = pd.concat([gdf, bounds], axis=1, copy=False)
        return knut.to_table(gdf, exec_context)


############################################
# Total bounds node
############################################
@knext.node(
    name="Total Bounds",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "TotalBounds.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to compute the total bounds for.",
)
@knext.output_table(
    name="Table with the total bounds",
    description="Table with the total bounds for all geometric input objects",
)
@knut.geo_node_description(
    short_description="Computes the total bounds for all given geometric objects.",
    description="This node computes the total bounds for all given geometric input objects.",
    references={
        "Total bounds": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.total_bounds.html"
    },
)
class TotalBoundsNode:
    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the total bounds for."
    )

    _COL_MIN_X = "minx"
    _COL_MIN_Y = "miny"
    _COL_MAX_X = "maxx"
    _COL_MAX_Y = "maxy"

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return knext.Schema.from_columns(
            [
                knext.Column(knext.double(), self._COL_MIN_X),
                knext.Column(knext.double(), self._COL_MIN_Y),
                knext.Column(knext.double(), self._COL_MAX_X),
                knext.Column(knext.double(), self._COL_MAX_Y),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        bounds = gdf.total_bounds

        import pandas as pd

        result = pd.DataFrame(
            [bounds],
            columns=[
                self._COL_MIN_X,
                self._COL_MIN_Y,
                self._COL_MAX_X,
                self._COL_MAX_Y,
            ],
        )
        return knut.to_table(result, exec_context)


############################################
# Unary Union
############################################


@knext.node(
    name="Unary Union",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "UnaryUnion.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed Geo input table",
)
@knut.geo_node_description(
    short_description="This node returns a geometry containing the union of all geometries in the input column.",
    description="This node returns a geometry containing the union of all geometries in the input column.",
    references={
        "Unary union": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.unary_union.html",
    },
)
class UnaryUnionNode:
    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the unary union."
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = gp.GeoDataFrame(input.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting transformation..."
        )
        gdf_union = gdf.unary_union
        gdfunion = gp.GeoDataFrame(geometry=gp.GeoSeries(gdf_union), crs=gdf.crs)
        exec_context.set_progress(0.1, "Transformation done")
        return knext.Table.from_pandas(gdfunion)


############################################
# MinimumBoundingCircle
############################################


@knext.node(
    name="Bounding Circle",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "BoundCircle.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with the geometry column to compute the bounding circle for",
)
@knext.output_table(
    name="Geo table with bounding circle",
    description="Input table with the computed bounding circle",
)
@knut.geo_node_description(
    short_description="For each input geometry this node generates the minimum bounding circle that contains all its points.",
    description="""For each input geometry this node generates the minimum bounding circle that contains all its 
    points. The minimum bounding circle is determined by finding the two points that are farthest apart, and the 
    center and radius of the circle are calculated from these points.""",
    references={
        "Minimum bounding circle": "https://en.wikipedia.org/wiki/Smallest-circle_problem",
        "Scipy distance matrix": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance_matrix.html",
        "Shapely convex hull": "https://shapely.readthedocs.io/en/stable/manual.html#object.convex_hull",
    },
)
class BoundCircleNode:
    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the minimum bounding circle."
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name="Circle",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            knut.TYPE_POLYGON,
        )

    def execute(self, exec_context: knext.ExecutionContext, input):
        from shapely.geometry import Point, Polygon
        import numpy as np
        from scipy.spatial import distance_matrix

        def minimum_bounding_circle(points):
            # Compute the pairwise Euclidean distances
            D = distance_matrix(points, points)

            # Find the pair of points with the maximum distance
            i, j = np.unravel_index(np.argmax(D), D.shape)

            # Compute the center and radius of the minimum bounding circle
            lon1, lat1 = points[i]
            lon2, lat2 = points[j]
            center_lon = (lon1 + lon2) / 2
            center_lat = (lat1 + lat2) / 2
            center = Point(center_lon, center_lat)
            radius = D[i, j] / 2
            if radius == 0:
                radius = 1e-6  # set a small radius to avoid divide-by-zero errors
            return center, radius

        def compute_circle(geom):
            # Get the points of the convex hull
            hull = geom.convex_hull
            points = np.array(hull.boundary.coords)
            # Compute the minimum bounding circle
            center, radius = minimum_bounding_circle(points)
            return center.buffer(radius)

        return self.result_settings.get_computed_result_table(
            exec_context, input, self.geo_col, compute_circle
        )


############################################
# Line Endpoints
############################################


@knext.node(
    name="Line Endpoints",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path=__NODE_ICON_PATH + "LineEndpoints.png",
    category=__category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with the geometry column to compute the end points for",
)
@knext.output_table(
    name="Geo table with endpoints",
    description="Input table with the computed endpoints",
)
@knut.geo_node_description(
    short_description="This node returns two geometry columns representing the start and end points of each line in the input data.",
    description="""This node returns two geometry columns representing the start and end points of each line in the 
input data. Any none LineString geometry types will return a missing value for the start and endpoint.""",
    references={
        "LineString": "https://shapely.readthedocs.io/en/stable/reference/shapely.LineString.html#shapely.LineString",
        "Point": "https://shapely.readthedocs.io/en/stable/reference/shapely.Point.html",
    },
)
class LineEndpointNode:
    geo_col = knut.geo_col_parameter(
        description="Select the geometry column to compute the endpoints.",
    )
    _STARTP = "Start Point"
    _ENDP = "End Point"

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        from shapely.geometry import Point

        return input_schema.append(
            [
                knext.Column(
                    ktype=knext.logical(Point),
                    name=knut.get_unique_column_name(self._STARTP, input_schema),
                ),
                knext.Column(
                    ktype=knext.logical(Point),
                    name=knut.get_unique_column_name(self._ENDP, input_schema),
                ),
            ]
        )

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = knut.load_geo_data_frame(input, self.geo_col, exec_context)

        from shapely.geometry import Point, LineString

        startp = knut.get_unique_column_name(self._STARTP, input.schema)
        endp = knut.get_unique_column_name(self._ENDP, input.schema)

        def get_line_endpoints(row):
            geom = row[self.geo_col]
            geom_type = geom.geom_type
            if geom_type == "LineString":
                line = LineString(geom)
                start_point = Point(line.coords[0])
                end_point = Point(line.coords[-1])
            else:
                start_point = None
                end_point = None
            return (start_point, end_point)

        gdf[[startp, endp]] = gdf.apply(
            get_line_endpoints, axis=1, result_type="expand"
        )
        gdf[startp] = gp.GeoSeries(gdf[startp], crs=gdf.crs)
        gdf[endp] = gp.GeoSeries(gdf[endp], crs=gdf.crs)

        gdf.reset_index(drop=True, inplace=True)

        return knut.to_table(gdf, exec_context)
