# lingbo
import geopandas as gp
import logging
import knime_extension as knext
import util.knime_utils as knut
import util.projection as kproj

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/community/geo",
    level_id="transform",
    name="Spatial Transformation",
    description="Nodes that transform, decompose, and generate new geometric entities from single geometric objects.",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/icon/TransformationCategory.png",
    after="spatialtool",
)

############################################
# CRS Transformer
############################################


@knext.node(
    name="Projection",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/Projection.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column to transform.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Transformed geo output table.",
)
@knut.geo_node_description(
    short_description="Projection Transformation",
    description="""This node transforms the 
    [Coordinate reference system (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) of the selected 
    geometry column to the entered new coordinate reference system. The node will transform the points in all 
    objects individually. It has no notion of projecting entire geometries. All segments joining points are assumed 
    to be lines in the current projection, not geodesics. Objects crossing the dateline (or other projection boundary) 
    will have undesirable behavior.
    """,
    references={
        "Map projection (Wikipedia)": "https://en.wikipedia.org/wiki/Map_projection",
        "Comparison of map projection": "https://map-projections.net/index.php",
        "Collection of common map projections and their properties": "https://www.icsm.gov.au/sites/default/files/projections.pdf",
        "Projection wizard that helps to find a good projection": "https://projectionwizard.org/",
        "Coordinate Reference System (CRS) EPSG:4326": "https://epsg.io/4326",
        "geopandas.GeoSeries.to_crs()": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.to_crs.html",
    },
)
class CrsTransformerNode:
    """
    This node projects the data from its original CRS to the entered CRS.
    """

    geo_col = knut.geo_col_parameter()

    new_crs = knext.StringParameter(
        "New CRS", kproj.DEF_CRS_DESCRIPTION, kproj.DEFAULT_CRS
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Projected",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )
        # use the data type of the selected column as result type
        result_type = input_schema[self.geo_col].ktype
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            result_type,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # check whether multiple crs exist
        try:
            gdf = knut.load_geo_data_frame(input_table, self.geo_col, exec_context)
            if self.result_settings.mode == knut.ResultSettingsMode.APPEND.name:
                result_col = knut.get_unique_column_name(
                    self.result_settings.new_column_name, input_table.schema
                )
                gdf[result_col] = gdf[self.geo_col]
                gdf.set_geometry(result_col, inplace=True)
            gdf.to_crs(self.new_crs, inplace=True)
            return knut.to_table(gdf, exec_context)
        except:
            import pyarrow as pa
            import shapely
            import pyproj
            from shapely.ops import transform

            exec_context.set_progress(0.1, "Different input CRS detected")

            t1 = input_table.to_pyarrow()

            # extract the geometry column to project
            t_dict = t1.select([self.geo_col]).to_pydict()
            geometry_col = t_dict[self.geo_col]

            # The function that performs the projection
            def transform_geometry(geom, original_crs, target_crs):
                project = pyproj.Transformer.from_crs(
                    original_crs, target_crs, always_xy=True
                ).transform
                return transform(project, geom)

            transformed_geometries = []
            round_count = 0
            n_loop = len(geometry_col)
            for geo_value in geometry_col:
                round_count += 1
                exec_context.set_progress(
                    0.8 * round_count / n_loop,
                    f"Row {round_count} of {n_loop} processed",
                )
                knut.check_canceled(exec_context)

                if geo_value is None or not geo_value.wkb:
                    transformed_geometries.append(None)
                else:
                    # get WKB from GeoValue, transform to shapely objects
                    shapely_geom = shapely.wkb.loads(bytes(geo_value.wkb))
                    original_crs = geo_value.crs  # CRS
                    transformed_geom = transform_geometry(
                        shapely_geom, original_crs, self.new_crs
                    )
                    # save as WKT
                    transformed_geometries.append(transformed_geom.wkt)

            # Create geospatial column from the projected WKTs
            projected_geo = gp.GeoSeries.from_wkt(
                transformed_geometries, crs=self.new_crs
            )

            # Create a KNIME table with the projected geospatial column
            result_table = knext.Table.from_pandas(
                projected_geo.to_frame(name="projected")
            )
            # Convert the KNIME table to PyArrow and get the projected column
            projected_column = result_table.to_pyarrow().column("projected")

            if self.result_settings.mode == knut.ResultSettingsMode.APPEND.name:
                result_col = knut.get_unique_column_name(
                    self.result_settings.new_column_name, input_table.schema
                )
                t2 = t1.append_column(result_col, projected_column)
            else:
                # get the selected geometry  column index
                geo_col_idx = t1.column_names.index(self.geo_col)
                # and replace it with the projected column
                t2 = t1.set_column(geo_col_idx, self.geo_col, projected_column)
            exec_context.set_progress(1, "Projection finished")

            return knext.Table.from_pyarrow(t2)


############################################
# Geometry To Point
############################################


@knext.node(
    name="Geometry To Point",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/FeatureToPoint.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column.",
)
@knut.geo_node_description(
    short_description="Returns a GeoSeries of points representing each geometry.",
    description="""This node returns a GeoSeries of points representing each geometry.
    There are two types of points, centroids and representative points. 
    The latter is guaranteed to be within each geometry.
    """,
    references={
        "GeoSeries.representative_point": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.representative_point.html",
        "GeoSeries.centroid": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.centroid.html",
    },
)
class GeometryToPointNode:
    """
    This node returns a GeoSeries of points representing each geometry.
    """

    geo_col = knut.geo_col_parameter()

    pointtype = knext.StringParameter(
        "Point Type Selection",
        "The point type to choose from.",
        "centroid",
        enum=["centroid", "representative_point"],
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Point",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo
        )

        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            knut.TYPE_POINT,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        if self.pointtype == "centroid":
            func = lambda l: l.centroid
        else:
            func = lambda l: l.representative_point()

        return self.result_settings.get_computed_result_table(
            exec_context, input_1, self.geo_col, func
        )


############################################
# Multipart to Singlepart
############################################


@knext.node(
    name="Multipart To Singlepart",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/Explode.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column.",
)
@knut.geo_node_description(
    short_description="Explode multi-part geometries into multiple single geometries.",
    description="""This node explodes multi-part geometries into multiple single geometries.
    Each row containing a multi-part geometry will be split into multiple rows with single geometries,
    thereby increasing the number rows in the output table.
    """,
    references={
        "GeoDataFrame.explode": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.explode.html",
    },
)
class ExplodeNode:
    """
    This node dismantles the multiparts into single parts.
    """

    geo_col = knut.geo_col_parameter()

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting explosion...")
        exploded = gdf.explode(ignore_index=True)
        # gdf[self.geo_col] = exploded.geometry
        exec_context.set_progress(0.1, "Explosion done")
        LOGGER.debug("Feature geometry " + self.geo_col + "exploded")
        return knext.Table.from_pandas(exploded)


############################################
# Polygon To Line
############################################


@knext.node(
    name="Polygon To Line",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/PolygonToLine.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column. ",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column. ",
)
@knut.geo_node_description(
    short_description="Returns the boundaries of each polygon.",
    description="""This node return the boundaries of each polygon with geopandas.GeoSeries.boundary,
    which Returns a GeoSeries of lower dimensional objects representing each geometry’s set-theoretic boundary.
    """,
    references={
        "GeoSeries.boundary": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.boundary.html",
    },
)
class PolygonToLineNode:
    """
    This node generate lines from the boundaries of polygons.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.boolean_or(knut.is_geo_polygon, knut.is_geo_multi_polygon),
        include_row_key=False,
        include_none_column=False,
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Line",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context,
            self.geo_col,
            input_schema_1,
            knut.boolean_or(knut.is_geo_polygon, knut.is_geo_multi_polygon),
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema_1,
            self.geo_col,
            knut.TYPE_LINE,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # extract the boundary for each geometry
        return self.result_settings.get_computed_result_table(
            exec_context, input_table, self.geo_col, lambda l: l.boundary
        )


############################################
# Points To Line
############################################


@knext.node(
    name="Points To Line",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/PointToLine.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column.",
)
@knut.geo_node_description(
    short_description="This node generate lines from points according to group id and serial label.",
    description="""This node generate lines from points according to group id and serial label. The result
    table contains one LineString object per group value. Each constructed LineString object represents one or 
    more connected linear splines between the given points of a group ordered by the serial column. 
    Repeated points in the ordered sequence are allowed, but may incur performance penalties and should be avoided. 
    A LineString may cross itself. A LineString has zero area and non-zero length.
    """,
    references={
        "Shapely.LineStrings": "https://shapely.readthedocs.io/en/stable/manual.html#linestrings",
    },
)
class PointsToLineNode:
    """
    This node generate lines from points according to group id and serial label.
    """

    geo_col = knut.geo_point_col_parameter()

    group_col = knext.ColumnParameter(
        "Group column",
        "Select the group column (string) as group id for points.",
        # Allow only string columns
        column_filter=knut.is_int_or_string,
        include_row_key=False,
        include_none_column=False,
    )

    seiral_col = knext.ColumnParameter(
        "Serial column",
        "Select the serial column (numeric) for each group .",
        # Allow only string columns
        column_filter=knut.is_numeric,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo_point
        )
        self.group_col = knut.column_exists_or_preset(
            configure_context, self.group_col, input_schema, knut.is_int_or_string
        )
        self.seiral_col = knut.column_exists_or_preset(
            configure_context, self.seiral_col, input_schema, knut.is_numeric
        )
        return None

    def execute(self, exec_context: knext.ExecutionContext, input):
        gdf = gp.GeoDataFrame(input.to_pandas(), geometry=self.geo_col)
        if self.geo_col != "geometry":
            gdf.rename_geometry("geometry", inplace=True)
        exec_context.set_progress(0.3, "Geo data frame loaded. Starting grouping...")
        from shapely.geometry import MultiPoint, LineString

        line_gdf = (
            gdf.sort_values(by=[self.seiral_col])
            .groupby([self.group_col], as_index=False)["geometry"]
            .apply(lambda x: LineString(x.tolist()))
        )
        line_gdf = gp.GeoDataFrame(line_gdf, geometry="geometry", crs=gdf.crs)
        exec_context.set_progress(0.1, "PointsToLine done")
        return knext.Table.from_pandas(line_gdf)


############################################
# Line To MultiPoint
############################################


@knext.node(
    name="Line To MultiPoint",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/LinePolygonToPoints.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column.",
)
@knut.geo_node_description(
    short_description="This node generate points from the lines.",
    description="""This node generate points from the lines.
    The list of coordinates that describe a geometry are represented as the CoordinateSequence object in Shapely 
    which is the dependence of GeoPandas. 
    """,
    references={
        "Coordinate sequences": "https://shapely.readthedocs.io/en/stable/manual.html",
    },
)
class GeometryToMultiPointNode:
    """
    This node generate points from the lines.
    """

    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo_line,
        include_row_key=False,
        include_none_column=False,
    )

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.REPLACE.name,
        new_name="Multipoint",
    )

    def configure(self, configure_context, input_schema):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema, knut.is_geo_line
        )
        return self.result_settings.get_result_schema(
            configure_context,
            input_schema,
            self.geo_col,
            knut.TYPE_MULTI_POINT,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_table):
        # extract coordinates of each geometry into a new MultiPoint geometry
        from shapely.geometry import MultiPoint

        return self.result_settings.get_computed_result_table(
            exec_context, input_table, self.geo_col, lambda l: MultiPoint(l.coords)
        )


############################################
# Create points in polygon
############################################


@knext.node(
    name="Create Random Points",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/RandomPoint.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Table with geometry column.",
)
@knext.output_table(
    name="Transformed geo table",
    description="Table with transformed geometry column.",
)
@knut.geo_node_description(
    short_description="This node generates random points in polygons.",
    description="""This node generates random points from a uniform distribution within (or along) each input geometry.
For polygons, the points will be sampled within the area of the polygon. For lines, they will be sampled along
the length of the LineString. For multi-part geometries, the weights of each part are selected according to their
relevant attribute (area for Polygons, length for LineStrings), and then points are sampled from each part.
Any other geometry type (e.g., Point, MultiPoint, GeometryCollection) is ignored, and an empty MultiPoint geometry
is returned.
    
The numerical value column is used to determine the number of points to be generated 
for the geometry of the same row. Additionally, you need to provide an ID column that will be used to 
identify the original input row.
    
The node will create a new MultiPoint geometry that includes the random set of points for each input geometry, 
which can be exploded into individual points using the 
[Multipart To Singlepart node.](https://hub.knime.com/center%20for%20geographic%20analysis%20at%20harvard%20university/extensions/sdl.harvard.features.geospatial/latest/org.knime.python3.nodes.extension.ExtensionNodeSetFactory$DynamicExtensionNodeFactory:55ec235c/)
    """,
    references={
        "Sampling points user guide": "https://geopandas.org/en/stable/docs/user_guide/sampling.html",
        "sample_points method": "https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.sample_points.html",
    },
)
class RandomPointNode:
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False,
    )
    num_col = knext.ColumnParameter(
        "Number of points column",
        "Select the column for the number of points to draw.",
        column_filter=knut.is_long,
        include_row_key=False,
        include_none_column=False,
    )
    use_seed = knext.BoolParameter(
        "Use random seed",
        "I selected you may enter a fixed seed here to get reproducible results upon re-execution. "
        + "If you do not specify a seed, a new random seed is taken for each execution.",
        False,
    )
    seed = knext.IntParameter(
        "Seed",
        "A seed to initialize the random number generator.",
        1234,
    ).rule(knext.OneOf(use_seed, [True]), knext.Effect.SHOW)

    result_settings = knut.ResultSettings(
        mode=knut.ResultSettingsMode.APPEND.name,
        new_name="Points",
    )

    def configure(self, configure_context, input_schema_1):
        self.geo_col = knut.column_exists_or_preset(
            configure_context, self.geo_col, input_schema_1, knut.is_geo
        )

        knut.column_exists(self.num_col, input_schema_1, knut.is_long)

        return self.result_settings.get_result_schema(
            configure_context,
            input_schema_1,
            self.geo_col,
            knut.TYPE_GEO,
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.geo_col, exec_context)
        if self.use_seed:
            seed = self.seed
        else:
            seed = None
        if self.result_settings.mode == knut.ResultSettingsMode.APPEND.name:
            result_col = knut.get_unique_column_name(
                self.result_settings.new_column_name, input_1.schema
            )
        else:
            result_col = self.geo_col
        gdf[result_col] = gdf[self.geo_col].sample_points(
            size=gdf[self.num_col], method="uniform", seed=seed
        )
        return knut.to_table(gdf, exec_context)


############################################
# Directed Bezier Curve
############################################


@knext.node(
    name="Directed Bezier Curve",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/icon/GeometryTransformation/ODtoCurve.png",
    category=category,
    after="",
)
@knext.input_table(
    name="Geo table",
    description="Input table containing geometry columns representing origin and destination points.",
)
@knext.output_table(
    name="Table with directed Bézier curve",
    description="Output table with geometry columns representing Bézier curves connecting origin and destination points.",
)
@knut.geo_node_description(
    short_description="Generate Bézier curves between origin and destination points.",
    description="""This node generates a GeoSeries containing geometries representing smooth Bézier curves between origin and destination points. 
    The Bézier curves are created based on user-defined parameters,including the number of points, height scaling, and curve angle.
    This transformation is particularly useful for visualizing flows or movements in a more intuitive manner compared to straight lines.    
    """,
    references={
        "Bézier curve": "https://en.wikipedia.org/wiki/B%C3%A9zier_curve",
    },
)
class ODtoCurveNode:
    """
    This node generates geometries representing connections between origin and destination points as Bézier curves.
    """

    o_geo_col = knext.ColumnParameter(
        "Origin geopoint column",
        "Select the geometry column that describes the origins.",
        port_index=0,
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

    d_geo_col = knext.ColumnParameter(
        "Destination geopoint column",
        "Select the geometry column that describes the destination.",
        port_index=0,
        column_filter=knut.is_geo_point,
        include_row_key=False,
        include_none_column=False,
    )

    col_name = knext.StringParameter(
        "Column name for result geometry",
        "Specify the column name to define the Bézier curve.",
        "Curve",
    )

    num_points = knext.IntParameter(
        "Number of points ",
        "Specify the number of points to define the Bézier curve.",
        default_value=100,
        is_advanced=True,
    )

    height_scale = knext.DoubleParameter(
        "Height scale ",
        "Set the scale factor for the curve's height, controlling the distance of the curve from the straight line.",
        default_value=0.3,
        min_value=0,
        max_value=1,
        is_advanced=True,
    )

    angle_degrees = knext.IntParameter(
        "Curve angle (degrees) ",
        "Define the angle that controls the shape and direction of the Bézier curve.",
        default_value=60,
        min_value=0,
        max_value=360,
        is_advanced=True,
    )

    def configure(self, configure_context, input_schema):
        self.o_geo_col = knut.column_exists_or_preset(
            configure_context, self.o_geo_col, input_schema, knut.is_geo
        )
        self.d_geo_col = knut.column_exists_or_preset(
            configure_context, self.d_geo_col, input_schema, knut.is_geo
        )
        result_col = knut.get_unique_column_name(self.col_name, input_schema)
        return input_schema.append(
            knext.Column(
                knut.TYPE_LINE,
                result_col,
            )
        )

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = knut.load_geo_data_frame(input_1, self.o_geo_col, exec_context)

        import numpy as np
        from shapely.geometry import LineString

        def create_line_geometry(row):
            p0 = row[self.o_geo_col].coords[0]
            p2 = row[self.d_geo_col].coords[0]
            return LineString([p0, p2])

        def bezier_curve(p0, cp, p2, num_points=100):
            t_values = np.linspace(0, 1, num_points)
            curve = (
                np.outer((1 - t_values) ** 2, p0)
                + np.outer(2 * (1 - t_values) * t_values, cp)
                + np.outer(t_values**2, p2)
            )
            return curve

        def calculate_control_point(p0, p2, height_scale=0.3, angle_degrees=60):
            midpoint = (p0 + p2) / 2
            distance = np.linalg.norm(p2 - p0)
            offset_distance = distance * height_scale
            vec = p2 - p0
            angle_radians = np.radians(angle_degrees)
            rotation_matrix = np.array(
                [
                    [np.cos(angle_radians), np.sin(angle_radians)],
                    [-np.sin(angle_radians), np.cos(angle_radians)],
                ]
            )
            vec_perpendicular = np.dot(rotation_matrix, vec)
            vec_perpendicular /= np.linalg.norm(vec_perpendicular)
            vec_perpendicular *= offset_distance
            control_point = midpoint + vec_perpendicular
            return control_point

        def create_bezier_line_geometry(
            row,
            height_scale=self.height_scale,
            angle_degrees=self.angle_degrees,
            num_points=self.num_points,
        ):
            p0 = np.array(row[self.o_geo_col].coords[0])
            p2 = np.array(row[self.d_geo_col].coords[0])

            # Calculate the control point
            cp = calculate_control_point(
                p0, p2, height_scale=self.height_scale, angle_degrees=self.angle_degrees
            )

            # Calculate the points on the Bezier curve
            curve_points = bezier_curve(p0, cp, p2, num_points=self.num_points)

            # Convert the curve points to a LineString geometry
            line_geom = LineString(curve_points)

            return line_geom

        result_col = knut.get_unique_column_name(self.col_name, input_1.schema)

        gdf[result_col] = gdf.apply(create_bezier_line_geometry, axis=1)

        return knut.to_table(gdf)
