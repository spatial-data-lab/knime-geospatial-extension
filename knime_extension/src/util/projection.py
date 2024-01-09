import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


DEFAULT_CRS = "epsg:4326"
"""Default coordinate reference system."""

DEF_CRS_DESCRIPTION = """Enter the 
        [Coordinate reference system (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) to use. 
        The input field supports the following types:
        
        - An authority string (i.e. 'epsg:4326')
        - An EPSG code (i.e. 4326)
        - [CRS WKT string](https://www.ogc.org/standards/wkt-crs)

        Common [EPSG codes](https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset) that can be universally 
        used for mapping coordinates everywhere in the world are:

        - [epsg:4326 (WGS 84, Unit: degree)](https://epsg.io/4326): Latitude/longitude coordinate system based 
        on the Earth's center of mass;  Used by the Global Positioning System among others. 
        This is also the default projection that the geospatial nodes use if not otherwise specified.
        - [epsg:3857 (Unit: meter)](https://epsg.io/3857): Web Mercator projection used by many web-based mapping tools,
        including Google Maps and OpenStreetMap.
        
        There are EPSG codes for specific regions that provide a higher accuracy in the corresponding regions:

        - North America:
            - [epsg:4269 (NAD83, Unit: degree)](https://epsg.io/4269) 
            - [epsg:26918 (NAD83 18N, Unit: meter)](https://epsg.io/26918)
        - China
            - [epsg:4490 (CGCS2000, Unit: degree)](https://epsg.io/4490) 
            - [epsg:4479 (CGCS2000, Unit: meter)](https://epsg.io/4479)

        For a selection of projections that preserve different properties see this 
        [Wikipedia article.](https://en.wikipedia.org/wiki/Map_projection#Projections_by_preservation_of_a_metric_property)
        Once you have found the appropriate projection name or coordinate reference system you can search for its 
        EPSG code at [https://epsg.io/.](https://epsg.io/) To do so simply type the projection name into the search 
        field (e.g. [Pseudo-Mercator).](https://epsg.io/?q=Pseudo-Mercator) The result page will show you the EPSG 
        code that you can enter in this field (e.g. EPSG:3857) but also the distance unit e.g. meter or degree and 
        the area of use where the projection works best.

        If you are looking for a projection for a specific area, try out the 
        [Projection wizard page](https://projectionwizard.org/) which suggests projection with specific properties 
        for a defined area on the globe. Once you have found the appropriate projection simply click on the 
        [WKT](https://www.ogc.org/standards/wkt-crs) link next to the suggested projection name, 
        copy it to your clipboard and paste it into this field.
        """

DEFAULT_DISTANCE_COLUMN_NAME = "Distance"
"""Default column name that contains a distance."""


def is_geographic(crs):
    """
    This checks if the CRS is geographic. It will check if it has a geographic CRS in the sub CRS if it is a
    compound CRS and will check if the source CRS is geographic if it is a bound CRS.

    Returns
    bool:
        True if the CRS is in geographic (lon/lat) coordinates.
    """
    from pyproj import CRS  # For CRS Units check

    crs_input = CRS.from_user_input(crs)
    return crs_input.is_geographic


def is_projected(crs):
    """
    This checks if the CRS is projected. It will check if it has a projected CRS in the sub CRS if it is a compound
    CRS and will check if the source CRS is projected if it is a bound CRS.

    Returns
    bool:
        True if CRS is projected.
    """
    from pyproj import CRS  # For CRS Units check

    crs_input = CRS.from_user_input(crs)
    return crs_input.is_projected


class Distance:
    """
    Helper class that takes care of projection to different CRS if required by the selected distance unit.
    """

    @staticmethod
    def get_distance_parameter(
        label: str = "Distance",
        description: str = "The buffer distance for the input geometry.",
        default_value: float = 1000.0,
        min_distance: float = 0,
        max_distance: float = None,
    ):
        "Double parameter that stores the distance value. Usually named 'distance'."
        return knext.DoubleParameter(
            label=label,
            description=description,
            default_value=default_value,
            min_value=min_distance,
            max_value=max_distance,
        )

    @staticmethod
    def get_unit_parameter(
        label: str = "Distance unit",
        description: str = "Choose the distance unit to use.",
        default_value: str = "INPUT",
    ):
        "Distance unit parameter. Usually named 'unit'"
        return knext.EnumParameter(
            label=label,
            description=description,
            default_value=default_value,
            enum=Distance.Unit,
            since_version="1.1.0",
        )

    @staticmethod
    def get_keep_input_crs_parameter(
        label: str = "Keep CRS from input table",
        description: str = "If checked the CRS of the input table is retained even if a re-projection was necessary "
        + "for the selected distance unit.",
        default_value: bool = False,
    ):
        "Boolean parameter that indicates if the input CRS should be retained. Usually named 'keep_input_crs'."
        return knext.BoolParameter(
            label=label,
            description=description,
            default_value=default_value,
            since_version="1.1.0",
        )

    class Unit(knext.EnumParameterOptions):
        """
        Distance units supported by this class.
        """

        INPUT = (
            "Use unit from input CRS",
            "Do not re-project but use the unit of the input CRS.",
        )
        METER = (
            "Meter",
            """Use the default projection with meter as unit. The projection is the UTM CRS for the UTM zone in 
            which the centroid of the input geometries lies. This projection works well for most latitudes but may 
            not work for some extreme northern locations like Svalbard or far northern Norway in which case you
            might want to use your own projection prior using this node.""",
        )
        KILOMETER = (
            "Kilometer",
            """Use the default projection with kilometer as unit. The projection is the UTM CRS for the UTM zone in 
            which the centroid of the input geometries lies. This projection works well for most latitudes but may 
            not work for some extreme northern locations like Svalbard or far northern Norway in which case you
            might want to use your own projection prior using this node.""",
        )
        MILES = (
            "Miles",
            """Use the default projection with miles as unit. The projection is the UTM CRS for the UTM zone in 
            which the centroid of the input geometries lies. This projection works well for most latitudes but may 
            not work for some extreme northern locations like Svalbard or far northern Norway in which case you
            might want to use your own projection prior using this node.""",
        )
        DEGREE = (
            "Degree",
            "Use the [epsg:4326](https://epsg.io/4326) with degree as unit.",
        )

        @classmethod
        def get_default(cls):
            # update default value of get_unit_parameter() if this changes!
            return cls.INPUT

    orig_crs = None

    def __init__(self, unit: str, keep_orig_crs: bool):
        self.unit = unit
        self.keep_input_crs = keep_orig_crs

    def pre_processing(
        self,
        exec_context: knext.ExecutionContext,
        gdf: gp.GeoDataFrame,
        in_place: bool = False,
    ) -> gp.GeoDataFrame:
        self.orig_crs = gdf.crs
        unit = self.unit
        if unit == Distance.Unit.INPUT.name:
            return gdf

        new_crs = None
        if unit == Distance.Unit.DEGREE.name:
            if is_projected(gdf.crs):
                new_crs = DEFAULT_CRS
            else:
                return gdf
        if (
            unit == Distance.Unit.METER.name
            or unit == Distance.Unit.KILOMETER.name
            or unit == Distance.Unit.MILES.name
        ):
            if is_projected(gdf.crs):
                knut.check_canceled(exec_context)
                exec_context.set_progress(0.4, "Preparing projection to new CRS")
                un_projected_gdf = gdf.to_crs(DEFAULT_CRS, inplace=False)
            else:
                un_projected_gdf = gdf

            # from OSMNX https://github.com/gboeing/osmnx/blob/main/osmnx/projection.py#L104
            # calculate approximate longitude of centroid of union of all geometries in gdf
            avg_lng = un_projected_gdf.representative_point().x.mean()
            import numpy as np

            # calculate UTM zone from avg longitude to define CRS to project to
            utm_zone = int(np.floor((avg_lng + 180) / 6) + 1)
            new_crs_proj = f"+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

            from pyproj import CRS  # For CRS Units check

            new_crs = CRS.from_user_input(new_crs_proj).to_epsg()
            if new_crs is None:
                # fallback to WKT which shouldn't be required
                new_crs = CRS.from_user_input(new_crs_proj).to_wkt()

        if new_crs is None:
            # this should not happen
            raise ValueError(f"Invalid distance unit: {unit}")

        knut.check_canceled(exec_context)
        exec_context.set_progress(0.3, "Projection to new CRS for distance computation")
        if in_place:
            gdf.to_crs(new_crs, inplace=True)
            return gdf
        else:
            return gdf.to_crs(new_crs, inplace=False)

    def get_distance_factor(self) -> float:
        "Returns the factor to use for the given unit."
        unit = self.unit
        factor = None
        if (
            unit == Distance.Unit.INPUT.name
            or unit == Distance.Unit.DEGREE.name
            or unit == Distance.Unit.METER.name
        ):
            factor = 1
        elif unit == Distance.Unit.KILOMETER.name:
            factor = 1000
        elif unit == Distance.Unit.MILES.name:
            factor = 1609.34
        else:
            raise ValueError(f"Invalid distance unit: {unit}")
        return factor

    def convert_input_distance(self, distance: float) -> float:
        "Returns the given distance converted to the given unit."
        return distance * self.get_distance_factor()

    def convert_result_distance(self, distance: float) -> float:
        return distance / self.get_distance_factor()

    def post_processing(
        self,
        exec_context: knext.ExecutionContext,
        gdf: gp.GeoDataFrame,
        in_place: bool = False,
    ) -> gp.GeoDataFrame:
        if not self.keep_input_crs:
            return gdf
        knut.check_canceled(exec_context)
        exec_context.set_progress(0.9, "Projecting back to input CRS")
        if in_place:
            gdf.to_crs(self.orig_crs, inplace=True)
            return gdf
        else:
            return gdf.to_crs(self.orig_crs, inplace=False)


def string_distances_parser(value: str, separator: str = ","):
    """Parses the value as a comma separated list of double values."""
    try:
        distances = [float(i) for i in value.split(separator)]
        return distances
    except ValueError:
        raise ValueError(
            f"'{value}' is not a valid buffer distance. Enter numbers separated by a comma."
        )
