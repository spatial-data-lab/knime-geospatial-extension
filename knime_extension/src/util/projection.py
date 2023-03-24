import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut


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
    def get_distance_parameter():
        "Double parameter that stores the distance value. Usually named 'distance'."
        return knext.DoubleParameter(
            label="Distance",
            description="The buffer distance for geometry.",
            default_value=1000.0,
        )

    @staticmethod
    def get_unit_parameter():
        "Distance unit parameter. Usually named 'unit'"
        return knext.EnumParameter(
            label="Distance unit",
            description="Choose distance unit to use.",
            default_value=Distance.Unit.get_default().name,
            enum=Distance.Unit,
            since_version="1.1.0",
        )

    @staticmethod
    def get_keep_input_crs_parameter():
        "Boolean parameter that indicates if the input CRS should be retained. Usually named 'keep_input_crs'."
        return knext.BoolParameter(
            label="Keep input CRS",
            description="If checked the CRS of the input table is retained even if a re-projection was necessary "
            + "for the selected distance unit.",
            default_value=False,
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
            "Use the default projection and convert distances to meter.",
        )
        KILOMETER = (
            "Kilometer",
            "Use the default projection and convert distances to kilometer.",
        )
        MILES = (
            "Miles",
            "Use the default projection and convert distances to miles.",
        )
        DEGREE = (
            "Degree",
            "Use the epsg:4326 which uses degree as distance unit.",
        )

        @classmethod
        def get_default(cls):
            return cls.INPUT

    orig_crs = None

    def __init__(self, distance: float, unit: str, keep_orig_crs: bool):
        self.distance = distance
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
                new_crs = knut.DEFAULT_CRS
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
                un_projected_gdf = gdf.to_crs(knut.DEFAULT_CRS, inplace=False)
            else:
                un_projected_gdf = gdf

            # from OSMNX https://github.com/gboeing/osmnx/blob/main/osmnx/projection.py#L104
            # calculate approximate longitude of centroid of union of all geometries in gdf
            avg_lng = un_projected_gdf.representative_point().x.mean()
            import numpy as np

            # calculate UTM zone from avg longitude to define CRS to project to
            utm_zone = int(np.floor((avg_lng + 180) / 6) + 1)
            new_crs = f"+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
            # TODO: workaround until Java part supports proj strings AP-20328
            from pyproj import CRS  # For CRS Units check

            new_crs = CRS.from_user_input(new_crs).to_wkt()

        if new_crs is None:
            # this should not happen
            raise ValueError(f"Invalid distance unit: {unit}")

        knut.check_canceled(exec_context)
        exec_context.set_progress(0.5, "Projection to new CRS for distance computation")
        if in_place:
            gdf.to_crs(new_crs, inplace=True)
            return gdf
        else:
            return gdf.to_crs(new_crs, inplace=False)

    def get_distance(self) -> float:
        unit = self.unit
        if (
            unit == Distance.Unit.INPUT.name
            or unit == Distance.Unit.DEGREE.name
            or unit == Distance.Unit.METER.name
        ):
            return self.distance

        factor = None
        if unit == Distance.Unit.KILOMETER.name:
            factor = 1000

        if unit == Distance.Unit.MILES.name:
            factor = 1609.34

        if factor is None:
            # this should not happen
            raise ValueError(f"Invalid distance unit: {unit}")
        return self.distance * factor

    def post_processing(
        self, exec_context: knext.ExecutionContext, gdf: gp.GeoDataFrame
    ) -> gp.GeoDataFrame:
        if not self.keep_input_crs:
            return gdf
        knut.check_canceled(exec_context)
        exec_context.set_progress(0.8, "Projecting back to input CRS")
        gdf.to_crs(self.orig_crs, inplace=True)
        return gdf
