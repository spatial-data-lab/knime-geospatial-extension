import logging

import geopandas as gp
import knime_extension as knext
import util.knime_utils as knut

LOGGER = logging.getLogger(__name__)


category = knext.category(
    path="/geo",
    level_id="transform",
    name="Project & Transformation",
    description="Geospatial transformation nodes",
    # starting at the root folder of the extension_module parameter in the knime.yml file
    icon="icons/projection.png",
)


@knext.node(name="CRS-Transformer", node_type=knext.NodeType.MANIPULATOR, icon_path="icons/icon.png", category=category)
@knext.input_table(name="Geo table", description="Table with geometry column to transform")
@knext.output_table(name="Transformed geo table", description="Transformed Geo input table")
class CrsTransformerNode:
    """
    This node projects the data from its original CRS to the entered CRS.
    """
    geo_col = knext.ColumnParameter(
        "Geometry column",
        "Select the geometry column to transform.",
        # Allow only GeoValue compatible columns
        column_filter=knut.is_geo,
        include_row_key=False,
        include_none_column=False
    )

    new_crs = knext.StringParameter("New CRS", "The new CRS system to use", "EPSG:4326")    
    @new_crs.validator
    def crs_validator(value):
        if not value:
            raise ValueError("New CRS must not be empty")

    def configure(self, configure_context, input_schema_1):
        knut.column_exists(self.geo_col, input_schema_1)
        return input_schema_1

    def execute(self, exec_context: knext.ExecutionContext, input_1):
        gdf = gp.GeoDataFrame(input_1.to_pandas(), geometry=self.geo_col)
        exec_context.set_progress(
            0.3, "Geo data frame loaded. Starting projection...")
        gdf = gdf.to_crs(self.new_crs)
        exec_context.set_progress(1.0, "Projection done")
        LOGGER.debug("CRS converted to " + self.new_crs)
        return knext.Table.from_pandas(gdf)
