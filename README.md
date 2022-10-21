# Geospatial Analytics Extension for KNIME

This repository is the home of the Geospatial Analytics Extension for [KNIME Analytics Platform](https://www.knime.com/knime-analytics-platform). The extension is developed by the [Center for Geographic Analysis](https://gis.harvard.edu/) at [Harvard University](https://www.harvard.edu/) and [KNIME](https://www.knime.com/) as part of a two-year project of the [Spatiotemporal Innovation Center](https://www.stcenter.net/). The goal of the collaboration is to develop KNIME Analytics Platform extensions and best-practice workflows to provide a consistent and compatible platform for spatial data analysis across disciplines.

## Package Organization

* `knime_extension`: This folder contains all files of the KNIME Geospatial extension such as the source code of each node. The folder itself is structured as suggested in the Python Best Practices file.
* `docs`: Additional material to get you started with the development of this extension such as development setup instructions or best practices.
* `tests`: Test data and workflows used to test the node functionality.
* `config.yml`: Example `config.yml` file that should point to the `knime_extension` folder on your local hard drive during local development and debugging as described [here](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#tutorial-writing-first-py-node).


## Installation
The extension will be is installable via the [KNIME Hub](https://hub.knime.com/) once it is publicly available. 
If you want to test the latest version you can follow the instructions in the Setup section of the [Contribution guide](https://github.com/spatial-data-lab/knime-geospatial-extension/blob/main/CONTRIBUTING.md#setup).


## Contribute
We are very happy about every contributor. Please read the [Contributors' Guide](https://github.com/spatial-data-lab/knime-geospatial-extension/blob/main/CONTRIBUTING.md) for more details.


## License
The repository is released under the [MIT License](https://opensource.org/licenses/MIT).

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)