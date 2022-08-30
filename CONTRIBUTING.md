# How to contribute
Thanks a lot for taking the time to read this little guide and contributing to this project
Everybody is welcome to contribute to this Geospatial extension. The repository contains a Python based KNIME extension for geospatial analysis. To get started make yourself familiar with the general information on how to implement a Python based node extension for KNIME Analytics Platform that is available [here](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#introduction).

## Submitting changes

Please send a [GitHub Pull Request to KNIME Geospatial Extension](https://github.com/spatial-data-lab/knime-python-spatial-statistic-nodes/pull/new/master) with a clear list of what you've done (read more about [pull requests](http://help.github.com/pull-requests/)). When you send a pull request, please also include an example or testing workflow that demonstrates or tests the new functionality. We can always use more test coverage. Please follow our coding conventions (below) and make sure all of your commits are atomic (one feature per commit).


## Coding conventions
We use the [PEP-8 Style guide](https://peps.python.org/pep-0008/) as a guideline for our source code.


## Setup

This is a short step by step guide on how to setup your development environment to get started with the node development.

1. [Download](https://www.knime.com/nightly-build-downloads) and install the latest KNIME nightly build

2.  Install KNIME extension: Open the downloaded KNIME Analytics Platform nightly build and click KNIME -> File-> Install KNIME Extensions, Search for 'Python' and select `KNIME Python Extension Development (Labs)`,`KNIME Python Integration (Labs) `, Disable the “Group items by category” option and search for 'geo' and select  `KNIME Geospatial Nodes`.

3. Follow the instruction in the Tutorials section of the [Create a New Python based KNIME Extension guide](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#_tutorials). When building the new Python environment for Geospatial node development also install the `geopandas` package as following: 
```bash
conda create -n my_python_env python=3.9 knime-python-base knime-extension geopandas -c knime -c conda-forge 

conda activate my_python_env

conda install libpysal scipy # if you use these packages, please install here too

conda info
# Record the env location path such as D:\ProgramData\Anaconda3\envs\my_python_env 
```

4. Configure the following KNIME Python settings:
   * KNIME -> File-> Preference> KNIME> Conda  Choose anaconda directory
   * KNIME -> File-> Preference> KNIME> Python(Labs) Choose my_python_env

5. Clone this repository
```bash
git clone  https://github.com/spatial-data-lab/knime-python-spatial-statistic-nodes.git

```

6. Adapt the `config.yml` file to point to the `knime_extension` folder on your local hard drive as described [here](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#tutorial-writing-first-py-node). The file tells the KNIME Analytics Platform where to search for the Python extension on your local hard drive during development and debugging.

7. Add the following line to your `knime.ini`, located in the KNIME nightly root folder, to point to your adjusted `config.yml` file e.g. `-Dknime.python.extension.config=D:\Software\knime_470\knimespace\geo\config.yml`

8.  Reopen the KNIME Analytics Platform and you should see the Geospatial extension in your KNIME node repository.

9. If you do not see the extension have a look at the KNIME log file via View -> Open KNIME log. One common problem is that not all required Python packages are install in you Python environment (for example `libpysal` or `geopandas`)


## Testing
Every KNIME node that is part of this extension should be accompanied by a KNIME workflow that tests its functionality and also known corner cases. To write the workflows KNIME provides the [KNIME Testing Framework UI extension](https://kni.me/e/ufBEiCcvH9QIFePn) that contains several helper nodes e.g. the [Testflow Configuration node](https://kni.me/n/SrlKL_mJ63P7BVXh) and the [Table Difference Checker node](https://kni.me/n/dWyH_vs7JoIWPRsJ) to create test data and to compare results.

[Test guideline ](https://docs.google.com/document/d/1XrJFvqVCreyBGRcP9M-M9afq0l5gFGIMxk_SUpno1yc/edit?usp=sharing)


## Build
When testing locally you do not need to build anything. All you need to do is adjust the `config.yml` to your local setup e.g. adjust the `src` path to the location of the `knime_extension` folder on your hard drive and the `conda_env_path` path to the path of the local conda environment that contains all dependencies required by the extension e.g. geopandas and add the path to the adjusted `config.yml` file to your `knime.ini` file as described above.
However if you want to build the extension locally follow the instructions [here](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#extension-bundling).