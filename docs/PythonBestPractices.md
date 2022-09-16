# KNIME Python Development  Best Practices

## Project Setup

This section extends the existing [Create a New Python based KNIME Extension](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#extension-bundling) guide with useful information. 

### General suggestions
- Create one repository per Python node extension
- If the extension contains several nodes create a file for each node group e.g., IO, analysis, views

### Proposed File Structure
The following shows the proposed file structure for a KNIME Extension with several nodes grouped into different categories and written in Python:
```
my_extension
+-- icons
    +-- my_node.svg
+-- src
    +--  my_extension.py
    +-- nodes
        +-- category1.py
        +-- category2.py
    +-- util
        +-- general_helper.py
+-- test
    +-- test_my_extension.py
+-- knime.yml
+-- LICENSE.TXT
+-- my_conda_env.yml
```

## Working with file paths
### Importing Nodes
All nodes in the nodes folder need to be imported into the *my_extension.py* file like the following: 
```python
import nodes.category1 etc.
```

### Referencing Icons
When referencing icons the root path is *my_extension* (the folder of the _extension\_module_ parameter in the _knime.yml_ file). So to reference *node1\_icon.png* the path is *icons/node1_icon.png*.

## Local development
For local development create a _config.yml_ file in the parent folder of the  *my_extension* folder. In the file use the following lines:
```
my.group: # {group\_id}.{name} from the knime.yml_
  src: C:\DEV\Python\my\_git\_repo # Path to folder containing knime.yml file_
  conda_env_path: C:\Users\XYZ\.conda\envs\my\_python\_env # Path to the local Python env to use_
  debug_mode: true # Optional line, if set to true, it will always use the latest changes of execute/configure, when that method is used within the KNIME Analytics Platform_
```

The _config.yml_ file can contain list several repositories below each other.

To add the Python extension in  _my\_extension_ to your KNIME Analytics Platform add the following line to the _knime.ini_ file:

_-Dknime.python.extension.config=C:\DEV\Python\config.yml_

## Local Debugging

To debug your Python code in Visual Studio Code (VSCode) even though it is executed from KNIME, perform the following steps:

1. In the Python environment that you use for running the tests, install **debugpy** using
_pip install debugpy_ or Anaconda or similar
2. Open the Python file you want to edit, e.g., _my\_node__.py_ (by pressing _Ctrl+P_ and start typing the file name)
3. Inject the following lines of code in your Python file from where you want to start debugging e.g. in the _execute_ method of your node:
    ```python
    import debugpy
    import logging
    LOGGER = logging.getLogger(__name__)
    debugpy.listen(5678)
    LOGGER.error("Waiting for debugger attach")
    debugpy.wait_for_client()
    debugpy.breakpoint()
    ```
4. Start KNIME Analytics Platform and execute the Python node you want to debug
5. Once you see "Waiting for debugger attach" in the KNIME Analytics Platform console, go to VSCode and open the **Run & Debug** view (*Ctrl+Shift+D*), click  **Run and Debug** , choose  **Remote**  **Attach** , localhost (should be default), and Port 5678 (default).
 (If you do not have a Python file open while clicking here, it will not show the Python debugging options!)