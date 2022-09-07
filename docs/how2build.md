
# How to build and share the extension


See [here](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#_accessing_flow_variables)


### Setup

To ensure that the users you have shared your extension with are able to utilise its functionality fully and error-free, we bundle the source files together with the required packages using `conda` as the bundling channel.

The `knime.yml` file (refer to the [Python Node Extension Setup](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#python-node-extension-setup) section for an example of this configuration file) contains the information required to bundle your extension, including:

-   `extension_module`: the name of the `.py` file containing the node definitions of your extension.
    
-   `env_yml_path`: the path to the `.yml` file containing the configuration of the `conda` environment that is used with your extension.
    

The YAML file containing the `conda` environment definition needs to contain all the dependencies necessary for your nodes to work.

<table><tbody><tr><td><i title=""></i></td><td>Since the bundled extension needs to be operational on all operating systems supported by KNIME Analytics Platform, it is important to <em>not</em> strictly enforce package versions (unless you really have to) when generating the YAML file for the environment (see <a href="https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#environment-yml"><code>environment.yml</code></a> for an example). You can generate a YAML file for the environment by activating the environment with the <code>conda activate &lt;env_name&gt;</code> command, and then running the <code>conda env export --from-history &gt; &lt;env_yml_filename.yml&gt;</code> command.</td></tr></tbody></table>

#### `environment.yml`:

```
name: knime-python-scripting
channels:
- conda-forge
- knime
dependencies:
- python=3.9        # base dependency
- knime-python-base # base dependency
- knime-extension   # base dependency
...
```

Lastly, a new extension needs a `LICENSE.TXT` that will be displayed during the installation process.


### Option 1: Bundling a Python extension to share a zipped update site

Once you have finished implementing your Python extension, you can bundle it, together with the appropriate `conda` environment, into a local update site. This allows other users to install your extension in the KNIME Analytics Platform.

Follow the steps of [`extension setup`](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#extension-setup). Once you have prepared the YAML configuration file for the environment used by your extension, and have set up the `knime.yml` file, you can proceed to generating the local update site.

We provide a special `conda` package, `knime-extension-bundling`, which contains the necessary tools to automatically build your extension. Run the following commands in your terminal (Linux/macOS) or Anaconda Prompt (Windows). They will setup a `conda` environment, which gives the tools to bundle extensions. Then the extension will be bundled.

1.  Create a fresh environment prepopulated with the `knime-extension-bundling` package:
    
    ```
    conda create -n knime-ext-bundling -c knime -c conda-forge knime-extension-bundling
    ```
    
2.  Activate the environment:
    
    ```
    conda activate knime-ext-bundling
    ```
    
3.  With the environment activated, run the following command to bundle your Python extension:
    
    -   macOS/Linux:
        
        ```
        build_python_extension.py <path/to/directoryof/myextension/> <path/to/directoryof/output>
        ```
        
    -   Windows:
        
        ```
        build_python_extension.bat <path/to/directoryof/myextension/> <path/to/directoryof/output>
        ```
        
        where `<path/to/directoryof/myextension/>` is the path to the directory containing your `.py` extension module and the `knime.yml` file, and `<path/to/directoryof/output>` is the path to the directory where the bundled extension **repository** will be stored.
        
        <table><tbody><tr><td><i title=""></i></td><td>The bundling process can take several minutes to complete.</td></tr></tbody></table>
        
    
4.  Add the generated **repository** folder to KNIME AP as a Software Site in _File → Preferences → Install/Update → Available Software Sites_
    
5.  Install it via _File → Install KNIME Extensions_
    

The generated repository can now be shared with and installed by other users.




### Option 2: Publish your extension on KNIME Hub

Once you have finished implementing your Python extension, you can share it, together with the appropriate `conda` environment, to KNIME Hub.

#### Provide the extension

Follow the steps of [`extension setup`](https://docs.knime.com/latest/pure_python_node_extensions_guide/index.html#extension-setup) to prepare the `environment.yml` or some other `yml` defining your Python environment and the `knime.yml`.

Upload your extension into a Git repository, where the `knime.yml` is found top-level. A `config.yml` is not needed.

Some recommended project structure:

```
https://github.com/user/my_knime_extension
├── icons
│   └── my_node_icon.png
├── knime.yml
├── LICENSE.txt
├── environment.yml
└── my_extension.py
```