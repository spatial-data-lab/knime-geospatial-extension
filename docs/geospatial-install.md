# Setup Custom Executor Docker Image with the Geospatial Extension for KNIME Business Hub

Running a workflow on KNIME Business Hub requires an Execution Context (EC). An EC runs one or multiple pods. A pod holds Docker containers in turn, which run Docker images in the end. In general, one can simply use the official Docker images from the public [KNIME Artifact Registry](https://registry.hub.knime.com/) to create an EC.

However, the official Docker images hold KNIME Extensions only; no Community (trusted or experimental) or Partner Extensions included. Since the "Geospatial Analytics Extension for KNIME" is a Trusted Community Extension, it needs to be integrated into the Docker image, which requires to enhance the official image and build your own custom Docker image.

## Requirements
To be able to build your custom Docker image it is required to run [Docker](https://www.docker.com/), more precisely the [Docker Engine](https://docs.docker.com/engine/install/).
Please make yourself familiar with the install instructions according to your operating system.
In general, you can build the image locally. It is not required that the image is built somewhere on the KNIME Business Hub.

## Dockerfile
The basis of a Docker image is a Dockerfile. It defines
- the source image (in this case the official Docker image from KNIME; check out all available images [here](https://docs.knime.com/latest/business_hub_admin_guide/index.html#docker-executor-images)),
- the libraries/packages to be installed, and subsequently 
- which extensions should be installed into the executor.

Make sure to save that file with the name ```Dockerfile``` and add no file extension.

Within the Dockerfile, we define two variables at the bottom that are important while installing the extension(s):

```KNIME_UPDATE_SITES``` - the URL to the necessary Update Site(s)

```KNIME_FEATURES``` - the Feature ID of the extension to be installed

Check out [this video](https://www.youtube.com/watch?v=-dO79Id3VAo&t=143s) from the KNIMETV channel that helps you to find the right update site URL and the Feature ID.

    # Define the base image
    FROM registry.hub.knime.com/knime/knime-full:r-5.2.1-369

    # Change to root user to be able to install system packages
    USER root

    # Update/upgrade package manager and install ca-certificates to enable ca certificates that micromamba (for python) is asking for
    RUN apt-get update && \
        apt-get upgrade -yq && \
        apt-get install -yq \
            ca-certificates && \
        # cleanup
        rm -rf /var/lib/apt/lists/*


    # Change to knime user to handle extensions
    USER knime

    # Define the list of update sites and features
    ENV KNIME_UPDATE_SITES="https://update.knime.com/analytics-platform/5.2,https://update.knime.com/community-contributions/trusted/5.2"

    # Install a feature from the Community Trusted update site
    ENV KNIME_FEATURES=sdl.harvard.features.geospatial.feature.group

    # Execute extension installation script
    RUN ./install-extensions.sh

## Build and push the Docker image
Now, the Docker image can be built from the Dockerfile and made available to the KNIME Business Hub with the following commands. Please make sure to adapt the following lines accordingly, when running them in the CLI of your laptop. 

    # 1) Build the Docker image
    docker build . -f .\Dockerfile -t knime-full-geo:5.2.1 --no-cache

    # 2) Retag the image to make it readable for KNIME Business Hub
    docker tag knime-full-geo:5.2.1 registry.<base-url>/knime-full-geo:5.2.1

    # 3) Login to the embedded registry of the KNIME Business Hub (if configured)
    # Authenticate with the registry:
    # You will need to provide a password. In case you do not have one you can manage this from the KOTS Admin Console > Config > Embedded Registry.
    docker login --username <username> registry.<base-url>

    # 4) Push the image to registry
    docker push registry.<base-url>/knime-full-geo:5.2.1

In case, you are not using the embedded registry but an external one, you can configure your KNIME Business Hub to use that within the KOTS Admin Console. Using an external docker registry, requires adapting the steps above accordingly.

![registry-settings](./imgs/registry-settings.png)

Please also refer to the official documentation on [Connecting to an External Registry](https://docs.replicated.com/vendor/packaging-private-images).

## Check if the Docker image is available on the registry
Once the upload is finished you can double check if that image can be found on the registry by searching the Docker registry API. Therefore, it is basically enough to use your web browser.

    # Check which image names are available
    $ registry.<hub-url>/v2/_catalog
    {
        "repositories":["knime-full-geo"]
    }

    # Check which tags are available for a dedicated image name
    $ registry.<hub-url>/v2/<image-name>/tags/list
    {
        "name":"knime-full-geo",
        "tags":["5.2.1"]
    }

## Creating an Execution Context
Now you can apply the Docker image and create an EC from it by passing the registry string under the "Docker image" section while creating an EC. Please read the [KNIME docs](https://docs.knime.com/latest/business_hub_user_guide/index.html#execution_contexts) to find out more about how to setup a new execution context.

![create-execution-context](./imgs/create-execution-context.png)
