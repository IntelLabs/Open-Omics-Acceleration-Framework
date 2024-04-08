#!/bin/bash
set -e



#Pre-req: Docker installation
# Check if Docker is installed
if [ "$(command -v docker)" ]; then
    # Docker is installed
    echo "Docker is installed."
    
    # Check Docker version
    docker_version=$("docker" --version | awk '{print $3}')
    echo "Docker version: $docker_version"
else
    # Docker is not installed
    echo "Docker is not installed on this system."
	exit 1
fi


# Build docker

# This will save deepvariant images
cd ../../../../applications/deepvariant
sudo docker build -t deepvariant .

# check the built and print the image ID

sudo docker images | grep "deepvariant:latest"

#save image(~7 GB) to tar file if you are using multiple nodes.

echo "Saving deepvariant:latest image as deepvariant.tar..."
cd - # Move to pipelines/deepvariant
sudo docker save -o deepvariant.tar deepvariant:latest
