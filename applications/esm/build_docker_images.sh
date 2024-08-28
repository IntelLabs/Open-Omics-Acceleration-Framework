
#!/bin/bash

#constants
BASE_IMAGE=esm_base_image
ESM_IMAGE=esm_image
ESMFOLD_IMAGE=esmfold_image

# Function to check if a command exists
command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# Function to check if a Docker image exists
image_exists() {
  local image_name="$1"
  $runtime images -q "$image_name" | grep -q .
}


# Prompt for container runtime
echo "Select container runtime (docker or podman):"
read -r runtime

# Check if the selected runtime is available
if [[ "$runtime" != "docker" && "$runtime" != "podman" ]]; then
  echo "Invalid container runtime selected. Please choose either 'docker' or 'podman'."
  exit 1
fi

if ! command_exists "$runtime"; then
  echo "$runtime is not installed on your system. Please install it first."
  exit 1
fi

build_esm=false
build_esmfold=false

# Prompt for tasks to build
echo "Which images do you want to build?"
echo "1 esm"
echo "2 esm_fold"
echo "3 Both esm and esm_fold"
read -r task_option

case $task_option in
    1)

        build_esm=true
        ;;
    2)
        build_esmfold=true
        ;;
    3)
        build_esm=true
        build_esmfold=true
        ;;
    *)
        echo "Invalid option selected. Please choose 1, 2, or 3."
        exit 1
        ;;
esac

echo "build_esm = ${build_esm}"
echo "build_esmfold = ${build_esmfold}"

if ! image_exists "$BASE_IMAGE"; then
    echo "Building base image..."
    $runtime build -f Dockerfile.base -t $BASE_IMAGE .
fi

# Build image for esm
if $build_esm && ! image_exists "$ESM_IMAGE"; then
    echo "Building image for esm..."
    $runtime build --build-arg BASE_IMAGE=$BASE_IMAGE -f Dockerfile.esm -t $ESM_IMAGE .
fi

# Build image for esm_fold
if $build_esmfold && ! image_exists "$ESMFOLD_IMAGE"; then
    echo "Building image for esm_fold..."
    $runtime build --build-arg BASE_IMAGE=$BASE_IMAGE -f Dockerfile.esmfold -t $ESMFOLD_IMAGE .
fi

echo "Build process completed."
