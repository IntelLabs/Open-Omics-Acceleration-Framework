#!/bin/bash

# Constants
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

# Parse optional proxy args from command-line
for arg in "$@"; do
  case $arg in
    --http_proxy=*)
      http_proxy="${arg#*=}"
      ;;
    --https_proxy=*)
      https_proxy="${arg#*=}"
      ;;
    --no_proxy=*)
      no_proxy="${arg#*=}"
      ;;
    *)
      echo "Unknown option: $arg"
      exit 1
      ;;
  esac
done

runtime=docker
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
    1) build_esm=true ;;
    2) build_esmfold=true ;;
    3) build_esm=true; build_esmfold=true ;;
    *) echo "Invalid option selected. Please choose 1, 2, or 3."; exit 1 ;;
esac

echo "build_esm = ${build_esm}"
echo "build_esmfold = ${build_esmfold}"

# Function to build Docker image
build_image() {
  local image_name="$1"
  local dockerfile="$2"
  local args=(--build-arg BASE_IMAGE=$BASE_IMAGE)

  [[ -n "$http_proxy" ]] && args+=(--build-arg http_proxy=$http_proxy)
  [[ -n "$https_proxy" ]] && args+=(--build-arg https_proxy=$https_proxy)
  [[ -n "$no_proxy" ]] && args+=(--build-arg no_proxy=$no_proxy)

  $runtime build "${args[@]}" -f "$dockerfile" -t "$image_name" .
}

# Build base image
if ! image_exists "$BASE_IMAGE"; then
  echo "Building base image..."
  build_image "$BASE_IMAGE" "Dockerfile.base"
fi

# Build esm image
if $build_esm && ! image_exists "$ESM_IMAGE"; then
  echo "Building image for esm..."
  build_image "$ESM_IMAGE" "Dockerfile.esm"
fi

# Build esm_fold image
if $build_esmfold && ! image_exists "$ESMFOLD_IMAGE"; then
  echo "Building image for esm_fold..."
  build_image "$ESMFOLD_IMAGE" "Dockerfile.esmfold"
fi

echo "Build process completed."
