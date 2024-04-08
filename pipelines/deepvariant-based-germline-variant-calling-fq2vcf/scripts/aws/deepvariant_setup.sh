#!/bin/bash
set -e

echo "Step 1: Basic installation.." 
bash basic_setup.sh

echo "Step 2: Building applications.."
bash build_tools.sh

echo "Step 3: Building and saving Deepvaraint image.."
bash build_deepvariant_docker_image.sh

echo "Setup done!!"

echo "Next step, modify the \"config\" file according to the reference sequence and the read datasets, and run \"bash create_reference_index.sh \""







