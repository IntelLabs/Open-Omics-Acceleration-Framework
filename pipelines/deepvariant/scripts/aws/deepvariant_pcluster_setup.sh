

echo "Step 1: Basic installation.." 
bash basic_setup.sh

echo "Step 2: Building applications.."
bash build_tools.sh

echo "Step 3: Building and saving Deepvaraint image.."
bash build_deepvariant_docker_image.sh

echo "Setup done!!"


----- Login node 


1. Allocate nodes

2. configure hostfile


1 and 2 are running inside the script
run pcluster_compute_nodes_setup.sh







