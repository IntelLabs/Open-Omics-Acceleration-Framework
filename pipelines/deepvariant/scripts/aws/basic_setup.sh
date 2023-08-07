
# All basic dev tools for Ubuntu 22.04

sudo apt update 

#sudo apt -y upgrade 

sudo apt -y install make 

sudo apt -y install autoconf 

sudo apt -y install numactl 

sudo apt -y install build-essential 

sudo apt -y install zlib1g-dev 

sudo apt -y install libncurses5-dev 

sudo apt -y update 

#sudo apt -y upgrade 

sudo apt -y install libbz2-dev 

sudo apt -y install liblzma-dev 

sudo apt-get -qq -y update
sudo apt-get -qq -y install wget

# All dependencies for bcftools Docker
echo "Installing Docker"
sudo apt-get -qq -y install apt-transport-https ca-certificates curl gnupg-agent software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -


sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

sudo apt-get -qq -y update
sudo apt-get -qq -y install docker-ce
sudo systemctl start docker

sudo docker --version

echo "Running Docker installation hello world!! test"
sudo docker run hello-world

#echo "Creating and activating a conda environment"
#source setup_env.sh deepvaraint_env
