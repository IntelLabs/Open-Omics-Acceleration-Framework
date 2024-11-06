export DATA_DIR=/home/ubuntu/data/singlecell/input/
export OUTPUT_DIR=/home/ubuntu/data/singlecell/output/
if [ -f $DATA_DIR/1M_brain_cells_10X.sparse.h5ad ]
then
	echo "input file 1M_brain_cells_10X.sparse.h5ad exists at $DATA_DIR, skipping the download"
else	
wget -P  $DATA_DIR https://rapids-single-cell-examples.s3.us-east-2.amazonaws.com/1M_brain_cells_10X.sparse.h5ad 
fi

docker run -it -p 8888:8888 -v $DATA_DIR:/data -it --cpuset-cpus="0-47" scanpy
