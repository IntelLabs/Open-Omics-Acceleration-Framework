import os
import sys
# for linux env.
sys.path.insert(0,'..')
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import time
from data.data_frame_parser import DataFrameParser
from data.data_loader import NumpyTupleDataset
from data.smile_to_graph import GGNNPreprocessor


def parse():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data_name', type=str, default='qm9',
                        choices=['qm9', 'zinc250k'],
                        help='dataset to be downloaded')
    parser.add_argument("--data_dir", type=str, default='../data')
    parser.add_argument('--data_type', type=str, default='relgcn',
                        choices=['gcn', 'relgcn'],)
    args = parser.parse_args()
    return args

def ensure_directory_exists(data_dir):
    data_dir = args.data_dir.strip()

    # Validate if the path is provided
    if not data_dir:
        raise ValueError("Directory path must be provided")

    # Check if the path exists
    if not os.path.exists(data_dir):
        print(f"The directory {data_dir} does not exist. Creating it now...")
        os.makedirs(data_dir, exist_ok=True)  # Creates the directory, including parents if necessary
        print(f"Directory {data_dir} created.")

    # Ensure it's a directory, not a file
    if not os.path.isdir(data_dir):
        raise ValueError(f"The path {data_dir} is not a directory.")

    # Check for read and write permissions
    if not os.access(data_dir, os.R_OK):    
        raise ValueError(f"Read permission denied for {data_dir}.")
    if not os.access(data_dir, os.W_OK):
        raise ValueError(f"Write permission denied for {data_dir}.")

    # Normalize the path
    data_dir = os.path.normpath(data_dir)

    return data_dir

start_time = time.time()
args = parse()
data_name = args.data_name
data_type = args.data_type
print('args', vars(args))

if data_name == 'qm9':
    max_atoms = 9
elif data_name == 'zinc250k':
    max_atoms = 38
else:
    raise ValueError("[ERROR] Unexpected value data_name={}".format(data_name))


if data_type == 'relgcn':
    # preprocessor = GGNNPreprocessor(out_size=max_atoms, kekulize=True, return_is_real_node=False)
    preprocessor = GGNNPreprocessor(out_size=max_atoms, kekulize=True)
else:
    raise ValueError("[ERROR] Unexpected value data_type={}".format(data_type))

data_dir = ensure_directory_exists(args.data_dir)

if data_name == 'qm9':
    print('Preprocessing qm9 data:')
    df_qm9 = pd.read_csv('qm9.csv', index_col=0)
    labels = ['A', 'B', 'C', 'mu', 'alpha', 'homo', 'lumo', 'gap', 'r2',
              'zpve', 'U0', 'U', 'H', 'G', 'Cv']
    parser = DataFrameParser(preprocessor, labels=labels, smiles_col='SMILES1')
    result = parser.parse(df_qm9, return_smiles=True)
    dataset = result['dataset']
    smiles = result['smiles']
elif data_name == 'zinc250k':
    print('Preprocessing zinc250k data')
    # dataset = datasets.get_zinc250k(preprocessor)
    df_zinc250k = pd.read_csv('zinc250k.csv', index_col=0)
    # Caution: Not reasonable but used in used in chain_chemistry\datasets\zinc.py:
    # 'smiles' column contains '\n', need to remove it.
    # Here we do not remove \n, because it represents atom N with single bond
    labels = ['logP', 'qed', 'SAS']
    parser = DataFrameParser(preprocessor, labels=labels, smiles_col='smiles')
    result = parser.parse(df_zinc250k, return_smiles=True)
    dataset = result['dataset']
    smiles = result['smiles']
else:
    raise ValueError("[ERROR] Unexpected value data_name={}".format(data_name))

NumpyTupleDataset.save(os.path.join(data_dir, '{}_{}_kekulized_ggnp.npz'.format(data_name, data_type)), dataset)
print('Total time:', time.strftime("%H:%M:%S", time.gmtime(time.time()-start_time)) )
