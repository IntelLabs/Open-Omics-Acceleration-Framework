# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..","..")))

import socket
import time
from typing import Optional
from pydantic import BaseModel
from fastapi import HTTPException


import torch
import base64
from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry

from pathlib import Path
import logging
import typing as T
import json

import pandas as pd
import numpy as np


from data import transform_qm9, transform_zinc250k
from data.transform_zinc250k import zinc250_atomic_num_list
from mflow.models.hyperparams import Hyperparameters
from mflow.models.utils import check_validity, adj_to_smiles, check_novelty
from mflow.utils.model_utils import load_model
from mflow.models.model import rescale_adj
import mflow.utils.environment as env
from mflow.generate import generate_mols,visualize_interpolation_between_2_points,visualize_interpolation

# from IPython.display import SVG, display
from data.data_loader import NumpyTupleDataset

import functools
print = functools.partial(print, flush=True)
import subprocess



logger = logging.getLogger()
logger.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%y/%m/%d %H:%M:%S",
)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


PathLike = T.Union[str, Path]

from fastapi import HTTPException

def validate_client_inference_input(random_generation, reconstruct,intgrid,inter_times,n_experiments):
    try:
        # Check if all parameters are None or False/0
        if not any([random_generation, reconstruct, intgrid, inter_times, n_experiments]):
            logger.error(
                f"Invalid input received: random_generation={random_generation}, "
                f"reconstruct={reconstruct}, intgrid={intgrid}, "
                f"inter_times={inter_times}, n_experiments={n_experiments}"
            )
            raise HTTPException(
                status_code=400,
                detail=(
                    f"Missing required input: at least one of "
                    f"[random_generation, reconstruct, intgrid, inter_times, n_experiments] "
                    f"must be provided. Received: "
                    f"random_generation={random_generation}, "
                    f"reconstruct={reconstruct}, intgrid={intgrid}, "
                    f"inter_times={inter_times}, n_experiments={n_experiments}"
                )
            )

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in MoFlow microservice validation: {e}")
        raise HTTPException(status_code=500, detail=str(e))


WORKLOAD_NAME = "moflow_generate"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")

model_path = None
model = None
alphabet = None

class MoflowInput(BaseModel):
    batch_size:  Optional[int] = 100
    delta: Optional[float] = 0.1
    n_experiments:  Optional[int] = 1
    save_fig: Optional[bool] = False
    save_score: Optional[bool] = False
    temperature: Optional[float] = 1.0
    additive_transformations: Optional[bool] = False
    random_generation: Optional[bool] = False
    reconstruct: Optional[bool] = False
    int2point: Optional[bool] = False
    intgrid: Optional[bool] = False
    inter_times:  Optional[int] = 5
    correct_validity: Optional[bool] = False
    timing: Optional[bool] = False


class MoflowOutput(BaseModel):
    status: str
    result: Optional[str] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_MOFLOW")
class Opea_Moflow(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)

        logger.info("Loading model")
        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")
        self.model_dir = config["model_dir"]
        self.snapshot_path = config["snapshot_path"]
        self.hyperparams_path = config["hyperparams_path"]
        self.data_name = config["data_name"]
        self.data_dir = config["data_dir"]
        self.seed = config["seed"]
        self.debug = config["debug"]
        model_load_time = time.time()
        self.model = self._load_model(self.model_dir,self.snapshot_path,self.hyperparams_path,
                                      self.data_name,self.data_dir,self.seed,self.debug)
        print("Model load time",time.time() - model_load_time)

    def _load_model(self,model_dir,snapshot_path,hyperparams_path,
                                      data_name,data_dir,seed,debug):
        # chainer.config.train = False
        # Check if model_dir is empty or doesn't exist
        if not os.path.exists(model_dir) or not os.listdir(model_dir):
            print(f"{model_dir} is empty. Downloading model...")
            subprocess.run(['gdown', 'https://drive.google.com/uc?export=download&id=1GpXHrfP1vyzKu97aCReygLT7lg47pMi7'], check=True)
            results_dir = '/'+model_dir.split("/")[-2]
            os.makedirs(results_dir, exist_ok=True)
            subprocess.run(['unzip', 'results.zip', '-d', results_dir], check=True)
        snapshot_path = os.path.join(model_dir, snapshot_path)
        hyperparams_path = os.path.join(model_dir, hyperparams_path)
        print("loading hyperparamaters from {}".format(hyperparams_path))
        model_params = Hyperparameters(path=hyperparams_path)
        model = load_model(snapshot_path, model_params, debug=True)
        if len(model.ln_var) == 1:
            print('model.ln_var: {:.2f}'.format(model.ln_var.item()))
        elif len(model.ln_var) == 2:
            print('model.ln_var[0]: {:.2f}, model.ln_var[1]: {:.2f}'.format(model.ln_var[0].item(), model.ln_var[1].item()))
        gpu=-1
        noipex=True
        bf16=False
        if gpu >= 0:
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        model.to(self.device)
        model.eval()  # Set model for evaluation
        if not noipex:
            dtype = torch.bfloat16 if bf16 else torch.float32
            import intel_extension_for_pytorch as ipex
            model = ipex.optimize(model, dtype=dtype)
        if noipex and bf16:
            model.bfloat16()
        # true_data = NumpyTupleDataset.load(os.path.join(data_dir, molecule_file))
        self.enable_autocast = bf16
        self.device_type ="cuda" if gpu >=0 else "cpu"

        with torch.amp.autocast(device_type=self.device_type , enabled=self.enable_autocast):
            if data_name == 'qm9':
                self.atomic_num_list = [6, 7, 8, 9, 0]
                transform_fn = transform_qm9.transform_fn
                # true_data = TransformDataset(true_data, transform_qm9.transform_fn)
                file_path = 'data/valid_idx_qm9.json'
                valid_idx = transform_qm9.get_val_ids(file_path)
                # valid_idx = transform_qm9.get_val_ids()

                molecule_file = 'qm9_relgcn_kekulized_ggnp.npz'
            elif data_name == 'zinc250k':
                self.atomic_num_list = zinc250_atomic_num_list
                # transform_fn = transform_qm9.transform_fn
                transform_fn = transform_zinc250k.transform_fn_zinc250k
                # true_data = TransformDataset(true_data, transform_fn_zinc250k)
                file_path = 'data/valid_idx_zinc.json'
                valid_idx = transform_zinc250k.get_val_ids(file_path)
                molecule_file = 'zinc250k_relgcn_kekulized_ggnp.npz'
            file_path = os.path.join(data_dir, molecule_file)
            # Project root (adjust ".." depending on your script's depth)
            BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",".."))
            # Run preprocessing if the file doesn't exist
            if not os.path.exists(file_path):
                print(f"{molecule_file} not found. Running data_preprocess.py...")
                preprocessing_script_path = os.path.join(BASE_DIR, "data", "data_preprocess.py")

                if os.path.exists(preprocessing_script_path):
                    subprocess.run(['python', preprocessing_script_path, '--data_name', data_name,'--data_dir',data_dir], cwd=os.path.dirname(preprocessing_script_path))
                else:
                    print(f"Error: {preprocessing_script_path} not found."); exit(1)

        self.dataset = NumpyTupleDataset.load(os.path.join(data_dir,molecule_file), transform=transform_fn)

        assert len(valid_idx) > 0
        train_idx = [t for t in range(len(self.dataset)) if t not in valid_idx]  # 120803 = 133885-13082
        n_train = len(train_idx)  # 120803
        self.train = torch.utils.data.Subset(self.dataset, train_idx)  # 120803
        self.test = torch.utils.data.Subset(self.dataset, valid_idx)  # 13082  not used for generation
        return model


    def _run_inference(self,batch_size,delta,save_fig,save_score,n_experiments,temperature,additive_transformations,random_generation,reconstruct,int2point,intgrid,inter_times,correct_validity):

        batch_size = batch_size
        print('{} in total, {}  training data, {}  testing data, {} batchsize, train/batchsize {}'.format(
            len(self.dataset ),
            len(self.train),
            len(self.test),
            batch_size,
            len(self.train)/batch_size)
        )

        import time
        start = time.time()
        print("Start at Time: {}".format(time.ctime()))
        # 1. Reconstruction
        if reconstruct:

            train_dataloader = torch.utils.data.DataLoader(self.train, batch_size=batch_size)
            results = []
            reconstruction_rate_list = []
            max_iter = len(train_dataloader)
            for i, batch in enumerate(train_dataloader):
                x = batch[0].to(self.device)  # (256,9,5)
                adj = batch[1].to(self.device)  # (256,4,9, 9)
                adj_normalized = rescale_adj(adj).to(self.device)
                z, sum_log_det_jacs = self.model(adj, x, adj_normalized)
                z0 = z[0].reshape(z[0].shape[0], -1)
                z1 = z[1].reshape(z[1].shape[0], -1)
                adj_rev, x_rev = self.model.reverse(torch.cat([z0,z1], dim=1))
                # val_res = check_validity(adj_rev, x_rev, self.atomic_num_list)
                adj_rev=adj_rev.to(dtype=torch.float32)
                reverse_smiles = adj_to_smiles(adj_rev.cpu(), x_rev.cpu(), self.atomic_num_list)
                train_smiles = adj_to_smiles(adj.cpu(), x.cpu(), self.atomic_num_list)
                lb = np.array([int(a!=b) for a, b in zip(train_smiles, reverse_smiles)])
                idx = np.where(lb)[0]
                if len(idx) > 0:
                    for k in idx:
                        print(i*batch_size+k, 'train: ', train_smiles[k], ' reverse: ', reverse_smiles[k])
                reconstruction_rate = 1.0 - lb.mean()
                reconstruction_rate_list.append(reconstruction_rate)
                results.append({
                    "iter": i,
                    "total": max_iter,
                    "reconstruction_rate": reconstruction_rate
                })
                # novel_r, abs_novel_r = check_novelty(val_res['valid_smiles'], train_smiles, x.shape[0])
                print("iter/total: {}/{}, reconstruction_rate:{}".format(i, max_iter, reconstruction_rate))
            reconstruction_rate_total = np.array(reconstruction_rate_list).mean()
            results_summary = {
                        "reconstruction_rate_total": reconstruction_rate_total,
                        "num_samples": len(self.train)
                    }
            print("reconstruction_rate for all the train data:{} in {}".format(reconstruction_rate_total, len(self.train)))
            return {
                    "iterations": results,
                    "summary": results_summary
                }


        def _file_to_base64(file_path):
            with open(file_path, "rb") as f:
                return base64.b64encode(f.read()).decode("utf-8")
        ### interpolate 2 points
        if int2point:
            mol_smiles = None
            gen_dir = os.path.join(self.model_dir, 'generated')
            print('Dump figure in {}'.format(gen_dir))
            if not os.path.exists(gen_dir):
                os.makedirs(gen_dir)
            results = []
            for seed in range(inter_times):
                filepath = os.path.join(gen_dir, '2points_interpolation-2point_molecules_seed{}'.format(seed))
                visualize_interpolation_between_2_points(filepath, self.model, mol_smiles=mol_smiles, mols_per_row=15,
                                                        n_interpolation=50,
                                                        atomic_num_list=self.atomic_num_list, seed=seed, true_data=self.train,
                                                        device=self.device, data_name=self.data_name)
                # Save results as dict
                results.append({
                    "seed": seed,
                    "image_base64": _file_to_base64(filepath+".png"),
                    "pdf_base64": _file_to_base64(filepath+".pdf"),
                })

            return results

        ### interpolate in 2d grid
        if intgrid:
            mol_smiles = None
            gen_dir = os.path.join(self.model_dir, 'generated')
            print('Dump figure in {}'.format(gen_dir))
            if not os.path.exists(gen_dir):
                os.makedirs(gen_dir)
            results = []
            for seed in range(inter_times):
                filepath = os.path.join(gen_dir, 'generated_interpolation-grid_molecules_seed{}'.format(seed))
                print('saving {}'.format(filepath))
                visualize_interpolation(filepath, self.model, mol_smiles=mol_smiles, mols_per_row=9, delta=delta,
                                        atomic_num_list=self.atomic_num_list, seed=seed, true_data=self.train,
                                        device=self.device, data_name=self.data_name, keep_duplicate=True)

                filepath_unique = os.path.join(gen_dir, 'generated_interpolation-grid_molecules_seed{}_unique'.format(seed))
                visualize_interpolation(filepath_unique, self.model, mol_smiles=mol_smiles, mols_per_row=9, delta=delta,
                                        atomic_num_list=self.atomic_num_list, seed=seed, true_data=self.train,
                                        device=self.device, data_name=self.data_name, keep_duplicate=False)
                results.append({
                    "seed": seed,
                    "image_base64_1": _file_to_base64(filepath+".png"),
                    "pdf_base64_1": _file_to_base64(filepath+".pdf"),
                    "image_base64_2": _file_to_base64(filepath_unique+".png"),
                    "pdf_base64_2": _file_to_base64(filepath_unique+".pdf"),
                })
            return results

        # 3. Random generation
        if random_generation:

            train_x = [a[0] for a in self.train]
            train_adj = [a[1] for a in self.train]
            train_smiles = adj_to_smiles(train_adj, train_x, self.atomic_num_list)
            print('Load trained model and data done! Time {:.2f} seconds'.format(time.time() - start))

            # save_fig = save_fig
            valid_ratio = []
            unique_ratio = []
            novel_ratio = []
            abs_unique_ratio = []
            abs_novel_ratio = []

            results = []
            for i in range(n_experiments):
                # 1. Random generation
                adj, x = generate_mols(self.model, batch_size=batch_size, true_adj=None, temp=temperature,
                                    device=self.device)
                adj=adj.to(dtype=torch.float32)
                val_res = check_validity(adj, x, self.atomic_num_list, correct_validity=correct_validity)
                novel_r, abs_novel_r = check_novelty(val_res['valid_smiles'], train_smiles, x.shape[0])
                novel_ratio.append(novel_r)
                abs_novel_ratio.append(abs_novel_r)

                unique_ratio.append(val_res['unique_ratio'])
                abs_unique_ratio.append(val_res['abs_unique_ratio'])
                valid_ratio.append(val_res['valid_ratio'])
                n_valid = len(val_res['valid_mols'])

                if save_score:
                    assert len(val_res['valid_smiles']) == len(val_res['valid_mols'])
                    smiles_qed_plogp = [(sm, env.qed(mol), env.penalized_logp(mol))
                                        for sm, mol in zip(val_res['valid_smiles'], val_res['valid_mols'])]
                    smiles_qed_plogp.sort(key=lambda tup: tup[2], reverse=True)
                    gen_dir = os.path.join(self.model_dir, 'generated')
                    os.makedirs(gen_dir, exist_ok=True)
                    filepath = os.path.join(gen_dir, 'smiles_qed_plogp_{}_RankedByPlogp.csv'.format(i))
                    df = pd.DataFrame(smiles_qed_plogp, columns =['Smiles', 'QED', 'Penalized_logp'])
                    df.to_csv(filepath, index=None, header=True)

                    smiles_qed_plogp.sort(key=lambda tup: tup[1], reverse=True)
                    filepath2 = os.path.join(gen_dir, 'smiles_qed_plogp_{}_RankedByQED.csv'.format(i))
                    df2 = pd.DataFrame(smiles_qed_plogp, columns=['Smiles', 'QED', 'Penalized_logp'])
                    df2.to_csv(filepath2, index=None, header=True)
                    results.append({
                        "iteration": i,
                        "RankedByPlogp": df.to_dict(orient="records"),
                        "RankedByQED": df2.to_dict(orient="records")
                    })

                # saves a png image of all generated molecules
                if save_fig:
                    gen_dir = os.path.join(self.model_dir, 'generated')
                    os.makedirs(gen_dir, exist_ok=True)
                    filepath = os.path.join(gen_dir, 'generated_mols_{}.png'.format(i))
                    img = Draw.MolsToGridImage(val_res['valid_mols'], legends=val_res['valid_smiles'],
                                            molsPerRow=20, subImgSize=(300, 300))  # , useSVG=True
                    img.save(filepath)

            print("validity: mean={:.2f}%, sd={:.2f}%, vals={}".format(np.mean(valid_ratio), np.std(valid_ratio), valid_ratio))
            print("novelty: mean={:.2f}%, sd={:.2f}%, vals={}".format(np.mean(novel_ratio), np.std(novel_ratio), novel_ratio))
            print("uniqueness: mean={:.2f}%, sd={:.2f}%, vals={}".format(np.mean(unique_ratio), np.std(unique_ratio),
                                                                        unique_ratio))
            print("abs_novelty: mean={:.2f}%, sd={:.2f}%, vals={}".
                format(np.mean(abs_novel_ratio), np.std(abs_novel_ratio), abs_novel_ratio))
            print("abs_uniqueness: mean={:.2f}%, sd={:.2f}%, vals={}".
                format(np.mean(abs_unique_ratio), np.std(abs_unique_ratio),
                                                                        abs_unique_ratio))
            print('Task random generation done! Time {:.2f} seconds, Data: {}'.format(time.time() - start, time.ctime()))

            metrics = {
                "validity": {"mean": float(np.mean(valid_ratio)), "sd": float(np.std(valid_ratio)), "vals": valid_ratio},
                "novelty": {"mean": float(np.mean(novel_ratio)), "sd": float(np.std(novel_ratio)), "vals": novel_ratio},
                "uniqueness": {"mean": float(np.mean(unique_ratio)), "sd": float(np.std(unique_ratio)), "vals": unique_ratio},
                "abs_novelty": {"mean": float(np.mean(abs_novel_ratio)), "sd": float(np.std(abs_novel_ratio)), "vals": abs_novel_ratio},
                "abs_uniqueness": {"mean": float(np.mean(abs_unique_ratio)), "sd": float(np.std(abs_unique_ratio)), "vals": abs_unique_ratio},
                # "runtime": {"time_seconds": round(time.time() - start_time, 2), "finished_at": time.ctime()}
            }

            return {"results": results, "metrics": metrics}


    async def invoke(self, input: MoflowInput) -> MoflowOutput:

        try:
            start_time = time.time()

            validate_client_inference_input(input.random_generation,input.reconstruct,input.intgrid,input.inter_times,input.n_experiments)

            results = self._run_inference(input.batch_size,input.delta,input.save_fig,input.save_score,input.n_experiments,input.temperature,input.additive_transformations,input.random_generation,input.reconstruct,input.int2point,input.intgrid,input.inter_times,input.correct_validity)
            # assuming results is dict/list
            results_json = json.dumps(results)
            results_bytes = results_json.encode("utf-8")
            results_b64 = base64.b64encode(results_bytes).decode("utf-8")

            print("total run time",time.time() - start_time)
            return MoflowOutput(
                status="success",
                result=results_b64,
                message="Inference completed successfully."
            )

        except HTTPException as http_exc:
            raise http_exc



    def check_health(self) -> bool:
        """Checks the health of the Moflow service."""
        if self.model is None:
            logger.error("Moflow model is not initialized.")
            return False

        return True
