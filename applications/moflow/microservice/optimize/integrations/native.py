# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import sys
import os
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..","..")))

import time
from typing import  Optional
from pydantic import BaseModel
from fastapi import HTTPException


import torch
import base64
from comps import  OpeaComponent, OpeaComponentRegistry

from pathlib import Path
import sys
import logging
import typing as T
from timeit import default_timer as timer
import json



from data.data_loader import NumpyTupleDataset
import torch
from data import transform_qm9, transform_zinc250k
from data.transform_zinc250k import zinc250_atomic_num_list
from mflow.models.hyperparams import Hyperparameters
from mflow.utils.model_utils import load_model

import time
import functools
print = functools.partial(print, flush=True)
import subprocess

from mflow.optimize_property import MoFlowProp as _MoFlowProp
import __main__

if not hasattr(__main__, "MoFlowProp"):
    setattr(__main__, "MoFlowProp", _MoFlowProp)


from mflow.optimize_property import optimize_mol,fit_model,constrain_optimization_smiles,find_top_score_smiles,load_property_csv

print("all import is done")


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


def validate_client_inference_input(topk, sim_cutoff, topscore, consopt):
    try:
        if topk is None:
            raise HTTPException(
                status_code=400,
                detail="Missing required input: topk"
            )
        if topscore is None and consopt is None:
            raise HTTPException(
                status_code=400,
                detail="Either 'topscore' or 'consopt' must be provided"
            )

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in Moflow microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))



class MoflowOptimizeInput(BaseModel):
    delta: Optional[float] = 0.1
    temperature: Optional[float] = 1.0
    additive_transformations: Optional[bool] = False
    topk:  Optional[int] = 5
    sim_cutoff: Optional[float] = 1.0
    topscore: Optional[bool] = False
    consopt: Optional[bool] = False

class MoflowOptimizeOutput(BaseModel):
    status: str
    result: Optional[str] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_MOFLOW_OPTIMIZE")
class Opea_Moflow_Optimize(OpeaComponent):

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
        self.property_model_path = config["property_model_path"]
        self.property_name=config["property_name"]
        self.lr_decay = config["lr_decay"]
        self.weight_decay = config["weight_decay"]
        self.hidden = config["hidden"]
        self.max_epochs = config["max_epochs"]
        self.debug = config["debug"]
        self.batch_size = config["batch_size"]
        model_load_time = time.time()
        self.model = self._load_model(self.model_dir,self.snapshot_path,self.hyperparams_path,
                                      self.data_name,self.data_dir,self.property_model_path,self.lr_decay,self.weight_decay,self.hidden,self.max_epochs,self.debug,self.property_name,self.batch_size)
        print("Model load time",time.time() - model_load_time)


    def _load_model(self,model_dir,snapshot_path,hyperparams_path,
                                      data_name,data_dir,property_model_path,lr_decay,weight_decay,hidden,max_epochs,debug,property_name,batch_size):

        # Device configuration
        start = time.time()
        print("Start at Time: {}".format(time.ctime()))
        device = -1
        gpu=-1
        noipex=True
        bf16=False
        if gpu >= 0:
            # device = gpu
            #device = torch.device('cuda:' + str(gpu) if torch.cuda.is_available() else 'cpu')
            device = torch.device('cuda')
        else:
            device = torch.device('cpu')

        property_name = property_name.lower()
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
        model_params = Hyperparameters(path=hyperparams_path)
        model = load_model(snapshot_path, model_params, debug=True)  # Load moflow model

        if hidden in ('', ','):
            hidden = []
        else:
            hidden = [int(d) for d in hidden.strip(',').split(',')]
        self.property_model = _MoFlowProp(model, hidden)
        # model.eval()  # Set model for evaluation

        if data_name == 'qm9':
            self.atomic_num_list = [6, 7, 8, 9, 0]
            transform_fn = transform_qm9.transform_fn
            file_path = 'data/valid_idx_qm9.json'
            valid_idx = transform_qm9.get_val_ids(file_path)
            molecule_file = 'qm9_relgcn_kekulized_ggnp.npz'
            # smile_cvs_to_property('qm9')
        elif data_name == 'zinc250k':
            self.atomic_num_list = zinc250_atomic_num_list
            transform_fn = transform_zinc250k.transform_fn_zinc250k
            file_path = 'data/valid_idx_zinc.json'
            valid_idx = transform_zinc250k.get_val_ids(file_path)
            molecule_file = 'zinc250k_relgcn_kekulized_ggnp.npz'
            # smile_cvs_to_property('zinc250k')
        else:
            raise ValueError("Wrong data_name{}".format(data_name))
        file_path = os.path.join(data_dir, molecule_file)
        # Run preprocessing if the file doesn't exist
        # Project root (adjust ".." depending on your script's depth)
        BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",".."))

        # Absolute path to data_preprocess.py


        # print("Preprocessing script:", preprocessing_script_path)
        if not os.path.exists(file_path):
            print(f"{molecule_file} not found. Running data_preprocess.py...")
            preprocessing_script_path = os.path.join(BASE_DIR, "data", "data_preprocess.py")

            if os.path.exists(preprocessing_script_path):
                subprocess.run(['python', preprocessing_script_path, '--data_name', data_name,'--data_dir',data_dir], cwd=os.path.dirname(preprocessing_script_path))
            else:
                print(f"Error: {preprocessing_script_path} not found."); exit(1)
        # dataset = NumpyTupleDataset(os.path.join(data_dir, molecule_file), transform=transform_fn)  # 133885
        dataset = NumpyTupleDataset.load(os.path.join(data_dir, molecule_file), transform=transform_fn)

        print('Load {} done, length: {}'.format(os.path.join(data_dir, molecule_file), len(dataset)))
        assert len(valid_idx) > 0
        train_idx = [t for t in range(len(dataset)) if t not in valid_idx]  # 224568 = 249455 - 24887
        n_train = len(train_idx)  # 120803 zinc: 224568
        train = torch.utils.data.Subset(dataset, train_idx)  # 120803
        test = torch.utils.data.Subset(dataset, valid_idx)  # 13082  not used for generation

        train_dataloader = torch.utils.data.DataLoader(train, batch_size=batch_size)

        # print("loading hyperparamaters from {}".format(hyperparams_path))

        if property_model_path is None:
            print("Training regression model over molecular embedding:")
            prop_list = load_property_csv(data_name, normalize=True)
            train_prop = [prop_list[i] for i in train_idx]
            test_prop = [prop_list[i] for i in valid_idx]
            print('Prepare data done! Time {:.2f} seconds'.format(time.time() - start))
            property_model_path = os.path.join(model_dir, '{}_model.pt'.format(property_name))
            property_model = fit_model(property_model, self.atomic_num_list, train_dataloader, train_prop, device,
                                    property_name=property_name, max_epochs=max_epochs,
                                    learning_rate=learning_rate, weight_decay=weight_decay)
            print("saving {} regression model to: {}".format(property_name, property_model_path))
            torch.save(property_model, property_model_path)
            print('Train and save model done! Time {:.2f} seconds'.format(time.time() - start))
        else:
            print("Loading trained regression model for optimization")
            prop_list = load_property_csv("",data_name, normalize=False)
            self.train_prop = [prop_list[i] for i in train_idx]
            test_prop = [prop_list[i] for i in valid_idx]
            print('Prepare data done! Time {:.2f} seconds'.format(time.time() - start))
            property_model_path = os.path.join(model_dir, property_model_path)

            print("loading {} regression model from: {}".format(property_name, property_model_path))
            self.device = torch.device('cpu')
            self.property_model = torch.load(property_model_path, map_location=device)
            print('Load model done! Time {:.2f} seconds'.format(time.time() - start))

            self.property_model.to(self.device )
            self.property_model.eval()

            model.to(self.device )
            model.eval()
            if not noipex:
                dtype = torch.bfloat16 if bf16 else torch.float32
                import intel_extension_for_pytorch as ipex
                model = ipex.optimize(model, dtype=dtype)
            if noipex and bf16:
                model.bfloat16()

        return model


    def _run_inference(self,delta,temperature,additive_transformations,topk,sim_cutoff,topscore,consopt):
        gpu=-1
        enable_autocast = self.use_bf16
        device_type ="cpu" if gpu >=0 else "cuda"
        infer_time= time.time()
        with torch.amp.autocast(device_type=device_type , enabled=enable_autocast):
            if topscore:
                print('Finding top score:')
                result= find_top_score_smiles(self.model,self.property_model ,self.device , self.data_name, self.property_name, self.train_prop, topk, self.atomic_num_list, self.debug)

            if consopt:
                print('Constrained optimization:')
                result=constrain_optimization_smiles(self.model,self.property_model, self.device, self.data_name, self.property_name, self.train_prop, topk,   # train_prop
                                        self.atomic_num_list, self.debug, sim_cutoff=sim_cutoff)

        print('Total Time {:.2f} seconds'.format(time.time() - infer_time))
        return result

    async def invoke(self, input: MoflowOptimizeInput) -> MoflowOptimizeOutput:

        try:
            start_time = time.time()

            validate_client_inference_input(input.topk,input.sim_cutoff,input.topscore,input.consopt)

            results = self._run_inference(input.delta,input.temperature,input.additive_transformations,input.topk,input.sim_cutoff,input.topscore,input.consopt)
            print("results.....",results)
            results_json = json.dumps(results)
            results_bytes = results_json.encode("utf-8")
            results_b64 = base64.b64encode(results_bytes).decode("utf-8")

            print("total run time",time.time() - start_time)
            return MoflowOptimizeOutput(
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
