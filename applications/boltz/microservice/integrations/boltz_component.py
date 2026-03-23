import os
import shutil
import uuid
import requests
from pathlib import Path
from typing import Optional, Tuple
from comps import OpeaComponent, OpeaComponentRegistry, CustomLogger
from pydantic import BaseModel, root_validator
import torch
from dataclasses import asdict
from pytorch_lightning import Trainer

# Boltz imports
from boltz.main import (
    download_boltz2, process_inputs, Boltz2, Boltz2DiffusionParams,
    PairformerArgsV2, MSAModuleArgs, BoltzSteeringParams
)
from boltz.data.module.inferencev2 import Boltz2InferenceDataModule
from boltz.data.write.writer import BoltzWriter
from boltz.data.types import Manifest

logger = CustomLogger("Opea_Boltz_Component")

class BoltzInput(BaseModel):
    yaml_content: Optional[str] = None
    protein_sequence: Optional[str] = None
    ligand_smiles: Optional[str] = None

    @root_validator(pre=True)
    def check_inputs_provided(cls, values):
        if values.get('yaml_content') is None and values.get('protein_sequence') is None:
            raise ValueError('Either "yaml_content" or "protein_sequence" must be provided.')
        return values

class BoltzOutput(BaseModel):
    status: str
    result: Optional[str] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_BOLTZ")
class Opea_Boltz(OpeaComponent):
    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)
        
        # 1. Setup Dirs
        self.cache_dir = Path(os.getenv("BOLTZ_CACHE", "~/.boltz")).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.output_base_dir = Path(os.getenv("BOLTZ_OUTPUT_DIR", "boltz_all_outputs")).resolve()
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        
        # 2. Download Weights
        logger.info("Downloading/Checking Boltz-2 weights...")
        download_boltz2(self.cache_dir)
        
        torch.set_grad_enabled(False)
        torch.set_float32_matmul_precision("highest")

        # 3. Load Model
        checkpoint_path = self.cache_dir / "boltz2_conf.ckpt"
        
        predict_args = {
            "recycling_steps": 3, "sampling_steps": 200, "diffusion_samples": 1, 
            "max_parallel_samples": 1, "write_confidence_summary": True, 
            "write_full_pae": False, "write_full_pde": False
        }

        logger.info(f"Loading model from {checkpoint_path}...")
        
        self.model_module = Boltz2.load_from_checkpoint(
            checkpoint_path, 
            strict=False,            
            predict_affinity=False,  # Try to set via init
            predict_args=predict_args, 
            map_location="cpu",
            diffusion_process_args=asdict(Boltz2DiffusionParams(step_scale=1.5)), 
            ema=False,
            pairformer_args=asdict(PairformerArgsV2()), 
            msa_args=asdict(MSAModuleArgs(use_paired_feature=True)),
            steering_args=asdict(BoltzSteeringParams())
        )
        
        # --- CRITICAL FIX: FORCE DISABLE AFFINITY ---
        # Some versions of Lightning/Boltz might ignore the init kwarg if reloading from ckpt params.
        self.model_module.predict_affinity = False
        # --------------------------------------------

        self.model_module.eval()
        
        self.accelerator = "gpu" if torch.cuda.is_available() else "cpu"
        if self.accelerator == "gpu": self.model_module.to("cuda")
        logger.info(f"Model loaded. Accelerator: {self.accelerator}")

    async def invoke(self, input: BoltzInput) -> Tuple[Path, Path]:
        request_id = str(uuid.uuid4())
        temp_dir = self.output_base_dir / request_id
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Processing request {request_id} in {temp_dir}")

        try:
            boltz_prediction_dir = temp_dir / "boltz_prediction"
            boltz_prediction_dir.mkdir(exist_ok=True)
            yaml_file = temp_dir / "input.yaml"

            # Create Input File
            if input.yaml_content:
                with open(yaml_file, "w") as f: f.write(input.yaml_content)
            else:
                with open(temp_dir / "protein.fasta", "w") as f: 
                    f.write(f">protein\n{input.protein_sequence}\n")
                with open(temp_dir / "ligand.smi", "w") as f: 
                    f.write(f"{input.ligand_smiles or 'C'}\tligand_name\n")
                with open(yaml_file, "w") as f:
                    f.write(f"sequences:\n- protein:\n    file: {temp_dir}/protein.fasta\n- ligand:\n    file: {temp_dir}/ligand.smi\n")

            # Process Inputs
            process_inputs(
                data=[yaml_file], out_dir=boltz_prediction_dir, 
                ccd_path=self.cache_dir / "ccd.pkl", mol_dir=self.cache_dir / "mols", 
                use_msa_server=True, msa_server_url="https://api.colabfold.com", 
                msa_pairing_strategy="greedy", boltz2=True, preprocessing_threads=1
            )

            processed_dir = boltz_prediction_dir / "processed"
            manifest = Manifest.load(processed_dir / "manifest.json")
            
            # Write Predictions
            pred_writer = BoltzWriter(
                data_dir=processed_dir / "structures", 
                output_dir=boltz_prediction_dir / "predictions",
                output_format="pdb", boltz2=True
            )

            trainer = Trainer(
                default_root_dir=boltz_prediction_dir, callbacks=[pred_writer],
                accelerator=self.accelerator, devices=1, precision="32" 
                # Note: using precision 32 for CPU stability
            )

            data_module = Boltz2InferenceDataModule(
                manifest=manifest, target_dir=processed_dir / "structures",
                msa_dir=processed_dir / "msa", mol_dir=self.cache_dir / "mols",
                constraints_dir=processed_dir / "constraints", template_dir=processed_dir / "templates",
                extra_mols_dir=processed_dir / "mols", num_workers=0,
            )

            trainer.predict(self.model_module, datamodule=data_module)

            # Zip and Clean
            zip_path = shutil.make_archive(
                str(temp_dir / "boltz_prediction"), 'zip', 
                root_dir=str(temp_dir), base_dir="boltz_prediction"
            )
            shutil.rmtree(boltz_prediction_dir)
            
            return Path(zip_path), temp_dir

        except Exception as e:
            logger.error(f"Error processing request: {e}")
            raise e

    def check_health(self) -> bool:
        return self.model_module is not None