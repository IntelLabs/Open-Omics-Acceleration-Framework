# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License
import sys
import os
# sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
)



import argparse
import os
import time

from comps import (
    CustomLogger,
    OpeaComponentLoader,
    SDInputs,
    SDOutputs,
    ServiceType,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
# from integrations.native import Opea_esmlmdesign
from integrations.native import ESM_foldInput,ESM_foldOutput
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup

from pydantic import BaseModel
from typing import Dict, Any, Optional
from integrations.native import Opea_ESM_fold
from pathlib import Path


logger = CustomLogger("opea_service_omics_esmfold_microservice")

@register_statistics(names=["opea_service@omics_esmfold"])
async def esm_fold_microservice(input: ESM_foldInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_esmfold"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS ESM Fold invocation: {e}")
        raise


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9002, help="Port to run the microservice on")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    parser.add_argument("-m", "--model_dir", help="Parent path to Pretrained ESM data directory. ", type=Path, default=None)
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=None,
        help="Chunks axial attention computation to reduce memory usage from O(L^2) to O(L). "
        "Equivalent to running a for loop over chunks of of each dimension. Lower values will "
        "result in lower memory usage at the cost of speed. Recommended values: 128, 64, 32. "
        "Default: None.",
    )
    parser.add_argument("--cpu_offload", help="Enable CPU offloading", action="store_true")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_esmfold_component_name = os.getenv("OPEA_OMICS_ESMFOLD_COMPONENT_NAME", "OPEA_OMICS_ESMFOLD")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_esmfold_component_name,
        description=f"OPEA OMICS esmfold Component: {opea_omics_esmfold_component_name}",
        config=args.__dict__,
    )
    
    machine_ip = get_machine_ip()
    logger.info(f"Starting esmfold microservice at http://{machine_ip}:{args.port}/v1/esmfold")
    
    register_microservice(
        name="opea_service@omics_esmfold",
        endpoint="/v1/esmfold",
        host="0.0.0.0",
        port=args.port,
        input_datatype=ESM_foldInput,
        output_datatype=ESM_foldOutput,
    )(esm_fold_microservice)
    setup_cleanup(args.port,workload_name="esmfold_inputs")
    logger.info("OPEA_OMICS_esmfold server started.")
    opea_microservices["opea_service@omics_esmfold"].start()