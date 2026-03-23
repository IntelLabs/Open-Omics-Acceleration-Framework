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
from integrations.native import ESM_lmdesignInput,ESM_lmdesignOutput
# from integrations.native import get_machine_ip,is_port_in_use,setup_cleanup
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup

from pydantic import BaseModel
from typing import Dict, Any, Optional
from integrations.native import Opea_ESM_lmdesign

logger = CustomLogger("opea_service_omics_esmlmdesign_microservice")

@register_statistics(names=["opea_service@omics_esmlmdesign"])
async def esm_lmdesign_microservice(input: ESM_lmdesignInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_esmlmdesign"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS ESM LMdesign invocation: {e}")
        raise


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9005, help="Port to run the microservice on")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_esmlmdesign_component_name = os.getenv("OPEA_OMICS_ESMEMBDEDDING_COMPONENT_NAME", "OPEA_OMICS_esmlmdesign")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_esmlmdesign_component_name,
        description=f"OPEA OMICS esmlmdesign Component: {opea_omics_esmlmdesign_component_name}",
        config=args.__dict__,
    )
    
    machine_ip = get_machine_ip()
    logger.info(f"Starting esmlmdesign microservice at http://{machine_ip}:{args.port}/v1/esmlmdesign")
    register_microservice(
        name="opea_service@omics_esmlmdesign",
        endpoint="/v1/esmlmdesign",
        host="0.0.0.0",
        port=args.port,
        input_datatype=ESM_lmdesignInput,
        output_datatype=ESM_lmdesignOutput,
    )(esm_lmdesign_microservice)
    setup_cleanup(args.port,workload_name="inversefold_inputs")
    logger.info("OPEA_OMICS_esmlmdesign server started.")
    opea_microservices["opea_service@omics_esmlmdesign"].start()