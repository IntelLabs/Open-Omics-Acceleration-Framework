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
from integrations.native import ESM_embeddingInput,ESM_embeddingOutput
# from integrations.native import get_machine_ip,is_port_in_use,setup_cleanup

from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup
from pydantic import BaseModel
from typing import Dict, Any, Optional
from integrations.native import Opea_ESM_embedding

logger = CustomLogger("opea_service_omics_esmembedding_microservice")

@register_statistics(names=["opea_service@omics_esmembedding"])
async def esm_embedding_microservice(input: ESM_embeddingInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_esmembedding"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS ESM Embedding invocation: {e}")
        raise


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9001, help="Port to run the microservice on")
    parser.add_argument("--model_name", type=str, default="esm2_t33_650M_UR50D", help="Model to use on the server")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_esmembedding_component_name = os.getenv("OPEA_OMICS_ESMEMBDEDDING_COMPONENT_NAME", "OPEA_OMICS_esmembedding")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_esmembedding_component_name,
        description=f"OPEA OMICS esmembedding Component: {opea_omics_esmembedding_component_name}",
        config=args.__dict__,
    )
    
    machine_ip = get_machine_ip()
    logger.info(f"Starting esmembedding microservice at http://{machine_ip}:{args.port}/v1/esmembedding")
    register_microservice(
        name="opea_service@omics_esmembedding",
        endpoint="/v1/esmembedding",
        host="0.0.0.0",
        port=args.port,
        input_datatype=ESM_embeddingInput,
        output_datatype=ESM_embeddingOutput,
    )(esm_embedding_microservice)
    setup_cleanup(args.port,workload_name="embedding_inputs")
    logger.info("OPEA_OMICS_esmembedding server started.")
    opea_microservices["opea_service@omics_esmembedding"].start()