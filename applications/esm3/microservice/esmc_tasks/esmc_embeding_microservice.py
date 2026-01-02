# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import os
import time

from comps import (
    CustomLogger,
    OpeaComponentLoader,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
from integrations.native import ESMC_embeddingInput,ESMC_embeddingOutput
import sys
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
)
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup

logger = CustomLogger("opea_service_omics_esmcembedding_microservice")

@register_statistics(names=["opea_service@omics_esmcembedding"])
async def esmc_embedding_microservice(input: ESMC_embeddingInput):
    start = time.time()
    try:
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_esmcembedding"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS ESM Embedding invocation: {e}")
        raise


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9009, help="Port to run the microservice on")
    parser.add_argument("--model_name", type=str, default="esmc_600m", help="Name of the model for ESMC.")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_esmcembedding_component_name = os.getenv("OPEA_OMICS_ESMCEMBDEDDING_COMPONENT_NAME", "OPEA_OMICS_esmcembedding")
    loader = OpeaComponentLoader(
        opea_omics_esmcembedding_component_name,
        description=f"OPEA OMICS esmcembedding Component: {opea_omics_esmcembedding_component_name}",
        config=args.__dict__,
    )

    machine_ip = get_machine_ip()
    endpoint="/v1/esmcembedding"
    logger.info(f"Starting esmcembedding microservice at http://{machine_ip}:{args.port}{endpoint}")
    register_microservice(
        name="opea_service@omics_esmcembedding",
        endpoint=endpoint,
        host="0.0.0.0",
        port=args.port,
        input_datatype=ESMC_embeddingInput,
        output_datatype=ESMC_embeddingOutput,
    )(esmc_embedding_microservice)
    setup_cleanup(args.port, workload_name="esmc_embedding")
    logger.info("OPEA_OMICS_esmcembedding server started.")
    opea_microservices["opea_service@omics_esmcembedding"].start()
