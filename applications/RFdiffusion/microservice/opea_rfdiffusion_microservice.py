# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import os
import sys
import time

from comps import (
    CustomLogger,
    OpeaComponentLoader,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
from integrations.native import RFDiffusionInput,RFDiffusionOutput
import argparse

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
)
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup

logger = CustomLogger("opea_service_omics_rfdiffusion_microservice")

@register_statistics(names=["opea_service@omics_rfdiffusion"])
async def rfdiffusion_service(input: RFDiffusionInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_rfdiffusion"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during text2image invocation: {e}")
        raise


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9000, help="Port to run the microservice on")
    parser.add_argument("--bfloat16", action="store_true", help="Use bfloat16")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_rfdiffusion_component_name = os.getenv("OPEA_OMICS_RFDIFFUSION_COMPONENT_NAME", "OPEA_OMICS_RFDIFFUSION")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_rfdiffusion_component_name,
        description=f"OPEA OMICS RFDIFFUSION Component: {opea_omics_rfdiffusion_component_name}",
        config=args.__dict__,
    )

    machine_ip = get_machine_ip()
    logger.info(f"Starting RFdiffusion microservice at http://{machine_ip}:{args.port}/v1/rfdiffusion")

    register_microservice(
        name="opea_service@omics_rfdiffusion",
        endpoint="/v1/rfdiffusion",
        host="0.0.0.0",
        port=args.port,
        input_datatype=RFDiffusionInput,
        output_datatype=RFDiffusionOutput,
    )(rfdiffusion_service)
    setup_cleanup(args.port, workload_name="rfdiffusion")
    logger.info("OPEA_OMICS_RFDIFFUSION server started.")
    opea_microservices["opea_service@omics_rfdiffusion"].start()
