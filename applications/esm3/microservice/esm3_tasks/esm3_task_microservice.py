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
from integrations.native import ESM3_taskInput,ESM3_taskOutput
import sys

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
)
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup

logger = CustomLogger("opea_service_omics_esm3task_microservice")

@register_statistics(names=["opea_service@omics_esm3task"])
async def esm3_task_microservice(input: ESM3_taskInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_esm3task"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS ESM Task invocation: {e}")
        raise


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9008, help="Port to run the microservice on")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_esm3task_component_name = os.getenv("OPEA_OMICS_ESM3EMBDEDDING_COMPONENT_NAME", "OPEA_OMICS_esm3task")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_esm3task_component_name,
        description=f"OPEA OMICS esm3task Component: {opea_omics_esm3task_component_name}",
        config=args.__dict__,
    )
    endpoint = "/v1/esm3task"
    machine_ip = get_machine_ip()
    logger.info(f"Starting esm3task microservice at http://{machine_ip}:{args.port}{endpoint}")
    register_microservice(
        name="opea_service@omics_esm3task",
        endpoint=endpoint,
        host="0.0.0.0",
        port=args.port,
        input_datatype=ESM3_taskInput,
        output_datatype=ESM3_taskOutput,
    )(esm3_task_microservice)
    setup_cleanup(args.port, workload_name="esm3_tasks")
    logger.info("OPEA_OMICS_esm3task server started.")
    opea_microservices["opea_service@omics_esm3task"].start()
