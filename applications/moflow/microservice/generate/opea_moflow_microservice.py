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
from integrations.native import MoflowInput,MoflowOutput

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
)
from common.utils.microservice_utils import get_machine_ip,is_port_in_use

logger = CustomLogger("opea_service_omics_moflow_microservice")

@register_statistics(names=["opea_service@omics_moflow"])
async def moflow_microservice(input: MoflowInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_moflow"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS Moflow invocation: {e}")
        raise


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9006, help="Port to run the microservice on")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    parser.add_argument("--model_dir", type=str, default='./results')
    parser.add_argument("--snapshot_path", "-snapshot", type=str, required=True)
    parser.add_argument("--hyperparams_path", type=str, default='moflow-params.json', required=True)
    parser.add_argument('--data_name', type=str, default='qm9', choices=['qm9', 'zinc250k'], help='dataset name')
    parser.add_argument("--data_dir", type=str, default='../data')
    parser.add_argument('--seed', type=int, default=123)
    parser.add_argument('--debug', action='store_true', default=False)


    parser.add_argument("--timing", action="store_true", help="Enable timing for inference")

    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_moflow_component_name = os.getenv("OPEA_OMICS_MOFLOW_COMPONENT_NAME", "OPEA_OMICS_MOFLOW")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_moflow_component_name,
        description=f"OPEA OMICS Molfow Component: {opea_omics_moflow_component_name}",
        config=args.__dict__,
    )

    machine_ip = get_machine_ip()
    logger.info(f"Starting moflow microservice at http://{machine_ip}:{args.port}/v1/moflow")

    register_microservice(
        name="opea_service@omics_moflow",
        endpoint="/v1/moflow",
        host="0.0.0.0",
        port=args.port,
        input_datatype=MoflowInput,
        output_datatype=MoflowOutput,
    )(moflow_microservice)

    logger.info("OPEA_OMICS_moflow server started.")
    opea_microservices["opea_service@omics_moflow"].start()
