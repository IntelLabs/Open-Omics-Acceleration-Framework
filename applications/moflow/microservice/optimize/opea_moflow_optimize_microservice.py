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
from integrations.native import MoflowOptimizeInput,MoflowOptimizeOutput
from distutils.util import strtobool

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
)
from common.utils.microservice_utils import get_machine_ip,is_port_in_use

logger = CustomLogger("opea_service_omics_moflow_optimize_microservice")

@register_statistics(names=["opea_service@omics_moflow_optimize"])
async def moflow_microservice(input: MoflowOptimizeInput):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_moflow_optimize"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during OMICS Moflow invocation: {e}")
        raise


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9007, help="Port to run the microservice on")
    parser.add_argument("--bf16", action="store_true", help="Use bfloat16")
    parser.add_argument("--model_dir", type=str, default='./results', required=True)
    parser.add_argument("--data_dir", type=str, default='../data')
    parser.add_argument('--data_name', type=str, default='qm9', choices=['qm9', 'zinc250k'],
                        help='dataset name')
    parser.add_argument("--snapshot_path", "-snapshot", type=str, required=True)
    parser.add_argument("--hyperparams_path", type=str, default='moflow-params.json', required=True)
    parser.add_argument("--property_model_path", type=str, default=None)
    parser.add_argument("--property_name", type=str, default='qed', choices=['qed', 'plogp'])
    # parser.add_argument('--molecule_file', type=str, default='qm9_relgcn_kekulized_ggnp.npz',
    #                     help='path to molecule dataset')parser.add_argument('-l', '--learning_rate', type=float, default=0.001, help='Base learning rate')
    parser.add_argument('-e', '--lr_decay', type=float, default=0.999995,
                        help='Learning rate decay, applied every step of the optimization')
    parser.add_argument('-w', '--weight_decay', type=float, default=1e-5,
                        help='L2 norm for the parameters')
    parser.add_argument('--hidden', type=str, default="",
                        help='Hidden dimension list for output regression')
    parser.add_argument('-x', '--max_epochs', type=int, default=5, help='How many epochs to run in total?')
    parser.add_argument('--debug', type=strtobool, default='true', help='To run optimization with more information')
    parser.add_argument("--batch_size", type=int, default=100)

    parser.add_argument("--timing", action="store_true", help="Enable timing for inference")

    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_moflow_optimize_component_name = os.getenv("OPEA_OMICS_MOFLOW_OPTIMIZE_COMPONENT_NAME", "OPEA_OMICS_MOFLOW_OPTIMIZE")
    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_moflow_optimize_component_name,
        description=f"OPEA OMICS Molfow Optimize Component: {opea_omics_moflow_optimize_component_name}",
        config=args.__dict__,
    )

    machine_ip = get_machine_ip()
    logger.info(f"Starting moflow microservice at http://{machine_ip}:{args.port}/v1/moflow_optimize")

    register_microservice(
        name="opea_service@omics_moflow_optimize",
        endpoint="/v1/moflow_optimize",
        host="0.0.0.0",
        port=args.port,
        input_datatype=MoflowOptimizeInput,
        output_datatype=MoflowOptimizeOutput,
    )(moflow_microservice)

    logger.info("OPEA_OMICS_moflow server started.")
    opea_microservices["opea_service@omics_moflow_optimize"].start()
