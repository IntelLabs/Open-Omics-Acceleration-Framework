#!/usr/bin/env python
# encoding: utf-8

# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import os
import time
import socket

from comps import (
    CustomLogger,
    OpeaComponentLoader,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
from integrations.native import Opea_ProteinMPNN
from integrations.native import ProteinMPNNInput, ProteinMPNNOutput

import sys
# sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
)
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup

logger = CustomLogger("opea_service_protein_mpnn_microservice")


@register_statistics(names=["opea_service@protein_mpnn"])
async def protein_mpnn_service(input: ProteinMPNNInput):
    start = time.time()
    try:
        results = await loader.invoke(input)
        statistics_dict["opea_service@protein_mpnn"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during ProteinMPNN invocation: {e}")
        raise


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--port", type=int, default=9100, help="Port to run the ProteinMPNN microservice on")
    argparser.add_argument("--ca_only", action="store_true", default=False, help="Use CA-only models")
    argparser.add_argument("--path_to_model_weights", type=str, default="", help="Path to model weights folder")
    argparser.add_argument("--model_name", type=str, default="v_48_020", help="ProteinMPNN model name")
    argparser.add_argument("--use_soluble_model", action="store_true", default=False, help="Use soluble model weights")

    args = argparser.parse_args()

    # Check if port is in use
    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please choose a different port.")
        exit(1)

    # Get component name from environment or default
    opea_protein_mpnn_component_name = os.getenv("OPEA_PROTEINMPNN_COMPONENT_NAME", "OPEA_OMICS_PROTEINMPNN")
    logger.info(f"Using component name: {opea_protein_mpnn_component_name}")

    # Initialize the component loader with provided arguments as config
    loader = OpeaComponentLoader(
        component_name=opea_protein_mpnn_component_name,
        description=f"OPEA ProteinMPNN Component: {opea_protein_mpnn_component_name}",
        config=args.__dict__  # pass CLI arguments as configuration
    )

    # Get machine IP and log service start
    machine_ip = get_machine_ip()
    logger.info(f"Starting ProteinMPNN microservice at http://{machine_ip}:{args.port}/v1/protein_mpnn")

    # Register the microservice endpoint
    register_microservice(
        name="opea_service@protein_mpnn",
        endpoint="/v1/protein_mpnn",
        host="0.0.0.0",
        port=args.port,
        input_datatype=ProteinMPNNInput,
        output_datatype=ProteinMPNNOutput,
    )(protein_mpnn_service)

    logger.info("ProteinMPNN microservice server started.")
    opea_microservices["opea_service@protein_mpnn"].start()

