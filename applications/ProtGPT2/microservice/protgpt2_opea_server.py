# Copyright (C) 2024 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

import argparse
import os
import time

import sys
import os
# sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
)

from comps import (
    CustomLogger,
    OpeaComponentLoader,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
from integrations.native import Opea_ProtGPT2
from integrations.native import ProtGPT2Input, ProtGPT2Output
from common.utils.microservice_utils import get_machine_ip,is_port_in_use,setup_cleanup
from pydantic import BaseModel

logger = CustomLogger("opea_service_omics_protgpt2_microservice")

@register_statistics(names=["opea_service@omics_protgpt2"])
async def protgpt2_service(input: ProtGPT2Input):
    start = time.time()
    try:
        # Use the loader to invoke the active component
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_protgpt2"].append_latency(time.time() - start, None)
        return results
    except Exception as e:
        logger.error(f"Error during ProtGPT2 invocation: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=9010, help="Port to run the ProtGPT2 microservice on")
    parser.add_argument("--fp32", action="store_true", help="Use float32")
    args = parser.parse_args()

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    opea_omics_protgpt2_component_name = os.getenv("OPEA_OMICS_PROTGPT2_COMPONENT_NAME", "OPEA_OMICS_PROTGPT2")

    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        opea_omics_protgpt2_component_name,
        description=f"OPEA OMICS ProtGPT2 Component: {opea_omics_protgpt2_component_name}",
        config=args.__dict__,
    )

    machine_ip = get_machine_ip()
    logger.info(f"Starting ProtGPT2 microservice at http://{machine_ip}:{args.port}/v1/protgpt2")

    register_microservice(
        name="opea_service@omics_protgpt2",
        endpoint="/v1/protgpt2",
        host="0.0.0.0",
        port=args.port,
        input_datatype=ProtGPT2Input,
        output_datatype=ProtGPT2Output,
    )(protgpt2_service)

    logger.info("OPEA_OMICS_PROTGPT2 server started.")
    opea_microservices["opea_service@omics_protgpt2"].start()
