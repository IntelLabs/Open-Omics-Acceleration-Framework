#!/usr/bin/env python

import argparse
import time
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from common.utils.microservice_utils import (
    get_machine_ip,
    is_port_in_use,
    setup_cleanup,
)


from comps import (
    CustomLogger,
    OpeaComponentLoader,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
from integrations.native import (
    OpeaRelion,
    RelionInput,
    RelionOutput,
)

logger = CustomLogger("opea_service_omics_relion")

# Loader will be set in __main__
loader = None  # type: ignore


@register_statistics(names=["opea_service@omics_relion"])
async def relion_service(input: RelionInput):
    """
    Async wrapper that forwards the request to the active OPEA component (OpeaRelion).
    """
    start = time.time()
    try:
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_relion"].append_latency(
            time.time() - start, None
        )
        return results
    except Exception as e:
        logger.error(f"Error during RELION invocation: {e}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="OPEA OMICS RELION microservice (classification + auto-refine)."
    )
    parser.add_argument(
        "--port",
        type=int,
        default=9000,
        help="Port to run the microservice on.",
    )
    args = parser.parse_args()
    server_id = os.getenv("HOSTNAME", "local")
    _ = setup_cleanup(server_id=server_id, workload_name="relion")

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        exit(1)

    # Component name from env or default
    relion_component_name = os.getenv(
        "OPEA_OMICS_RELION_COMPONENT_NAME", "OPEA_OMICS_RELION"
    )

    # Initialize OpeaComponentLoader
    loader = OpeaComponentLoader(
        relion_component_name,
        description=f"OPEA OMICS RELION Component: {relion_component_name}",
        config=args.__dict__,
    )

    machine_ip = get_machine_ip()
    logger.info(
        f"Starting RELION microservice at http://{machine_ip}:{args.port}/v1/relion"
    )

    # Register FastAPI endpoint
    register_microservice(
        name="opea_service@omics_relion",
        endpoint="/v1/relion",
        host="0.0.0.0",
        port=args.port,
        input_datatype=RelionInput,
        output_datatype=RelionOutput,
    )(relion_service)

    logger.info("OPEA OMICS RELION server started.")
    opea_microservices["opea_service@omics_relion"].start()

