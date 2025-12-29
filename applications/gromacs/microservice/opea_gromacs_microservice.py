#!/usr/bin/env python3
import argparse
import socket
import time
import os, sys
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from common.utils.microservice_utils import (
    get_machine_ip,
    is_port_in_use,
    setup_cleanup,
)

from fastapi import HTTPException

from comps import (
    CustomLogger,
    OpeaComponentLoader,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)

from integrations.native import (
    OpeaGromacs,
    GromacsInput,
    GromacsOutput,
)

logger = CustomLogger("opea_service_openomics_gromacs")

_component_loader = None

@register_statistics(names=["opea_service@openomics_gromacs"])
async def gromacs_service(input: GromacsInput):
    """
    Thin service wrapper that delegates to the OpeaGromacs component loader.
    Also records latency into statistics_dict.

    This service supports ONLY workflow mode with:
      - pdb_filename
      - pdb_b64
    and returns only metadata_<job_id>.json as artifact.
    """
    global _component_loader
    start_time = time.time()

    if _component_loader is None:
        raise HTTPException(
            status_code=503,
            detail="Component loader not initialized",
        )

    try:
        res: GromacsOutput = await _component_loader.invoke(input)
        return res

    except HTTPException:
        # propagate HTTP errors (400 / 500 with metadata, etc.)
        raise
    except Exception as e:
        logger.error(f"Error during GROMACS invocation: {e}")
        raise HTTPException(
            status_code=500,
            detail=str(e),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="OPENOMICS GROMACS microservice (workflow-only)."
    )
    parser.add_argument(
        "--port", type=int, default=9000, help="Port to run the microservice on"
    )
    args = parser.parse_args()
    server_id = os.getenv("HOSTNAME", "local")
    _ = setup_cleanup(server_id=server_id, workload_name="openomics_gromacs")

    if is_port_in_use(args.port):
        logger.error(f"Port {args.port} is already in use. Please use a different port.")
        raise SystemExit(1)

    gromacs_component_name = os.getenv(
        "OPENOMICS_GROMACS_COMPONENT_NAME", "OPENOMICS_GROMACS"
    )

    machine_ip = get_machine_ip()
    logger.info(
        f"Starting GROMACS microservice at http://{machine_ip}:{args.port}/v1/gromacs "
        f"(component={gromacs_component_name})"
    )

    # Register FastAPI endpoint
    register_microservice(
        name="opea_service@openomics_gromacs",
        endpoint="/v1/gromacs",
        host="0.0.0.0",
        port=args.port,
        input_datatype=GromacsInput,
        output_datatype=GromacsOutput,
    )(gromacs_service)

    try:
        # Create the component loader so it can build OpeaGromacs
        _component_loader = OpeaComponentLoader(
            gromacs_component_name,
            description=f"OPENOMICS_GROMACS Component: {gromacs_component_name}",
            config=args.__dict__,
        )
    except Exception as e:
        logger.error(f"Failed to initialize components: {e}")
        raise SystemExit(1)

    # Start the microservice (uvicorn under the hood via comps)
    opea_microservices["opea_service@openomics_gromacs"].start()

