import argparse
import os
import time
import sys
from typing import Any, Dict, Optional
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".")))
from fastapi import HTTPException
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
    Opea_AutoDock,
    AutoDockInput,
    AutoDockOutput,
)

logger = CustomLogger("opea_service_omics_autodock_microservice")
loader = None

def _api_error_server(
    status_code: int,
    message: str,
    error_type: str = "ServerError",
    details: Optional[Dict[str, Any]] = None,
) -> HTTPException:
    """
    Same schema as ApiError in native.py, but local helper for server-level failures.
    """
    payload = {
        "error_type": error_type,
        "message": message,
        "details": details,
    }
    return SystemExit(1)


@register_statistics(names=["opea_service@omics_autodock"])
async def autodock_service(input: AutoDockInput) -> AutoDockOutput:
    start = time.time()
    if loader is None:
        raise _api_error_server(500, "OpeaComponentLoader not initialised.",
                details={"reason": "Server startup incomplete or loader missing"},
                )
    try:
        results = await loader.invoke(input)
        statistics_dict["opea_service@omics_autodock"].append_latency(time.time() - start, None)
        return results
    except HTTPException as http_exc:
        # This already follows our standard schema from native.py or _api_error_server
        logger.error(f"HTTPException in autodock_service: {http_exc.detail}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in autodock_service wrapper: {e}")
        raise _api_error_server(
            500,
            message="Unexpected error in AutoDock microservice wrapper",
            details={"reason": str(e)},
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--port",
        type=int,
        default=None,
        help="Port to run the AutoDock-GPU microservice on",
    )
    args = parser.parse_args()

    env_port = os.getenv("AUTODOCK_PORT")
    if args.port:
        port = args.port
    elif env_port:
        port = int(env_port)
    else:
        port = 9012  # default
    server_id = os.getenv("HOSTNAME", "local")
    _ = setup_cleanup(server_id=server_id, workload_name="autodock")

    if is_port_in_use(port):
        logger.error(f"Port {port} is already in use. Please use a different port.")
        # Use standard error schema even here
        raise _api_error_server(
            500,
            message=f"Port {port} is already in use.",
            details={"port": port},
        )

    autodock_component_name = os.getenv(
        "OPEA_OMICS_AUTODOCK_COMPONENT_NAME", "OPEA_OMICS_AUTODOCK"
    )

    loader = OpeaComponentLoader(
        autodock_component_name,
        description=f"OPEA OMICS AutoDock Component: {autodock_component_name}",
        config={"port": port},
    )

    machine_ip = get_machine_ip()
    logger.info(
        f"Starting AutoDock-GPU microservice at http://{machine_ip}:{port}/v1/autodock"
    )

    register_microservice(
        name="opea_service@omics_autodock",
        endpoint="/v1/autodock",
        host="0.0.0.0",
        port=port,
        input_datatype=AutoDockInput,
        output_datatype=AutoDockOutput,
    )(autodock_service)

    logger.info("OPEA_OMICS_AUTODOCK server started.")
    opea_microservices["opea_service@omics_autodock"].start()

