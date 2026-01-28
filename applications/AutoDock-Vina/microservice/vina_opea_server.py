import argparse
import os,sys
from comps import (
    CustomLogger,
    opea_microservices,
    register_microservice,
    register_statistics,
)
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from common.utils.microservice_utils import (
    get_machine_ip,
    is_port_in_use,
    setup_cleanup,
)
from integrations.native import VinaInput, VinaOutput, init_vina_component, autodock_vina_service
logger = CustomLogger("autodock_vina_microservice_server")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="OPEA microservice for AutoDock Vina docking."
    )

    parser.add_argument(
        "--port",
        type=int,
        default=9000,
        help="Port to run the Vina microservice on (default: 9000).",
    )

    parser.add_argument(
        "--workspace-root",
        type=str,
        default="/workspaces",
        help="Root directory for per-user/per-job workspaces (default: /workspaces).",
    )

    parser.add_argument(
        "--max-tar-size-mb",
        type=int,
        default=500,
        help=(
            "Maximum allowed size (in MB) for the decoded dataset tar.gz. "
            "Requests exceeding this limit will be rejected."
        ),
    )

    args = parser.parse_args()
    server_id = os.getenv("HOSTNAME", "local")
    _ = setup_cleanup(server_id=server_id, workload_name="autodock_vina")

    if is_port_in_use(args.port):
        logger.error(
            f"Port {args.port} is already in use. Please choose a different port."
        )
        raise SystemExit(1)

    machine_ip = get_machine_ip()
    logger.info(
        f"Starting AutoDock Vina microservice at "
        f"http://{machine_ip}:{args.port}/v1/vina"
    )

    # ---- Build component config and init global component ----
    component_config = {
        "workspace_root": args.workspace_root,
        "max_tar_bytes": int(args.max_tar_size_mb) * 1024 * 1024,
        "port": args.port,
    }

    # Initialize global component in native.py
    init_vina_component(component_config)

    register_microservice(
        name="opea_service@autodock_vina",
        endpoint="/v1/vina",
        host="0.0.0.0",
        port=args.port,
        input_datatype=VinaInput,
        output_datatype=VinaOutput,
    )(autodock_vina_service)

    # ---- Start the FastAPI server ----
    opea_microservices["opea_service@autodock_vina"].start()


if __name__ == "__main__":
    main()

