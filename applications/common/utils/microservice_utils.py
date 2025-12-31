import sys
import os

import socket

import re

import shutil
import signal
import subprocess
from comps import CustomLogger


def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def get_machine_ip():
    result = subprocess.run(
        ["ip", "-o", "-4", "addr", "show", "scope", "global"],
        capture_output=True, text=True
    )

    preferred_interfaces = ["eth", "eno", "ens", "enp"]  # common ethernet prefixes

    for line in result.stdout.splitlines():
        if any(line.split()[1].startswith(prefix) for prefix in preferred_interfaces):
            match = re.search(r'inet (\d+\.\d+\.\d+\.\d+)', line)
            if match:
                return match.group(1)

        return "127.0.0.1"  # Fallback

def setup_cleanup(server_id: str, workload_name: str):
    logger = CustomLogger(f"{workload_name}_microservice")
    work_dir = f"/tmp/{workload_name}/server_{server_id}"

    def cleanup(signum, frame):
        try:
            if os.path.exists(work_dir):
                shutil.rmtree(work_dir)
                logger.info(f"Cleaned up {work_dir}")
        except Exception as e:
            logger.error(f"Error cleaning {work_dir}: {e}")
        sys.exit(0)  # exit gracefully

    # Register handlers for Ctrl+C and kill
    signal.signal(signal.SIGINT, cleanup)
    signal.signal(signal.SIGTERM, cleanup)

    # Create folder when starting
    os.makedirs(work_dir, exist_ok=True)
    logger.info(f"Work dir created: {work_dir}")

    return work_dir
