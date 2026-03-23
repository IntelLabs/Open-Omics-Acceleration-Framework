import argparse
import os
import time
import subprocess
from pathlib import Path
import base64
from fastapi.responses import FileResponse

from comps import (
    CustomLogger, OpeaComponentLoader, opea_microservices,
    register_microservice, register_statistics, statistics_dict
)

from integrations.boltz_component import BoltzInput, BoltzOutput

logger = CustomLogger("opea_service_omics_boltz")

# --------------------------------------------------------------------
# UTILITY: CALCULATE TIME BASED ON FILE SIZE
# --------------------------------------------------------------------
def hours_to_keep_based_on_file_size(zip_path: Path) -> int:
    try:
        size_gb = zip_path.stat().st_size / (1024 ** 3)
        # E.g. 0.1 GB -> 2.4 hours -> min 1 hour. 2GB -> 48 hours.
        hours = max(1, int(size_gb * 24))  
        return hours
    except Exception:
        return 1 # Fallback

# --------------------------------------------------------------------
# API ENDPOINT
# --------------------------------------------------------------------
@register_microservice(
    name="opea_service@omics_boltz",
    endpoint="/v1/boltz",
    host="0.0.0.0",
    port=8000,
    input_datatype=BoltzInput,
)
@register_statistics(names=["opea_service@omics_boltz"])
async def boltz_service(input: BoltzInput):
    start_time = time.time()

    # Invoke returns (path_to_zip, path_to_parent_temp_dir)
    zip_file_path, temp_dir = await loader.invoke(input)

    latency = time.time() - start_time
    statistics_dict["opea_service@omics_boltz"].append_latency(latency, None)

    # ---- Calculate dynamic cleanup time ----
    hours_to_keep = hours_to_keep_based_on_file_size(zip_file_path)
    seconds_to_sleep = hours_to_keep * 3600
    
    logger.info(f"🕒 Scheduling cleanup for {temp_dir} in {hours_to_keep} hours.")

    # ---- Create Cron Job / Background Cleanup Script ----
    # We create a script inside the temp dir and execute it in the background
    cleanup_script_path = temp_dir / "cleanup.sh"
    
    # Creates a robust script that sleeps, then deletes the specific parent temp dir
    script_content = f"""#!/bin/bash
echo "Cleanup script started for: {temp_dir}"
sleep {seconds_to_sleep}
rm -rf "{temp_dir}"
"""
    
    with open(cleanup_script_path, "w") as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(cleanup_script_path, 0o755)

    # Execute detached (nohup) so it survives if the main server process hiccups (optional but safer)
    # or just run in background.
    subprocess.Popen(
        ["nohup", "bash", str(cleanup_script_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True
    )
    with open(zip_file_path, "rb") as f:
            # Read binary -> Base64 Bytes -> Decode to UTF-8 String
            encoded_content = base64.b64encode(f.read()).decode('utf-8')

    # ---- Return file directly ----
    # FastAPI FileResponse will stream the file from disk.
    return BoltzOutput(
            status="success",
            result=encoded_content,  # <--- The Base64 string
            message="Prediction successful. Result contains Base64 encoded ZIP."
        )


# MAIN ENTRY
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OPEA Boltz Microservice")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()

    loader = OpeaComponentLoader(
        "OPEA_OMICS_BOLTZ",
        description="OPEA OMICS Boltz Component",
        config=args.__dict__,
    )

    opea_microservices["opea_service@omics_boltz"].port = args.port
    logger.info(f"Starting on http://0.0.0.0:{args.port}/v1/boltz")

    opea_microservices["opea_service@omics_boltz"].start()
