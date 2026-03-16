# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import requests
import time
import base64  # <--- ADDED for decoding
from pathlib import Path

def main():
    """
    Client script to send requests to the Boltz OPEA microservice.
    """
    parser = argparse.ArgumentParser(description="Boltz OPEA Client")
    
    # --- Server Connection Arguments ---
    parser.add_argument("--host", type=str, default="localhost", help="Host of the Boltz microservice.")
    parser.add_argument("--port", type=int, default=8000, help="Port of the Boltz microservice.")

    # --- Boltz Model Input Arguments ---
    # Only YAML file is supported now
    parser.add_argument(
        "--yaml_file",
        type=Path,
        required=True,
        help="Path to a YAML file containing the full prediction input."
    )

    args = parser.parse_args()

    # Construct the URL
    url = f"http://{args.host}:{args.port}/v1/boltz"
    payload = {}

    # --- 1. Prepare Payload ---
    print(f"📄 Reading YAML input from: {args.yaml_file}")
    if not args.yaml_file.exists():
        print(f"❌ Error: YAML file not found at '{args.yaml_file}'")
        return

    with open(args.yaml_file, 'r') as f:
        payload["yaml_content"] = f.read()

    print(f"📡 Sending request to Boltz microservice at {url}...")
    
    start_time = time.time()
    try:
        # We don't need stream=True anymore because we are getting a JSON object
        response = requests.post(url, json=payload)
        
        # --- 2. Handle Response ---
        if response.status_code == 200:
            try:
                data = response.json()
                
                # Check the application-level status from BoltzOutput
                if data.get("status") == "success":
                    encoded_content = data.get("result")
                    message = data.get("message", "")
                    output_filename = "boltz_prediction_results.zip"

                    print("\n✅ Prediction Successful!")
                    print(f"💬 Server Message: {message}")
                    print(f"🔓 Decoding Base64 content...")

                    # --- 3. Decode and Save ---
                    # Decode the Base64 string back to binary
                    binary_data = base64.b64decode(encoded_content)

                    print(f"💾 Saving results to '{output_filename}'...")
                    with open(output_filename, "wb") as f:
                        f.write(binary_data)
                    
                    print(f"🎉 Done. You can now unzip '{output_filename}' to see the results.")

                else:
                    # Application level error (e.g., Boltz failed inside the container)
                    print(f"\n❌ Service Error: {data.get('message')}")

            except requests.exceptions.JSONDecodeError:
                print("\n❌ Error: Received invalid JSON from server.")
                print(f"Raw Response: {response.text[:200]}...") # Print first 200 chars

        else:
            # HTTP level error (e.g., 404, 500)
            print(f"\n❌ HTTP Error {response.status_code}: Request failed.")
            print(f"   - Response: {response.text}")

    except requests.exceptions.ConnectionError:
        print(f"\n❌ Connection Error: Could not connect to the server at {url}.")
        print(f"   - Please ensure the Boltz microservice is running.")

    except Exception as e:
        print(f"\n❌ An unexpected error occurred: {e}")

    finally:
        total_time = time.time() - start_time
        print(f"\n⏱️  Total time: {total_time:.2f} seconds")


if __name__ == "__main__":
    main()