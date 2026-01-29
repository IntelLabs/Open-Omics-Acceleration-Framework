# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import requests
import time
def main():
    parser = argparse.ArgumentParser(description="ProtGPT2 OPEA Client")
    parser.add_argument("--host", type=str, required=True, help="Server host")
    parser.add_argument("--port", type=int, required=True, help="Server port")
    parser.add_argument("--max_length", type=int, default=10)
    parser.add_argument("--num_return_sequences", type=int, default=5)
    parser.add_argument("--iterations", type=int, default=1, help="Number of iterations to run")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility")
    args = parser.parse_args()

    start_time = time.time()
    url = f"http://{args.host}:{args.port}/v1/protgpt2"
    payload = {
        "max_length": args.max_length,
        "num_return_sequences": args.num_return_sequences,
        "iterations": args.iterations,
        "seed":args.seed,
    }

    print(f"📡 Sending request to {url} ...")
    response = requests.post(url, json=payload)

    if response.status_code == 200:
        sequences = response.json().get("sequences", [])
        print("✅ Generated Sequences:")
        for seq in sequences:
            print(seq)

        # Save to file
        with open("sequence_output.txt", "w") as f:
            for seq in sequences:
                f.write(seq + "\n")
        print(f"💾 Output saved to sequence_output.txt")

    else:
        print(f"❌ Error {response.status_code}: {response.text}")
    #print(f"total time:, {time.time() - start_tgme}")
    print(f"Total time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
