# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import requests
import json
from distutils.util import strtobool
import pandas as pd
import base64

def main():
    parser = argparse.ArgumentParser(description="Client for Moflow microservice")
    parser.add_argument("--host", default="0.0.0.0", help="Microservice host")
    parser.add_argument("--port", default="9007", help="Microservice port")
    parser.add_argument('--delta', type=float, default=0.1)
    parser.add_argument('--temperature', type=float, default=1.0,
                        help='temperature of the gaussian distribution')
    parser.add_argument('--additive_transformations', type=strtobool, default='false',
                        help='apply only additive coupling layers')
    parser.add_argument('--topk', type=int, default=500, help='Top k smiles as seeds')
    parser.add_argument("--sim_cutoff", type=float, default=0.00)
    parser.add_argument('--topscore', action='store_true', default=False, help='To find top score')
    parser.add_argument('--consopt', action='store_true', default=False, help='To do constrained optimization')

    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}/v1/moflow_optimize"

    payload = {
        "delta": args.delta,
        "temperature": args.temperature,
        "additive_transformations": args.additive_transformations,
        "topk":args.topk,
        "sim_cutoff": args.sim_cutoff,
        "topscore": args.topscore,
        "consopt": args.consopt,
    }
    print("Sending POST to:", url)

    try:
        response = requests.post(url, json=payload)
        print(f"Status Code: {response.status_code}")

        if response.status_code == 200:
            if args.topscore:
                data = response.json()
                decoded_results = json.loads(base64.b64decode(data["result"]).decode("utf-8"))
                df = pd.DataFrame(decoded_results)
                df.to_csv("discovered_sorted.csv", index=False)
                print("CSV saved as discovered_sorted.csv")
            else: 
                data = response.json()
                decoded_results = json.loads(base64.b64decode(data["result"]).decode("utf-8"))
                print("decoded_results-----",decoded_results)
                df_res = pd.DataFrame(decoded_results["results"])
                df_res.to_csv("constrain_optimization.csv", index=False)
                df_stat = pd.DataFrame(decoded_results["stats"])
                df_stat.to_csv("constrain_optimization_stats.csv", index=False)
                print("summary",decoded_results["summary"])
                print("constrain_optimization.csv saved!")
                 
        else:
            print("Server Error:")
            print(response.text)

    except Exception as e:
        print(f"Request failed: {e}")

if __name__ == "__main__":
    main()