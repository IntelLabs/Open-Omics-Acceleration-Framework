# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import requests
import base64
import json
from distutils.util import strtobool
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Client for Moflow microservice")
    parser.add_argument("--host", default="0.0.0.0", help="Microservice host")
    parser.add_argument("--port", default="9006", help="Microservice port")
    parser.add_argument("--batch_size", type=int, default=100)
    parser.add_argument('--delta', type=float, default=0.1)
    parser.add_argument('--n_experiments', type=int, default=1, help='number of times generation to be run')
    parser.add_argument('--save_fig', type=strtobool, default='true')
    parser.add_argument('--save_score', type=strtobool, default='true')
    parser.add_argument('--temperature', type=float, default=1.0,
                        help='temperature of the gaussian distribution')
    parser.add_argument('--additive_transformations', type=strtobool, default='false',
                        help='apply only additive coupling layers')
    parser.add_argument('--random_generation', action='store_true', default=False)
    parser.add_argument('--reconstruct', action='store_true', default=False)
    parser.add_argument('--int2point', action='store_true', default=False)
    parser.add_argument('--intgrid', action='store_true', default=False)

    parser.add_argument('--inter_times', type=int, default=5)
    parser.add_argument('--correct_validity', type=strtobool, default='true',
                        help='if apply validity correction after the generation')
    
    args = parser.parse_args()

    # Construct URL
    url = f"http://{args.host}:{args.port}/v1/moflow"

    payload = {
        "batch_size": args.batch_size,
        "delta": args.delta,
        "n_experiments": args.n_experiments,
        "save_fig":args.save_fig,
        "save_score":args.save_score,
        "temperature": args.temperature,
        "additive_transformations": args.additive_transformations,
        "random_generation":args.random_generation,
        "reconstruct": args.reconstruct,
        "int2point": args.int2point,
        "intgrid": args.intgrid,
        "inter_times": args.inter_times,
        "correct_validity": args.correct_validity,
    }
    print("Sending POST to:", url)

    try:
        response = requests.post(url, json=payload)
        print(f"Status Code: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            if args.reconstruct:
                if "result" in data:
                    b64_str = data["result"]
                    decoded_str = base64.b64decode(b64_str).decode("utf-8")
                    decoded_dict = json.loads(decoded_str)
                    with open("results.txt", "w") as f:
                        f.write(str(decoded_dict))
                    
            elif args.random_generation:
                b64_str = data["result"]
                decoded_str = base64.b64decode(b64_str).decode("utf-8")
                result = json.loads(decoded_str)
                if isinstance(result, dict) and "results" in result:
                    results_list = result["results"]
                elif isinstance(result, list):
                    results_list = result
                else:
                    raise ValueError(f"Unexpected result format: {type(result)}, keys: {list(result.keys()) if isinstance(result, dict) else None}")

                for res in results_list:
                    iteration = res["iteration"]

                    df_plogp = pd.DataFrame(res["RankedByPlogp"])
                    df_plogp.to_csv(f"smiles_qed_plogp_{iteration}_RankedByPlogp.csv", index=False)
                    
                    df_qed = pd.DataFrame(res["RankedByQED"])
                    df_qed.to_csv(f"smiles_qed_plogp_{iteration}_RankedByQED.csv", index=False)
                with open("metrics_score.txt", "w") as f:
                    f.write(f"Metrics Score: {result['metrics']}")
                    
            elif args.int2point:
                data = response.json()
                b64_str = data["result"]
                decoded_str = base64.b64decode(b64_str).decode("utf-8")
                results = json.loads(decoded_str)
                for i, item in enumerate(results):
                    image_b64 = item["image_base64"]
                    image_path = f"2points_interpolation-2point_molecules_seed0_{i}.png"
                    with open(image_path, "wb") as f:
                        f.write(base64.b64decode(image_b64))

                    pdf_b64 = item["pdf_base64"]
                    pdf_path = f"2points_interpolation-2point_molecules_seed_{i}.pdf"
                    with open(pdf_path, "wb") as f:
                        f.write(base64.b64decode(pdf_b64))
                    
            elif args.intgrid:
                data = response.json()
                b64_str = data["result"]
                decoded_str = base64.b64decode(b64_str).decode("utf-8")
                results = json.loads(decoded_str)

                for i, item in enumerate(results):

                    image_b64 = item["image_base64_1"]
                    image_path = f"generated_interpolation-grid_molecules_seed{i}.png"
                    with open(image_path, "wb") as f:
                        f.write(base64.b64decode(image_b64))

                    pdf_b64 = item["pdf_base64_1"]
                    pdf_path = f"generated_interpolation-grid_molecules_seed{i}.pdf"
                    with open(pdf_path, "wb") as f:
                        f.write(base64.b64decode(pdf_b64))
                        
                    image_b64 = item["image_base64_2"]
                    image_path = f"generated_interpolation-grid_molecules_seed{i}_unique.png"
                    with open(image_path, "wb") as f:
                        f.write(base64.b64decode(image_b64))

                    pdf_b64 = item["pdf_base64_2"]
                    pdf_path = f"generated_interpolation-grid_molecules_seed{i}_unique.pdf"
                    with open(pdf_path, "wb") as f:
                        f.write(base64.b64decode(pdf_b64))
            else:
                print("Key 'results' not found in response. Available keys:", list(data.keys()))
        else:
            print("Server Error:")
            print(response.text)
    except Exception as e:
        print(f"Request failed: {e}")

if __name__ == "__main__":
    main()