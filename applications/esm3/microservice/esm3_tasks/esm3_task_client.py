import argparse
from email.mime import base
import requests
import json
import torch
import base64
import csv
import sys

def write_fasta(filename, sequences):
    """Writes generated sequences to a FASTA file."""
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences, 1):
            f.write(f">sampled_seq_{i}\n{seq}\n")

def main():
    # Root parser without auto help to avoid conflicts
    parser = argparse.ArgumentParser(description="Client for ESM Microservice",add_help=False)
    parser.add_argument("--host", default="0.0.0.0", help="Microservice host")
    parser.add_argument("--port", default="9008", help="Microservice port")

    sub = parser.add_subparsers(dest="task", required=True, help="Task to run")

    # === Logits Embedding Parser ===
    logits_embedding_parser = sub.add_parser("logits_embedding")
    logits_embedding_parser.add_argument("--fasta_file", type=str, required=True, help="Path to the input FASTA file.")
    logits_embedding_parser.add_argument("--structure", action="store_true", help="Enable structure logits")
    logits_embedding_parser.add_argument("--secondary_structure", action="store_true", help="Enable secondary structure logits")
    logits_embedding_parser.add_argument("--sasa", action="store_true", help="Enable SASA logits")
    logits_embedding_parser.add_argument("--function", action="store_true", help="Enable function logits")
    logits_embedding_parser.add_argument("--residue_annotations", action="store_true", help="Enable residue annotations")
    logits_embedding_parser.add_argument("--return_hidden_states", action="store_true", help="Return hidden states")
    logits_embedding_parser.add_argument("--ith_hidden_layer", type=int, default=-1, help="Specify which hidden layer to return")
    logits_embedding_parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes")

    # === Fold Parser ===
    fold_parser = sub.add_parser("fold")
    fold_parser.add_argument("--fasta_file", type=str, required=True, help="Path to the input FASTA file.")
    fold_parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine", help="Schedule type (cosine or linear).")
    fold_parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy", help="Unmasking strategy (random or entropy).")
    fold_parser.add_argument("--num_steps", type=int, default=1, help="Number of steps for generation.")
    fold_parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature.")
    fold_parser.add_argument("--temperature_annealing", action="store_true", help="Enable temperature annealing.")
    fold_parser.add_argument("--top_p", type=float, default=1.0, help="Top-p sampling value.")
    fold_parser.add_argument("--condition_on_coordinates_only", action="store_true", help="Condition only on coordinates.")
    fold_parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes.")

    # === InverseFold Parser ===
    inversefold_parser = sub.add_parser("inversefold")
    inversefold_parser.add_argument("--pdb_file", type=str, required=True, help="Path to the input PDB file.")
    inversefold_parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine", help="Schedule type (cosine or linear).")
    inversefold_parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy", help="Unmasking strategy (random or entropy).")
    inversefold_parser.add_argument("--num_steps", type=int, default=1, help="Number of steps for generation.")
    inversefold_parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature.")
    inversefold_parser.add_argument("--temperature_annealing", action="store_true", help="Enable temperature annealing.")
    inversefold_parser.add_argument("--top_p", type=float, default=1.0, help="Top-p sampling value.")
    inversefold_parser.add_argument("--condition_on_coordinates_only", action="store_true", help="Condition only on coordinates.")
    inversefold_parser.add_argument("--num_sequences", type=int, default=8, help="Number of sequences to generate in batch mode.")
    inversefold_parser.add_argument("--batch_run", action="store_true", help="Run in batch mode.")
    inversefold_parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes.")

    # === Chain of Thought Parser ===
    chain_of_thought_parser = sub.add_parser("chain_of_thought")
    chain_of_thought_parser.add_argument("--csv_file", type=str, required=True, help="Path to input CSV file.")
    chain_of_thought_parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine")
    chain_of_thought_parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy")
    chain_of_thought_parser.add_argument("--num_steps", type=int, default=1)
    chain_of_thought_parser.add_argument("--temperature", type=float, default=1.0)
    chain_of_thought_parser.add_argument("--temperature_annealing", action="store_true")
    chain_of_thought_parser.add_argument("--top_p", type=float, default=1.0)
    chain_of_thought_parser.add_argument("--condition_on_coordinates_only", action="store_true")
    chain_of_thought_parser.add_argument("--sequence_length", type=int)
    chain_of_thought_parser.add_argument("--protein_complex", action="store_true")

    # === Function Prediction Parser ===
    function_prediction_parser = sub.add_parser("function_prediction")
    function_prediction_parser.add_argument("--pdb_file", type=str, required=True)
    function_prediction_parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine")
    function_prediction_parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy")
    function_prediction_parser.add_argument("--num_steps", type=int, default=1)
    function_prediction_parser.add_argument("--temperature", type=float, default=1.0)
    function_prediction_parser.add_argument("--temperature_annealing", action="store_true")
    function_prediction_parser.add_argument("--top_p", type=float, default=1.0)
    function_prediction_parser.add_argument("--condition_on_coordinates_only", action="store_true")
    function_prediction_parser.add_argument("--protein_complex", action="store_true")

    # === Prompt Task Parser ===
    prompt_task_parser = sub.add_parser("prompt_task")
    prompt_task_parser.add_argument("--fasta_file", type=str, required=True)
    prompt_task_parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine")
    prompt_task_parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy")
    prompt_task_parser.add_argument("--num_steps", type=int, default=1)
    prompt_task_parser.add_argument("--temperature", type=float, default=1.0)
    prompt_task_parser.add_argument("--temperature_annealing", action="store_true")
    prompt_task_parser.add_argument("--top_p", type=float, default=1.0)
    prompt_task_parser.add_argument("--condition_on_coordinates_only", action="store_true")
    prompt_task_parser.add_argument("--protein_complex", action="store_true")

    # ---------- Parse ----------
    args = parser.parse_args()

    # --- Task-specific parsers ---
    if args.task == "logits_embedding":

        with open(args.fasta_file, "rb") as f:
            fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

        payload = {
            "fasta_base64": fasta_b64,
            "structure": args.structure,
            "secondary_structure": args.secondary_structure,
            "sasa": args.sasa,
            "function": args.function,
            "residue_annotations": args.residue_annotations,
            "return_hidden_states": args.return_hidden_states,
            "ith_hidden_layer": args.ith_hidden_layer,
            "protein_complex": args.protein_complex,
            "task": args.task,
            "port": args.port
        }

    # --- Fold task ---
    elif args.task == "fold":

        with open(args.fasta_file, "rb") as f:
            fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

        payload = {
            "fasta_base64": fasta_b64,
            "schedule": args.schedule,
            "strategy": args.strategy,
            "num_steps": args.num_steps,
            "temperature": args.temperature,
            "temperature_annealing": args.temperature_annealing,
            "top_p": args.top_p,
            "condition_on_coordinates_only": args.condition_on_coordinates_only,
            "protein_complex": args.protein_complex,
            "task": args.task,
            "port": args.port
        }

    # --- InverseFold task ---
    elif args.task == "inversefold":

        with open(args.pdb_file, "rb") as f:
            pdb_b64 = base64.b64encode(f.read()).decode("utf-8")

        payload = {
            "fasta_base64": pdb_b64,
            "schedule": args.schedule,
            "strategy": args.strategy,
            "num_steps": args.num_steps,
            "temperature": args.temperature,
            "temperature_annealing": args.temperature_annealing,
            "top_p": args.top_p,
            "condition_on_coordinates_only": args.condition_on_coordinates_only,
            "protein_complex": args.protein_complex,
            "task": args.task,
            "num_sequences": args.num_sequences,
            "batch_run": args.batch_run,
            "port": args.port
        }

    # --- Chain of Thought ---
    elif args.task == "chain_of_thought":

        with open(args.csv_file, "rb") as f:
            csv_b64 = base64.b64encode(f.read()).decode("utf-8")

        payload = {
            "fasta_base64": csv_b64,
            "schedule": args.schedule,
            "strategy": args.strategy,
            "num_steps": args.num_steps,
            "temperature": args.temperature,
            "temperature_annealing": args.temperature_annealing,
            "top_p": args.top_p,
            "condition_on_coordinates_only": args.condition_on_coordinates_only,
            "protein_complex": args.protein_complex,
            "task": args.task,
            "port": args.port
        }

    # --- Function Prediction ---
    elif args.task == "function_prediction":

        with open(args.pdb_file, "rb") as f:
            pdb_b64 = base64.b64encode(f.read()).decode("utf-8")

        payload = {
            "fasta_base64": pdb_b64,
            "schedule": args.schedule,
            "strategy": args.strategy,
            "num_steps": args.num_steps,
            "temperature": args.temperature,
            "temperature_annealing": args.temperature_annealing,
            "top_p": args.top_p,
            "condition_on_coordinates_only": args.condition_on_coordinates_only,
            "protein_complex": args.protein_complex,
            "task": args.task,
            "port": args.port
        }

    # --- Prompt Task ---
    elif args.task == "prompt_task":

        with open(args.fasta_file, "rb") as f:
            fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

        payload = {
            "fasta_base64": fasta_b64,
            "schedule": args.schedule,
            "strategy": args.strategy,
            "num_steps": args.num_steps,
            "temperature": args.temperature,
            "temperature_annealing": args.temperature_annealing,
            "top_p": args.top_p,
            "condition_on_coordinates_only": args.condition_on_coordinates_only,
            "protein_complex": args.protein_complex,
            "task": args.task,
            "port": args.port
        }

    else:
        raise ValueError(f"Unknown task: {args.task}")

    # --- Send request to server ---
    url = f"http://{args.host}:{args.port}/v1/esm3task"

    try:
        response = requests.post(url, json=payload)
        print(f"Status Code: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            encoded_result = data.get("result")
            decoded_json = base64.b64decode(encoded_result).decode("utf-8")
            decoded_data = json.loads(decoded_json)

            if args.task == "logits_embedding":
                for item in decoded_data:
                    for name, content in item.items():
                        filename = f"{name}.pt"
                        torch.save(content, filename)
                        print(f"Saved {filename}")

            elif args.task == "fold":
                for pdb_entry in decoded_data:
                    for key, content in pdb_entry.items():
                        filename = f"{key}.pdb"
                        with open(filename, "w") as f:
                            f.write(content.strip() + "\n")
                        print(f"Saved {filename}")

            elif args.task == "inversefold":
                filename = "inversefold_output.fasta"
                write_fasta(filename, decoded_data)

            elif args.task == "chain_of_thought":
                filename = "chain_of_thought_output.pdb"
                with open(filename, "w") as f:
                    f.write(decoded_data.strip() + "\n")
                print(f"Saved {filename}")

            elif args.task == "function_prediction":
                with open("function_prediction_output.csv", "w", newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(["Label", "Start", "End"])
                    for annotation in decoded_data:
                        writer.writerow([annotation["label"], annotation["start"], annotation["end"]])

            elif args.task == "prompt_task":
                filename = "prompt_output.fasta"
                write_fasta(filename, decoded_data)

        else:
            print("Server Error:", response.text)

    except Exception as e:
        print(f"Request failed: {e}")


if __name__ == "__main__":
    main()
