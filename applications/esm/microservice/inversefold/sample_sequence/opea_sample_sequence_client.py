import argparse
import requests
import base64
import json
from pathlib import Path


def read_pdb_file(pdb_path: Path) -> str:
    """Extract ATOM and HETATM lines from a PDB file."""
    if not pdb_path.exists():
        raise FileNotFoundError(f"File not found: {pdb_path}")
    return "\n".join(
        line for line in pdb_path.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )


def main():
    parser = argparse.ArgumentParser(description="Client for ESM Inversefold microservice")
    parser.add_argument("--host", default="0.0.0.0", help="Microservice host")
    parser.add_argument("--port", default="9003", help="Microservice port")
    parser.add_argument('--pdb_file', type=str, required=True, help='Input PDB filepath')
    parser.add_argument('--chain', type=str, default=None, help='Chain ID')
    parser.add_argument('--temperature', type=float, default=1.0, help='Sampling temperature')
    parser.add_argument('--num_samples', type=int, default=1, help='Number of sequences to sample')
    parser.set_defaults(multichain_backbone=False)
    parser.add_argument('--multichain_backbone', action='store_true', help='Use all chains for conditioning')
    parser.add_argument('--singlechain_backbone', dest='multichain_backbone', action='store_false',
                        help='Use only target chain for conditioning')

    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}/v1/esminversefold"
    pdb_str = read_pdb_file(Path(args.pdb_file))

    payload = {
        "pdb_str": pdb_str,
        "chain": args.chain,
        "temperature": args.temperature,
        "num_samples": args.num_samples,
        "multichain_backbone": args.multichain_backbone,
        "port": args.port
    }

    print(f" Sending POST request to {url} ...")

    try:
        response = requests.post(url, json=payload)
        response.raise_for_status()
        data = response.json()
        if response.status_code == 200:
            decoded_bytes = base64.b64decode(data["result"])
            print("decoded_bytes,...........", decoded_bytes)

            decoded_json = json.loads(decoded_bytes.decode("utf-8"))

            # FIX: check decoded_json, not data
            if isinstance(decoded_json, list):
                fasta_path = "output_sequences.fasta"
                with open(fasta_path, "w") as fasta_file:
                    for item in decoded_json:     # FIX: write decoded_json
                        seq_id = item.get("id", "unknown")
                        sequence = item.get("sequence", "")
                        fasta_file.write(f">{seq_id}\n{sequence}\n")

                print(f"FASTA file saved: {fasta_path}")

        else:
            print("Server returned failure:", data.get("message"))

    except Exception as e:
        print(f" Unexpected error: {e}")


if __name__ == "__main__":
    main()
