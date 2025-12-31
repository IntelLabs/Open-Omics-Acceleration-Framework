import argparse
import requests
import base64
import json
from pathlib import Path
import csv

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
    parser.add_argument("--port", default="9004", help="Microservice port")
    parser.add_argument('--pdbfile', type=str,help='input filepath, either .pdb or .cif')
    parser.add_argument('--seqfile', type=str,help='input filepath for variant sequences in a .fasta file')
    parser.add_argument('--chain', type=str, default=None, help='Chain ID')
    parser.set_defaults(multichain_backbone=False)
    parser.add_argument('--multichain_backbone', action='store_true', help='Use all chains for conditioning')
    parser.add_argument('--singlechain_backbone', dest='multichain_backbone', action='store_false',
                        help='Use only target chain for conditioning')

    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}/v1/esminversefoldscore"
    pdb_str = read_pdb_file(Path(args.pdbfile))
    # Build payload with base64-encoded FASTA
    with open(args.seqfile, "rb") as f:
        fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

    payload = {
        "pdb_str": pdb_str,
        "fasta_str":fasta_b64,
        "chain": args.chain,
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
            decoded_json = json.loads(decoded_bytes.decode("utf-8"))
            if isinstance(decoded_json, list) and all(isinstance(item, list) and len(item) == 2 for item in decoded_json):


                # Save results
                with open("score_ouput.csv", "w") as fout:
                    fout.write("seqid,log_likelihood\n")
                    for seqid, ll in decoded_json:
                        fout.write(f"{seqid},{ll}\n")

                print(f"Results saved to score_ouput.csv")

        else:
            print(" Server returned failure:", parsed.get("message"))

    except Exception as e:
        print(f" Unexpected error: {e}")


if __name__ == "__main__":
    main()
