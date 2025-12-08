# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: --- License
# Author: Rahamathullah [shaikx.rahamathullah@intel.com]
import os

#from some_pdb_library import parse_pdb, renum_pdb_str  # Replace with your actual functions
#import jax
#import jax.numpy as jnp
import numpy as np

#from colabdesign.af.alphafold.common import residue_constants
from string import ascii_uppercase, ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)

def renum_pdb_str(pdb_str, Ls=None, renum=True, offset=1):
  new_chain = None
  if Ls is not None:
    L_init = 0
    new_chain = {}
    for L,c in zip(Ls, alphabet_list):
      new_chain.update({i:c for i in range(L_init,L_init+L)})
      L_init += L

  n,num,pdb_out = 0,offset,[]
  resnum_ = None
  chain_ = None
  new_chain_ = new_chain[0] if new_chain is not None else None #new_chain_ = new_chain[0]
  for line in pdb_str.split("\n"):
    if line[:4] == "ATOM":
      chain = line[21:22]
      resnum = int(line[22:22+5])
      if resnum_ is None: resnum_ = resnum
      if chain_ is None: chain_ = chain
      if resnum != resnum_ or chain != chain_:
        num += (resnum - resnum_)
        n += 1
        resnum_,chain_ = resnum,chain
      if Ls is not None:
        if new_chain[n] != new_chain_:
          num = offset
          new_chain_ = new_chain[n]
      N = num if renum else resnum
      if Ls is None: pdb_out.append("%s%4i%s" % (line[:22],N,line[26:]))
      else: pdb_out.append("%s%s%4i%s" % (line[:21],new_chain[n],N,line[26:]))
  return "\n".join(pdb_out)

#################################################################################

def main():
    pdb_path = os.environ.get("PDB_PATH")
    lengths_path = os.environ.get("LENGTHS_PATH")
    output_path = os.environ.get("OUTPUT_PDB", "final_output.pdb")

    with open(pdb_path, 'r') as f:
        pdb_str = f.read()

    length_data = np.load(lengths_path)
    lengths = length_data["lengths"]
    print("✅ Loaded lengths:", lengths)

    # Remove MODEL/ENDMDL if present
    pdb_lines = [line for line in pdb_str.splitlines() if "MODEL" not in line and "ENDMDL" not in line]
    pdb_str_clean = "\n".join(pdb_lines)

    # Apply renumbering
    pdb_str_renum = renum_pdb_str(pdb_str_clean, lengths)
    pdb_final = f"MODEL        1\n{pdb_str_renum}\nENDMDL\n"

    with open(output_path, 'w') as out:
        out.write(pdb_final)
    print(f"✅ Final PDB written to: {output_path}")

if __name__ == "__main__":
    main()

