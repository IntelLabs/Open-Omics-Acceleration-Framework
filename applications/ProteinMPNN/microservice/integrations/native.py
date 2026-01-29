#!/usr/bin/env python
# encoding: utf-8

# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import socket
import logging
import os
from pathlib import Path
import sys
import json
import time
import copy
import random
import subprocess
import warnings
import numpy as np
import torch
import netifaces
import base64
import tempfile
from torch.utils.data import DataLoader
from pydantic import BaseModel
from typing import Optional
from fastapi import HTTPException

from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry
# Go two levels up to reach the project root (ProteinMPNN/)

project_root = Path(__file__).resolve().parents[3]

# Add it to Python path
if str(project_root) not in sys.path:
    sys.path.append(str(project_root))

from ProteinMPNN.protein_mpnn_utils import (
    loss_nll, loss_smoothed, gather_edges, gather_nodes, gather_nodes_t,
    cat_neighbors_nodes, _scores, _S_to_seq, tied_featurize, parse_PDB,
    parse_fasta, StructureDataset, StructureDatasetPDB, ProteinMPNN
)

logger = CustomLogger("ProteinMPNN")

def load_model(model_weights_path=None, model_name='v_48_020', ca_only=False,
               use_soluble_model=False, backbone_noise=0.0, device='cpu',
               hidden_dim=128, num_layers=3, seed=None):
    """Load ProteinMPNN model"""
    if model_weights_path:
        model_folder_path = model_weights_path
        if not model_folder_path.endswith('/'):
            model_folder_path += '/'
    else:
        file_path = os.path.realpath(__file__)
        base_dir = os.path.abspath(os.path.join(os.path.dirname(file_path), "../../"))

        if ca_only:
            print("Using CA-ProteinMPNN!")
            model_folder_path = os.path.join(base_dir, 'ca_model_weights') #file_path[:k] + '/ca_model_weights/'
            if use_soluble_model:
                print("WARNING: CA-SolubleMPNN is not available yet")
                sys.exit()
        else:
            if use_soluble_model:
                print("Using ProteinMPNN trained on soluble proteins only!")
                model_folder_path = os.path.join(base_dir, 'soluble_model_weights') #file_path[:k] + '/soluble_model_weights/'
            else:
                model_folder_path = os.path.join(base_dir, 'vanilla_model_weights') #file_path[:k] + '/vanilla_model_weights/'
    checkpoint_path = os.path.join(model_folder_path, f'{model_name}.pt')
    checkpoint = torch.load(checkpoint_path, map_location=device)
    noise_level_print = checkpoint['noise_level']

    model = ProteinMPNN(
        ca_only=ca_only,
        num_letters=21,
        node_features=hidden_dim,
        edge_features=hidden_dim,
        hidden_dim=hidden_dim,
        num_encoder_layers=num_layers,
        num_decoder_layers=num_layers,
        augment_eps=backbone_noise,
        k_neighbors=checkpoint['num_edges']
    )

    model.to(device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    model.model_name = model_name

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)

    return model

class ProteinMPNNInput(BaseModel):
    suppress_print: int
    seed: int
    save_score: bool
    save_probs: bool
    score_only: bool
    path_to_fasta: Optional[str] = None
    conditional_probs_only: bool
    conditional_probs_only_backbone: bool
    unconditional_probs_only: bool
    backbone_noise: float
    num_seq_per_target: int
    batch_size: int
    max_length: int
    sampling_temp: str
    pdb_path: Optional[str] = None
    pdb_path_chains: Optional[str] = None
    jsonl_path: Optional[str] = None
    chain_id_jsonl: str
    fixed_positions_jsonl: str
    omit_AAs: str
    bias_AA_jsonl: str
    bias_by_res_jsonl: str
    omit_AA_jsonl: str
    pssm_jsonl: str
    pssm_multi: float
    pssm_threshold: float
    pssm_log_odds_flag: int
    pssm_bias_flag: int
    tied_positions_jsonl: str

class ProteinMPNNOutput(BaseModel):
    status: str
    message: Optional[str] = None
    results: dict
###############################################################
### it has been taken from the ProteinMPNN/protein_mpnn_run.py#
###############################################################

def run_inference(
    suppress_print,
    model,
    model_name,
    ca_only,
    seed,
    save_score,
    save_probs,
    score_only,
    path_to_fasta,
    conditional_probs_only,
    conditional_probs_only_backbone,
    unconditional_probs_only,
    backbone_noise,
    num_seq_per_target,
    batch_size,
    max_length,
    sampling_temp,
    pdb_path,
    pdb_path_chains,
    jsonl_path,
    chain_id_jsonl,
    fixed_positions_jsonl,
    omit_AAs,
    bias_AA_jsonl,
    bias_by_res_jsonl,
    omit_AA_jsonl,
    pssm_jsonl,
    pssm_multi,
    pssm_threshold,
    pssm_log_odds_flag,
    pssm_bias_flag,
    tied_positions_jsonl):

    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)
        torch.manual_seed(seed)
        print(f"[INFO] CPU random seed set to {seed}")

    hidden_dim = 128
    num_layers = 3

    out_folder = "/tmp/"
    folder_for_outputs = out_folder

    NUM_BATCHES = num_seq_per_target//batch_size
    BATCH_COPIES = batch_size
    temperatures = [float(item) for item in sampling_temp.split()]
    omit_AAs_list = omit_AAs
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
    alphabet_dict = dict(zip(alphabet, range(21)))
    print_all = suppress_print == 0
    omit_AAs_np = np.array([AA in omit_AAs_list for AA in alphabet]).astype(np.float32)
    device = torch.device("cpu")

    # Helper function to load and decode Base64 JSON strings
    def load_json_base64(file_path, name, print_all=False):
        if file_path:
            try:
                decoded_str = base64.b64decode(file_path).decode('utf-8')
                return [json.loads(line) for line in decoded_str.strip().splitlines() if line.strip()]
            except Exception as e:
                if print_all:
                    print(40 * '-')
                    print(f'Failed to decode {name}: {e}')
                return None
        else:
            if print_all:
                print(40 * '-')
                print(f'{name} is NOT loaded (no base64 string)')
            return None

    # Initialize dicts
    chain_id_dict = None
    fixed_positions_dict = None
    pssm_dict = None
    omit_AA_dict = None
    bias_AA_dict = None
    tied_positions_dict = None
    bias_by_res_dict = None

    # Decode and load base64 JSON content
    json_list = load_json_base64(chain_id_jsonl, 'chain_id_jsonl', print_all)
    if json_list:
        for d in json_list:
            chain_id_dict = d

    json_list = load_json_base64(fixed_positions_jsonl, 'fixed_positions_jsonl', print_all)
    if json_list:
        for d in json_list:
            fixed_positions_dict = d

    json_list = load_json_base64(pssm_jsonl, 'pssm_jsonl', print_all)
    if json_list:
        pssm_dict = {}
        for d in json_list:
            pssm_dict.update(d)

    json_list = load_json_base64(omit_AA_jsonl, 'omit_AA_jsonl', print_all)
    if json_list:
        for d in json_list:
            omit_AA_dict = d

    json_list = load_json_base64(bias_AA_jsonl, 'bias_AA_jsonl', print_all)
    if json_list:
        for d in json_list:
            bias_AA_dict = d

    json_list = load_json_base64(tied_positions_jsonl, 'tied_positions_jsonl', print_all)
    if json_list:
        for d in json_list:
            tied_positions_dict = d

    json_list = load_json_base64(bias_by_res_jsonl, 'bias_by_res_jsonl', print_all)
    if json_list:
        for d in json_list:
            bias_by_res_dict = d
        if print_all:
            print('bias by residue dictionary is loaded')

    # Initialize bias array
    bias_AAs_np = np.zeros(len(alphabet))
    if bias_AA_dict:
        for n, AA in enumerate(alphabet):
            if AA in bias_AA_dict:
                bias_AAs_np[n] = bias_AA_dict[AA]

    if pdb_path:
        pdb_base64_decoded = base64.b64decode(pdb_path).decode('utf-8')
        pdb_path = "decoded_input.pdb"
        with open(pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_base64_decoded.strip() + "\n")
        pdb_dict_list = parse_PDB(pdb_path, ca_only=ca_only)
        dataset_valid = StructureDatasetPDB(pdb_dict_list, truncate=None, max_length=max_length)
        all_chain_list = [item[-1:] for item in list(pdb_dict_list[0]) if item[:9]=='seq_chain'] #['A','B', 'C',...]
        if pdb_path_chains:
            designed_chain_list = [str(item) for item in pdb_path_chains.split()]
        else:
            designed_chain_list = all_chain_list
        fixed_chain_list = [letter for letter in all_chain_list if letter not in designed_chain_list]
        chain_id_dict = {}
        chain_id_dict[pdb_dict_list[0]['name']]= (designed_chain_list, fixed_chain_list)
    else:
        jsonl_base64 = base64.b64decode(jsonl_path).decode('utf-8')
        jsonl_path = "decoded_output.json"
        with open(jsonl_path, "w", encoding="utf-8") as f:
            f.write(jsonl_base64.strip() + "\n")
        dataset_valid = StructureDataset(jsonl_path, truncate=None, max_length=max_length, verbose=print_all)
    base_folder = folder_for_outputs
    for sub in ["seqs", "scores", "score_only", "conditional_probs_only", "unconditional_probs_only", "probs"]:
        os.makedirs(os.path.join(base_folder, sub), exist_ok=True)
    results = {
        'seqs': [],
        'score': [],
        'score_only': [],
        'conditional_probs_only': [],
        'unconditional_probs_only': [],
        'probs': []
    }

    # Timing
    start_time = time.time()
    total_residues = 0
    protein_list = []
    total_step = 0
    # Validation epoch
    with torch.no_grad():
        test_sum, test_weights = 0., 0.
        for ix, protein in enumerate(dataset_valid):
            score_list = []
            global_score_list = []
            all_probs_list = []
            all_log_probs_list = []
            S_sample_list = []
            batch_clones = [copy.deepcopy(protein) for i in range(BATCH_COPIES)]
            X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list, visible_list_list, masked_list_list, masked_chain_length_list_list, chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask, tied_pos_list_of_lists_list, pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res_all, tied_beta = tied_featurize(batch_clones, device, chain_id_dict, fixed_positions_dict, omit_AA_dict, tied_positions_dict, pssm_dict, bias_by_res_dict, ca_only=ca_only)
            pssm_log_odds_mask = (pssm_log_odds_all > pssm_threshold).float() #1.0 for true, 0.0 for false
            name_ = batch_clones[0]['name']
            if score_only:
                loop_c = 0
                if path_to_fasta:
                    fasta_text = base64.b64decode(path_to_fasta).decode("utf-8")
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
                        tmp.write(fasta_text.encode("utf-8"))
                        tmp_fasta_path = tmp.name
                    fasta_names, fasta_seqs = parse_fasta(tmp_fasta_path,omit=["/"])
                    loop_c = len(fasta_seqs)
                for fc in range(1+loop_c):
                    if fc == 0:
                        structure_sequence_score_file = base_folder + '/score_only/' + batch_clones[0]['name'] + f'_pdb'
                    else:
                        structure_sequence_score_file = base_folder + '/score_only/' + batch_clones[0]['name'] + f'_fasta_{fc}'
                    native_score_list = []
                    global_native_score_list = []
                    if fc > 0:
                        input_seq_length = len(fasta_seqs[fc-1])
                        S_input = torch.tensor([alphabet_dict[AA] for AA in fasta_seqs[fc-1]], device=device)[None,:].repeat(X.shape[0], 1)
                        S[:,:input_seq_length] = S_input #assumes that S and S_input are alphabetically sorted for masked_chains
                    for j in range(NUM_BATCHES):
                        randn_1 = torch.randn(chain_M.shape, device=X.device)
                        log_probs = model(X, S, mask, chain_M*chain_M_pos, residue_idx, chain_encoding_all, randn_1)
                        mask_for_loss = mask*chain_M*chain_M_pos
                        scores = _scores(S, log_probs, mask_for_loss)
                        native_score = scores.cpu().data.numpy()
                        native_score_list.append(native_score)
                        global_scores = _scores(S, log_probs, mask)
                        global_native_score = global_scores.cpu().data.numpy()
                        global_native_score_list.append(global_native_score)
                    native_score = np.concatenate(native_score_list, 0)
                    global_native_score = np.concatenate(global_native_score_list, 0)
                    ns_mean = native_score.mean()
                    ns_mean_print = np.format_float_positional(np.float32(ns_mean), unique=False, precision=4)
                    ns_std = native_score.std()
                    ns_std_print = np.format_float_positional(np.float32(ns_std), unique=False, precision=4)

                    global_ns_mean = global_native_score.mean()
                    global_ns_mean_print = np.format_float_positional(np.float32(global_ns_mean), unique=False, precision=4)
                    global_ns_std = global_native_score.std()
                    global_ns_std_print = np.format_float_positional(np.float32(global_ns_std), unique=False, precision=4)

                    ns_sample_size = native_score.shape[0]
                    seq_str = _S_to_seq(S[0,], chain_M[0,])
                    np.savez(structure_sequence_score_file, score=native_score, global_score=global_native_score, S=S[0,].cpu().numpy(), seq_str=seq_str)
                    npz_path = f"{structure_sequence_score_file}.npz"
                    with open(npz_path, "rb") as f:
                        encoded_content = base64.b64encode(f.read()).decode("utf-8")
                        if results['score_only'] is None:
                            results['score_only'] = []
                        elif not isinstance(results['score_only'], list):
                            results['score_only'] = [results['score_only']]
                        results['score_only'].append(encoded_content)

                    if print_all:
                        if fc == 0:
                            print(f'Score for {name_} from PDB, mean: {ns_mean_print}, std: {ns_std_print}, sample size: {ns_sample_size},  global score, mean: {global_ns_mean_print}, std: {global_ns_std_print}, sample size: {ns_sample_size}')
                        else:
                            print(f'Score for {name_}_{fc} from FASTA, mean: {ns_mean_print}, std: {ns_std_print}, sample size: {ns_sample_size},  global score, mean: {global_ns_mean_print}, std: {global_ns_std_print}, sample size: {ns_sample_size}')
            elif conditional_probs_only:
                if print_all:
                    print(f'Calculating conditional probabilities for {name_}')
                conditional_probs_only_file = base_folder + '/conditional_probs_only/' + batch_clones[0]['name']
                log_conditional_probs_list = []
                for j in range(NUM_BATCHES):
                    randn_1 = torch.randn(chain_M.shape, device=X.device)
                    log_conditional_probs = model.conditional_probs(X, S, mask, chain_M*chain_M_pos, residue_idx, chain_encoding_all, randn_1, conditional_probs_only_backbone)
                    log_conditional_probs_list.append(log_conditional_probs.cpu().numpy())
                concat_log_p = np.concatenate(log_conditional_probs_list, 0) #[B, L, 21]
                mask_out = (chain_M*chain_M_pos*mask)[0,].cpu().numpy()
                np.savez(conditional_probs_only_file, log_p=concat_log_p, S=S[0,].cpu().numpy(), mask=mask[0,].cpu().numpy(), design_mask=mask_out)
                npz_path = f"{conditional_probs_only_file}.npz"
                with open(npz_path, "rb") as f:
                    encoded_content = base64.b64encode(f.read()).decode("utf-8")
                    if results['conditional_probs_only'] is None:
                        results['conditional_probs_only'] = []
                    elif not isinstance(results['conditional_probs_only'], list):
                        results['conditional_probs_only'] = [results['conditional_probs_only']]
                    results['conditional_probs_only'].append(encoded_content)
            elif unconditional_probs_only:
                if print_all:
                    print(f'Calculating sequence unconditional probabilities for {name_}')
                unconditional_probs_only_file = base_folder + '/unconditional_probs_only/' + batch_clones[0]['name']
                log_unconditional_probs_list = []
                for j in range(NUM_BATCHES):
                    log_unconditional_probs = model.unconditional_probs(X, mask, residue_idx, chain_encoding_all)
                    log_unconditional_probs_list.append(log_unconditional_probs.cpu().numpy())
                concat_log_p = np.concatenate(log_unconditional_probs_list, 0) #[B, L, 21]
                mask_out = (chain_M*chain_M_pos*mask)[0,].cpu().numpy()
                np.savez(unconditional_probs_only_file, log_p=concat_log_p, S=S[0,].cpu().numpy(), mask=mask[0,].cpu().numpy(), design_mask=mask_out)
                npz_path = f"{unconditional_probs_only_file}.npz"
                with open(npz_path, "rb") as f:
                    encoded_content = base64.b64encode(f.read()).decode("utf-8")
                    if results['unconditional_probs_only'] is None:
                        results['unconditional_probs_only'] = []
                    elif not isinstance(results['unconditional_probs_only'], list):
                        results['unconditional_probs_only'] = [results['unconditional_probs_only']]
                    results['unconditional_probs_only'].append(encoded_content)
            else:
                randn_1 = torch.randn(chain_M.shape, device=X.device)
                log_probs = model(X, S, mask, chain_M*chain_M_pos, residue_idx, chain_encoding_all, randn_1)
                mask_for_loss = mask*chain_M*chain_M_pos
                scores = _scores(S, log_probs, mask_for_loss) #score only the redesigned part
                native_score = scores.cpu().data.numpy()
                global_scores = _scores(S, log_probs, mask) #score the whole structure-sequence
                global_native_score = global_scores.cpu().data.numpy()
                # Generate some sequences
                ali_file = base_folder + '/seqs/' + batch_clones[0]['name'] + '.fa'
                score_file = base_folder + '/scores/' + batch_clones[0]['name'] + '.npz'
                probs_file = base_folder + '/probs/' + batch_clones[0]['name'] + '.npz'
                if print_all:
                    print(f'Generating sequences for: {name_}')
                t0 = time.time()
                with open(ali_file, 'w') as f:
                    for temp in temperatures:
                        for j in range(NUM_BATCHES):
                            randn_2 = torch.randn(chain_M.shape, device=X.device)
                            if tied_positions_dict == None:
                                sample_dict = model.sample(X, randn_2, S, chain_M, chain_encoding_all, residue_idx, mask=mask, temperature=temp, omit_AAs_np=omit_AAs_np, bias_AAs_np=bias_AAs_np, chain_M_pos=chain_M_pos, omit_AA_mask=omit_AA_mask, pssm_coef=pssm_coef, pssm_bias=pssm_bias, pssm_multi=pssm_multi, pssm_log_odds_flag=bool(pssm_log_odds_flag), pssm_log_odds_mask=pssm_log_odds_mask, pssm_bias_flag=bool(pssm_bias_flag), bias_by_res=bias_by_res_all)
                                S_sample = sample_dict["S"]
                            else:
                                sample_dict = model.tied_sample(X, randn_2, S, chain_M, chain_encoding_all, residue_idx, mask=mask, temperature=temp, omit_AAs_np=omit_AAs_np, bias_AAs_np=bias_AAs_np, chain_M_pos=chain_M_pos, omit_AA_mask=omit_AA_mask, pssm_coef=pssm_coef, pssm_bias=pssm_bias, pssm_multi=pssm_multi, pssm_log_odds_flag=bool(pssm_log_odds_flag), pssm_log_odds_mask=pssm_log_odds_mask, pssm_bias_flag=bool(pssm_bias_flag), tied_pos=tied_pos_list_of_lists_list[0], tied_beta=tied_beta, bias_by_res=bias_by_res_all)
                            # Compute scores
                                S_sample = sample_dict["S"]
                            log_probs = model(X, S_sample, mask, chain_M*chain_M_pos, residue_idx, chain_encoding_all, randn_2, use_input_decoding_order=True, decoding_order=sample_dict["decoding_order"])
                            mask_for_loss = mask*chain_M*chain_M_pos
                            scores = _scores(S_sample, log_probs, mask_for_loss)
                            scores = scores.cpu().data.numpy()

                            global_scores = _scores(S_sample, log_probs, mask) #score the whole structure-sequence
                            global_scores = global_scores.cpu().data.numpy()

                            all_probs_list.append(sample_dict["probs"].cpu().data.numpy())
                            all_log_probs_list.append(log_probs.cpu().data.numpy())
                            S_sample_list.append(S_sample.cpu().data.numpy())
                            for b_ix in range(BATCH_COPIES):
                                masked_chain_length_list = masked_chain_length_list_list[b_ix]
                                masked_list = masked_list_list[b_ix]
                                seq_recovery_rate = torch.sum(torch.sum(torch.nn.functional.one_hot(S[b_ix], 21)*torch.nn.functional.one_hot(S_sample[b_ix], 21),axis=-1)*mask_for_loss[b_ix])/torch.sum(mask_for_loss[b_ix])
                                seq = _S_to_seq(S_sample[b_ix], chain_M[b_ix])
                                score = scores[b_ix]
                                score_list.append(score)
                                global_score = global_scores[b_ix]
                                global_score_list.append(global_score)
                                native_seq = _S_to_seq(S[b_ix], chain_M[b_ix])
                                if b_ix == 0 and j==0 and temp==temperatures[0]:
                                    start = 0
                                    end = 0
                                    list_of_AAs = []
                                    for mask_l in masked_chain_length_list:
                                        end += mask_l
                                        list_of_AAs.append(native_seq[start:end])
                                        start = end
                                    native_seq = "".join(list(np.array(list_of_AAs)[np.argsort(masked_list)]))
                                    l0 = 0
                                    for mc_length in list(np.array(masked_chain_length_list)[np.argsort(masked_list)])[:-1]:
                                        l0 += mc_length
                                        native_seq = native_seq[:l0] + '/' + native_seq[l0:]
                                        l0 += 1
                                    sorted_masked_chain_letters = np.argsort(masked_list_list[0])
                                    print_masked_chains = [masked_list_list[0][i] for i in sorted_masked_chain_letters]
                                    sorted_visible_chain_letters = np.argsort(visible_list_list[0])
                                    print_visible_chains = [visible_list_list[0][i] for i in sorted_visible_chain_letters]
                                    native_score_print = np.format_float_positional(np.float32(native_score.mean()), unique=False, precision=4)
                                    global_native_score_print = np.format_float_positional(np.float32(global_native_score.mean()), unique=False, precision=4)
                                    script_dir = os.path.dirname(os.path.realpath(__file__))
                                    try:
                                        commit_str = subprocess.check_output(f'git --git-dir {script_dir}/.git rev-parse HEAD', shell=True, stderr=subprocess.DEVNULL).decode().strip()
                                    except subprocess.CalledProcessError:
                                        commit_str = 'unknown'
                                    if ca_only:
                                        print_model_name = 'CA_model_name'
                                    else:
                                        print_model_name = 'model_name'
                                    f.write('>{}, score={}, global_score={}, fixed_chains={}, designed_chains={}, {}={}, git_hash={}, seed={}\n{}\n'.format(name_, native_score_print, global_native_score_print, print_visible_chains, print_masked_chains, print_model_name, model_name, commit_str, seed, native_seq)) #write the native sequence
                                start = 0
                                end = 0
                                list_of_AAs = []
                                for mask_l in masked_chain_length_list:
                                    end += mask_l
                                    list_of_AAs.append(seq[start:end])
                                    start = end

                                seq = "".join(list(np.array(list_of_AAs)[np.argsort(masked_list)]))
                                l0 = 0
                                for mc_length in list(np.array(masked_chain_length_list)[np.argsort(masked_list)])[:-1]:
                                    l0 += mc_length
                                    seq = seq[:l0] + '/' + seq[l0:]
                                    l0 += 1
                                score_print = np.format_float_positional(np.float32(score), unique=False, precision=4)
                                global_score_print = np.format_float_positional(np.float32(global_score), unique=False, precision=4)
                                seq_rec_print = np.format_float_positional(np.float32(seq_recovery_rate.detach().cpu().numpy()), unique=False, precision=4)
                                sample_number = j*BATCH_COPIES+b_ix+1
                                f.write('>T={}, sample={}, score={}, global_score={}, seq_recovery={}\n{}\n'.format(temp,sample_number,score_print,global_score_print,seq_rec_print,seq)) #write generated sequence
                with open(ali_file, "rb") as f:
                    encoded_content = base64.b64encode(f.read()).decode("utf-8")
                    if results['seqs'] is None:
                        results['seqs'] = []
                    elif not isinstance(results['seqs'], list):
                        results['seqs'] = [results['seqs']]
                    results['seqs'].append(encoded_content)
                if save_score:
                    np.savez(score_file, score=np.array(score_list, np.float32), global_score=np.array(global_score_list, np.float32))
                    with open(score_file, "rb") as f:
                        encoded_content = base64.b64encode(f.read()).decode("utf-8")
                        if results['score'] is None:
                            results['score'] = []
                        elif not isinstance(results['score'], list):
                            results['score'] = [results['score']]
                        results['score'].append(encoded_content)
                if save_probs:
                    all_probs_concat = np.concatenate(all_probs_list)
                    all_log_probs_concat = np.concatenate(all_log_probs_list)
                    S_sample_concat = np.concatenate(S_sample_list)
                    np.savez(probs_file, probs=np.array(all_probs_concat, np.float32), log_probs=np.array(all_log_probs_concat, np.float32), S=np.array(S_sample_concat, np.int32), mask=mask_for_loss.cpu().data.numpy(), chain_order=chain_list_list)
                    with open(probs_file, "rb") as f:
                        encoded_content = base64.b64encode(f.read()).decode("utf-8")
                        if results['probs'] is None:
                            results['probs'] = []
                        elif not isinstance(results['probs'], list):
                            results['probs'] = [results['probs']]
                        results['probs'].append(encoded_content)
                t1 = time.time()
                dt = round(float(t1-t0), 4)
                num_seqs = len(temperatures)*NUM_BATCHES*BATCH_COPIES
                total_length = X.shape[1]
                if print_all:
                    print(f'{num_seqs} sequences of length {total_length} generated in {dt} seconds')
    return results

@OpeaComponentRegistry.register("OPEA_OMICS_PROTEINMPNN")
class Opea_ProteinMPNN(OpeaComponent):
    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)
        self.model_name_or_path = (os.getenv("MODEL") or (config.get("model_name") if config else None))
        # 🔸 If nothing provided, use a default model
        if not self.model_name_or_path:
            model_name_or_path = "v_48_020"
        # 🔹 Other optional parameters
        self.ca_only = config.get("ca_only", False) if config else False
        self.use_soluble_model = config.get("use_soluble_model", False) if config else False
        self.path_to_model_weights = config.get("path_to_model_weights", "") if config else ""
        self.device = "cpu"

        # 🔹 Load model dynamically
        model_load_time = time.time()
        self.model = load_model(
                model_weights_path=self.path_to_model_weights,
                model_name=self.model_name_or_path,
                ca_only=self.ca_only,
                use_soluble_model=self.use_soluble_model,
                device=self.device
        )
        print(f"✅ Model loaded in {time.time() - model_load_time:.2f} seconds")

    def check_health(self) -> bool:
        return self.model is not None
        print(f"🧩 Server initialized with model: {model_name_or_path}")
    async def invoke(self, input: ProteinMPNNInput) -> ProteinMPNNOutput:
        try:
            results = run_inference(
                        suppress_print=input.suppress_print,
                        model = self.model,
                        model_name=self.model_name_or_path,
                        ca_only=self.ca_only,
                        seed=input.seed,
                        save_score=input.save_score,
                        save_probs=input.save_probs,
                        score_only=input.score_only,
                        path_to_fasta=input.path_to_fasta,
                        conditional_probs_only=input.conditional_probs_only,
                        conditional_probs_only_backbone=input.conditional_probs_only_backbone,
                        unconditional_probs_only=input.unconditional_probs_only,
                        backbone_noise=input.backbone_noise,
                        num_seq_per_target=input.num_seq_per_target,
                        batch_size=input.batch_size,
                        max_length=input.max_length,
                        sampling_temp=input.sampling_temp,
                        pdb_path=input.pdb_path,
                        pdb_path_chains=input.pdb_path_chains,
                        jsonl_path=input.jsonl_path,
                        chain_id_jsonl=input.chain_id_jsonl,
                        fixed_positions_jsonl=input.fixed_positions_jsonl,
                        omit_AAs=input.omit_AAs,
                        bias_AA_jsonl=input.bias_AA_jsonl,
                        bias_by_res_jsonl=input.bias_by_res_jsonl,
                        omit_AA_jsonl=input.omit_AA_jsonl,
                        pssm_jsonl=input.pssm_jsonl,
                        pssm_multi=input.pssm_multi,
                        pssm_threshold=input.pssm_threshold,
                        pssm_log_odds_flag=input.pssm_log_odds_flag,
                        pssm_bias_flag=input.pssm_bias_flag,
                        tied_positions_jsonl=input.tied_positions_jsonl
                        )
            return ProteinMPNNOutput(
                status="success",
                message="Inference completed successfully.",
                results=results
            )
        except HTTPException as http_exc:
            raise http_exc
        except Exception as e:
            logger.error(f"Unexpected error in ProteinMPNN inference: {e}")
            raise HTTPException(status_code=500, detail=str(e))
