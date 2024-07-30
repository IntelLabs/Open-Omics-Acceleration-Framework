## ProtGPT2
ProtGPT2 (peer-reviewed paper) is a language model that speaks the protein language and can be used for de novo protein design and engineering. ProtGPT2 generated sequences conserve natural proteins' critical features (amino acid propensities, secondary structural content, and globularity) while exploring unseen regions of the protein space.

## Model description
ProtGPT2 is based on the GPT2 Transformer architecture and contains 36 layers with a model dimensionality of 1280, totalling 738 million parameters.

ProtGPT2 is a decoder-only transformer model pre-trained on the protein space, database UniRef50 (version 2021_04). The pre-training was done on the raw sequences without FASTA headers. Details of training and datasets can be found here: https://huggingface.co/datasets/nferruz/UR50_2021_04

ProtGPT2 was trained in a self-supervised fashion, i.e., the raw sequence data was used during training without including the annotation of sequences. In particular, ProtGPT2 was trained using a causal modelling objective, in which the model is trained to predict the next token (or, in this case, oligomer) in the sequence. By doing so, the model learns an internal representation of proteins and is able to speak the protein language.

# Downloading the Model

If you want to download the model, use the following option:

--model_source download

By default, `--model_source` is set to `local`, and we have provided the model in the project.

## How to use ProtGPT2

```bash
cd ~/TransOmics.OpenOmicsInternal/applications/protgpt2
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t protgpt2 . 
```
# Run a docker container

```bash
export DATA_DIR=<path-to-database-directory>
export OUTPUT_DIR=<path-to-output-directory>

# outside of the docker
docker run -v $DATA_DIR:/model 
     -v $OUTPUT:/app/output protgpt2:latest 
     python protgpt.py --max_length 150 --do_sample False --model_source local/download
     --top_k 900 --repetition_penalty 1.5 --num_return_sequences 5 
     --eos_token_id 1  --dtype float32/float16 --iterations 5

# inside of the docker

docker run -v $DATA_DIR:/model -v $OUTPUT:/app/output -it sequence:latest bash

python protgpt.py --max_length 150 --do_sample False --model_source local/download
--top_k 900 --repetition_penalty 1.5 --num_return_sequences 5 
--eos_token_id 1  --dtype float32/bfloat16 --iterations 5

```

# Example 1: Generating de novo proteins in a zero-shot fashion

In the example below, ProtGPT2 generates sequences that follow the amino acid 'M'. Any other amino acid, oligomer, fragment, or protein of choice can be selected instead. The model will generate the most probable sequences that follow the input. Alternatively, the input field can also be left empty and it will choose the starting tokens.


```bash
>>> from transformers import pipeline
>>> protgpt2 = pipeline('text-generation', model="nferruz/ProtGPT2")
# length is expressed in tokens, where each token has an average length of 4 amino acids.
>>> sequences = protgpt2("<|endoftext|>", max_length=100, do_sample=True, top_k=950, repetition_penalty=1.2, num_return_sequences=10, eos_token_id=0)
>>> for seq in sequences:
        print(seq):
 {'generated_text': 'MINDLLDISRIISGKMTLDRAEVNLTAIARQVVEEQRQAAEAKSIQLLCSTPDTNHYVFG\nDFDRLKQTLWNLLSNAVKFTPSGGTVELELGYNAEGMEVYVKDSGIGIDPAFLPYVFDRF\nRQSDAADSRNYGGLGLGLAIVKHLLDLHEGNVSAQSEGFGKGATFTVLLPLKPLKRELAA\nVNRHTAVQQSAPLNDNLAGMKILIVEDRPDTNEMVSYILEEAGAIVETAESGAAALTSLK\nSYSPDLVLSDIGMPMMDGYEMIEYIREWKTTKGG'}
{'generated_text': 'MQGDSSISSSNRMFT\nLCKPLTVANETSTLSTTRNSKSNKRVSKQRVNLAESPERNAPSPASIKTNETEEFSTIKT\nTNNEVLGYEPNYVSYDFVPMEKCNLCNENCSIELASLNEETFVKKTICCHECRKKAIENA\nENNNTKGSAVSNNSVTSSSGRKKIIVSGSQILRNLDSLTSSKSNISTLLNPNHLAKLAKN\nGNLSSLSSLQSSASSISKSSSTSSTPTTSPKVSSPTNSPSSSPINSPTP'}
{'generated_text': 'M\nSTHVSLENTLASLQATFFSLEARHTALETQLLSTRTELAATKQELVRVQAEISRADAQAQ\nDLKAQILTLKEKADQAEVEAAAATQRAEESQAALEAQTAELAQLRLEKQAPQHVAEEGDP\nQPAAPTTQAQSPVTSAAAAASSAASAEPSKPELTFPAYTKRKPPTITHAPKAPTKVALNP\nSTLSTSGSGGGAKADPTPTTPVPSSSAGLIPKALRLPPPVTPAASGAKPAPSARSKLRGP\nDAPLSPSTQS'}
{'generated_text': 'MVLLSTGPLPILFLGPSLAELNQKYQVVSDTLLRFTNTV\nTFNTLKFLGSDS\n'}
{'generated_text': 'M\nNNDEQPFIMSTSGYAGNTTSSMNSTSDFNTNNKSNTWSNRFSNFIAYFSGVGWFIGAISV\nIFFIIYVIVFLSRKTKPSGQKQYSRTERNNRDVDSIKRANYYG\n'}
{'generated_text': 'M\nEAVYSFTITETGTGTVEVTPLDRTISGADIVYPPDTACVPLTVQPVINANGTWTLGSGCT\nGHFSVDTTGHVNCLTGGFGAAGVHTVIYTVETPYSGNSFAVIDVNVTEPSGPGDGGNGNG\nDRGDGPDNGGGNNPGPDPDPSTPPPPGDCSSPLPVVCSDRDCADFDTQAQVQIYLDRYGG\nTCDLDGNHDGTPCENLPNNSGGQSSDSGNGGGNPGTGSTHQVVTGDCLWNIASRNNGQGG\nQAWPALLAANNESITNP'}
{'generated_text': 'M\nGLTTSGGARGFCSLAVLQELVPRPELLFVIDRAFHSGKHAVDMQVVDQEGLGDGVATLLY\nAHQGLYTCLLQAEARLLGREWAAVPALEPNFMESPLIALPRQLLEGLEQNILSAYGSEWS\nQDVAEPQGDTPAALLATALGLHEPQQVAQRRRQLFEAAEAALQAIRASA\n'}
{'generated_text': 'M\nGAAGYTGSLILAALKQNPDIAVYALNRNDEKLKDVCGQYSNLKGQVCDLSNESQVEALLS\nGPRKTVVNLVGPYSFYGSRVLNACIEANCHYIDLTGEVYWIPQMIKQYHHKAVQSGARIV\nPAVGFDSTPAELGSFFAYQQCREKLKKAHLKIKAYTGQSGGASGGTILTMIQHGIENGKI\nLREIRSMANPREPQSDFKHYKEKTFQDGSASFWGVPFVMKGINTPVVQRSASLLKKLYQP\nFDYKQCFSFSTLLNSLFSYIFNAI'}
{'generated_text': 'M\nKFPSLLLDSYLLVFFIFCSLGLYFSPKEFLSKSYTLLTFFGSLLFIVLVAFPYQSAISAS\nKYYYFPFPIQFFDIGLAENKSNFVTSTTILIFCFILFKRQKYISLLLLTVVLIPIISKGN\nYLFIILILNLAVYFFLFKKLYKKGFCISLFLVFSCIFIFIVSKIMYSSGIEGIYKELIFT\nGDNDGRFLIIKSFLEYWKDNLFFGLGPSSVNLFSGAVSGSFHNTYFFIFFQSGILGAFIF\nLLPFVYFFISFFKDNSSFMKLF'}
{'generated_text': 'M\nRRAVGNADLGMEAARYEPSGAYQASEGDGAHGKPHSLPFVALERWQQLGPEERTLAEAVR\nAVLASGQYLLGEAVRRFETAVAAWLGVPFALGVASGTAALTLALRAYGVGPGDEVIVPAI\nTFIATSNAITAAGARPVLVDIDPSTWNMSVASLAARLTPKTKAILAVHLWGQPVDMHPLL\nDIAAQANLAVIEDCAQALGASIAGTKVGTFGDAAAFSFYPTKNMTTGEGGMLVTNARDLA\nQAARMLRSHGQDPPTAYMHSQVGFN'}
```
