## Open-Omics-Protgpt2
Open-Omics-Acceleration ProtGPT2 is an optimized tool for protein design and engineering. It uses Intel Extension for PyTorch (IPEX) to run efficiently on modern CPUs and performs calculations in bfloat16 (bf16) for faster performance. ProtGPT2 generates protein sequences that mimic natural proteins' key features, like structure and amino acid composition, while also exploring new possibilities in protein design.

## Downloading the Model
```bash
git clone https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal.git
cd TransOmics.OpenOmicsInternal/applications/protgpt2
bash model_script.sh
```
## Run a Protgpt2 Standalone 

```bash 
conda env create -f env.yml

conda activate protgpt2
```

```bash
export OUTPUT_FOLDER=$PWD
chmod a+w $OUTPUT_FOLDER
python protgpt2.py --model_dir ./model_dir --max_length <Integer> --do_sample True --top_k 950 --repetition_penalty 1.2 --num_return_sequences <Integer> --eos_token_id 0  --dtype <float32/bfloat16> --iterations <Integer> --output_file <output_seq.txt>
```

```bash
python protgpt2.py --model_dir ./model_dir --max_length 100 --do_sample True --top_k 950 --repetition_penalty 1.2 --num_return_sequences 10 --eos_token_id 0  --dtype float32 --iterations 5 --output_file output_seq.txt
```

## How to use Docker

## Model Management
The Docker setup provides two ways to handle the required model files:

## Mount Pre-Downloaded Models:
If you already have the model files downloaded on your local machine, you can mount the directory into the Docker container.

## Download Models During Build:
The models can be downloaded directly inside the Docker container during the image build process.


```bash
git clone https://github.com/intel-sandbox/TransOmics.OpenOmicsInternal.git
cd TransOmics.OpenOmicsInternal/applications/protgpt2
```
## Downloading the model
```bash
bash model_script.sh
```
## build the docker image
```bash
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t protgpt2 . 
```
## Run a docker container

## Running
### Information on flags
Performance optimization with bfloat16 and Intel Extension for PyTorch

 `--bf16` flag accelerates performance by utilizing bfloat16 precision, which enhances computational efficiency without compromising accuracy.

When you run the Docker command, the models will automatically be downloaded during the first run. The user will need to wait for the download to complete, but from the second run onward, the inference will run directly without downloading the models again.

```bash
mkdir -p model_dir

export MODELS=$PWD/model_dir
export OUTPUT_FOLDER=$PWD

docker run -it -v $OUTPUT_FOLDER:/app protgpt2:latest python protgpt2.py  --model_dir ./model_dir --max_length <Integer> --do_sample <True> --top_k <950> --repetition_penalty <1.2> --num_return_sequences <Integer> --eos_token_id <0>  --dtype <float32/bfloat16> --iterations <Integer> --output_file <output_seq.txt>


```

## Pre-Downloaded Models docker commnad:

```bash
export OUTPUT_FOLDER=$PWD
chmod a+w $OUTPUT_FOLDER

docker run -it -v $OUTPUT_FOLDER:/app protgpt2_1:latest python protgpt2.py  --model_dir ./model_dir --max_length 100 --do_sample True --top_k 950 --repetition_penalty 1.2 --num_return_sequences 10 --eos_token_id 0  --dtype float32 --iterations 5 --output_file /app/output_seq.txt
```

## Download Models During Build

```bash
mkdir -p model_dir
export MODELS=$PWD/model_dir
export OUTPUT_FOLDER=$PWD
chmod a+w $MODELS $OUTPUT_FOLDER

docker run -it -v $MODELS:/model_dir -v $OUTPUT_FOLDER:/app protgpt2_1:latest python protgpt2.py --max_length 100 --do_sample True --top_k 950 --repetition_penalty 1.2 --num_return_sequences 5 --eos_token_id 0  --dtype float32 --iterations 5 --output_file /app/output_seq.txt
```

# All steps are ended here for optimized Protgpt2.

# The following lines are stock information of the Original Protgpt2 repo:

# **ProtGPT2**

ProtGPT2 ([peer-reviewed paper](https://www.nature.com/articles/s41467-022-32007-7)) is a language model that speaks the protein language and can be used for de novo protein design and engineering. ProtGPT2 generated sequences conserve natural proteins' critical features (amino acid propensities, secondary structural content, and globularity) while exploring unseen regions of the protein space.



## **Model description**
ProtGPT2 is based on the GPT2 Transformer architecture and contains 36 layers with a model dimensionality of 1280, totalling 738 million parameters.

ProtGPT2 is a decoder-only transformer model pre-trained on the protein space, database UniRef50 (version 2021_04). The pre-training was done on the raw sequences without FASTA headers. Details of training and datasets can be found here: https://huggingface.co/datasets/nferruz/UR50_2021_04

ProtGPT2 was trained in a self-supervised fashion, i.e., the raw sequence data was used during training without including the annotation of sequences. In particular, ProtGPT2 was trained using a causal modelling objective, in which the model is trained to predict the next token (or, in this case, oligomer) in the sequence.
 By doing so, the model learns an internal representation of proteins and is able to <em>speak</em> the protein language.

### **How to use ProtGPT2**
ProtGPT2 can be used with the HuggingFace transformer python package. Detailed installation instructions can be found here: https://huggingface.co/docs/transformers/installation

Since ProtGPT2 has been trained on the classical language model objective, it excels at generating protein sequences. It can be used to generate sequences in a zero-shot fashion or to generate sequences of a particular type after finetuning on a user-defined dataset.

**Example 1: Generating _de novo_ proteins in a zero-shot fashion**

In the example below, ProtGPT2 generates sequences that follow the amino acid 'M'. Any other amino acid, oligomer, fragment, or protein of choice can be selected instead. The model will generate the most probable sequences that follow the input. Alternatively, the input field can also be left empty and it will choose the starting tokens.

```
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

**Example 2: Finetuning on a set of user-defined sequences**

This alternative option to the zero-shot generation permits introducing direction in the generation process. User-defined training and validation files containing the sequences of interest are provided to the model. After a short update of the model's weights, ProtGPT2 will generate sequences that follow the input properties.

To create the validation and training file, it is necessary to (1) substitute the FASTA headers for each sequence with the expression "<|endoftext|>" and (2) split the originating dataset into training and validation files (this is often done with the ratio 90/10, 80/20 or 95/5). Then, to finetune the model to the input sequences, we can use the example below. Here we show a learning rate of 1e-06, but ideally, the learning rate should be optimised in separate runs. After training, the finetuned model will be stored in the ./output folder. Lastly, ProtGPT2 can generate the tailored sequences as shown in Example 1:

```
python run_clm.py --model_name_or_path nferruz/ProtGPT2 --train_file training.txt --validation_file validation.txt --tokenizer_name nferruz/ProtGPT2
 --do_train --do_eval --output_dir output --learning_rate 1e-06

```
The HuggingFace script run_clm.py can be found here: https://github.com/huggingface/transformers/blob/master/examples/pytorch/language-modeling/run_clm.py

### **How to select the best sequences**
We've observed that perplexity values correlate with AlphaFold2's plddt.
We recommend computing perplexity for each sequence as follows:

```
sequence='MGEAMGLTQPAVSRAVARLEERVGIRIFNRTARAITLTDEGRRFYEAVAPLLAGIEMHGYR\nVNVEGVAQLLELYARDILAEGRLVQLLPEWAD'

#Convert the sequence to a string like this
#(note we have to introduce new line characters every 60 amino acids,
#following the FASTA file format).

sequence = "<|endoftext|>\nMGEAMGLTQPAVSRAVARLEERVGIRIFNRTARAITLTDEGRRFYEAVAPLLAGIEMHGY\nRVNVEGVAQLLELYARDILAEGRLVQLLPEWAD\n<|endoftext|>"

# ppl function
def calculatePerplexity(sequence, model, tokenizer):
    input_ids = torch.tensor(tokenizer.encode(sequence)).unsqueeze(0)
    input_ids = input_ids.to(device)
    with torch.no_grad():
        outputs = model(input_ids, labels=input_ids)
    loss, logits = outputs[:2]
    return math.exp(loss)

#And hence:
ppl = calculatePerplexity(sequence, model, tokenizer)

```

Where `ppl` is a value with the perplexity for that sequence.
We do not yet have a threshold as to what perplexity value gives a 'good' or 'bad' sequence, but given the fast inference times, the best is to sample many sequences, order them by perplexity, and select those with the lower values (the lower the better).


### **Training specs**
The model was trained on 128 NVIDIA A100 GPUs for 50 epochs, using a block size of 512 and a total batch size of 1024 (65,536 tokens per batch). The optimizer used was Adam (beta1 = 0.9, beta2 = 0.999) with a learning rate of 1e-3.
