STAR 2.7.10a --- 2022/01/14 ::: New features, behavior changes and bug fixes
============================================================================
See [CHANGES.md](CHANGES.md)

STAR 2.7.9a --- 2021/05/05 ::: STARsolo updates
=====================================================

* [**Counting *multi-gene* (multimapping) reads**](#multi-gene-reads)
* STARsolo uses [SIMDe](https://github.com/simd-everywhere/simde) package which support different types of SIMD extensions. For processors that do not support AVX extensions, specify the target SIMD architecture, e.g.
```
make STAR CXXFLAGS_SIMD=sse
```


STAR 2.7.8a --- 2021/02/20
===========================
**Major STARsolo updates and many bug fixes**

* [**Cell calling (filtering) similar to CellRanger:**](docs/STARsolo.md#emptydrop-like-filtering)
    * ```--soloCellFilter EmptyDrops_CR``` option for cell filtering (calling) nearly identical to that of CellRanger 3 and 4
    * ```--runMode soloCellFiltering``` option for cell filtering (calling) of the raw count matrix, without re-mapping
* [**Input from BAM files for STARsolo:**](docs/STARsolo.md#input-reads-from-bam-files)
    * Input from unmapped or mapped SAM/BAM for STARsolo, with options ```--soloInputSAMattrBarcodeSeq``` and ```--soloInputSAMattrBarcodeQual``` to specify SAM tags for the barcode read sequence and qualities
* [**Read trimming similar to CellRanger4:**](docs/STARsolo.md#matching-cellranger-4xx-and-5xx-results)
    * ```--clipAdapterType CellRanger4``` option for 5' TSO adapter and 3' polyA-tail clipping of the reads to better match CellRanger >= 4.0.0 mapping results
* [**Support for barcodes embedded in mates (such as 10X 5' protocol):**](docs/STARsolo.md#barcode-and-cdna-on-the-same-mate)
    * ```--soloBarcodeMate``` to support scRNA-seq protocols in which one of the paired-end mates contains both barcode sequence and cDNA (e.g. 10X 5' protocol)

STAR 2.7.7a --- 2020/12/28
==========================
**Major new feature:
STARconsensus: mapping RNA-seq reads to consensus genome.**

* Provide the VCF file with consensus SNVs and InDels at the genome generation stage with ```--genomeTransformVCF Variants.vcf  --genomeTransformType Haploid```.
The alternative alleles in this VCF will be inserted to the reference genome to create a "transformed" genome.
Both the genome sequence and transcript/gene annotations are transformed.

* At the mapping stage, the reads will be mapped to the transformed (consensus) genome.
The quantification in the transformed annotations can be performed with standard ```--quantMode TranscriptomeSAM and/or GeneCounts``` options.
If desired, alignments (SAM/BAM) and spliced junctions (SJ.out.tab) can be transformed back to the original (reference) coordinates with ```--genomeTransformOutput SAM and/or SJ```.
This is useful if downstream processing relies on reference coordinates.

STAR 2.7.6a --- 2020/09/19
==========================
**Major new feature:**
Output multimapping chimeric alignments in BAM format using
```
--chimMultimapNmax N>1 --chimOutType WithinBAM --outSAMtype BAM Unsorted [and/or] SortedByCoordinate
```
Many thanks to Sebastian @suhrig who implemented this feature!
More detailed description from Sebastian in PR #802.

STAR 2.7.5a 2020/06/16
======================
**Major new features:  
~ support for Plate-based (Smart-seq) scRNA-seq  
~ manifest file to list the input reads FASTQ files**

* Typical STAR command for mapping and quantification of plate-based (Smart-seq) scRNA-seq  will look like:
```
 --soloType SmartSeq --readFilesManifest /path/to/manifest.tsv --soloUMIdedup Exact --soloStrand Unstranded
```
For detailed description, see [Plate-based (Smart-seq) scRNA-seq](docs/STARsolo.md#plate-based-Smart-seq-scRNA-seq)

* The convenient way to list a large number of reads FASTQ files and their IDs is to create a file manifest and supply it in `--readFilesManifest /path/to/manifest.tsv`. The manifest file should contain 3 tab-separated columns. For paired-end reads:
```
Read1-file-name \t Read2-file-name \t File-id
```
For single-end reads, the 2nd column should contain the dash - :
```
Read1-file-name \t - \t File-id
```
File-id can be any string without spaces. File-id will be added as ReadGroup tag (*RG:Z:*) for each read in the SAM/BAM output. If File-id starts with *ID:*, it can contain several fields separated by tab, and all the fields will be copied verbatim into SAM *@RG* header line.


STAR 2.7.4a 2020/06/01
======================
This release fixes multiple bugs and issues.  
The biggest issue fixed was a seg-fault for small genome which previously required scaling down `--genomeSAindexNbases`. Such scaling is still recommended but is no longer required.  
**This release requires re-generation of the genome indexes**

STAR 2.7.3a 2019/10/08
======================
Major new features in STARsolo
------------------------------
* **Output enhancements:**
    * Summary.csv statistics output for raw and filtered cells useful for quick run quality assessment.
    * --soloCellFilter option for basic filtering of the cells, similar to the methods used by CellRanger 2.2.x.
* [**Better compatibility with CellRanger 3.x.x:**](docs/STARsolo.md#matching-cellranger-3xx-results)
    * --soloUMIfiltering MultiGeneUMI option introduced in CellRanger 3.x.x for filtering UMI collisions between different genes.
    * --soloCBmatchWLtype 1MM_multi_pseudocounts option, introduced in CellRanger 3.x.x, which slightly changes the posterior probability calculation for CB with 1 mismatch.
* [**Velocyto spliced/unspliced/ambiguous quantification:**](docs/STARsolo.md#velocyto-splicedunsplicedambiguous-quantification)
    * --soloFeatures Velocyto option to produce Spliced, Unspliced, and Ambiguous counts similar to the [velocyto.py](http://velocyto.org/) tool developed by [LaManno et al](https://doi.org/10.1038/s41586-018-0414-6). This option is under active development and the results may change in the future versions.
* [**Support for complex barcodes, e.g. inDrop:**](docs/STARsolo.md#barcode-geometry)
    * Complex barcodes in STARsolo with --soloType CB_UMI_Complex, --soloCBmatchWLtype --soloAdapterSequence, --soloAdapterMismatchesNmax, --soloCBposition,--soloUMIposition
* [**BAM tags:**](#bam-tags)
    * CB/UB for corrected CellBarcode/UMI
    * GX/GN for gene ID/name
* STARsolo most up-to-date [documentation](docs/STARsolo.md).

STAR 2.7.2a 2019/08/13
======================

Chimeric read alignment reporting updates
-----------------------------------------
The chimeric read alignment reporting has been changed to improve on the specificity of chimeric read detection.
Only those chimeric read aligments with alignment scores that exceed the corresponding score of the non-chimeric
alignment to the reference genome are now reported.
The Chimeric.junction.out file formatting has been updated to include column headers and the alignment scores for
both chimeric alignments and the non-chimeric alternative. See the latest STAR manual for full details.



STAR 2.7.0c 2019/02/05
======================

STARsolo: mapping, demultiplexing and gene quantification for single cell RNA-seq
---------------------------------------------------------------------------------
STARsolo is a turnkey solution for analyzing droplet single cell RNA sequencing data (e.g. 10X Genomics Chromium System) built directly into STAR code.
STARsolo inputs the raw FASTQ reads files, and performs the following operations
* error correction and demultiplexing of cell barcodes using user-input whitelist
* mapping the reads to the reference genome using the standard STAR spliced read alignment algorithm
* error correction and collapsing (deduplication) of Unique Molecular Identifiers (UMIa)
* quantification of per-cell gene expression by counting the number of reads per gene

STARsolo output is designed to be a drop-in replacement for 10X CellRanger gene quantification output.
It follows CellRanger logic for cell barcode whitelisting and UMI deduplication, and produces nearly identical gene counts in the same format.
At the same time STARsolo is ~10 times faster than the CellRanger.

The STAR solo algorithm is turned on with:
```
--soloType Droplet
```

Presently, the cell barcode whitelist has to be provided with:
```
--soloCBwhitelist /path/to/cell/barcode/whitelist
```

The 10X Chromium whitelist file can be found inside the CellRanger distribution,
e.g. [10X-whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-).
Please make sure that the whitelist is compatible with the specific version of the 10X chemistry (V1,V2,V3 etc).

Importantly, in the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e.
```
--readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz
```

Important: the genome index has to be re-generated with the latest 2.7.0x release.

Other parameters that control STARsolo output are listed below. Note that default parameters are compatible with 10X Chromium V2 protocol.

```
soloCBstart                 1
    int>0: cell barcode start base

soloCBlen                   16
    int>0: cell barcode length

soloUMIstart                17
    int>0: UMI start base

soloUMIlen                  10
    int>0: UMI length

soloStrand                  Forward
    string: strandedness of the solo libraries:
                            Unstranded  ... no strand information
                            Forward     ... read strand same as the original RNA molecule
                            Reverse     ... read strand opposite to the original RNA molecule

soloFeatures                Gene
    string(s)               genomic features for which the UMI counts per Cell Barcode are collected
                            Gene            ... genes: reads match the gene transcript
                            SJ              ... splice junctions: reported in SJ.out.tab

soloUMIdedup                1MM_All
    string(s)               type of UMI deduplication (collapsing) algorithm
                            1MM_All             ... all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)
                            1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).
                            1MM_NotCollapsed      ... UMIs with 1 mismatch distance to others are not collapsed (i.e. all counted)

soloOutFileNames            Solo.out/ genes.tsv barcodes.tsv matrix.mtx matrixSJ.mtx
    string(s)               file names for STARsolo output
                            1st word    ... file name prefix
                            2nd word    ... barcode sequences
                            3rd word    ... gene IDs and names
                            4th word    ... cell/gene counts matrix
                            5th word    ... cell/splice junction counts matrix
```




STAR 2.6.0a 2018/04/23
======================

Major new features:
-------------------

**1. Merging and mapping of overlapping paired-end reads.**  
This feature improves mapping accuracy  for paired-end libraries with short insert sizes, where many reads have overlapping mates. Importantly, it allows detection of chimeric junction in the overlap region.
STAR will search for an overlap between mates larger or equal to --peOverlapNbasesMin bases with proportion of mismatches in the overlap area not exceeding --peOverlapMMp .
If the overlap is found, STAR will map merge the mates and attempt to map the resulting (single-end) sequence.
If requested, the chimeric detection will be performed on the merged-mate sequence, thus allowing chimeric detection in the overlap region.
If the score of this alignment higher than the original one, or if a chimeric alignment is found, STAR will report the merged-mate aligment instead of the original one.
In the output, the merged-mate aligment will be converted back to paired-end format.  
The developmment of this algorithm was supported by Illumina, Inc.
Many thanks to June Snedecor, Xiao Chen, and Felix Schlesinger for their extensive help in developing this feature.


**2. Detection of personal variants overlapping alignments.**  
Option --varVCFfile /path/to/vcf/file is used to input VCF file with personal variants. Only single nucleotide variants (SNVs) are supported at the moment.
Each variant is expected to have a genotype with two alleles.
To output variants that overlap alignments, vG and vA have to be added to --outSAMattributes list.
SAM attribute vG outputs the genomic coordinate of the variant, allowing for identification of the variant.
SAM attribute vA outputs which allele is detected in the read: 1 or 2 match one of the genotype alleles, 3 - no match to genotype.

**3. WASP filtering of allele specific alignments.**  
This is re-implementation of the original WASP algorithm by Bryce van de Geijn, Graham McVicker, Yoav Gilad & Jonathan K Pritchard. Please cite the original [WASP paper: Nature Methods 12, 1061–1063 (2015)   ](https://www.nature.com/articles/nmeth.3582).
WASP filtering is activated with --waspOutputMode SAMtag, which will add vW tag to the SAM output:
vW:i:1 means alignment passed WASP filtering, while all other values mean it did not pass.  
Many thanks to Bryce van de Geijn for fruitful discussions.

**4. Detection of multimapping chimeras.**  
Previous STAR chimeric detection algorithm only detected uniquely mapping chimeras, which reduced its sensitivity in some cases.
The new algorithm can detect and output multimapping chimeras. Presently, the only output into Chimeric.out.junction is supported.
This algorithm is activated with >0 value in --chimMultimapNmax, which defines the maximum number of chimeric multi-alignments.
The --chimMultimapScoreRange (=1 by default) parameter defines the score range for multi-mapping chimeras below the best chimeric score, similar to the --outFilterMultimapScoreRange parameter for normal alignments.
The --chimNonchimScoreDropMin (=20 by default) defines the threshold triggering chimeric detection: the drop in the best non-chimeric alignment score with respect to the read length has to be smaller than this value.  
Many thanks to Brian Haas for testing and feedback.


Minor new features:
-------------------
* --outSAMtlen 1/2 option to select the calculation method for the TLEN field in the SAM/BAM files:
              1 ... leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate
              2 ... leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from 1 for overlapping mates with protruding ends
* --alignInsertionFlush option which defines how to flush ambiguous insertion positions: None: old method, insertions are not flushed; Right: insertions are flushed to the right.
* --outBAMsortingBinsN option to control the number of sorting bins. Increasing this number reduces the amount of RAM required for sorting.


STAR 2.5.0a 2015/11/06
======================

**STAR now uses essential c++11 features. Compiling from sources requires gcc 4.7.0 or later.**

Major new features:
-------------------
1. It is now possible to add extra sequences to the reference genome ont the fly (without re-generating the genome) by specifying
_--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2_ at the mapping stage.

2. By default, the order of the multi-mapping alignments for each read is not truly random.
The _--outMultimapperOrder Random_ option outputs multiple alignments for each read in random order,
and also also randomizes the choice of the primary alignment from the highest scoring alignments.
Parameter _--runRNGseed_ can be used to set the random generator seed.
With this option, the ordering of multi-mapping alignments of each read,
and the choice of the primary alignment will vary from run to run, unless only one thread is used and the seed is kept constant.

3. The --outSAMmultNmax parameter limits the number of output alignments (SAM lines) for multimappers.
For instance, _--outSAMmultNmax 1_ will output exactly one SAM line for each mapped read.


STAR 2.4.2a 2015/06/19
======================

New features:

Counting reads per gene while mapping with --quantMode GeneCounts option.
A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the paired-end read are checked for overlaps.
The counts coincide with those produced by htseq-count with default parameters.

Requires annotations (GTF or GFF with --sjdbGTFfile option) used at the genome generation step, or at the mapping step.

Outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:
column 1: gene ID
column 2: counts for unstranded RNA-seq
column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
Select the output according to the strandedness of your data.
Note, that if you have stranded data and choose one of the columns 3 or 4, the other column (4 or 3) will give you the count of antisense reads.

With --quantMode TranscriptomeSAM GeneCounts, and get both the Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab outputs.



STAR 2.4.1a 2015/04/17
======================

New features:

1. The annotations can now be included on the fly at the mapping step, without including them at the genome generation step.
   At the mapping step, specify --sjdbGTFfile /path/to/ann.gtf and/or --sjdbFileChrStartEnd /path/to/sj.tab, as well as --sjdbOverhang, and any other --sjdb* options.
   The genome indices can be generated with or  without another set of annotations/junctions. In the latter case the new junctions will added to the old ones.
   STAR will insert the junctions into genome indices on the fly before mapping, which takes 1~2 minutes.
   The on the fly genome indices can be saved (for reuse) with "--sjdbInsertSave All", into _STARgenome directory inside the current run directory.
   Default --sjdbOverhang is now set at 100, and does not have to be specified unless you need to change this value.

   The "all-sample" 2-pass method can be simplified using this on the fly junction insertion option:
   (i) run the 1st pass for all samples as usual, with or without annotations
   (ii) run 2nd pass for all samples, listing SJ.out.tab files from all samples in --sjdbFileChrStartEnd /path/to/sj1.tab /path/to/sj2.tab ...

2. New option to activate on the fly "per sample" 2-pass method: "--twopassMode Basic".
   Default --twopass1readsN is now -1, i.e. using all reads in the 1st pass.
   2-pass mode can now be used with annotations, which can be included either at the run-time (see #1), or at the genome generation step.
   Annotated junctions will be included in both the 1st and 2nd passes.

3. Included link (submodule) to Brian Haas' STAR-Fusion code for detecting fusion transcript from STAR chimeric output:
   https://github.com/STAR-Fusion/STAR-Fusion

4. Included Gery Vessere's shared memory implementation for POSIX and SysV.
   To compile STAR with POSIX shared memory, use `make POSIXSHARED`

5. New option "--chimOutType WithinBAM" to include chimeric alignments together with normal alignments in the main (sorted or unsorted) BAM file(s).
   Formatting of chimeric alignments follows the latest SAM/BAM specifications. Thanks to Felix Schlesinger for thorough testing of this option.

6. New option "--quantTranscriptomeBan Singleend" allows insertions, deletions ans soft-clips in the transcriptomic alignments, which can be used by some expression quantification software (e.g. eXpress).

7. New option "--alignEndsTypeExtension Extend5pOfRead1" to enforce full extension of the 5p of the read1, while all other ends undergo local alignment and may be soft-clipped.
