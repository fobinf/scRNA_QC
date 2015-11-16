# scRNA_QC
Quality control of single cell RNA sequencing data


A pipeline for mapping and quality assessment single cell RNA-seq data

# See usage and parameters by tryping:
python pipeline_master.py -h 

# What you need to run the program:
Directory of single cell sequencing data (parameter -i)
Mapping tool (parameter -m) - at least one of: bowtie, bowtie2, STAR, BWA, gsnap, salmon
Reference genome - either an index that fits the mapping tool (parameter -g) or as fasta and gtf (use -f and -gtf)
Quantification tool - at least one of: Tophat (run with bowtie or bowtie2), htseq-count or cufflinks 
Ulitity tools: samtools, gtf_splicesites, iit_store, picard-tools, bam2fastq, iit_store
R (with access to the internet for library intstallation and access to biomarRt)

## list is to be continued - provide links to 3rd party tools ##



# Running the program
From the provided data path (input dir relative to provided root - see config file) the pipeline expects a directory named "raw" containing all the files in the experiment.
Files should be named with a hash and a number, that identifies the cell. e.g. for cell number 45. 
SLX-999.GXNNN4.N17-N19#45_1.fq
The two latter "_1" and ".fq" are identifiers for paired read number and fastq-file extention, respectively; and can be set in the config file.

For the most simple running of the pipeline create a directory (eg. /data1/scRNAdata/MouseKO/) for running of the pipeline. Store fastq files in a directory called "raw" (/data1/scRNAdata/MouseKO/raw). Store a fasta file and a gtf files of the reference gennome (e.g. /data1/reference_genomes/mm9.fasta, /data1/reference_genomes/mm9.gtf). Set appropriate paths in the config file (see config_example.txt)


# The config file explained
Please refer to the config_example.txt in order to understand the structure. In order avoid path problems provide all {path}s as full path (starting and ending with "/")

[DIRECTORIES]
ROOT_DIR={path}                     # prefix for input directory. Could be a single cell project directory. Following fastq files must be in {ROOT_DIR}/{INPUT DIRECTORY}/raw/
TEMP_DIR={path}                     # all temporary files will be stored here - could be a fast non-backed up drive 
REF_DIR={path}                      # path to reference genomes should contain directories named according to /{species_name}/{mapper_name} (-g and -m parameters). this can be build uding the. The reference genomes can be build by the pipeline, just requires setting -f parameter
TOOLS_MAPPING =                     # path to mapping tools - will be added to path and must thus contain executables in child directories
TOOLS_QUANTIFICATION = {path}       # prefix for quantification tools.
MAPPING_ROOT = mapped               # prefix for mapping tools
QUANTIFICATION_ROOT = counts        # directory name for count data
PREPROCESSED_ROOT = preprocessed    # directory name for preprocessed files
SORTED_BAM_ROOT = sorted_bam        # directory name for sorted bam files.
SAM_ROOT = sam                      # directory name for samfiles


[EXTENSIONS]
EXT_BAM = .bam
EXT_FASTQ= .fq
EXT_FASTQ_GZ = .fq.gz
EXT_SAM = .sam
EXT_SORTED = .sorted
EXT_COUNTS = .counts
EXT_METRICS = .metrics
EXT_MERGED = .merged
EXT_MARKED_DUPL= .marked.duplicates
EXT_FASTA = .fa
EXT_LOG = .log
EXT_GTF= .gtf
EXT_GSNAP_SPLICE = .splice_sites
EXT_GSNAP_IIT = .iit
EXT_DICT = .dict
EXT_SUMMARY = .stat
FORWARD_STRAND = _1                 # Identifier used if analysing paired end reads.
REVERSE_STRAND = _2                 

[SOFTWARE]
GSNAP=/bin/gsnap                    # give path relative to MAPPING_ROOT
GSNAP_BUILD=/bin/gmap_build         # -
BOWTIE1=/bowtie                     # -
BOWTIE1_BUILD=/bowtie-build         # -
BOWTIE2=/bowtie2                    # -
BOWTIE2_BUILD=/bowtie2-build        # -
BWA=/bwa                            # -
BWA_BUILD=/bwa                      # -
STAR = /source/STAR                 # -
STAR_BUILD=/source/STAR             # -
SALMON = /bin/salmon                # -
SALMON_BUILD = /bin/salmon          # -
HTSEQTOOL=/scripts/htseq-count                                      # give path relative to QUANTIFICATION_ROOT
TOPHAT=/tophat2                                                     # -
CUFFLINKS=/cufflinks                                                # -
PICARD_TOOL = /nfs/research2/teichmann/tools/picard-tools-1.113     # -
BAM_2_FASTQ=/homes/ti1/tools/bam2fastq-1.1.0/bam2fastq              # -
SAMTOOLS=/homes/ti1/tools/samtools-0.1.19/samtools                  # -
GTF_SPLICE_TOOL = /bin/gtf_splicesites                              # -
GTF_IIT_TOOL = /bin/iit_store                                       # -


Link to all 3rd party programs:
### DOTO ####

