# Targeted Sequencing Analyis of HIV Integration Sites
#
# Author : Christopher Nobles, Ph.D.

import os
import sys
import re
import yaml
import configparser
from pathlib import Path
from chivalib import import_sample_info, choose_sequence_data, igv_render
sys.path.append('/home/kevin/anaconda3/envs/qiime1/lib/python2.7/site-packages/')
sys.path.append('/home/kevin/dev/sunbeam')
#from snakemake.utils import R
#from subprocess import DEVNULL, STDOUT, check_call

if not config:
    raise SystemExit("No config file specified.")

#for key, value in config.items() :
#    print(key, value)

# Import sampleInfo
if ".csv" in config["Sample_Info"]:
    delim = ","
elif ".tsv" in config["Sample_Info"]:
    delim = "\t"
else:
    raise SystemExit("Sample Info file must contain extention '.csv' or '.tsv'.")

# Slashes should always be interpreted as directory separators.
wildcard_constraints:
    sample="[^/]+"

## Sample information
print("Attempting to read sample_info file: {}".format(config["Sample_Info"]))
sampleInfo = import_sample_info(
    config["Sample_Info"], config["Sample_Name_Column"], delim)

SAMPLES=sampleInfo[config["Sample_Name_Column"]]
SAMPLES_U3={x: SAMPLES[x] for x in SAMPLES if sampleInfo[config["Unique_Region_Column"]][x]=="U3"}
SAMPLES_U5={x: SAMPLES[x] for x in SAMPLES if sampleInfo[config["Unique_Region_Column"]][x]=="U5"}
TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]

# Trimming data references
R1_LEAD = choose_sequence_data(config["R1_Leading_Trim"], sampleInfo)
R1_OVER = choose_sequence_data(config["R1_Overreading_Trim"], sampleInfo)
R2_LEAD_PRIMER = choose_sequence_data(config["R2_Primer_Seq"], sampleInfo)
R2_LTR = choose_sequence_data(config["R2_LTRBit_Seq"], sampleInfo)
TERM_CA = config["Require_TermCA"]
R2_OVER = config["R2_Overreading_Trim"]

# Default params if not included in config
if not "maxNcount" in config:
    config["maxNcount"] = 1

if not "demultiCores" in config: 
    demulti_cores = snakemake.utils.available_cpu_count()
else:
    demulti_cores = min(
        config["demultiCores"], snakemake.utils.available_cpu_count()
    )


# Working paths
ROOT_DIR = ""
try:
    ROOT_DIR = os.environ["CHIVA_DIR"]
except KeyError:
    raise SystemExit(
        "CHIVA_DIR environment variable not defined. Are you sure you "
        "activated the chiva conda environment?")

if not Path(ROOT_DIR).exists():
    raise SystemExit(
        "Cannot identify Install Directory for cHIVa. Check CHIVA_DIR variable"
        "and make sure you've activated the 'chiva' conda environment."
    )

print("Using ROOT_DIR: {}".format(ROOT_DIR))

RUN = config["Run_Name"]

print("Analyzing RUN: {}".format(RUN))

CODE_DIR = Path(ROOT_DIR) / "tools/rscripts"
if not CODE_DIR.exists():
    raise SystemExit(
        "Cannot identify code directory for cHIVa. Check CHIVA_DIR variable "
        "and make sure you've activated the 'chiva' conda environment. "
        "Additionally, make sure your install of cHIVa is up-to-date."
    )

print("Using CODE_DIR: {}".format(CODE_DIR))

if "Processing_Path" in config:
    PROC_DIR = Path(config["Processing_Path"]) / RUN
    print("Trying PROC_DIR: {}".format(PROC_DIR))
    if not PROC_DIR.exists():
        abs_PROC_DIR = Path(ROOT_DIR) / str(PROC_DIR)
        print("Trying PROC_DIR: {}".format(abs_PROC_DIR))
        if not abs_PROC_DIR.exists():
            raise SystemExit(
                "Cannot locate processing directory: {}".format(abs_PROC_DIR)
            )
        else:
            PROC_DIR = abs_PROC_DIR
else:
    PROC_DIR = Path(ROOT_DIR) / "analysis" / RUN

if not PROC_DIR.exists():
    raise SystemExit(
        "Cannot identify processing directory for cHIVa run {}. Make sure "
        "you've set up your processing directory before trying to 'run' cHIVa. "
        "Also double check your config file.".format(
            config["Run_Name"]
        )
    )
    
if "Viral_Genomes" in config:
    VIRAL_DIR = Path(config["Viral_Genomes"])
    if not VIRAL_DIR.exists():
        abs_VIRAL_DIR = Path(ROOT_DIR) / str(VIRAL_DIR)
        if not abs_VIRAL_DIR.exists():
            raise SystemExit(
                "Cannot locate viral genome directory: {}".format(config["Viral_Genomes"])
            )
        else:
            VIRAL_DIR = abs_VIRAL_DIR
else:
    VIRAL_DIR = Path(ROOT_DIR) / "genomes/viral_genomes"

if not VIRAL_DIR.exists():
    raise SystemExit(
        "Cannot identify viral genome directory for cHIVa run {}. "
        "If you don't have any viral sequences, feel free to use the HXB2 "
        "reference provided with cHIVa. If the {}/genomes/viral_genomes "
        "path does not exist, then you may need to clone cHIVa again from GitHub."
        "Double check your config file.".format(
            config["Run_Name"], ROOT_DIR
        )
    )

if "IGV_Genome" in config:
    IGV_GENOME = Path(config["IGV_Genome"])
    if not IGV_GENOME.exists():
        abs_IGV_GENOME = Path(os.path.join(ROOT_DIR, IGV_GENOME))
        if not abs_IGV_GENOME.exists():
            raise SystemExit(
                "Cannot locate viral genome (for IGV & Bowtie) file: {} (tried {})".format(config["IGV_Genome"], abs_IGV_GENOME)
            )
        else:
            IGV_GENOME = abs_IGV_GENOME
else:
    IGV_GENOME = Path(ROOT_DIR) / "genomes/viral_genomes/HIV_HXB2.fasta"

print("Using IGV_Genome: {}".format(IGV_GENOME))

if "IGV_Genome_Name" in config:
    IGV_GENOME_NAME = config["IGV_Genome_Name"]
else:
    IGV_GENOME_NAME = "K03455|HIVHXB2CG"

print("Using IGV_Genome: {}".format(IGV_GENOME))

if "IGV_GFF3" in config:
    IGV_GFF3 = Path(config["IGV_GFF3"])
    if not IGV_GFF3.exists():
        abs_IGV_GFF3 = Path(os.path.join(ROOT_DIR, IGV_GFF3))
        if not abs_IGV_GFF3.exists():
            raise SystemExit(
                "Cannot locate viral GFF3 (for IGV) file: {} (tried {})".format(config["IGV_GFF3"], abs_IGV_GFF3)
            )
        else:
            IGV_GFF3 = abs_IGV_GFF3
else:
    IGV_GFF3 = Path(ROOT_DIR) / "genomes/viral_genomes/HIV_HXB2.gff3"

print("Using IGV_GFF3: {}".format(IGV_GFF3))

if "BLAST_DB" in config:
    BLAST_DB = config["BLAST_DB"]
else:
    BLAST_DB = ""

print("Using BLAST_DB: {}".format(BLAST_DB))

# Change to strings
ROOT_DIR = str(ROOT_DIR)
CODE_DIR = str(CODE_DIR)
PROC_DIR = str(PROC_DIR) 
VIRAL_DIR = str(VIRAL_DIR)
IGV_GENOME = str(IGV_GENOME)
IGV_BOWTIE_INDEX = IGV_GENOME[:-6]
IGV_GFF3 = str(IGV_GFF3)
BLAST_DB = str(BLAST_DB)

# Check for input files
R1_SEQ_INPUT = Path(config["Seq_Path"]) / config["R1"]
print(R1_SEQ_INPUT)
if not R1_SEQ_INPUT.exists():
    proc_R1_SEQ_INPUT = Path(PROC_DIR) / config["R1"]
    if not proc_R1_SEQ_INPUT.exists():
        skmk_R1_SEQ_INPUT = Path(ROOT_DIR) / config["R1"]
        if not skmk_R1_SEQ_INPUT:
            raise SystemExit(
                "Cannot find sequencing files: {}".format(config["R1"])
            )
        else:
            R1_SEQ_INPUT = skmk_R1_SEQ_INPUT
            R2_SEQ_INPUT = Path(ROOT_DIR) / config["R2"]
            I1_SEQ_INPUT = Path(ROOT_DIR) / config["I1"]
    else:
        R1_SEQ_INPUT = proc_R1_SEQ_INPUT
        R2_SEQ_INPUT = Path(PROC_DIR) / config["R2"]
        I1_SEQ_INPUT = Path(PROC_DIR) / config["I1"]
else:
    R2_SEQ_INPUT = Path(config["Seq_Path"]) / config["R2"]
    I1_SEQ_INPUT = Path(config["Seq_Path"]) / config["I1"]

print(R1_SEQ_INPUT)
R1_SEQ_INPUT = str(R1_SEQ_INPUT)
R2_SEQ_INPUT = str(R2_SEQ_INPUT)
I1_SEQ_INPUT = str(I1_SEQ_INPUT)

## Memory param defaults
if not "demultiMB" in config:
    config["demultiMB"] = 16000

if not "trimMB" in config:
    config["trimMB"] = 4000

if not "filtMB" in config:
    config["filtMB"] = 4000

if not "consolMB" in config:
    config["consolMB"] = 4000

if not "alignMB" in config:
    config["alignMB"] = 4000

if not "coupleMB" in config:
    config["coupleMB"] = 4000

if not "processMB" in config:
    config["processMB"] = 4000

if not "reportMB" in config:
    config["reportMB"] = 4000


# Target Rules
rule all:
  input: 
    #IFR1=expand(PROC_DIR + "/analysis_data/internalfrags/{sample}.R1.internalfrags.fastq.gz", sample=SAMPLES),
    #IFR2=expand(PROC_DIR + "/analysis_data/internalfrags/{sample}.R2.internalfrags.fastq.gz", sample=SAMPLES),
#   demultiFiles=expand(PROC_DIR + "/analysis_data/{sample}.{type}.fastq.gz", sample=SAMPLES, type=TYPES),
    condSites=PROC_DIR + "/output_data/condensed_sites." + RUN + ".csv",
    fragMat=PROC_DIR + "/output_data/fragment_site_matrix." + RUN + ".csv",
    readMat=PROC_DIR + "/output_data/read_site_matrix." + RUN + ".csv",
    report=PROC_DIR + "/output_data/report." + RUN + "." + config["reportFormat"],
    stdSites=PROC_DIR + "/output_data/standardized_uniq_sites." + RUN + ".rds",
    sumTbl=PROC_DIR + "/output_data/summary_table." + RUN + ".csv",
    xofilSites=PROC_DIR + "/output_data/xofil_condensed_sites." + RUN + ".csv",
    sampleSites=expand(PROC_DIR + "/analysis_data/uniqSites/{sample}.uniq.csv", sample=SAMPLES)
#    rawAlign1=expand(PROC_DIR + "/analysis_data/{sample}.raw.R1.psl", sample=SAMPLES),
#    rawAlign2=expand(PROC_DIR + "/analysis_data/{sample}.raw.R2.psl", sample=SAMPLES)
#    LTR_all_metadata=expand(PROC_DIR + "/metadata/all_combined.csv", sample=SAMPLES),
#    LTR_metadata=expand(PROC_DIR + "/metadata/samples/{sample}.csv", sample=SAMPLES),
#    QC=PROC_DIR + "/output_data/QC_table." + RUN + ".csv",
#    QC_report=PROC_DIR + "/output_data/QC_table." + RUN + ".pdf"
#    BLAST_summary=expand(PROC_DIR + "/analysis_data/{sample}.blast.summary.csv", sample=SAMPLES),
    #sampleFasta=expand(PROC_DIR + "/analysis_data/{sample}.fasta", sample=SAMPLES)

# Processing Rules
include: "rules/demulti.rules"
include: "rules/trim.rules"
include: "rules/filter.rules"
include: "rules/consol.rules"
include: "rules/align.blat.rules"
include: "rules/process.rules"
include: "rules/ifseqs.rules"
include: "rules/igv.rules"
#include: "rules/QC.rules"
include: "rules/BLAST.rules"
#include: "rules/ref.map.rules"
