"""
The snakefile that runs the pipeline.
# HPC
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
# STAR can only use 16 threads on hpc for some reason
# https://github.com/alexdobin/STAR/issues/1074
MediumJobCpu = config["MediumJobCpu"]
MediumJobMem = config["MediumJobMem"]


### DIRECTORIES
include: "rules/directories.smk"

READS = config['Reads']
OUTPUT = config['Output']
INPUT = config['Input']
ReadLength = config['ReadLength']


# Parse the samples and read files
include: "rules/samples.smk"

if INPUT == 'fastq':
    sampleReads = parseSamplesFastq(READS)
elif INPUT == 'bam':
    sampleReads = parseSamplesBam(READS)

# samples
SAMPLES = sampleReads.keys()

# Import rules and functions
include: "rules/targets.smk"

if INPUT == 'fastq':
    include: "rules/fastq_parse_fastq.smk"
elif INPUT == 'bam':
    include: "rules/fastq_parse_bam.smk"

include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/salmon_map.smk"
include: "rules/kallisto_map.smk"
include: "rules/receptor.smk"

rule all:
    input:
        PreprocessingFiles,
        AlignFiles,
        FeatureCountFiles,
        SalmonFiles,
        KallistoFiles
