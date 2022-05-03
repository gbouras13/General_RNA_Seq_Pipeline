"""
All target output files are declared here
"""


# if bam is the input


# Preprocessing files
PreprocessingFiles = [
    os.path.join(LOGS, "fastq_parse.txt"),
    os.path.join(MULTIQC,"aggr_qc.txt")
]

if ReadLength == 150:
    AlignFiles = os.path.join(LOGS, "star_150_align.txt")
    FeatureCountFiles = [
        os.path.join(RESULTS,"geneCounts_150.out"),
        os.path.join(RESULTS,"geneCounts_150.txt")
    ]
else:
    AlignFiles = os.path.join(LOGS, "star_200_align.txt")
    FeatureCountFiles = [
        os.path.join(RESULTS,"geneCounts_200.out"),
        os.path.join(RESULTS,"geneCounts_200.txt")
    ]


SalmonFiles = os.path.join(LOGS, "salmon_agr.txt")
KallistoFiles = os.path.join(LOGS, "kallisto_agr.txt")

