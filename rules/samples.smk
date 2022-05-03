"""
Function for parsing the 'Reads' config and identifying samples and read files
"""

from itertools import chain

def samplesFromDirectoryBam(dir):
    """Parse samples from a directory"""
    outDict = {}
    # https://stackoverflow.com/questions/11860476/how-to-unnest-a-nested-list
    samples= glob_wildcards(os.path.join(dir,'{sample}.bam'))
    samples2 = chain(*samples)
    for sample in samples2:
        outDict[sample] = {}
        bam = os.path.join(dir,f'{sample}.bam')
        if os.path.isfile(bam):
            outDict[sample]['bam'] = bam
        else:
            sys.stderr.write("\n"
                             "    FATAL: Error globbing files."
                             f"    {bam} \n"
                             "    does not exist. Ensure consistent formatting and file extensions."
                             "\n")
            sys.exit(1)
    return outDict

def parseSamplesBam(readFileDir):
    """Parse samples from a directory"""
    if os.path.isdir(readFileDir):
        sampleDict = samplesFromDirectoryBam(readFileDir)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {readFileDir} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any bam samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict

def writeSamplesTsvBam(dict, outfh):
    """Write the samples to a TSV file"""
    with open(outfh, 'w') as out:
        for sample in dict.keys():
            out.write(f'{sample}\t{dict[sample]["bam"]}\n')
    return None


#### fastq


def samplesFromDirectoryFastq(dir):
    """Parse samples from a directory"""
    outDict = {}
    # https://stackoverflow.com/questions/11860476/how-to-unnest-a-nested-list
    samples= glob_wildcards(os.path.join(dir,'{sample}_1.fastq.gz'))
    samples2 = chain(*samples)
    for sample in samples2:
        outDict[sample] = {}
        fastq = os.path.join(dir,f'{sample}_1.fastq.gz')
        if os.path.isfile(fastq):
            outDict[sample]['fastq'] = fastq
        else:
            sys.stderr.write("\n"
                             "    FATAL: Error globbing files."
                             f"    {fastq} \n"
                             "    does not exist. Ensure consistent formatting and file extensions."
                             "\n")
            sys.exit(1)
    return outDict

def parseSamplesFastq(readFileDir):
    """Parse samples from a directory"""
    if os.path.isdir(readFileDir):
        sampleDict = samplesFromDirectoryFastq(readFileDir)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {readFileDir} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any fastq samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict

def writeSamplesTsvFastq(dict, outfh):
    """Write the samples to a TSV file"""
    with open(outfh, 'w') as out:
        for sample in dict.keys():
            out.write(f'{sample}\t{dict[sample]["fastq"]}\n')
    return None

