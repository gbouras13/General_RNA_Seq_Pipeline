rule copy_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(READS, "{sample}_1.fastq.gz"),
        os.path.join(READS, "{sample}_1.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        cp {input[0]} > {output[0]}
        cp {input[1]} > {output[1]}
        """

#### aggregation rule

rule aggr_fastq_parse:
    """aggregate."""
    input:
        expand(os.path.join(TMP,"{sample}_R1.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "fastq_parse.txt")
    threads:
        1
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        touch {output[0]}
        """
