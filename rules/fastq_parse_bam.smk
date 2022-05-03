rule bam_to_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}.parse_bam.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools fastq -@ {threads} {input[0]} \
        -1 {output[0]} \
        -2 {output[1]} \
        -0 /dev/null -s /dev/null -n 2> {log}
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
