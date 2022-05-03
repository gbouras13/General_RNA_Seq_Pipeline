rule fastp_Trim:
    """remove adapters etc """
    input:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}_fastp.log")
    conda:
        os.path.join('..', 'envs','qc.yaml')
    threads:
        16
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w {threads}
        """


rule fastqc:
    """fastqc trimmed reads"""
    input:
        fwd = expand(os.path.join(TMP,"{sample}_trim_R1.fastq.gz"), sample = SAMPLES),
        rev = expand(os.path.join(TMP,"{sample}_trim_R2.fastq.gz"), sample = SAMPLES),
        dir = TMP
    output:
        os.path.join(MULTIQC,"multiqc_report.html")
    params:
        fastqc = FASTQC,
        multiqc = MULTIQC
    conda:
        os.path.join('..', 'envs','qc.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        fastqc -t {threads} -o {params.fastqc} {input.fwd}
        fastqc -t {threads} -o {params.fastqc} {input.rev}
        multiqc {params.fastqc} {input.dir} -o {params.multiqc}
        """



#### aggregation rule
rule aggr_qc:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(MULTIQC,"multiqc_report.html")
    output:
        os.path.join(LOGS, "aggr_qc.txt")
    threads:
        1
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        touch {output[0]}
        """
