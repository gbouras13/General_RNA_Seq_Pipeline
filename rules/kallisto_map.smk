rule map_kallisto:
    """map kallisto """
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(KALLISTO_OUTPUT,"{sample}", "abundance.tsv")
    log:
        os.path.join(LOGS,"{sample}_kallisto.log")
    params:
        os.path.join(Kallisto_dir, 'gencode.v39.idx'),
        os.path.join(KALLISTO_OUTPUT, "{sample}"),
    conda:
        os.path.join('..', 'envs','kallisto.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        "kallisto quant -i {params[0]} -o {params[1]} -t {threads} {input[0]}  {input[1]} "


rule aggr_kallisto:
    input:
        expand(os.path.join(KALLISTO_OUTPUT,"{sample}", "abundance.tsv"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "kallisto_agr.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """

