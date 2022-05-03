rule map_salmon:
    """map salmon """
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf")
    log:
        os.path.join(LOGS,"{sample}_salmon.log")
    params:
        os.path.join(Salmon_dir, 'gencode.v39_decoys_salmon'),
        os.path.join(SALMON_OUTPUT, "{sample}"),
    conda:
        os.path.join('..', 'envs','salmon.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        "salmon quant -i {params[0]} -l A -1 {input[0]} -2 {input[1]} --validateMappings -p {threads} -o {params[1]}"

rule aggr_salmon:
    input:
        expand(os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "salmon_agr.txt")
    threads:
        1
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        touch {output[0]}
        """



