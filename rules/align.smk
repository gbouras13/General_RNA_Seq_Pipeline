rule align_star_150:
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(STAR_BAMS,"{sample}_star_150_Aligned.sortedByCoord.out.bam")
    log:
        os.path.join(LOGS,"{sample}_star.log")
    params:
        os.path.join(HG38_dir, 'hg38_150'),
        os.path.join(STAR_BAMS,"{sample}_star_150_")
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        MediumJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params[0]} \
            --readFilesIn {input[0]} {input[1]} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params[1]} \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate
        """

rule align_star_200:
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(STAR_BAMS,"{sample}_star_200_Aligned.sortedByCoord.out.bam")
    log:
        os.path.join(LOGS,"{sample}_star.log")
    params:
        os.path.join(HG38_dir, 'hg38_200'),
        os.path.join(STAR_BAMS,"{sample}_star_200_")
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        MediumJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params[0]} \
            --readFilesIn {input[0]} {input[1]} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params[1]} \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate
        """


#### aggregation rule

rule aggr_align_150:
    input:
        expand(os.path.join(STAR_BAMS,"{sample}_star_150_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "star_150_align.txt")
    threads:
        1
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        touch {output[0]}
        """

rule aggr_align_200:
    input:
        expand(os.path.join(STAR_BAMS,"{sample}_star_200_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "star_200_align.txt")
    threads:
        1
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        touch {output[0]}
        """



rule feature_count_150:
    """feature_counts """
    input:
        expand(os.path.join(STAR_BAMS,"{sample}_star_150_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(RESULTS,"geneCounts_150.out")
    log:
        os.path.join(LOGS,"feature_count.log")
    params:
        os.path.join(HG38_dir, 'gencode.v39.primary_assembly.annotation.gtf')
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        featureCounts -Q 10 -s 0 -T {threads} -p -a {params[0]} -o {output[0]} {input}
        """


rule feature_count_200:
    """feature_counts """
    input:
        expand(os.path.join(STAR_BAMS,"{sample}_star_200_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(RESULTS,"geneCounts_200.out")
    log:
        os.path.join(LOGS,"feature_count.log")
    params:
        os.path.join(HG38_dir, 'gencode.v39.primary_assembly.annotation.gtf')
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=MediumJobMem
    shell:
        """
        featureCounts -Q 10 -s 0 -T {threads} -p -a {params[0]} -o {output[0]} {input}
        """



rule feature_count_cut_150:
    """feature_counts """
    input:
         os.path.join(RESULTS,"geneCounts_150.out")
    output:
        os.path.join(RESULTS,"geneCounts_150.txt")
    log:
        os.path.join(LOGS,"feature_count_cut.log")
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        1
    shell:
        """
        cut -f1,7- {input[0]} | sed 1d > {output[0]}
        """

rule feature_count_cut_200:
    """feature_counts """
    input:
         os.path.join(RESULTS,"geneCounts_200.out")
    output:
        os.path.join(RESULTS,"geneCounts_200.txt")
    log:
        os.path.join(LOGS,"feature_count_cut.log")
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        1
    shell:
        """
        cut -f1,7- {input[0]} | sed 1d > {output[0]}
        """







