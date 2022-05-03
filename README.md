# General_RNA_Seq_Pipeline
Snakemake Pipeline to Extract and Align RNA-Seq Data from fastq and bam files

#### This is a Work in Progress.

* input fastqc to be fixed
* Some code (namely, the sample parsing in samples.smk) and inspiration for the general structure has been borrowed from https://github.com/shandley/hecatomb
* Inputs required are the relevant Fastq/Bam files, which all must be placed in a directory - (the directory they are in must be specified with Reads={directory}).
* Only software requirement is that snakemake be in the $PATH.
* All snakemake commands are assumed to be run from the pipeline directory (such that for offline use the conda envs are available)

# Usage

1. The hg38 fasta and gtf files need to be downloaded (on login node) to the directory specified with HG38_dir:
```console
snakemake -c 1 -s rules/Download_hg38.smk --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes'
```

2. Download the conda envs for offline use:
```console
snakemake -c 1 -s rna_seq_runner.smk --use-conda --config Reads=Bams Output=test hg38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --conda-create-envs-only --conda-frontend conda
```

3. The STAR indices need to be built - 200 (HiSeq 2000) read length index will be created (run on the compute node through slurm):

* TCGA consists of 2x48 paired end reads (after trimming - av read length is 95 bp), EGA is 2x75 paired end reads
* PRJNA453457 has varying read lengths (due to 2 different platforms being used). So safest to use the EGA setting.
* PRJNA523380 has 200 read length (HiSeq 2000)

```console
snakemake -s rules/create_star_indices_hg38.smk --use-conda --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --profile tcga_wgs
```

5. Kallisto and Salmon indices need to be built if they haven't already:

```console
snakemake -s rules/create_salmon_indices.smk --use-conda --config HG38_dir="/hpcfs/users/a1667917/STAR_Ref_Genomes" Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes' --profile wgs_tcga
snakemake -s rules/create_kallisto_indices.smk --use-conda --config Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes' Kallisto_dir="/hpcfs/users/a1667917/Kallisto_Ref_Genomes" --profile wgs_tcga
```

5. Run the pipeline

```console
snakemake -s rna_seq_runner.smk --use-conda --config Input={specify `fastq` or `bam`} Reads={reads_directory} Output={Output_directory} HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' ReadLegnth=200 --profile wgs_tcga
```

* With a Slurm profile (see https://snakemake.readthedocs.io/en/stable/executing/cli.html https://github.com/Snakemake-Profiles/slurm https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated)
* You will need to cd to the pipeline directory in your jobscript before running if you want to run this offline (to use the premade conda envs)
* Use the ReadLength flag to set the STAR index

#### On slurm
```console
# creates the indices
sbatch rna_tcga_create_indices.sh
# runs the pipeline
sbatch rna_tcga_align.sh
```

#### Other notes
* The fastq format for PRJNA523380 is "{sample}_1/2.fastq.gz". Otherwise, this needs to be changed in the fastq_parse_fastq.smk rule
