"""
Snakefile for downloading STAR index for hg38
snakemake -c 1 -s rules/Download_hg38.smk --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes'
"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '..', 'config', 'config.yaml')

if config['HG38_dir'] is None:
    HG38_dir = "/hpcfs/users/a1667917/STAR_Ref_Genomes"
else:
    HG38_dir = config["HG38_dir"]

# needs to be created before (should exist)
if not os.path.exists(os.path.join(HG38_dir)):
  os.makedirs(os.path.join(HG38_dir))

rule all:
    input:
        os.path.join(HG38_dir, 'download_hg38_fasta.dlflag'),
        os.path.join(HG38_dir, 'download_hg38_gtf.dlflag'),
        os.path.join(HG38_dir, 'unzip.dlflag')

rule download_hg_38_fasta:
    """Rule to Download hg38 fasta."""
    output:
        os.path.join(HG38_dir,'download_hg38_fasta.dlflag'),
        os.path.join(HG38_dir, 'GRCh38.primary_assembly.genome.fa.gz')
    threads:
        1
    shell:
        """
        wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz" -O GRCh38.primary_assembly.genome.fa.gz
        mv "GRCh38.primary_assembly.genome.fa.gz" {output[1]}
        touch {output[0]}
        """

rule download_hg_38_gtf:
    """Rule to Download hg38 gtf."""
    output:
        os.path.join(HG38_dir,'download_hg38_gtf.dlflag'),
        os.path.join(HG38_dir, 'gencode.v39.primary_assembly.annotation.gtf.gz')
    threads:
        1
    shell:
        """
        wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz" -O  gencode.v39.primary_assembly.annotation.gtf.gz
        mv "gencode.v39.primary_assembly.annotation.gtf.gz" {output[1]}
        touch {output[0]}
        """

# no -C on the hpc

rule unzip:
    """gunzip files."""
    input:
        os.path.join(HG38_dir, 'GRCh38.primary_assembly.genome.fa.gz'),
        os.path.join(HG38_dir, 'gencode.v39.primary_assembly.annotation.gtf.gz')
    output:
        os.path.join(HG38_dir, 'GRCh38.primary_assembly.genome.fa'),
        os.path.join(HG38_dir, 'gencode.v39.primary_assembly.annotation.gtf'),
        os.path.join(HG38_dir, 'unzip.dlflag')
    threads:
        1
    shell:
        """
        gunzip {input[0]}
        gunzip {input[1]}
        touch {output[2]}
        """
