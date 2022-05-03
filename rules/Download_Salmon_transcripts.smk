"""
Snakefile for downloading files for Salmon indices
snakemake -c 1 -s rules/Download_Salmon_transcripts.smk --config Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes' 
"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '..', 'config', 'config.yaml')

if config['Salmon_dir'] is None:
    Salmon_dir = "/hpcfs/users/a1667917/Salmon_Ref_Genomes"
else:
    Salmon_dir = config["Salmon_dir"]

# needs to be created before (should exist)
if not os.path.exists(os.path.join(Salmon_dir)):
  os.makedirs(os.path.join(Salmon_dir))

rule all:
    input:
        os.path.join(Salmon_dir,'download_transcripts.dlflag'),
        os.path.join(Salmon_dir,'download_fasta.dlflag'),
        os.path.join(Salmon_dir, 'unzip.dlflag')

rule download_salmon_transcripts:
    """Rule to Download encode transcripts."""
    output:
        os.path.join(Salmon_dir,'download_transcripts.dlflag'),
        os.path.join(Salmon_dir, 'gencode.v39.transcripts.fa.gz')
    threads:
        1
    shell:
        """
        wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz" -O gencode.v39.transcripts.fa.gz
        mv "gencode.v39.transcripts.fa.gz" {output[1]}
        touch {output[0]}
        """

rule download_salmon_fasta:
    """Rule to Download encode fasta."""
    output:
        os.path.join(Salmon_dir,'download_fasta.dlflag'),
        os.path.join(Salmon_dir, 'GRCh38.primary_assembly.genome.fa.gz')
    threads:
        1
    shell:
        """
        wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz" -O GRCh38.primary_assembly.genome.fa.gz
        mv "GRCh38.primary_assembly.genome.fa.gz" {output[1]}
        touch {output[0]}
        """

rule unzip_salmon:
    """gunzip files."""
    input:
        os.path.join(Salmon_dir, 'gencode.v39.transcripts.fa.gz'),
        os.path.join(Salmon_dir, 'GRCh38.primary_assembly.genome.fa.gz')
    output:
        os.path.join(Salmon_dir, 'gencode.v39.transcripts.fa'),
        os.path.join(Salmon_dir, 'GRCh38.primary_assembly.genome.fa'),
        os.path.join(Salmon_dir, 'unzip.dlflag')
    threads:
        1
    shell:
        """
        gunzip {input[0]}
        gunzip {input[1]}
        touch {output[2]}
        """



##############################
# refgenie wont work on hpc for some reason (likely firewall) so build manually
###############################
# if it does workuse this

# rule all:
#     input:
#         os.path.join(Salmon_dir,'download_salmon_indices.dlflag'),
#         os.path.join(Salmon_dir, 'genome_config.yaml') 


# rule download_salmon_index:
#     """Rule to Download salmon transcripts."""
#     output:
#         os.path.join(Salmon_dir, 'genome_config.yaml'),
#         os.path.join(Salmon_dir,'download_salmon_indices.dlflag')
#     conda:
#         os.path.join('..', 'envs','refgenie.yaml')
#     threads:
#         1
#     script:
#         '../scripts/refgenie.py'


