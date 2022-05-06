#!/bin/bash -l

#SBATCH --job-name=rna_tcga_snkg
#SBATCH --mail-user=george.bouras@adelaide.edu.au
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --err="rna_tcga_indices_snk.err"
#SBATCH --output="rna_tcga_indices_snk.out"

# Resources allocation request parameters
#SBATCH -p batch
#SBATCH -N 1                                    
#SBATCH -n 1                                                   
#SBATCH --time=12:00:00                                   
#SBATCH --mem=1GB                                             


SNK_DIR="/hpcfs/users/a1667917/Kenny/RNA_Seq_SRA/General_RNA_Seq_Pipeline"
PROF_DIR="/hpcfs/users/a1667917/snakemake_slurm_profile"

cd $SNK_DIR

module load Anaconda3/2020.07
conda activate snakemake_clean_env

snakemake -s rules/create_star_indices_hg38.smk --use-conda --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --profile wgs_tcga

conda deactivate