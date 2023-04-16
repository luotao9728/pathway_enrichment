#!/bin/bash
#SBATCH -N 4
#SBATCH -p RM
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=128

# activate local batch environment
source ~/.bashrc
echo "*******************************************************************************************"
echo "Bash environment activated!"
echo "*******************************************************************************************"

# activate conda environment for the pipeline
conda activate annotation
echo "*******************************************************************************************"
echo "Conda environment activated!"
echo "*******************************************************************************************"

# change direction to data location
cd /ocean/projects/bio230007p/tluo1/data_processing_new

# Trimmomatic
echo "*******************************************************************************************"
echo "Start trimming illumina short RNA-seq!"
echo "*******************************************************************************************"
trimmomatic PE -phred33 liver_illumina_R1.fastq liver_illumina_R2.fastq trimmed_liver_illumina_R1.fastq output_forward_unpaired.fq.gz trimmed_liver_illumina_R2.fastq output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

rm output_forward_unpaired.fq.gz
rm output_reverse_unpaired.fq.gz
echo "*******************************************************************************************"
echo "Finish trimming illumina short RNA-seq!"
echo "*******************************************************************************************"

# LorDEC
## lordec-correct
echo "*******************************************************************************************"
echo "Start correcting PacBio long RNA-seq!"
echo "*******************************************************************************************"

lordec-correct -2 trimmed_liver_illumina_R1.fastq trimmed_liver_illumina_R2.fastq -k 19 -s 3 -i liver_pacbio.fastq -o liver_pacbio_lordec_corrected.fasta

echo "*******************************************************************************************"
echo "Finish correcting PacBio long RNA-seq!"
echo "*******************************************************************************************"

## lordec-trim
echo "*******************************************************************************************"
echo "Start trimming corrected PacBio long RNA-seq!"
echo "*******************************************************************************************"

lordec-trim -i liver_pacbio_lordec_corrected.fasta -o trim_liver_pacbio.fasta

echo "*******************************************************************************************"
echo "Finish trimming corrected PacBio long RNA-seq!"
echo "*******************************************************************************************"

## lordec-trim-split
echo "*******************************************************************************************"
echo "Start trim-splitting corrected PacBio long RNA-seq!"
echo "*******************************************************************************************"

lordec-trim-split -i liver_pacbio_lordec_corrected.fasta -o trim_split_liver_pacbio.fasta

echo "*******************************************************************************************"
echo "Finish trim-splitting corrected PacBio long RNA-seq!"
echo "*******************************************************************************************"

# Magic-BLAST
echo "*******************************************************************************************"
echo "Start alignments! Please be patient, this might take few (or many) hours!"
echo "*******************************************************************************************"

mkdir alignments

makeblastdb -in ref_genome_horse.fna -out horse_reference -parse_seqids -dbtype nucl

magicblast -query trim_split_liver_pacbio.fasta -db horse_reference -out alignments/liver_align_pacbio.sam
magicblast -query trimmed_liver_illumina_R1.fastq -db horse_reference -infmt fastq -out alignments/liver_align_illumina_R1.sam
magicblast -query trimmed_liver_illumina_R2.fastq -db horse_reference -infmt fastq -out alignments/liver_align_illumina_R2.sam

echo "*******************************************************************************************"
echo "Finish alignments!"
echo "*******************************************************************************************"

# deactivate conda environment
conda deactivate