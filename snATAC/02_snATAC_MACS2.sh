$ module load purge && module load shared slurm python/3.8.x-anaconda
$ module load gcc && module load gsl && module load hdf5 && module load htslib
$ module load bedtools/2.29.2

$ bamToBed -i possorted_bam.bam > READS.bed

#!/bin/sh

# SBATCH --job-name=MACS.ATAC

# job submit as
# sbatch --cpus-per-task=12 macs.sh

# scontrol show job 384
# scancel 384
# sinfo
# Interactive job
# srun --pty bash

PROJ_DIR=$project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC
OUT_DIR=$project/Neuroinformatics_Core/Konopka_lab/s204365/ATAC/MACS

if [ ! -d ${OUT_DIR} ]; then
  mkdir -p ${OUT_DIR}
fi

#FILE_TRT=`ls -la ~/ATAC/BWA/NG01-ATAQ-Seq-DNA-119*/NG01-ATAQ-Seq-DNA-*.dedupe.bam | awk '{print $9}'`

#  "--f BAMPE " or "--nomodel --f BAM --shift -100 --extsize 200"
# DEFAULT: 0.05. -q, and -p are mutually exclusive.

# If you followed original protocol for ATAC-Seq, you should get Paired-End reads. 
# If so, I would suggest you just use "--format BAMPE" to let MACS2 pileup the whole fragments in general. 
# But if you want to focus on looking for where the 'cutting sites' are, then '--nomodel --shift -100 --extsize 200' should work.

# With -f BAMPE on, MACS2 read the left mate and the insertion length information from BAM file, and discard right mate. With -f BAM, MACS2 only keeps the left mate, so it's definitely not a feature you want for your paired-end data.

macs2 callpeak \
 -n bam_pe \
 --format BAMPE \
 -g mm -q 0.01 --keep-dup all \
 --control $FILE_CTRL \
 --treatment $FILE_TRT \
 -B --trackline --SPMR \
 --call-summits \
 --outdir $OUT_DIR

#--nomodel --f BAM --shift -100 --extsize 200
