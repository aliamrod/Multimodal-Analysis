$ wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
$ tar xzgf refdata-gex-mm10-2020-A.tar.gz

# AO1
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO1 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND26,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO1 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO2
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO2 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND26,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO2 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO3
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO3 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND26,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO3 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO4
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO4 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND26,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO4 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO5
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO5 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND28,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO5 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO6
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO6 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND28,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO6 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO7
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO7 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND28,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO7 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48

# AO8
$ module load cellranger/7.0.0
$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger

$ cellranger count --id=AO8 \
   --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND28,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/ROUND30 \
   --sample=AO8 \
   --transcriptome=/work/Neuroinformatics_Core/s204365/References/refdata-gex-mm10-2020-A
   --localcores=72 
   --localmem=48
