$ module load cellrangerATAC/2.0.0

$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/ATAC

# A1
$ cellranger-atac count --id=A1 \
                        --reference=/work/Neuroinformatics_Core/s204365/References/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                        --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/atac1,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/atac2 \
                        --sample=AOP0A1 \
                        --localcores=8 \
# A2
$ cellranger-atac count --id=A2 \
                        --reference=/work/Neuroinformatics_Core/s204365/References/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                        --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/atac1,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/atac2 \
                        --sample=AOP0A2 \
                        --localcores=8 \

# A3
$ cellranger-atac count --id=A3 \
                        --reference=/work/Neuroinformatics_Core/s204365/References/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                        --fastqs=/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/atac1,/project/Neuroinformatics_Core/Konopka_lab/shared/AO_shared/for_alia/atac2 \
                        --sample=AOP0A3 \
                        --localcores=8 \
