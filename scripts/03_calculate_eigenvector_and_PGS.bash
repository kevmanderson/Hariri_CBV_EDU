
#!/bin/bash

# PCA (plink 1.09)
genetics_dir=/autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/EUR
# Subset the GSP genetic data to only include those with acceptable TBV values
/autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/plink                           \
    --bfile ${genetics_dir}/gsp.eur.imputed                                                     \
    --keep-fam /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_usable_subs.txt            \
    --make-bed --out /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/gsp.n1735.eur.imputed

# Calculate top 10 PCA eigenvectors
/autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/plink                           \
    --bfile /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/gsp.n1735.eur.imputed --pca 10    \
    --keep-fam /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_usable_subs.txt            \
    --out /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/pca_n1735_imputed



# Calculate PGS with PRSice
sice_dir=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/software/PRSice_v1.25
plink_path=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/software/plink_1.09/plink
gsp_imputed=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/gsp.n1735.eur.imputed


# SSGAC_Edu_Main_2016
edu_gwas=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/summary_stats/SSGAC_Edu_Main_2016.txt
R --file=${sice_dir}/PRSice_v1.25.R --args \
    base ${edu_gwas} \
    target ${gsp_imputed} \
    slower 0 \
    supper 1 \
    sinc 0.1 \
    covary F \
    clump.snps F \
    plink ${plink_path} \
    figname /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/EduYears_Main \
    binary.target F


# GIANTHEIGHT2014
height_gwas=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/summary_stats/GIANTHEIGHT2014/GIANTHEIGHT2014.txt
R --file=${sice_dir}/PRSice_v1.25.R --args \
    base ${height_gwas} \
    target ${gsp_imputed} \
    slower 0 \
    supper 1 \
    sinc 0.1 \
    covary F \
    clump.snps F \
    plink ${plink_path} \
    figname /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GIANTHEIGHT2014_v2 \
    binary.target F

# GIANTHEIGHT2014
bmi_gwas=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/summary_stats/GIANTBMI2015.txt
R --file=${sice_dir}/PRSice_v1.25.R --args \
    base ${bmi_gwas} \
    target ${gsp_imputed} \
    slower 0 \
    supper 1 \
    sinc 0.1 \
    covary F \
    clump.snps F \
    plink ${plink_path} \
    figname /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GIANTBMI2015 \
    binary.target F
