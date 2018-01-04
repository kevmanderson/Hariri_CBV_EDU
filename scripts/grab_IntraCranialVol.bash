#!/bin/bash
base_dir=/cluster/nexus/3/users/holmes/Anderson/ArcGet_Output
cd $base_dir
#sub_list=`ls -d ?????????????`
sub_list=$(cat /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_subjects_postQC.txt | tr "\n" " ")
var=1
out_file=/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_ICV_vals.txt
rm ${out_file}
touch ${out_file}
echo 'subj,IntraCranialVol,EstimatedTotalIntraCranialVol,BrainSegVol,BrainSegNotVent,CorticalWhiteMatter,lhCortexVol,rhCortexVol,TotalGray' >> $out_file
for subj in ${sub_list}; do
    echo $subj
    sub_dir=${base_dir}/${subj}/T1w/${subj}/stats
    if [ -f ${sub_dir}/aseg.stats ]; then
        echo -n $subj','>> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep IntraCranialVol | awk -F, '{print $4}'`','                >> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep EstimatedTotalIntraCranialVol | awk -F, '{print $4}'`','  >> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep -w BrainSeg | awk -F, '{print $4}'`','                    >> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep -w BrainSegNotVent | awk -F, '{print $4}'`','             >> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep -w CorticalWhiteMatter | awk -F, '{print $4}'`','         >> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep lhCortexVol | awk -F, '{print $4}'`','                    >> $out_file
        echo -n `cat ${sub_dir}/aseg.stats | grep rhCortexVol | awk -F, '{print $4}'`','                    >> $out_file
        echo `cat ${sub_dir}/aseg.stats | grep TotalGray | awk -F, '{print $4}'` >> $out_file
    else
        echo 'No aseg for ' ${subj}
    fi
done



# R code
usable.brains <- read.csv('/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_ICV_vals.txt', header=TRUE)
ref.table     <- read.csv('/cluster/nexus/3/users/holmes/Anderson/aholmes_2_3_2017_12_32_16.csv', header=TRUE)
fam.id.out    <- as.character(ref.table$FamID[which(ref.table$Label %in% usable.brains$subj)])
write.table(x=fam.id.out, '/cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_usable_subs.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)



# PCA (plink 1.09)
cd /autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/EUR
# Subset the GSP genetic data ton only include those with acceptable TBV values
/autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/plink                           \
    --bfile ./gsp.eur.imputed                                                                   \
    --keep-fam /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_usable_subs.txt            \
    --make-bed --out /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/gsp.n1735.eur.imputed

cd /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI
/autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/plink                   \
    --bfile /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/gsp.n1735.eur.imputed --pca 10 \
    --keep-fam /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_usable_subs.txt    \
    --out /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/pca_n1735_imputed

/autofs/eris/sabuncu/cluster/con/3/users/mert/Tian/GSP/genotype/plink --noweb           \
    --bfile ./gsp.eur --pca 10                                                          \
    --keep-fam /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/GSP_usable_subs.txt    \
    --out /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/pca_n1735_NOTimputed


# PRSice
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
    figname /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/figures/EduYears_Main \
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
    figname /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/figures/GIANTHEIGHT2014 \
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
    figname /cluster/nexus/3/users/holmes/Anderson/FOR_HARIRI/figures/GIANTBMI2015 \
    binary.target F









