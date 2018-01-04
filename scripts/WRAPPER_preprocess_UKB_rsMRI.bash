#! /bin/bash


data_dir=/nexsan/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/UKB_MRI
cd ${data_dir}
ls -d ??????? > /nexsan/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/UKB_MRI/subject_list.txt
readarray -t subj_list < /nexsan/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/UKB_MRI/subject_list.txt
num_subs=`more /nexsan/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/UKB_MRI/subject_list.txt | wc -l`

sub_idx=0
while [ $sub_idx -le $num_subs ]; do
    # get the current subject ready for processing
    sub=${subj_list[$sub_idx]}
    
    bsub_file=${data_dir}/${sub}/qsub_${subj}.txt
    script_dir=/nexsan/holmes/Open_Data/DATA_UKBIOBANK/SCRIPTS/fMRI
    #${script_dir}/preprocess_UKB_rsMRI.csh -s ${sub} -low_f .01 -high_f .10
    my_jobs=`bjobs | wc -l`
    if [ "$my_jobs" -lt "100" ];
    then
        echo ${sub}
        echo ${script_dir}/preprocess_UKB_rsMRI.csh -s ${sub} -low_f .01 -high_f .10 > $bsub_file
        bsub -n 1 -R span[hosts=1] < $bsub_file
        let "sub_idx+=1"
        echo $sub_idx
    fi
done







