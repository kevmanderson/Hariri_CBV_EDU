#!/bin/bash

# project subject directory (on NMR server)
base_dir=/cluster/nexus/3/users/holmes/Anderson/ArcGet_Output_171019
cd ${base_dir}

# subjects that passed basic QC
sub_list=$(cat /cluster/nexus/3/users/holmes/Anderson/PROJECTS/Hariri_CBV_EDU/postQC_subject_list.txt | tr "\n" " ")

# file with relevant FS derived information
out_file=/cluster/nexus/3/users/holmes/Anderson/PROJECTS/Hariri_CBV_EDU/Hariri_CBV_EDU_GSP_vals.txt
rm ${out_file}
touch ${out_file}
echo 'subj,IntraCranialVol,EstimatedTotalIntraCranialVol,BrainSegVol,BrainSegNotVent,CorticalWhiteMatter,lhCortexVol,rhCortexVol,TotalGray' >> $out_file
for subj in ${sub_list}; do
    echo $subj
    sub_dir=${base_dir}/${subj}/${subj}/stats
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
    # catch, but very unlikely to have missing data at this stage
    else
        echo 'No aseg for ' ${subj}
    fi
done
