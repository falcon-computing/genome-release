#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $CURR_DIR/cloud-helper.sh
source $CURR_DIR/global.bash

WORKDIR=`pwd`

start_ts=$(date +%s)

echo -e "Testing feature in fcs-genome "
echo -e "=============================\n"
$CURR_DIR/bats/bats features_test/ >> test.log

echo -e "Regression Test In Progress"
echo -e "===========================\n"

chmod ag+wr test.log  nohup.out
if [ ! -d `pwd`/log ];then
   mkdir `pwd`/log
   chmod 777 -R `pwd`/log/
fi

echo -e "DNA samples"
echo -e "===========================\n"
array=(A15_sample CDMD1015_sample DSDEX72_sample)
#array=(NA12878)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $CURR_DIR/bats/bats results_test/  >> test.log
  done

rm -rf output.bam 

echo -e "Pair sample for Mutect2"
echo -e "===========================\n"
array=(mutect2_sample)
for id in ${array[@]}
  do
    echo "Processing $id"
    export id=$id
    $CURR_DIR/bats/bats mutect2_test/ >> test.log
  done

end_ts=$(date +%s)
echo "Time taken: $((end_ts - start_ts))s"  >> test.log

echo "Cleaning Folder:"
array=(A15_sample CDMD1015_sample DSDEX72_sample mutect2_sample)
for id in ${array[@]}
  do
    echo "rm -rf /local/work_dir/$id*"
          rm -rf /local/work_dir/$id*
  done

#aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --region us-east-1 --subject "Regression Test on Merlin3" --message file://test.log


