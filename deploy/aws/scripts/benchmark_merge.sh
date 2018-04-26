
# Defining Global Variables:
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
if [ ! -f "$CURR_DIR/globals.sh" ];then
   echo "$CURR_DIR/globals.sh does not exist"
   exit 1
else
   echo "Importing Variables defined in $CURR_DIR/globals.sh "
   source $CURR_DIR/globals.sh
fi

# FCS, BWA_BIN and GLOBAL VARIABLES defined in globals.sh                                                                                                                                                                                       
BWABIN_VERSION=$($BWA_BIN --version)
GATK_VERSION=$($FCS gatk --version)
echo "BWA VERSION $BWABIN_VERSION"
echo "GATK version $GATK_VERSION"

# Check if FASTQ folder is not empty:                                                                                                                                                                                                           
check_fastq_folder=`ls -1 ${fastq_dir}/*gz | wc -l`
if [ "${check_fastq_folder}" == "0" ];then
   echo "$fastq_dir has no data in it."
   exit 1
fi

array=(`ls -1 fastq/*_1.fastq.gz`)

INSTANCE=$1
CLOUD=$2
sample_id=$3
BACKUP=$4
tmp_dir=/local/${sample_id}
mkdir -p ${tmp_dir}

echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"$0 : Sample ${sample_id} running in instance ${INSTANCE} $CLOUD\" --message \"Sample ${sample_id} running in instance ${INSTANCE} $CLOUD\""
aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "$0 : Sample ${sample_id} running in instance ${INSTANCE} $CLOUD" --message "Sample ${sample_id} running in instance ${INSTANCE} $CLOUD"

function extra_info {
    echo "=========== " >> temp.log
    echo "Check Disk: " >> temp.log
    echo "=========== " >> temp.log
    df -h >> temp.log 
    echo "=================== " >> temp.log
    echo "Check Top Processes: " >> temp.log
    echo "==================== " >> temp.log
    ps aux | sort -nrk 3,3 | head -n 5 >> temp.log
    echo "============ " >> temp.log
    echo "Check DMESG: " >> temp.log
    echo "============ " >> temp.log
    dmesg | tail -n10  >> temp.log
    echo "============ " >> temp.log
}

count=1
for r1 in ${array[@]}; 
  do 
     r2=`echo $r1 | sed 's/_1/_2/g'`
     RG=`basename $r1 | sed 's/_/ /g' | awk '{print $1}'` 
     library=$RG
     logfile=$tmp_dir/${RG}_bwa.log
     echo "${FCS} align -1 $r1  -2 $r2  -o $tmp_dir/$sample_id  -r ${ref_genome} --align-only -S $sample_id -R $RG -L $library -P Illumina -f 2>>$logfile"
     ${FCS} align -1 $r1  -2 $r2  -o $tmp_dir/$sample_id  -r ${ref_genome} --align-only -S $sample_id -R $RG -L $library -P Illumina -f 2>>$logfile
     ALIGNMENT_TIME=`grep -e "bwa mem finishes" $logfile | awk '{print $(NF-1)}'`
     if [ "$count" == "1" ]; then
        echo -e "$RG\t$ALIGNMENT_TIME" > $tmp_dir/${sample_id}_bwa.log
     else
        echo -e "$RG\t$ALIGNMENT_TIME" >> $tmp_dir/${sample_id}_bwa.log
     fi 

     if [[ $? -ne 0 ]];then 
        cat $logfile > temp.log
        extra_info;
        aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "From ${INSTANCE} in ${CLOUD}" --message file://temp.log
        cat temp.log
     fi
     let count=$count+1
  done

status=0

logfile=$tmp_dir/${sample_id}_md.log
echo "${FCS} markDup  -i ${tmp_dir}/$sample_id  -o ${tmp_dir}/${sample_id}.bam  -f 2>>$logfile"
${FCS} markDup  -i ${tmp_dir}/$sample_id  -o ${tmp_dir}/${sample_id}.bam  -f 2>>$logfile 
if [[ $? -ne 0 ]];then
   cat $logfile > temp.log
   extra_info;
   aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "From ${INSTANCE} in ${CLOUD}" --message file://temp.log
   cat temp.log
else
   let status=$status+1
fi

logfile=$tmp_dir/${sample_id}_bqsr.log
echo "${FCS} bqsr --ref ${ref_genome} --input $tmp_dir/${sample_id}.bam -o $tmp_dir/${sample_id}_recal.bam -b  $tmp_dir/${sample_id}/recalibration_report.grp -K ${db138_SNPs} -K ${g1000_gold_standard_indels} -K ${g1000_indels} 2>>$logfile"
${FCS} bqsr --ref ${ref_genome} --input $tmp_dir/${sample_id}.bam -o $tmp_dir/${sample_id}_recal.bam -b  $tmp_dir/${sample_id}/recalibration_report.grp -K ${db138_SNPs} -K ${g1000_gold_standard_indels} -K ${g1000_indels} 2>>$logfile
if [[ $? -ne 0 ]];then
   cat $logfile > temp.log
   extra_info;
   aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "From ${INSTANCE} in ${CLOUD}" --message file://temp.log
   cat temp.log
else
   let status=$status+1
fi

logfile=$tmp_dir/${sample_id}_htc.log
echo "${FCS} htc -r ${ref_genome} --input $tmp_dir/Intel/${sample_id}_recal.bam --output $tmp_dir/${sample_id}_merged.vcf -v 2>>$logfile"
${FCS} htc -r ${ref_genome} --input $tmp_dir/${sample_id}_recal.bam --output $tmp_dir/${sample_id}_merged.vcf -v 2>>$logfile
if [[ $? -ne 0 ]];then
   cat $logfile > temp.log
   extra_info;
   aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "From ${INSTANCE} in ${CLOUD}" --message file://temp.log
   cat temp.log
else
   let status=$status+1
fi

echo "===========================" >> temp.log
echo "BWA VERSION $BWABIN_VERSION" >> temp.log
echo "GATK version $GATK_VERSION"  >> temp.log
echo "===========================" >> temp.log
echo "ALIGNMENT TIME" >> temp.log
echo "==============" >> temp.log
cat $tmp_dir/${sample_id}_bwa.log >> temp.log
echo "=====================" >> temp.log
echo "MARK DUPLICATES  TIME" >> temp.log
echo "=====================" >> temp.log
cat $tmp_dir/${sample_id}_md.log >> temp.log
echo "=====================" >> temp.log
echo "BQSR  TIME"            >> temp.log
echo "=====================" >> temp.log
cat $tmp_dir/${sample_id}_bqsr.log >> temp.log
echo "=====================" >> temp.log
echo "HTC  TIME"            >> temp.log
echo "=====================" >> temp.log
cat $tmp_dir/${sample_id}_htc.log >> temp.log
echo "==============" >>  temp.log
echo "CPU INFO      " >> temp.log
echo "==============" >> temp.log
lscpu >> temp.log
echo "==============" >  temp.log
echo "MEM INFO      " >> temp.log
echo "==============" >> temp.log
cat /proc/meminfo     >> temp.log

TAR=${sample_id}_${INSTANCE}_$(date +%Y%m%d_%H%M%S).tar
echo "Compressing Log Files:"
echo "tar -cvf $TAR  temp.log  ${tmp_dir}/*log  nohup.out"
tar -cvf $TAR  temp.log  ${tmp_dir}/*log  nohup.out  
echo "cp $TAR  $BACKUP"
cp $TAR  $BACKUP

# Cleaning Folder:
if [ "$status" == "3" ];then
   echo "rm -rf $tmp_dir/${sample_id}  temp/*"
   rm -rf $tmp_dir/${sample_id}  temp/*
fi

SUBJECT="$INSTANCE : Pipeline Complete for ${sample_id} on `date`"
echo "$SUBJECT"
echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"$SUBJECT\" --message file://temp.log"
aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "$SUBJECT" --message file://temp.log

echo "aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject \"$0 : Sample ${sample_id} completed in instance ${INSTANCE} $CLOUD\" --message \"Sample ${sample_id} completed in instance ${INSTANCE} $CLOUD\""
aws sns publish --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results --subject "$0 : Sample ${sample_id} completed in instance ${INSTANCE} $CLOUD" --message "Sample ${sample_id} completed in instance ${INSTANCE} $CLOUD"
