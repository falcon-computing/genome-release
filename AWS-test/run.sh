

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo -e "\n"
echo "====================================================="
echo "Testing Options Available in fcs-genome on AWS Server"
echo "====================================================="
echo -e "\n"

echo "Check if SMALL folder exist:"
if [ ! -s ${DIR}/SMALL ]; then
   echo "Test folder ${DIR}/SMALL is not available in the Test Folder"
   exit 1
else
   echo "${DIR}/SMALL is in Test Folder"
fi

echo "Check if the Input FASTQ Files exist:"
fastq_dir=/genome/fastq
fastq1=$fastq_dir/small_1.fastq.gz
fastq2=$fastq_dir/small_2.fastq.gz
if [[  ! -f ${fastq1} ]] && [[  ! -f ${fastq2} ]];then
   echo "$fastq1 and $fastq2 do not exist"
   exit 1
else
   echo "FASTQ Files exist in the path"
fi

echo "Check if Reference Genome exists:"
REF="/local/ref/human_g1k_v37.fasta"
if [ ! -f ${REF} ]; then
   echo "${REF} is not available"
   exit 1
else
   echo "Reference Genome ${REF} found"
fi

echo "Check if executable fcs-genome is defined:"
FCS="/curr/software/falcon-genome/latest/bin/fcs-genome"
if [ ! -f ${FCS} ]; then 
   echo "${FCS} is not available"   
   exit 1
else 
   echo -e "Executable ${FCS} found\n"
fi

echo -e "Begin Test\n"

for bat in $(ls $DIR/cases/*.bats); do
  ${DIR}/bats/bin/bats $bat
  if [ "$?" -ne 0 ]; then
    exit 1
  fi
done
