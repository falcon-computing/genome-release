
# Executable Location:
echo "Check if fcs-exome exists:"
#FCS_DIR="/local/auto/falcon"
FCS_DIR="/genome/disk1/benchmarking/auto/falcon"
#FCS_DIR="/usr/local/falcon"
export FCS="${FCS_DIR}/bin/fcs-genome"
if [ ! -f "${FCS}" ];then
   echo -e "${FCS} does not exist.\n"
   echo "Possible Solutions:"
   echo "1. Provide new path in the globals BASH script and re-run again"
   echo -e "2. Download the executables from AWS S3 and set the location properly\n"
   exit 1
else
   echo -e "${FCS} OK\n"
fi

echo "Verify License PATH: $LM_LICENSE_PATH"
export LM_LICENSE_PATH="${FCS_DIR}/license.lic"
if [ ! -f "$LM_LICENSE_PATH" ];then
   echo "License LM_LICENSE_PATH : $LM_LICENSE_PATH does not exist"
   exit 1
else 
   echo -e "$LM_LICENSE_PATH OK\n"
fi

echo "Check BWA_BIN and GATK in FCS:" 
export BWA_BIN="${FCS_DIR}/tools/bin/bwa-bin"
if [ ! -f "${BWA_BIN}" ];then
   echo "${BWA_BIN} does not exist."
   exit 1
else
   echo -e "${BWA_BIN} OK"
fi

export GATK="${FCS_DIR}/tools/package/GenomeAnalysisTK.jar"
if [ ! -f "${GATK}" ];then
   echo "${GATK} does not exist."
   exit 1
else
   echo -e "${GATK} OK\n"
fi

#Directories setup
export fastq_dir=$(pwd)/fastq
export output_dir=$(pwd)
export ref_dir=$(pwd)/ref

echo "Checking Global Variables"
# Check if the FASTQ folder exists:
if [ ! -d ${fastq_dir} ];then
   echo "${fastq_dir} does not exist. Creating it now"
   echo "mkdir ${fastq_dir}"
   mkdir ${fastq_dir}
   echo "Populating ${fastq_dir} with data from AWS S3 "
   echo "aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --include=\"NA12878-I47_S6_L002_R*_001.fastq.gz\""
   aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --exclude "*"  --include="NA12878-I47_S6_L002_R*_001.fastq.gz"
   echo "aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --include=\"NA12878-I33_S2_L002_R*_001.fastq.gz\"" 
   aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --exclude "*"  --include="NA12878-I33_S2_L002_R*_001.fastq.gz"
   echo "aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir} --recursive --exclude "*" --include=\"NA12878-Garvan-Vial1_R*.fastq.gz\" --exclude \"*\""
   aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --exclude "*"  --include="NA12878-Garvan-Vial1_R*.fastq.gz"
   if [ -d ${fastq_dir} ];then
      echo "${fastq_dir} OK"
   fi
else
   if [ -z ${fastq_dir} ];then
      echo "Populating ${fastq_dir} with data from AWS S3 "
      echo "aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --include=\"NA12878-I47_S6_L002_R*_001.fastq.gz\""
      aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --exclude "*"  --include="NA12878-I47_S6_L002_R*_001.fastq.gz"
      echo "aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --include=\"NA12878-I33_S2_L002_R*_001.fastq.gz\"" 
      aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --exclude "*"  --include="NA12878-I33_S2_L002_R*_001.fastq.gz"
      echo "aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir} --recursive --exclude "*" --include=\"NA12878-Garvan-Vial1_R*.fastq.gz\" --exclude \"*\""
      aws s3 cp s3://fcs-genome-data/fastq/wgs/ ${fastq_dir}  --recursive --exclude "*"  --include="NA12878-Garvan-Vial1_R*.fastq.gz"
   fi 
   echo "${fastq_dir} OK"
fi  

# Check if the Output folder exists:
if [ ! -d ${output_dir} ];then
   echo "${output_dir} does not exist"
   echo "Creating ${output_dir}:"
   mkdir ${output_dir}   
   if [ -d ${output_dir} ];then
      echo "${output_dir} OK"
   fi
else
   echo "${output_dir} OK"
fi  
# Check if the Reference folder exists:
if [ ! -d ${ref_dir} ];then
   echo "${ref_dir} does not exist"
   echo "Creating ${ref_dir}"
   echo "mkdir ${ref_dir}"   
   mkdir ${ref_dir}
   echo "Populating the ${ref_dir} with data from AWS S3:"
   echo "aws s3 cp --recursive s3://fcs-genome-data/ref/ ${ref_dir}"
   aws s3 cp --recursive s3://fcs-genome-data/ref/ ${ref_dir}    
   if [ -d ${ref_dir} ];then
      echo "${ref_dir} OK"
   fi
else
   if [ -z ${ref_dir} ];then
      echo "Populating the ${ref_dir} with data from AWS S3:"
      echo "aws s3 cp --recursive s3://fcs-genome-data/ref/ ${ref_dir}"
      aws s3 cp --recursive s3://fcs-genome-data/ref/ ${ref_dir}
   fi
   echo "${ref_dir} OK" 
fi
echo -e "\n"

# Reference Genome Info
export ref_genome=$ref_dir/human_g1k_v37.fasta
export db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
export g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
export g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
# Check if Key Files are available:
check_ref_files=`ls -1 $ref_genome $db138_SNPs $g1000_indels $g1000_gold_standard_indels | wc -l`
if [ "${check_ref_files}" == "4" ];then
   echo -e "Files\n"
   echo -e "$ref_genome OK \n$db138_SNPs OK \n$g1000_indels OK \n$g1000_gold_standard_indels OK\n\n"
else
   echo "Some Required Files in ref/ folder are missing"
   exit 1
fi

check_input() {
  filename=$1;
  if [ ! -f $filename ]; then
    echo "Input file $filename does not exist";
    exit 1;
  fi;
}

check_output() {
  filename=$1;
  if [ -f $filename ]; then
    echo "Output file $filename already exists";
    exit 1;
  fi;
  # Check if output folder writtable
  echo "1" > $filename;
  if [ ! -f $filename ]; then
    echo "Cannot write to folder $(dirname $filename)";
    exit 1;
  fi;
  rm $filename;
}

# check if output dir exists and is writable
check_output_dir() {
  dir=$1;
  if [[ ! -d "$dir" ]]; then
    echo "Output dir $dir does not exists";
    exit 1;
  fi;
  echo "1" > $dir/.test;
  if [ ! -f $dir/.test ]; then
    echo "Cannot write to folder $dir";
    exit 1;
  fi;
  rm $dir/.test;
}

# create output dir if it does not exist 
create_dir() {
  dir=$1
  if [[ ! -d "$dir" ]]; then
    mkdir $dir &> /dev/null;
    if [ $? -ne 0 ]; then
      echo "Cannot create dir $dir";
      exit 1;
    fi
  fi;
}
