
# Executable Location:
echo "Check if fcs-exome exists:"
FCS_DIR="/usr/local/falcon"
export FCS="${FCS_DIR}/bin/fcs-genome"
if [ ! -f "${FCS}" ];then
   echo -e "${FCS} does not exist.\n"
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
export fastq_dir=/local/fastq
export output_dir=/local/
export ref_dir=/local/ref

echo "Checking Global Variables"
# Check if the FASTQ folder exists:
if [ ! -d ${fastq_dir} ];then
   echo "${fastq_dir} does not exist. Creating it now"
   echo "mkdir ${fastq_dir}"
   mkdir ${fastq_dir}
   if [ -d ${fastq_dir} ];then
      echo "${fastq_dir} OK"
   fi
fi

# Check if the Output folder exists:
if [ ! -d ${output_dir} ];then
   echo "${output_dir} does not exist"
   echo "Creating ${output_dir}:"
   mkdir ${output_dir}   
   if [ -d ${output_dir} ];then
      echo "${output_dir} OK"
   fi
fi 
 
# Check if the Reference folder exists:
if [ ! -d ${ref_dir} ];then
   echo "${ref_dir} does not exist"
   echo "Creating ${ref_dir}"
   echo "mkdir ${ref_dir}"   
   mkdir ${ref_dir}
   if [ -d ${ref_dir} ];then
      echo "${ref_dir} OK"
   fi
fi
echo -e "\n"

# Reference Genome Info
export ref_genome=$ref_dir/human_g1k_v37.fasta
export db138_SNPs=$ref_dir/dbsnp_138.b37.vcf
export g1000_indels=$ref_dir/1000G_phase1.indels.b37.vcf
export g1000_gold_standard_indels=$ref_dir/Mills_and_1000G_gold_standard.indels.b37.vcf
export cosmicVCF=$ref_dir/b37_cosmic_v54_120711.vcf

# Check if Key Files are available:
check_ref_files=`ls -1 $ref_genome $db138_SNPs $g1000_indels $g1000_gold_standard_indels ${cosmicVCF} | wc -l`
if [ "${check_ref_files}" == "5" ];then
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
