#!/bin/bash

function get_fastq_files(){
   local sample_id=$1
   local tag=$2   

   if [[ "${sample_id}" == "NA12878" ]] && [[ "${tag}" == "WES" ]];then 
      export NA12878_r1=NA12878-Rep01_S1_L001_R1_001.fastq.gz
      export NA12878_r2=NA12878-Rep01_S1_L001_R2_001.fastq.gz
      export r1=${NA12878_r1}
      export r2=${NA12878_r2}
   fi

   if [[ "${sample_id}" == "NA12878" ]] && [[ "${tag}" == "WES" ]];then
      export NA12891_r1=NA12891-Rep01_S5_L001_R1_001.fastq.gz
      export NA12891_r2=NA12891-Rep01_S5_L001_R2_001.fastq.gz
      export r1=${NA12891_r1}
      export r2=${NA12891_r2}
   fi

   if [[ "${sample_id}" == "NA12878" ]] && [[ "${tag}" == "WES" ]];then
      export NA12892_r1=NA12892-Rep01_S9_L001_R1_001.fastq.gz
      export NA12892_r2=NA12892-Rep01_S9_L001_R2_001.fastq.gz
      export r1=${NA12892_r1}
      export r2=${NA12892_r2}
   fi

   if [[ "${sample_id}" == "NA12878-Garvan-Vial1" ]] && [[ "${tag}" == "WGS" ]];then
      export garvan_r1=NA12878-Garvan-Vial1_R1.fastq.gz
      export garvan_r2=NA12878-Garvan-Vial1_R2.fastq.gz
      export r1=${garvan_r1}
      export r2=${garvan_r2}
   fi

   if [[ "${sample_id}" == "NA12878-HG001" ]] && [[ "${tag}" == "WGS" ]];then
      export hg001_r1=HG001-NA12878-50x_1.fastq.gz
      export hg001_r2=HG001-NA12878-50x_2.fastq.gz
      export r1=${hg001_r1}
      export r2=${hg001_r2}
   fi

   if [[ "${sample_id}" == "TCRBOA1" ]] && [[ "${tag}" == "baylor" ]];then
      export normal_r1=TCRBOA1-N-WEX.read1.fastq.gz
      export normal_r2=TCRBOA1-N-WEX.read2.fastq.gz
      export tumor_r1=TCRBOA1-T-WEX.read1.fastq.gz
      export tumor_r2=TCRBOA1-T-WEX.read2.fastq.gz
   fi

}
