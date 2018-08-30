# ==============================================
# Processing Platinum Trio Genomes Samples      
# ==============================================
# ==============================================
# Processing NA12878 WES                         
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/NA12892-Rep01_S9_L001_R1_001.fastq.gz -2 /local/fastq/NA12892-Rep01_S9_L001_R2_001.fastq.gz -o /local/NA12878/NA12878_marked.bam --rg NA12878 --pl Illumina --lb LIB-NA12878 -f 2>/local/NA12878/NA12878_bwa_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Alignment Failed for NA12878"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878/NA12878_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12878/gatk3/NA12878_bqsr.report -o /local/NA12878/gatk3/NA12878.recal.bam -f   2>/local/NA12878/gatk3/NA12878_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878/NA12878_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12878/gatk4/NA12878_bqsr.report -o /local/NA12878/gatk4/NA12878.recal.bam -f    --gatk4 2>/local/NA12878/gatk4/NA12878_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for NA12878"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878/gatk3/NA12878.recal.bam -o /local/NA12878/gatk3/NA12878_htc.vcf -f   2>/local/NA12878/gatk3/NA12878_htc_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878/gatk4/NA12878.recal.bam -o /local/NA12878/gatk4/NA12878_htc.vcf -v -f  --gatk4 2>/local/NA12878/gatk4/NA12878_htc_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="HTC Failed for NA12878"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing NA12891 WES                         
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/NA12892-Rep01_S9_L001_R1_001.fastq.gz -2 /local/fastq/NA12892-Rep01_S9_L001_R2_001.fastq.gz -o /local/NA12891/NA12891_marked.bam --rg NA12891 --pl Illumina --lb LIB-NA12891 -f 2>/local/NA12891/NA12891_bwa_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Alignment Failed for NA12891"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12891/NA12891_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12891/gatk3/NA12891_bqsr.report -o /local/NA12891/gatk3/NA12891.recal.bam -f   2>/local/NA12891/gatk3/NA12891_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12891/NA12891_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12891/gatk4/NA12891_bqsr.report -o /local/NA12891/gatk4/NA12891.recal.bam -f    --gatk4 2>/local/NA12891/gatk4/NA12891_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for NA12891"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12891/gatk3/NA12891.recal.bam -o /local/NA12891/gatk3/NA12891_htc.vcf -f   2>/local/NA12891/gatk3/NA12891_htc_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12891/gatk4/NA12891.recal.bam -o /local/NA12891/gatk4/NA12891_htc.vcf -v -f  --gatk4 2>/local/NA12891/gatk4/NA12891_htc_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="HTC Failed for NA12891"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing NA12892 WES                         
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/NA12892-Rep01_S9_L001_R1_001.fastq.gz -2 /local/fastq/NA12892-Rep01_S9_L001_R2_001.fastq.gz -o /local/NA12892/NA12892_marked.bam --rg NA12892 --pl Illumina --lb LIB-NA12892 -f 2>/local/NA12892/NA12892_bwa_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Alignment Failed for NA12892"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12892/NA12892_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12892/gatk3/NA12892_bqsr.report -o /local/NA12892/gatk3/NA12892.recal.bam -f   2>/local/NA12892/gatk3/NA12892_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12892/NA12892_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12892/gatk4/NA12892_bqsr.report -o /local/NA12892/gatk4/NA12892.recal.bam -f    --gatk4 2>/local/NA12892/gatk4/NA12892_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for NA12892"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12892/gatk3/NA12892.recal.bam -o /local/NA12892/gatk3/NA12892_htc.vcf -f   2>/local/NA12892/gatk3/NA12892_htc_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12892/gatk4/NA12892.recal.bam -o /local/NA12892/gatk4/NA12892_htc.vcf -v -f  --gatk4 2>/local/NA12892/gatk4/NA12892_htc_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="HTC Failed for NA12892"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing NA12878-Garvan-Vial1 Samples       
# ==============================================
# ==============================================
# Processing NA12878-Garvan-Vial1 WGS                         
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/NA12878-Garvan-Vial1_R1.fastq.gz -2 /local/fastq/NA12878-Garvan-Vial1_R2.fastq.gz -o /local/NA12878-Garvan-Vial1/NA12878-Garvan-Vial1_marked.bam --rg NA12878-Garvan-Vial1 --pl Illumina --lb LIB-NA12878-Garvan-Vial1 -f 2>/local/NA12878-Garvan-Vial1/NA12878-Garvan-Vial1_bwa_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Alignment Failed for NA12878-Garvan-Vial1"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-Garvan-Vial1/NA12878-Garvan-Vial1_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12878-Garvan-Vial1/gatk3/NA12878-Garvan-Vial1_bqsr.report -o /local/NA12878-Garvan-Vial1/gatk3/NA12878-Garvan-Vial1.recal.bam -f   2>/local/NA12878-Garvan-Vial1/gatk3/NA12878-Garvan-Vial1_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-Garvan-Vial1/NA12878-Garvan-Vial1_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12878-Garvan-Vial1/gatk4/NA12878-Garvan-Vial1_bqsr.report -o /local/NA12878-Garvan-Vial1/gatk4/NA12878-Garvan-Vial1.recal.bam -f    --gatk4 2>/local/NA12878-Garvan-Vial1/gatk4/NA12878-Garvan-Vial1_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for NA12878-Garvan-Vial1"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-Garvan-Vial1/gatk3/NA12878-Garvan-Vial1.recal.bam -o /local/NA12878-Garvan-Vial1/gatk3/NA12878-Garvan-Vial1_htc.vcf -f   2>/local/NA12878-Garvan-Vial1/gatk3/NA12878-Garvan-Vial1_htc_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-Garvan-Vial1/gatk4/NA12878-Garvan-Vial1.recal.bam -o /local/NA12878-Garvan-Vial1/gatk4/NA12878-Garvan-Vial1_htc.vcf -v -f  --gatk4 2>/local/NA12878-Garvan-Vial1/gatk4/NA12878-Garvan-Vial1_htc_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="HTC Failed for NA12878-Garvan-Vial1"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing NA12878-HG001 WGS                         
# ==============================================
mkdir -p /local/NA12878-HG001/gatk3
mkdir -p /local/NA12878-HG001/gatk4
/usr/local/falcon-180828/bin/fcs-genome align -r /local/ref/human_g1k_v37.fasta -1 /local/fastq/HG001-NA12878-50x_1.fastq.gz -2 /local/fastq/HG001-NA12878-50x_2.fastq.gz -o /local/NA12878-HG001/NA12878-HG001_marked.bam --rg NA12878-HG001 --pl Illumina --lb LIB-NA12878-HG001 -f 2>/local/NA12878-HG001/NA12878-HG001_bwa_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Alignment Failed for NA12878-HG001"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-HG001/NA12878-HG001_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12878-HG001/gatk3/NA12878-HG001_bqsr.report -o /local/NA12878-HG001/gatk3/NA12878-HG001.recal.bam -f   2>/local/NA12878-HG001/gatk3/NA12878-HG001_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-HG001/NA12878-HG001_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/NA12878-HG001/gatk4/NA12878-HG001_bqsr.report -o /local/NA12878-HG001/gatk4/NA12878-HG001.recal.bam -f    --gatk4 2>/local/NA12878-HG001/gatk4/NA12878-HG001_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for NA12878-HG001"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-HG001/gatk3/NA12878-HG001.recal.bam -o /local/NA12878-HG001/gatk3/NA12878-HG001_htc.vcf -f   2>/local/NA12878-HG001/gatk3/NA12878-HG001_htc_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome htc -r /local/ref/human_g1k_v37.fasta -i /local/NA12878-HG001/gatk4/NA12878-HG001.recal.bam -o /local/NA12878-HG001/gatk4/NA12878-HG001_htc.vcf -v -f  --gatk4 2>/local/NA12878-HG001/gatk4/NA12878-HG001_htc_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="HTC Failed for NA12878-HG001"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing Pair Normal/Tumor for Mutect2      
# ==============================================
# ==============================================
# Processing TCRBOA1 for Mutect2                
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome align -r /local/ref/human_g1k_v37.fasta --sample_sheet /local/SampleSheetMutect2.csv -o /local/ 
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Alignment Failed for NA12878-HG001"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing TCRBOA1-N for Mutect2                 
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/TCRBOA1-N/TCRBOA1-N_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/TCRBOA1-N/gatk3/TCRBOA1-N_bqsr.report -o /local/TCRBOA1-N/gatk3/TCRBOA1-N.recal.bam -f   2>/local/TCRBOA1-N/gatk3/TCRBOA1-N_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/TCRBOA1-N/TCRBOA1-N_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/TCRBOA1-N/gatk4/TCRBOA1-N_bqsr.report -o /local/TCRBOA1-N/gatk4/TCRBOA1-N.recal.bam -f    --gatk4 2>/local/TCRBOA1-N/gatk4/TCRBOA1-N_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for TCRBOA1-N"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
# ==============================================
# Processing TCRBOA1-T for Mutect2                 
# ==============================================
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/TCRBOA1-T/TCRBOA1-T_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/TCRBOA1-T/gatk3/TCRBOA1-T_bqsr.report -o /local/TCRBOA1-T/gatk3/TCRBOA1-T.recal.bam -f   2>/local/TCRBOA1-T/gatk3/TCRBOA1-T_bqsr_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome bqsr -r /local/ref/human_g1k_v37.fasta -i /local/TCRBOA1-T/TCRBOA1-T_marked.bam -K /local/ref/dbsnp_138.b37.vcf -b /local/TCRBOA1-T/gatk4/TCRBOA1-T_bqsr.report -o /local/TCRBOA1-T/gatk4/TCRBOA1-T.recal.bam -f    --gatk4 2>/local/TCRBOA1-T/gatk4/TCRBOA1-T_bqsr_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="BQSR Failed for TCRBOA1-T"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
/usr/local/falcon-180828/bin/fcs-genome mutect2 -r /local/ref/human_g1k_v37.fasta  -n /local/TCRBOA1-T-N/gatk3/TCRBOA1-T-N.recal.bam -t /local/TCRBOA1-T-T/gatk3/TCRBOA1-T-T.recal.bam  -o /local/mutect2-TCRBOA1-T/gatk3/TCRBOA1-T_mutect2.vcf --dbsnp /local/ref/dbsnp_138.b37.vcf --cosmic /local/ref/b37_cosmic_v54_120711.vcf -f   2>/local/mutect2-TCRBOA1-T/gatk3/TCRBOA1-T_mutect2_201808301535666972.log
/usr/local/falcon-180828/bin/fcs-genome mutect2 -r /local/ref/human_g1k_v37.fasta -n /local/TCRBOA1-T/gatk4/TCRBOA1-T-N.recal.bam  --normal_name TCRBOA1-T-Normal   -t /local/TCRBOA1-T/gatk4/TCRBOA1-T-T.recal.bam --tumor_name  TCRBOA1-T-Tumor -o /local/mutect2_TCRBOA1-T/gatk4/TCRBOA1-T_mutect2.vcf -f   -p /local/mutect2_inputs/mutect_gatk4_pon.vcf  -m /local/mutect2_inputs/af-only-gnomad.raw.sites.b37.vcf.gz --gatk4 -f 2>/local/mutect2_TCRBOA1-T/gatk4/TCRBOA1-T_mutect2_201808301535666972.log
if [ $? -ne 0 ]; then
   ERROR_MESSAGE="Mutect2 Failed for TCRBOA1-T"
   echo $ERROR_MESSAGE
   echo $ERROR_MESSAGE >> /local/error.log
   echo "aws sns publish  --region us-east-1   --topic-arn arn:aws:sns:us-east-1:520870693817:Genomics_Pipeline_Results  --subject "ERROR : From CPU ID: 201808301535666972 running in merlin3 --message "file://error.log" > /local/sender.sh
   source /local/sender.sh
   return 1
fi
