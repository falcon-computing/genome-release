import re
import argparse
from datetime import datetime
import sys

input_file=sys.argv[1]
output_file=sys.argv[2]

#parser=argparse.ArgumentParser()
#parser.add_argument("-i","--id", help="display processes for particular run_id",type=str,default="NONE")
#args=parser.parse_args()

output=open(input_file,'r')
previous_line=""

#run_all=0
#if args.id=="NONE":
#  run_all=1
process=""
with open(output_file,'wb') as file:
  for line in output:
    output_match_fcs_version = re.search('fcs-genome version ([0-9a-zA-Z\.\-]*)',line)
    if output_match_fcs_version:
      fcs=output_match_fcs_version.group(1)
      file.write('\n')
      file.write('fcs-genome version '+fcs)
    output_match_release_version = re.search('release version ([0-9a-zA-Z\.\-]*)',line)
    if output_match_release_version:
      rel=output_match_release_version.group(1)
      file.write('\n')
      file.write('Release version '+rel)
    output_match_bwa_version = re.search('bwa version ([0-9a-zA-Z\.\-]*)',line)
    if output_match_bwa_version:
      bwa=output_match_bwa_version.group(1)
      file.write('\n')
      file.write('BWA version '+bwa)
    output_match_gatk_version = re.match('gatk version ([0-9a-zA-Z\.\-]+)',line)
    if output_match_gatk_version:
        gatk=output_match_gatk_version.group(1)
        file.write('\n')
        file.write('GATK version '+gatk+'\n'+'\n')
        file.write('Sample'+','+'BWA'+','+'MD'+','+'BQSR'+','+'PR'+','+'HTC'+','+'Total'+'\n') 
    output_match_run = re.search('RUN: (\w*)',line)
    if output_match_run:
      run=output_match_run.group(1)
      file.write('\n')
      file.write('RUN '+str(int(run)+1))
    output_match_id = re.search('ID: (\w*)',line)
    if output_match_id:
      id=output_match_id.group(1)
      file.write('\n')
      file.write(id+',')  
    output_match_process = re.search('INFO: Start doing (\w*\s\w*)',line)
    if output_match_process:
      process=output_match_process.group(1)
    output_match_perf = re.search('INFO: '+process+' finishes in (\w*) seconds',line)
    if output_match_perf:
      time=output_match_perf.group(1)
      file.write(time+',')
    output_match_fail_bwa = re.search('ERROR: bwa mem failed',line)
    if output_match_fail_bwa:
      file.write('Failed'+','+'Failed'+',')
    output_match_fail_bqsr = re.search('Failed base recalibration',line)
    if output_match_fail_bqsr:
      file.write('Failed'+',')
    output_match_fail_pr = re.search('Failed print reads',line)
    if output_match_fail_pr:
      file.write('Failed'+',')
    output_match_fail_htc = re.search('Failed haplotype caller',line)
    if output_match_fail_htc:
      file.write('Failed'+',')
        
