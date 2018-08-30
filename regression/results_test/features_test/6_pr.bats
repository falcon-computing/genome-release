#!/usr/bin/env bats

load ../global

@test "printReads without input arg" {
   run ${FCSBIN} printreads
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads without -bqsr specified" {
   run ${FCSBIN} printreads -r ${ref_genome}  i ${INPUT_BAM} -o out.bam
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads input bqsr does not exist" {
   run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -o outputdir -b report.grp
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]
}

@test "printReads without -i specified" {
   run ${FCSBIN} printreads -r ${ref_genome}  -o outputdir -b ${REPORT}
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome printreads"* ]]   
   run rm -rf log/  outputdir/
}

@test "printReads output file directory not writeable" {
   if [ `hostname` == "merlin3" ];then 
      run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -o fake/ -b ${REPORT}
      [ "$status" -ne 0 ]
      [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
   else
      skip
   fi
}
