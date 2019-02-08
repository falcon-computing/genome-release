#!/usr/bin/env bats

load ../../lib/common

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

@test "printReads output file in directory with no writing permission" {
   run mkdir temp_dir
   run chmod u-w temp_dir
   run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -o temp_dir/file.bam -b report.grp 
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   run rm -rf temp_dir
}

@test "printReads sample name not defined" {
   run ${FCSBIN} printreads -r ${ref_genome} -i ${INPUT_BAM} -o outputdir -b ${REPORT} --sample-name
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
}