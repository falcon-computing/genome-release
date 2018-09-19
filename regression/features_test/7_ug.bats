#!/usr/bin/env bats

load ../global

@test "UG without input arg" {
   run ${FCSBIN} ug
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG without Reference" {
   run ${FCSBIN} ug -i ${INPUT_BAM} -o ${OUTPUT_VCF}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG without -i specified" {
   run ${FCSBIN} ug -r ${ref_genome} -o ${OUTPUT_VCF}
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG without -o specified" {
   run ${FCSBIN} ug -r ${ref_genome}  -i ${INPUT_BAM}  
   [ "$status" -ne 0 ]
   [[ "${output}" == *"ERROR"* ]]
   [[ "${output}" == *"fcs-genome ug"* ]]
}

@test "UG output directory does not exist" {
   if [ `hostname` == "merlin3" ];then 
      run ${FCSBIN} ug -r ${ref_genome} -i ${INPUT_BAM} -o ${WORK_DIR}/fake/output.vcf
      [ "$status" -ne 0 ]
      [[ "${output}" == *"ERROR: Cannot write to output path"* ]]
   else
      skip
   fi
}

@test "UG input file does not exist" {
   if [ `hostname` == "merlin3" ];then
      run ${FCSBIN} ug -r ${ref_genome} -i ${WORK_DIR}/doesnotexist -o ${OUTPUT_VCF}
      [ "$status" -ne 0 ]
      [[ "${output}" == *"ERROR"* ]]
   else
     skip
   fi
}

@test "UG sample tag option set but not defined" {
   run ${FCSBIN} ug -r ${ref_genome} -i ${INPUT_BAM} -o ${OUTPUT_VCF} --sample-tag
   [ "$status" -ne 0 ]
   #[[ "${output}" == *"ERROR: Cannot find"* ]]
}


