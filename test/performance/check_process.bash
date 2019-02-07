
function check_process(){
    if [ $? -ne 0 ]; then
       ERROR_MESSAGE="$1 Failed for $2"
       echo $ERROR_MESSAGE
       echo $ERROR_MESSAGE >> /local/error.log
       SUBJECT_STRING="--subject \"ERROR : From "${INSTANCE_TYPE}" ID: "${INSTANCE_ID}" running "${include}" in "${CLOUD}"\" --message \"file://error.log\""
       echo "aws sns publish  ${REGION_STRING}   ${TOPIC}   ${SUBJECT_STRING}" > ${WORK_DIR}/sender.sh
       source ${WORK_DIR}/sender.sh
       return 1
    fi
}
