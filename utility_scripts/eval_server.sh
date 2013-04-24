#!/bin/bash

QUIP_DIR=~/QUIP/build.linux_x86_64_gfortran
export GFORTRAN_UNBUFFERED_ALL=Y

trap quit_eval SIGINT SIGTERM

function init_eval {
   if [[ $1 == "" || $1 == "--help" ]]
   then
      help_init_eval
      return
   else
      PARAM_FILE=$1
   fi

   INIT_ARGS=$2

   if [ ! -f ${PARAM_FILE} ]
   then
      echo "${PARAM_FILE} does not exist"
      echo "Exiting."
      return
   fi

   echo "eval server initialising..."
   EVAL_PORT=`( netstat  -atn | awk '{printf "%s\n%s\n", $4, $4}' | grep -oE '[0-9]*$'; seq 32768 61000 ) | sort -n | uniq -u | head -n 1`
   export EVAL_PORT

   TMP_INPUT=$(mktemp -u /tmp/temp.XXXXX)
   mkfifo $TMP_INPUT

   TMP_OUTPUT=$(mktemp /tmp/temp.XXXXX)

   tail -f ${TMP_INPUT} | grep --line-buffered '' | ${QUIP_DIR}/eval param_file=${PARAM_FILE} init_args={${INIT_ARGS}} e f relax relax_tol=1.0e-8  2>/dev/null | grep --line-buffered '' >${TMP_OUTPUT} &
   GREP1_PID=$!
   EVAL_PID=$(( ${GREP1_PID} - 1 ))
   GREP2_PID=$(( ${GREP1_PID} - 2 ))
   TAIL_PID=$(( ${GREP1_PID} - 3 ))
   EVAL_STARTED=$(nc -l ${EVAL_PORT})
   echo "done"
}


function run_eval {
   if [[ $1 == "" || $1 == "--help" ]]
   then
      help_run_eval
      return
   else
      INPUT_FILE=$1
   fi

   if [[ $2 == "" || $2 == "--help" ]]
   then
      help_run_eval
      return
   else
      OUTPUT_FILE=$2
   fi

   if [[ ! -e ${INPUT_FILE} ]]
   then
      echo "${INPUT_FILE} does not exist"
      echo "Exiting."
      return
   fi

   N_FRAMES=$(grep -c "Lattice" ${INPUT_FILE})

   if [[ ${N_FRAMES} != "1" ]]
   then
      echo "Input file ${INPUT_FILE} contains ${N_FRAMES} frames instead of a single frame."
      echo "Exiting."
      return
   fi
  
   cat ${INPUT_FILE} | grep --line-buffered '' >>${TMP_INPUT}
   N_ATOMS=$(nc -l ${EVAL_PORT})
   grep --line-buffered "^AT" ${TMP_OUTPUT} | tail -$((${N_ATOMS}+2)) | colrm 1 3 >>${OUTPUT_FILE}
}

function help_run_eval {
   echo "Usage: run_eval input.xyz output.xyz"
}

function help_init_eval {
   echo "Usage: init_eval param_file [init_args]"
}

function quit_eval {
   echo "eval server exiting..."
   trap - SIGINT SIGTERM
   kill -INT ${GREP1_PID} ${EVAL_PID} ${GREP2_PID} ${TAIL_PID} 
   wait ${GREP1_PID} ${EVAL_PID} ${GREP2_PID} ${TAIL_PID} 2>/dev/null
   rm -rf ${TMP_INPUT} ${TMP_OUTPUT}
   echo "done"
}
