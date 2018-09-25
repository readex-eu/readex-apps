#!/bin/sh

source ../../init.sh

job_name="reg_ptf"

job_ids=$(squeue -u readexci -n $job_name --format="%i")
job_ids_array=(${job_ids[@]})
job_id=${job_ids_array[1]}

job_out_file=slurm-${job_id}.out

# wait for job to finish
job_state_value="x"
while [ "$job_state_value" != "" ]; do
  sleep 60
  job_state=$(squeue -u readexci -n $job_name --format="%t")
  job_state_array=(${job_state[@]})
  job_state_value=${job_state_array[1]}
done

if [[ -s ${REL_TUNING_MODEL_FILE_NAME} ]]; then
  exit 1
else
  echo "check ${job_out_file} for output from PTF."
  exit 0
fi

