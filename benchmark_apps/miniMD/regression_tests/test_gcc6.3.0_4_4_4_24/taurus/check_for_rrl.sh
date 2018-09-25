#!/bin/sh

source ../../init.sh

job_name="reg_rrl"

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

if grep -q -e ERROR -e FATAL "${job_out_file}"; then
  echo "check ${job_out_file} for output from RRL."
  exit 0
else
  exit 1
fi

