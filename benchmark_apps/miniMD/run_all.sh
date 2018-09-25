#!/bin/sh

DT=$(date +"%d_%m_%G_%H_%M_%S")
mkdir /projects/p_readex/auto_eval/miniMD/results_$DT


sh compile_for_saf.sh
sh run_saf.sh
cp scorep.filt /projects/p_readex/auto_eval/miniMD/results_$DT
cp scorep_icc.filt /projects/p_readex/auto_eval/miniMD/results_$DT


sh compile_for_rdd_manual.sh
sh run_rdd.sh 1>rdd.out
mv rdd.out /projects/p_readex/auto_eval/miniMD/results_$DT


sh do_extend_readex_config.sh
cp ptf_readex_config.xml /projects/p_readex/auto_eval/miniMD/results_$DT

sh compile_for_ptf_manual.sh
sbatch run_ptf.sh
# wait for job to finish
job_name=e_ptf
job_ids=$(squeue -u $(whoami) -n $job_name --format="%i")
job_ids_array=(${job_ids[@]})
job_id=${job_ids_array[1]}
job_out_file=slurm-$job_id.out
job_state_value="x"
while [ "$job_state_value" != "" ]; do
  sleep 60
  job_state=$(squeue -u $(whoami) -n $job_name --format="%t")
  job_state_array=(${job_state[@]})
  job_state_value=${job_state_array[1]}
done
cp slurm-$job_id.out /projects/p_readex/auto_eval/miniMD/results_$DT/ptf.out
cp tuning_model.json /projects/p_readex/auto_eval/miniMD/results_$DT


sbatch run_sacct_rrl_plain.sh
# wait for job to finish
job_name=e_rrl
job_ids=$(squeue -u $(whoami) -n $job_name --format="%i")
job_ids_array=(${job_ids[@]})
job_id=${job_ids_array[1]}
job_out_file=slurm-$job_id.out
job_state_value="x"
while [ "$job_state_value" != "" ]; do
  sleep 60
  job_state=$(squeue -u $(whoami) -n $job_name --format="%t")
  job_state_array=(${job_state[@]})
  job_state_value=${job_state_array[1]}
done
cp slurm-$job_id.out /projects/p_readex/auto_eval/miniMD/results_$DT/rrl_plain.out

rm -rf scorep-* scorep_* slurm-*.out properties_tune*.psc
