#!/bin/sh

test_out_file_path=$1
test_name=$2
run_count=$3
nodes=$4
ntasks=$5
tasks_per_node=$6
cpus_per_task=$7
plain_app_command=$8
rrl_app_command=$9
tuning_model_file=${10}
app_args=${*:11}

run_file=${test_name}_r.sh
plain_test_out_file=${test_out_file_path}/${test_name}_plain_hdeem.out
rrl_test_out_file=${test_out_file_path}/${test_name}_rrl_hdeem.out
hdeem_out_file=${test_name}_temp_hdeem.out
host_names_file=${test_name}_host_names.out

###########################################
# generate run script
###########################################

printf \
"#!/bin/sh\n\
\n\
#SBATCH --time=4:00:00\n\
#SBATCH --nodes=$nodes\n\
#SBATCH --ntasks=$ntasks\n\
#SBATCH --ntasks-per-node=$tasks_per_node\n\
#SBATCH --cpus-per-task=$cpus_per_task\n\
#SBATCH --exclusive\n\
#SBATCH --partition=haswell\n\
#SBATCH --comment=\"cpufreqchown\"\n\
#SBATCH --mem-per-cpu=2500M\n\
#SBATCH -J \"${test_name}_r\"\n\
#SBATCH -A p_readex\n\
#SBATCH --output=${test_name}.out\n\
#SBATCH --error=${test_name}.out\n\
\n\
\n\
energy_label=\"Energy\"\n\
srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 hostname >> ${host_names_file}\n\
\n\
module purge\n\
source ./readex_env/set_env_plain.source\n\
source ./set_env_blasbench.source\n\
\n\
# start plain run
export SCOREP_ENABLE_PROFILING=\"false\"\n\
export SCOREP_ENABLE_TRACING=\"false\"\n\
export SCOREP_SUBSTRATE_PLUGINS=\"\"\n\
export SCOREP_RRL_PLUGINS=\"\"\n\
export SCOREP_RRL_TMM_PATH=\"\"\n\
export SCOREP_MPI_ENABLE_GROUPS=ENV\n\
\n\
for run_id in \$(seq 1 $run_count); do\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 clearHdeem\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 startHdeem\n\
  # run application\n\
  start_time=\$((\$(date +%%s%%N)/1000000))\n\
  srun --cpu_bind=verbose,sockets $plain_app_command $app_args\n\
  stop_time=\$((\$(date +%%s%%N)/1000000))\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 stopHdeem\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 sleep 5\n\
  exec < ${host_names_file}\n\
  while read host_name; do\n\
    srun -N 1 -n 1 --ntasks-per-node=1 -c 1 --nodelist=\$host_name checkHdeem >> $hdeem_out_file\n\
  done\n\
  \n\
  energy_total=0\n\
  if [ -e $hdeem_out_file ]; then\n\
    exec < $hdeem_out_file\n\
    while read max max_unit min min_unit average average_unit energy energy_unit; do\n\
      if [ \"\$energy\" == \"\$energy_label\" ]; then\n\
        read blade max_val min_val average_val energy_val\n\
        energy_total=\$(echo \"\${energy_total} + \${energy_val}\" | bc)\n\
      fi\n\
    done\n\
    time_total=\$(echo \"\${stop_time} - \${start_time}\" | bc)\n\
    echo \"$test_name \$run_id $nodes $ntasks $tasks_per_node $cpus_per_task \$time_total \$energy_total\" >> $plain_test_out_file\n\
    rm -rf $hdeem_out_file\n\
  fi\n\
  # end plain run\n\
done\n\
\n\
\n\
\n\
# start rrl-tuned run\n\
module purge\n\
source ./readex_env/set_env_rrl.source\n\
source ./set_env_blasbench.source\n\
export SCOREP_ENABLE_PROFILING=\"false\"\n\
export SCOREP_ENABLE_TRACING=\"false\"\n\
export SCOREP_SUBSTRATE_PLUGINS=\"rrl\"\n\
export SCOREP_RRL_PLUGINS=\"cpu_freq_plugin,uncore_freq_plugin,OpenMPTP\"\n\
export SCOREP_RRL_TMM_PATH=\"${tuning_model_file}\"\n\
export SCOREP_MPI_ENABLE_GROUPS=ENV\n\
export SCOREP_RRL_CHECK_IF_RESET=\"reset\"\n\
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib\n\
\n\
for run_id in \$(seq 1 $run_count); do\n\
  date_name=\$(date +%%G%%m%%d%%H%%M%%S%%N)\n\
  export SCOREP_EXPERIMENT_DIRECTORY=\"scorep-${test_name}_${nodes}_${ntasks}_${tasks_per_node}_${cpus_per_task}_\${date_name}\"\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 clearHdeem\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 startHdeem\n\
  # run application\n\
  start_time=\$((\$(date +%%s%%N)/1000000))\n\
  srun --cpu_bind=verbose,sockets $rrl_app_command $app_args\n\
  stop_time=\$((\$(date +%%s%%N)/1000000))\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 stopHdeem\n\
  srun -N $nodes -n $nodes --ntasks-per-node=1 -c 1 sleep 5\n\
  exec < ${host_names_file}\n\
  while read host_name; do\n\
    srun -N 1 -n 1 --ntasks-per-node=1 -c 1 --nodelist=\$host_name checkHdeem >> $hdeem_out_file\n\
  done\n\
  \n\
  energy_total=0\n\
  if [ -e $hdeem_out_file ]; then\n\
    exec < $hdeem_out_file\n\
    while read max max_unit min min_unit average average_unit energy energy_unit; do\n\
      if [ \"\$energy\" == \"\$energy_label\" ]; then\n\
        read blade max_val min_val average_val energy_val\n\
        energy_total=\$(echo \"\${energy_total} + \${energy_val}\" | bc)\n\
      fi\n\
    done\n\
    time_total=\$(echo \"\${stop_time} - \${start_time}\" | bc)\n\
    echo \"$test_name \$run_id $nodes $ntasks $tasks_per_node $cpus_per_task \$time_total \$energy_total\" >> $rrl_test_out_file\n\
    rm -rf $hdeem_out_file\n\
  fi\n\
  # end rrl-tuned run\n\
done\n\
\n\
rm -rf ${host_names_file}
\n\
" > $run_file
