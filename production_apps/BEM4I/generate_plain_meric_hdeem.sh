#!/bin/sh

# in create_plain_meric_hdeem.sh and generate_plain_meric_hdeem.sh, the term "test" refers to a production run of the application (with or without RRL)

tests_input_file=$1 # file with test configuration details;
                    # example file: meric_tests.txt
                    # each line: test_name nodes ntasks tasks_per_node cpus_per_task plain_app_executable meric_app_executable tuning_model_file app_args
                    # one batch script for each test lines, ${test_name}_r.sh, is created in the current directory by create_plain_meric_hdeem.sh

run_count=$2 # number of runs of each test in $test_input_file
             # this is used in ${test_name}_r.sh to repeat the untuned (plain) and RRL-tuned runs of the application multiple times

test_out_file_path=. # path to create output files from tests;
                     # two output files are created for each test in tests_input_file: ${test_name}_plain_hdeem.out and ${test_name}_meric_hdeem.out
                     # each output file is a multi-column data files for gnu plot;
                     # each row corresponds to a test run: test_name run_id nodes ntasks tasks_per_node cpus_per_task total_time total_energy

if [ -e $tests_input_file ]; then
  exec < $tests_input_file
  while read test_name nodes ntasks tasks_per_node cpus_per_task plain_app_command meric_app_command tuning_model_file app_args; do
    ./create_plain_meric_hdeem.sh $test_out_file_path $test_name $run_count $nodes $ntasks $tasks_per_node $cpus_per_task $plain_app_command $meric_app_command $tuning_model_file $app_args
#    sbatch ${test_name}_r.sh
  done
else
  echo "invalid tests input file"
fi
