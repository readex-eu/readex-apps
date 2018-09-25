module purge
module load Anaconda3/4.4.0

/scratch/work/user/bes0030/readex-radar/printFullReport.py -configFile ./config.py
cp "../amg2013_meric_dir/COUNTERS - RAPL:-AVG Energy summary [J]-"*.opts ../energy.opts
 
