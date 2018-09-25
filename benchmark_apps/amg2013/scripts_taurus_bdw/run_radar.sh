module --force purge
module load modenv/classic
module load python

/projects/p_readex/it4i/readex-radar-dev/printFullReport.py -configFile ./config.py
cp "../amg2013_meric_dir/COUNTERS - RAPL:-AVG Energy summary [J]-"*.opts ../energy.opts
 
