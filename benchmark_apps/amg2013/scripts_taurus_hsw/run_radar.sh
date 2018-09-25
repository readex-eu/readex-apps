module purge
module load python

/projects/p_readex/it4i/readex-radar-dev/printFullReport.py -configFile ./config.py
cp "../amg2013_meric_dir/Blade summary-Average energy consumption [J]-"taurusi*.opts ../energy.opts
 
