module purge
module load python
#module load python/3.6-anaconda4.3.1
#module load python/3.5.2-anaconda4.2.0
#module load Python/3.6.4-intel-2018a
/projects/p_readex/it4i/readex-radar-dev/printFullReport.py -configFile ./config.py
#cp "../amg2013_meric_dir/Blade summary-Average energy consumption [J]-"taurusi*.opts ../energy.opts
#cp "./lulesh_meric_dir/Blade summary-Average energy consumption [J]-"taurusi*.opts ./energy.opts

#cp "./lulesh_meric_dir/Blade summary-Average energy consumption [J]-"taurusi*.opts ./energy.opts
cp "../lulesh_meric_dir/COUNTERS - RAPL:-AVG Energy summary [J]-"taurusi*.opts ../energy.opts
