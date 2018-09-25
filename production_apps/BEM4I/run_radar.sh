module purge
module load python

/projects/p_readex/it4i/readex-radar-dev/printFullReport.py -configFile ./config.py
cp "./bem4i_meric_dir/Blade summary-Energy consumption [J] (Samples)-"taurusi*.opts energy.opts
