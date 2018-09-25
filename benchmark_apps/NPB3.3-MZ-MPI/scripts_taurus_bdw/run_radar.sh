module --force purge
module load modenv/classic
module load python

cd ..

/projects/p_readex/it4i/readex-radar-dev/printFullReport.py -configFile ./config.py
cp "/bt.C.x_meric_dir/Blade summary-Energy consumption [J]".opts /energy.opts

 
