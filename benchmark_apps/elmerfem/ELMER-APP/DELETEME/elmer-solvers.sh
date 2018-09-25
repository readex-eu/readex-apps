#!/bin/bash -l


# qsub -q qprod -A SERVICE -l select=5:ncpus=24:mpiprocs=24:ompthreads=1:accelerator=True -l x86_adapt=true -I -X -l walltime=6:00:00
# cd /scratch/work/user/lriha/elmer/bench/readex-appparam-1/ 

ELMERHOME=/scratch/work/user/lriha/elmer/

if [ "$1" = "clean" ]; then
	rm -rf tmp_* *.out TEST.PASSED_* *.log *.qsub STDIN.e* STDIN.o*
	rm -rf winkel.grd
fi

if [ "$1" = "cleanall" ]; then
        rm -rf tmp_* *.out TEST.PASSED_* *.log *.dat.* *.dat *.qsub STDIN.e* STDIN.o* 
fi


if [ "$1" = "backup" ]; then
		
		actualTime=$( date +%y%m%d_%H:%M:%S )
        backup_dir=BACKUP-$actualTime

        mkdir $backup_dir

        cp -PR tmp_* *.out TEST.PASSED_* *.log *.dat.* *.dat *.qsub STDIN.e* STDIN.o* $backup_dir/

        rm -rf tmp_* *.out TEST.PASSED_* *.log *.dat.* *.dat *.qsub STDIN.e* STDIN.o* 
fi


if [ "$1" = "cleanmesh" ]; then
        rm -rf winkel.grd winkel
fi

#if [ "$1" = "wslurm" ]; then
#	watch -n 1 tail -n 15 `ls slurm-* | tail -n 5`
#fi


if [ "$1" = "wres" ]; then
  watch -n 2 'ls winkel_*_*.dat | sort -t _ -k 2 -V | xargs tail -n 40'
fi

if [ "$1" = "wresen" ]; then
	watch "cat LOG-el_es-* | sort"
fi

if [ "$1" = "pres" ]; then
  if [ -z $2 ]; then
    ls winkel_*_*.dat | sort -t _ -k 2 -V | xargs tail -n 40 
  else 
    ls winkel_*_*.dat | sort -t _ -k 2 -V | xargs tail -n 40 | grep " $2 "
  fi 
fi


if [ "$1" = "presf" ]; then
  if [ -z $2 ]; then
    ls winkel_*_*.dat | sort -t _ -k 2 -V | xargs tail -n 40 | sort -k 7 -g | awk '!seen[$3]++' | sort -k 3 -g
  else
    ls winkel_*_*.dat | sort -t _ -k 2 -V | xargs tail -n 40 | grep " $2 " | sort -k 7 -g | awk '!seen[$3]++' | sort -k 3 -g
  fi
fi


if [ "$1" = "gen-weak" ]; then

	source $ELMERHOME/env.sh

	density=0.25
	actualTime=$( date +%y%m%d_%H:%M:%S )
	log_file=LOG-gen-$actualTime.log
	
	echo "ELMER generator log file" | tee $log_file

	for dec in 1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30
	do
		dens=`echo "$density $dec" | awk '{printf "%.10f\n", $1/$2}'`
		#echo $dens
		cp winkel-orig.grd winkel.grd
		sed -i -- "s/refdensity/$dens/g" winkel.grd
		#cat winkel.grd
		echo "----------------------------------------------------------------------------------"
		ElmerGrid 1 2 winkel.grd -partition $dec $dec $dec 2 | tee -a $log_file | grep -e "Saving mesh for" -e "part  elements   nodes" -e "ave" 
	done 
fi 

if [ "$1" = "gen-strong" ]; then

        source $ELMERHOME/env.sh

        density=0.1
	
        actualTime=$( date +%y%m%d_%H:%M:%S )
        log_file=LOG-gen-$actualTime.log

	dens=$density
      	cp winkel-orig.grd winkel.grd
      	sed -i -- "s/refdensity/$dens/g" winkel.grd

        echo "ELMER generator log file" | tee $log_file

        for dec in 4 6 8 10
        do
                echo "----------------------------------------------------------------------------------"
                ElmerGrid 1 2 winkel.grd -partition $dec $dec $dec 2 | tee -a $log_file | grep -e "Saving mesh for" -e "part  elements   nodes" -e "ave"
        done
fi


if [ "$1" = "run" ]; then

  qsub_command_0="#!/bin/bash;"
#  qsub_command_0+="ELMERHOME=$HOME/elmer;"
  qsub_command_0+="source $ELMERHOME/env.sh;"
  
	#qsub_command_0+="syslist=\"51 52\";"
	qsub_command_0+="syslist=\"50 51 52 53 54 55 56 43 14 1 2 3 4 5 6 7 8 9 10 11 12 13 16 17 18 19 20 21 22 23 24 25 26 30 31 32 33 40 41 42\";" 
	#qsub_command_0+="syslist=\"23 24 25 26 30 31 32 33 40 41 42 43 14\";"
	#not working 15 21 22
        
	# block strategies start from 40 + some good strategies from prev
	#qsub_command_0+="syslist=\"42 43 14 25 40 41 20 26 10 1 2 3 4 5 6 7 24 16 21 11\";"
	
	#LR selection
	#qsub_command_0+="syslist=\"43 16 10 21 24 26 33 40 41 \";"
	#qsub_command_0+="syslist=\"42\";"   # 60 61
	#qsub_command_0+="syslist=\"50 51 52 14\";"     #\"50 51 52 53 54 55 56 43 14\";"
	

	# complete list: scalar eqs
	#qsub_command_0+="syslist=\"25 10 1 2 3 4 5 6 7 16 20 24 21 26 14 11 8 9 10 12 13 15 17 18 19 22 23\";" 
	
	# to see the linear system method say: grep IterMethod linsys/*.sif
	# qsub_command_0+="syslist=\"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25\";"

	# modified list
	# qsub_command_0+="syslist=\"25 10 1 2 3 4 5 6 7 16 20 24 21 26 14 11\";"

	# qsub_command_0+="set -x;"
	
	actualTime=$( date +%y%m%d_%H:%M:%S )
	cases=($(ls winkel | grep partitioning | sed "s/partitioning.//g" | sort -n))
	
	for ((i=0; i < ${#cases[*]}; i++))
	# for ((i=0; i < 1; i++))
	do
		
    echo "----------------------------------------------------------------------------------"

	ranks=${cases[$i]}
	nodes=$(( (ranks/22)+1 ))
                
	echo "case# = " $((i))
	echo "ranks = " $((ranks))	
	echo "nodes = " $((nodes))	
    echo 

	log_file=LOG-el_es-$ranks-$actualTime.log
	log_file2=LOG-el_es-$ranks-$actualTime
    err_file=LOG-el_es-$ranks-$actualTime.err
	qsub_command=$qsub_command_0


	qsub_command+="for s in \$syslist;"
  
  
  	qsub_command+="do;"
	qsub_command+="cd $PWD;"

  	qsub_command+="mkdir tmp_${ranks}_\${s};"
  	qsub_command+="cd tmp_${ranks}_\${s};"

  	qsub_command+="cp -r ../*.sif .;"

 	qsub_command+="ln -s ../winkel .;"
 	qsub_command+="cp -r ../ELMERSOLVER_STARTINFO .;"
 	qsub_command+="cp -r ../linsys/*.xml .;"

  	qsub_command+="cp -r ../linsys/linsys\${s}.sif linsys.sif;"
  	qsub_command+="cp -r ../linsys/espreso\${s}.ecf espreso.ecf;"
  
  	qsub_command+="wait;"
		
	# qsub_command+="aprun -t 300 -n $((ranks)) ElmerSolver_mpi | tee ${log_file}_linsys\${s};"
  	qsub_command+="$ELMERHOME/staticMERICtool/multiNodeStaticMeasureStart.sh --rapl;"
  	#qsub_command+="I_MPI_JOB_TIMEOUT=600 mpirun -n $((ranks)) ElmerSolver_mpi > >(tee -a ${log_file}_linsys\${s})  2> >(tee -a ${err_file}_linsys\${s}) ;"  #    | tee ${log_file}_linsys\${s};"
   	qsub_command+="I_MPI_JOB_TIMEOUT=6000 mpirun -n $((ranks)) ElmerSolver_mpi > >(tee -a ${log_file}_linsys\${s})  2> >(tee -a ${err_file}_linsys\${s}) ;"  #    | tee ${log_file}_linsys\${s};"
  	qsub_command+="echo -e linsys'\t'\${s}'\t'ranks'\t'${ranks}'\t'nodes'\t'${nodes} | tr \"\\n\" \"\\t\" | tee -a $PWD/${log_file2}_energy.log;"
  	qsub_command+="$ELMERHOME/staticMERICtool/multiNodeStaticMeasureStop.sh --rapl | tee -a $PWD/${log_file2}_energy.log;"
  	qsub_command+="echo | tee -a $PWD/${log_file2}_energy.log;"
  
	qsub_command+="wait;"
	qsub_command+="cd ..;"
	qsub_command+="done;"

  	echo $qsub_command | tr ";" "\n" | tee run_file_$ranks-$actualTime.qsub

  	if [ "$2" = "mpi" ]; then
    	source run_file_$ranks-$actualTime.qsub
  	fi

  	if [ "$2" = "pbs" ]; then
   	 	echo $qsub_command | tr ";" "\n" | \
    	#qsub -q qmpp -A SERVICE -l select=$(( nodes )):ncpus=24:mpiprocs=22:ompthreads=1:accelerator=True -l x86_adapt=true -l walltime=01:59:00 
  		qsub -q qprod -A SERVICE -l select=$(( nodes )):ncpus=24:mpiprocs=22:ompthreads=1:accelerator=True -l x86_adapt=true -l walltime=47:59:00
  	fi

# qsub -q qprod -A SERVICE -l select=5:ncpus=24:mpiprocs=24:ompthreads=1:accelerator=True -l x86_adapt=true -I -X -l walltime=6:00:00



#		if [ $nodes -le 23 ]; then # 3 23
#			echo $qsub_command | tr ";" "\n" | \
#			sbatch -N $(( nodes ))  -p $queue -t 00:25:00 # --gid 2000190
#			echo "sbatch -N $(( nodes ))  -p $queue -t 00:25:00 " # --gid 2000190"
#		else 
#			echo $qsub_command | tr ";" "\n" | \
#			sbatch -N $(( nodes ))  -p $queue -t 00:25:00 # --gid 2000190
#			echo "sbatch -N $(( nodes ))  -p $queue -t 00:25:00 " # --gid 2000190"
#		fi


	done

fi


