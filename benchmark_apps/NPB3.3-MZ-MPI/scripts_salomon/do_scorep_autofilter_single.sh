#!/bin/sh

cp scorep.filt old_scorep.filt
if [ "$READEX_INTEL" == "1" ]; then
  cp scorep_icc.filt old_scorep_icc.filt
fi

if [ "$READEX_INTEL" == "1" ]; then
  scorep-autofilter -t $1 -f scorep -i scorep_icc ${SCOREP_EXPERIMENT_DIRECTORY}/profile.cubex 1>/dev/null
else
  scorep-autofilter -t $1 -f scorep ${SCOREP_EXPERIMENT_DIRECTORY}/profile.cubex 1>/dev/null
fi

rm -rf ${SCOREP_EXPERIMENT_DIRECTORY}-*

commString=$(comm -1 -2 scorep.filt old_scorep.filt)
if [ "$commString" == "" ]; then
  cp scorep.filt old_scorep.filt
  if [ "$READEX_INTEL" == "1" ]; then
    cp scorep_icc.filt old_scorep_icc.filt
  fi
  exit 1
else
  filtFileBeginString="SCOREP_REGION_NAMES_BEGIN"
  filtFileEndString="SCOREP_REGION_NAMES_END"
  filtFileExcludeString="EXCLUDE"        
  sort scorep.filt > sorted_scorep.filt
  sort old_scorep.filt > sorted_old_scorep.filt
  if [ "$READEX_INTEL" == "1" ]; then
    sort scorep_icc.filt > sorted_scorep_icc.filt
    sort old_scorep_icc.filt > sorted_old_scorep_icc.filt
  fi
  commString=$(comm -2 -3 sorted_scorep.filt sorted_old_scorep.filt)
  if [ "$commString" != "" ]; then
    echo $filtFileBeginString > new_scorep.filt
    echo $filtFileExcludeString >> new_scorep.filt
    comm -2 -3 sorted_scorep.filt sorted_old_scorep.filt > comm1.txt
    exec < "old_scorep.filt"
    while read -r line; do
      if [ "$line" != "$filtFileEndString" ] && [ "$line" != "$filtFileBeginString" ] && [ "$line" != "$filtFileExcludeString" ]; then
        echo $line >> new_scorep.filt
      fi
    done
    cat comm1.txt >> new_scorep.filt
    echo $filtFileEndString >> new_scorep.filt
    mv new_scorep.filt scorep.filt
    rm -f sorted_scorep.filt sorted_old_scorep.filt comm1.txt
    if [ "$READEX_INTEL" == "1" ]; then
      comm -2 -3 sorted_scorep_icc.filt sorted_old_scorep_icc.filt > comm1.txt
      exec < "old_scorep_icc.filt"
      while read -r line; do
        echo $line >> new_scorep_icc.filt
      done
      cat comm1.txt >> new_scorep_icc.filt
      mv new_scorep_icc.filt scorep_icc.filt
      rm -f sorted_scorep_icc.filt sorted_old_scorep_icc.filt comm1.txt
    fi
    exit 1
  else
    mv old_scorep.filt scorep.filt
    rm -f sorted_scorep.filt sorted_old_scorep.filt comm1.txt
    if [ "$READEX_INTEL" == "1" ]; then
      mv old_scorep_icc.filt scorep_icc.filt
      rm -f sorted_scorep_icc.filt sorted_old_scorep_icc.filt comm1.txt
    fi
    exit 0
  fi
fi
