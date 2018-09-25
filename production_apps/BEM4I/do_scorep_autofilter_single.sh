#!/bin/sh

cp scorep.filt old_scorep.filt

scorep-autofilter -t $1 -f scorep scorep-*/profile.cubex 1>/dev/null
echo "#############"
cat scorep.filt
echo "#############"

rm -rf scorep-*
commString=$(comm -1 -2 scorep.filt old_scorep.filt)
if [ "$commString" == "" ]; then
  cp scorep.filt old_scorep.filt
  exit 1
else
  filtFileBeginString="SCOREP_REGION_NAMES_BEGIN"
  filtFileEndString="SCOREP_REGION_NAMES_END"
  filtFileExcludeString="EXCLUDE"        
  sort scorep.filt > sorted_scorep.filt
  sort old_scorep.filt > sorted_old_scorep.filt
  commString=$(comm -2 -3 sorted_scorep.filt sorted_old_scorep.filt)
  if [ "$commString" != "" ]; then
    echo $filtFileBeginString > new_scorep.filt
    echo $filtFileExcludeString >> new_scorep.filt
    comm -2 -3 sorted_scorep.filt sorted_old_scorep.filt > comm1.txt
    exec < "old_scorep.filt"
    while read line; do
      if [ "$line" != "$filtFileEndString" ] && [ "$line" != "$filtFileBeginString" ] && [ "$line" != "$filtFileExcludeString" ]; then
        echo $line >> new_scorep.filt
      fi
    done
    cat comm1.txt >> new_scorep.filt
    echo $filtFileEndString >> new_scorep.filt
    mv new_scorep.filt scorep.filt
    rm -f sorted_scorep.filt sorted_old_scorep.filt comm1.txt
    echo "++++++++++++"
    echo $commString
    echo "++++++++++++"
    exit 1
  else
    mv old_scorep.filt scorep.filt
    rm -f sorted_scorep.filt sorted_old_scorep.filt comm1.txt
    echo "++++++++++++"
    echo $commString
    echo "++++++++++++"
    exit 0
  fi
fi
