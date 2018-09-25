#!/bin/sh

source ../../init.sh

if [[ -s ${REL_SAF_FILTER_FILE_NAME} ]]; then
  exit 1
else
  exit 0
fi
