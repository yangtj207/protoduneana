#!/bin/bash

input=$1
while IFS= read -r line
do
  nevts=$( filelisting.py protodune-sp_runset_${line}_michelremoving_merged_v09_09_01_v0 )
  echo "Run $line Total entries: $nevts"
done < "$input"
