#!/bin/bash

mkdir log

input=$1
while IFS= read -r line
do
  echo "Run $line"
  nevts=$( filelisting.py protodune-sp_runset_${line}_michelremoving_merged_v09_09_01_v0 )
  echo "Total entries: $nevts"
  if [ $nevts -ge 20000 ]; then
      make_yz_correction input_run${line}.txt 2 >& log/yz_correction_run${line}.out
  fi
  if [ $nevts -ge 1000 ]; then
      if [ ! -f "YZcalo_mich2_r${line}.root" ]; then
          YZcalo=`ls YZcalo_mich2_r*.root | tail -n 1`
          echo "Using $YZcalo for YZ correction"
          ln -s $YZcalo YZcalo_mich2_r${line}.root
      fi
      make_x_correction input_run${line}.txt 2 >& log/x_correction_run${line}.out
  fi
done < "$input"
