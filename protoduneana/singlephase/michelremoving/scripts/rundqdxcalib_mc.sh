#!/bin/bash

nevts=$( filelisting.py tjyang_test_mcmichelremoving )
echo "Total entries: $nevts"
if [ $nevts -ge 20000 ]; then
    make_yz_correction input_run0.txt 2
fi
if [ $nevts -ge 1000 ]; then
    if [ ! -f "YZcalo_mich2_r0.root" ]; then
        YZcalo=`ls YZcalo_mich2_r*.root | tail -n 1`
        echo "Using $YZcalo for YZ correction"
        ln -s $YZcalo YZcalo_mich2_r0.root
    fi
    make_x_correction input_run0.txt 2
fi
