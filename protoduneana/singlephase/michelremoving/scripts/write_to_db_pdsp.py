#!/usr/bin/env python

import ROOT as RT
import os, re
from glob import glob as ls
from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-p", type=str, help='Path to files', default='./')
parser.add_argument( "-u", type=str, help='Username', default='tjyang')
args = parser.parse_args()

#print(ls(args.p))

all_files = [i.split("/")[-1] for i in sorted(ls(args.p + "*"))]

for p in all_files:
    if '.csv' in p:
        print(p)
        run = "0"
        data_type = "mc"
        command = "print"
        if 'globalmedians' in p:
            run = p[p.find("run")+3:p.find(".csv")]
            print (run)
            if run != "0":
                data_type = "data"
            command = "python $CONDB_DIR/bin/write_data.py -h ifdbprod2.fnal.gov -p 5451 -U {} -c {} -d {} pdunesp_prod pdunesp.distcorrnorm norm,norm_err".format(args.u, p, data_type)
        else:
            pattern = "\_(.*?)\_"
            run = re.search(pattern,p).group(1)
            print (run)
            if run != "0":
                data_type = "data"
            if 'xcorr' in p:
                command = "python $CONDB_DIR/bin/write_data.py -h ifdbprod2.fnal.gov -p 5451 -U {} -c {} -d {} pdunesp_prod pdunesp.distcorrx x,dx,shape,shape_err".format(args.u, p, data_type)
            elif 'yzcorr' in p:
                command = "python $CONDB_DIR/bin/write_data.py -h ifdbprod2.fnal.gov -p 5451 -U {} -c {} -d {} pdunesp_prod pdunesp.distcorryz y,dy,z,dz,corr,corr_err".format(args.u, p, data_type)
        print (command)
        os.system(command)
