#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import samweb_cli
from ROOT import TFile, gDirectory

samweb = samweb_cli.SAMWebClient(experiment='dune')

files = samweb.listFiles("defname: %s " % (sys.argv[1]))

run = -1
totalevts = -1
#print files
for file in files:
    #print (file)
    loc = samweb.locateFile(file)
    if run == -1:
        metadata = samweb.getMetadata(file)
        #print (metadata)
        if metadata['file_type'] != 'mc':
            totalevts = 0
            run = metadata['runs'][0][0]
            #print (run)
        else:
            totalevts = 0
            run = 0
        f = open("input_run%s.txt" % run, "w")
    pnfs = loc[0]['full_path'][8:]
    stream = os.popen("pnfsToXRootD %s/%s" % (pnfs,file))
    xroot = stream.read()
    #print(xroot)
    tfile = TFile.Open(xroot.strip())
    mytree = gDirectory.Get('michelremoving2/Event')
    totalevts += mytree.GetEntries()
    #print(mytree.GetEntries())
    #print (xroot)
    f.write(xroot)
    #print("File name: %s, total events: %s" % (file, totalevts))
#f.close()

print(totalevts)
