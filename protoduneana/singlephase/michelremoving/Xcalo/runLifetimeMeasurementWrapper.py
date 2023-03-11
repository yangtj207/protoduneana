import os
from argparse import ArgumentParser as ap
runArray=[5770]
for i in runArray:
    os.system("cache_state.py  -p -d protodune-sp_runset_"+str(i)+"_michelremoving_merged_v09_09_01_v0")

for i in runArray:
    os.system("python filelisting.py -r "+str(i))
    os.system("make_yz_correction file_list_"+str(i)+".txt 1")
    os.system("make_x_correction file_list_"+str(i)+".txt 1 1 1")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",46)'")
