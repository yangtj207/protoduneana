# How to use it

## Create a sam definition
```bash
samweb create-definition ${USER}_michelremoving_10418 'run_type protodune-sp and file_type detector and run_number 10418 and data_tier root-tuple and file_name %michelremoving%.root'
```

## Create a file list
```bash
filelisting.py ${USER}_michelremoving_10418
```
It will creates a text file input_run####.txt, where #### is the run number

## Create a fcl file
There is an example run10418.fcl

## Run it
```bash
plotdeltaz -c run10418.fcl
```