## Calibration instructions

### List input files

```
filelisting.py PDSPProd4a_MC_1GeV_michelremoving_sce_datadriven_merged_v1
```

- The only argument is a sam dataset definition
- It will creates a text file `input_run####.txt`, where `####` is the run number (0 for MC sample), with xroot urls to all files in the sam definition.

### YZ calibration

```
make_yz_correction input_run0.txt 2
```

- The first argument is the input file list created in the first step.
- The second argument is the ID of the TTree. The number `2` takes the TTree `michelremoving2/Event`, which is after SCE and lifetime calibration.
- This will creates a root file `YZcalo_mich2_r####.root`, that includes YZ correction histograms.

### X calibration

```
make_x_correction input_run0.txt 2 1
```

- The first argument is the input file list created in the first step.
- The second argument is the ID of the TTree. The number `2` takes the TTree `michelremoving2/Event`, which is after SCE and lifetime calibration.
- The third argument is the flag for SCE. 1 means SCE on. 0 means SCE off.
- This should be run after YZ calibration and it takes the YZ calibration root file as an input.
- This will creates two root files `Xcalo_mich2_r####.root` and `globalmedians_cathanode_r####.root` and three txt files `global_median_{0,1,2}_r####.txt`.

### Calculate calibration constants

```
dEdX_calibration -i ../input_run0.txt -c protoDUNE_dEdx_calib.fcl -o output.root
```

- `-i` specifies the input text file.
- `-c` specifies a fcl file for configuration
- `-o` specifies the output root file name.

Here is an example fcl file:

```
FieldMap: "SCE_DataDriven_180kV_v4.root"

OutputFile: "dedx_cal_out.root"

YZCaloFile: "../YZcalo_mich2_r0.root"
XCaloFile:  "../Xcalo_mich2_r0.root"

NormFactors: [
  0,
  0,
  0
]

KE_Range: [
  [10,  0.70437],
  [14,  1.27937],
  [20,  2.37894],
  [30,  4.72636],
  [40,  7.5788],
  [80,  22.0917],
  [100,  30.4441],
  [140,  48.2235],
  [200,  76.1461],
  [300,  123.567],
  [400,  170.845],
  [800,  353.438],
  [1000,  441.476]
]

Plane0Start: 1.035e-3
Plane0End: 1.045e-3
Plane0Diff: .0001e-3

Plane1Start: 1.035e-3
Plane1End: 1.045e-3
Plane1Diff: .0001e-3

Plane2Start: 1.011e-3
Plane2End: 1.021e-3
Plane2Diff: .0001e-3

Method: Lite
```

```
calculate_calibration_constants output.root
```

### Validate calibration

```
validate_calibration -c data.fcl
```
An example fcl file:
```
infile: "./input_run5387.txt"
outfile: "Validate_mich2_r5387.root"
michelnumber: 2
recalib: true
sceon: true
corr_end: false
xbins: 144
xmin: -360
xmax: 360
rroffset: [0., 0., 0.]
```

### Landau Gaussian fit

```
fitlangaumpv Validate_mich2_r5387.root
```

### Convert constants to csv files

```
make_norm_csv.py
make_x_csv.py
make_yz_csv.py
```

### Upload csv files to database

```
setup dunetpc v09_00_01 -q e19:prof:py2
setup condb 2.3
# The following file needs to be copied from protoduneana repository
write_to_db_pdsp.py
```
