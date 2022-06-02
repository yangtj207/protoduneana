import ROOT as RT
import sys
from statistics import median
from argparse import ArgumentParser as ap
from array import array
import numpy as np

parser = ap()


parser.add_argument('-f', type=str, required=True)
parser.add_argument('-o', type=str, required=True)
parser.add_argument('--ny', type=int, default=120)
parser.add_argument('--nz', type=int, default=139)
parser.add_argument('--ymax', type=float, default=600.)
parser.add_argument('--zmax', type=float, default=695.)
parser.add_argument('--bin_max', type=int, default=10000)
parser.add_argument('-n', type=int, default=-1)
args = parser.parse_args()


ts = [RT.TChain('yz_tree_%i'%i) for i in range(3)]
with open(args.f, 'r') as f:
  files = [l.strip() for l in f.readlines() if '#' not in l]
print(files)
for f in files:
  for t in ts:
    t.AddFile(f)
    print(t.GetEntries())

neg_hists = [RT.TH2D('dqdx_ZvsY_negativeX_hist_%i'%i, '', args.nz, 0, args.zmax, args.ny, 0, args.ymax) for i in range(3)]
pos_hists = [RT.TH2D('dqdx_ZvsY_positiveX_hist_%i'%i, '', args.nz, 0, args.zmax, args.ny, 0, args.ymax) for i in range(3)]
coverage_neg_hists = [RT.TH2D('coverage_ZvsY_negativeX_hist_%i'%i, '', args.nz, 0, args.zmax, args.ny, 0, args.ymax) for i in range(3)]
coverage_pos_hists = [RT.TH2D('coverage_ZvsY_positiveX_hist_%i'%i, '', args.nz, 0, args.zmax, args.ny, 0, args.ymax) for i in range(3)]
correction_neg_hists = [RT.TH2D('correction_dqdx_ZvsY_negativeX_hist_%i'%i, '', args.nz, 0, args.zmax, args.ny, 0, args.ymax) for i in range(3)]
correction_pos_hists = [RT.TH2D('correction_dqdx_ZvsY_positiveX_hist_%i'%i, '', args.nz, 0, args.zmax, args.ny, 0, args.ymax) for i in range(3)]

for h in coverage_pos_hists: h.SetDirectory(0)
for h in coverage_neg_hists: h.SetDirectory(0)

for i in range(3):

  dqdx_pos = []
  dqdx_neg = []
  n_cell_pos = []
  n_cell_neg = []
  
  for j in range(args.nz):
    print(j)
    dqdx_pos.append([])
    dqdx_neg.append([])
    n_cell_pos.append([])
    n_cell_neg.append([])
    for k in range(args.ny):
      dqdx_pos[j].append(array('f', [0.])*args.bin_max)
      dqdx_neg[j].append(array('f', [0.])*args.bin_max)
      n_cell_pos[-1].append(0)
      n_cell_neg[-1].append(0)


  all_dqdx_pos = []
  all_dqdx_neg = []

  t = ts[i]
  a = 0
  nevents = t.GetEntries()
  n_full = 0
  for e in t:
    if not a % 100000: print('%i/%i'%(a, nevents), end='\r')
    a += 1
    if args.n > 0 and a > args.n: break
    if n_full >= 2*args.nz*args.ny: break
    x = e.hit_x
    y = e.hit_y
    z = e.hit_z
    hit_plane = e.hit_plane
    if i != hit_plane: print('ERROR SOMEHOW NOT THE RIGHT PLANE')
    y_bin = int(y/5.)
    z_bin = int(z/5.)
  
    if x < 0.:
      if y_bin >= args.ny or z_bin >= args.nz: continue
      coverage_neg_hists[hit_plane].Fill(z, y)
      if n_cell_neg[z_bin][y_bin] >= args.bin_max: continue
      dqdx_neg[z_bin][y_bin][n_cell_neg[z_bin][y_bin]] = e.dq_dx
      n_cell_neg[z_bin][y_bin] += 1
      if n_cell_neg[z_bin][y_bin] == args.bin_max:
        print('neg', hit_plane, z_bin, y_bin, 'full')
        n_full += 1
      all_dqdx_neg.append(e.dq_dx)
    else:
      if y_bin >= args.ny or z_bin >= args.nz: continue
      #coverage_pos_hists[hit_plane].AddBinContent(z_bin+1, y_bin+1)
      coverage_pos_hists[hit_plane].Fill(z, y)
      if n_cell_pos[z_bin][y_bin] >= args.bin_max: continue
      dqdx_pos[z_bin][y_bin][n_cell_pos[z_bin][y_bin]] = e.dq_dx
      n_cell_pos[z_bin][y_bin] += 1
      if n_cell_pos[z_bin][y_bin] == args.bin_max:
        print('pos', hit_plane, z_bin, y_bin, 'full')
        n_full += 1
      all_dqdx_pos.append(e.dq_dx)

  global_median_pos = median(all_dqdx_pos) if len(all_dqdx_pos) >= 5 else 1.
  global_median_neg = median(all_dqdx_neg) if len(all_dqdx_neg) >= 5 else 1.
  print(global_median_pos, global_median_neg)

  for j in range(args.nz):
    print(i, j, end='\r')
    for k in range(args.ny):
      nonzero_neg = [l for l in dqdx_neg[j][k] if l > 0.]
      nonzero_pos = [l for l in dqdx_pos[j][k] if l > 0.]

      f_neg = global_median_neg/median(nonzero_neg) if (n_cell_neg[j][k] >= 5) else 1.
      f_pos = global_median_pos/median(nonzero_pos) if (n_cell_pos[j][k] >= 5) else 1.

      correction_neg_hists[i].SetBinContent(j+1, k+1, f_neg)
      correction_pos_hists[i].SetBinContent(j+1, k+1, f_pos)
        
      #print(i, j, k, med_neg, med_pos, (global_median_neg[i]/med_neg), (global_median_pos[i]/med_pos))

  #print('Resetting')
  #for j in range(args.nz):
  #  for k in range(args.ny):
  #    for l in range(args.bin_max):
  #      dqdx_pos[j][k][l] = 0.
  #      dqdx_neg[j][k][l] = 0.
  #    n_cell_pos[j][k] = 0
  #    n_cell_neg[j][k] = 0

fOut = RT.TFile(args.o, 'recreate')
for h in correction_pos_hists: h.Write()
for h in correction_neg_hists: h.Write()
for hn, hp in zip(coverage_neg_hists, coverage_pos_hists):
  hn.Write()
  hp.Write()
fOut.Close()
