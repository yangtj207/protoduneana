import ROOT as RT
from array import array
from argparse import ArgumentParser as ap
from statistics import median
from math import sqrt
parser = ap()

parser.add_argument('-o', type=str, required=True)
parser.add_argument('-i', type=str, required=True)
parser.add_argument('-e', type=str, required=True)
parser.add_argument('-n', type=int, default=100)
parser.add_argument('--nskip', type=int, default=0)
parser.add_argument('--ymax', type=float, default=600.)
parser.add_argument('--zmax', type=float, default=695.)
parser.add_argument('--ywidth', type=float, default=5.)
parser.add_argument('--zwidth', type=float, default=5.)
args = parser.parse_args()

nbins_y = int(args.ymax/args.ywidth)
nbins_z = int(args.zmax/args.zwidth)

'''
dqdx_neg = []
dqdx_pos = []
all_dqdx_pos = []
all_dqdx_neg = []
for i in range(3):
  dqdx_neg.append([])
  dqdx_pos.append([])
  all_dqdx_neg.append([])
  all_dqdx_pos.append([])
  for j in range(nbins_z):
    dqdx_neg[-1].append([])
    dqdx_pos[-1].append([])
    for k in range(nbins_y):
      dqdx_neg[-1][-1].append(RT.vector('double')())
      dqdx_pos[-1][-1].append(RT.vector('double')())
'''


efield_file = RT.TFile(args.e, 'open')
pos_hists = [efield_file.Get('Reco_ElecField_%s_Pos'%i) for i in ['X', 'Y', 'Z']]
neg_hists = [efield_file.Get('Reco_ElecField_%s_Neg'%i) for i in ['X', 'Y', 'Z']]
E0 = 0.4867
def get_efield(x, y, z):
  if x > 0: hists = pos_hists
  else: hists = neg_hists

  ex = E0 + E0*hists[0].Interpolate(x, y, z)
  ey = E0*hists[1].Interpolate(x, y, z)
  ez = E0*hists[2].Interpolate(x, y, z)
  return sqrt(ex**2 + ey**2 + ez**2)

def filter_event_for_dedx(e, i):
  if e.trkstartx[i]*e.trkendx[i] > 0: return False 
  if (e.peakT_min[i] < 100 or e.peakT_max[i] > 5900 or e.trklen[i] < 100 or
      e.trklen[i] > 700 or (e.trkendz[i] > 226 and e.trkendz[i] < 236) or
      (e.trkstartz[i] > 226 and e.trkstartz[i] < 236) or
      (e.trkendz[i] > 456 and e.trkendz[i] < 472) or
      (e.trkstartz[i] > 456 and e.trkstartz[i] < 472)): return False #filter for plane 2
  if e.adjacent_hits[i] != 0 or e.dist_min[i] > 5: return False;
  return True

def filter_plane_for_dedx(e, i, hit_plane, track_index, hit_index, res, dq):
  if (hit_plane == 2 and (e.lastwire[i] <= 5 or e.lastwire[i] >= 475)):
    return False 
  if (hit_plane == 3 and (e.lastwire[i] <= 5 or e.lastwire[i] >= 795)):
    return False 

  ##Skipping flipped tracks because they have no effect anyways
  nhits = e.ntrkhits[i*3 + hit_plane]
  if (nhits < 0): return False
  #print(track_index, hit_index, nhits, len(e.trkhity))
  if ((e.trkhity[track_index + hit_index + 1] < e.trkhity[track_index + hit_index] and
       e.trkresrange[track_index + hit_index + 1] > e.trkresrange[track_index + hit_index]) or
      (e.trkhity[track_index + hit_index] < e.trkhity[track_index + hit_index + nhits - 1] and
       e.trkresrange[track_index + hit_index ] > e.trkresrange[track_index + hit_index + nhits - 1])):
    return False 

  if len(res) == 0: return False
  max_res = max(res)
  first_5_dq = [q for r, q in zip(res, dq) if r < 5]
  last_5_dq = [q for r, q in zip(res, dq) if r > max_res - 5]
  if len(first_5_dq) < 5: return False

  med1 = RT.TMath.Median(len(first_5_dq), array('d', first_5_dq))
  med2 = RT.TMath.Median(len(last_5_dq), array('d', last_5_dq))
  if med1/med2 < 1.4: return False
  return True

def filter_event_for_yz_corr(e, i):
  ###Check if the start or endpoint are in the FV. Both must be outside to enter sample
  if ((abs(e.trkstartx[i]) < 350 and e.trkstarty[i] > 40 and
       e.trkstarty[i] < 560 and e.trkstartz[i] > 40 and e.trkstartz[i] < 655) or
      (abs(e.trkendx[i]) < 350 and e.trkendy[i] > 40 and e.trkendy[i] < 560 and
       e.trkendz[i] > 40 and e.trkendz[i] < 655)):
    return False
  return True

with open(args.i, 'r') as f:
  input_files = [l.strip() for l in f.readlines()]
print(input_files)

tree = RT.TChain('michelremoving2/Event')
for i in input_files:
  tree.AddFile(i)
print(tree.GetEntries())

spline_range = [0.70437, 1.27937, 2.37894, 4.72636, 7.5788, 22.0917, 30.4441, 48.2235, 76.1461, 123.567, 170.845, 353.438, 441.476]
spline_ke = [10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000]
ke_spline = RT.TSpline3("Cubic Spline", array('d', spline_range), array('d', spline_ke), 13, "b2e2", 0, 0)

fOut = RT.TFile.Open(args.o, 'recreate')
out_tree = RT.TTree('tree', '')
corrected_dq_dx = array('d', [0])
dq_dx = array('d', [0])
hit_plane_out = array('i', [0])
hit_x = array('d', [0])
hit_y = array('d', [0])
hit_z = array('d', [0])
efield_out = array('d', [0.])
ke = array('d', [0.])
resrange_out = array('d', [0.])
range_bin = array('i', [0])

dedx_filter = array('i', [0])
yz_corr_filter = array('i', [0])
x_corr_filter = array('i', [0])
Cx = array('i', [0])
Cyz = array('i', [0])

out_tree.Branch('corrected_dq_dx', corrected_dq_dx, 'corrected_dq_dx/D')
out_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
out_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
out_tree.Branch('hit_x', hit_x, 'hit_x/D')
out_tree.Branch('hit_y', hit_y, 'hit_y/D')
out_tree.Branch('hit_z', hit_z, 'hit_z/D')
out_tree.Branch('efield', efield_out, 'efield/D')
out_tree.Branch('resrange', resrange_out, 'resrange/D')
out_tree.Branch('range_bin', range_bin, 'range_bin/I')
out_tree.Branch('ke', ke, 'ke/D')
out_tree.Branch('dedx_filter', dedx_filter, 'dedx_filter/I')
out_tree.Branch('yz_corr_filter', yz_corr_filter, 'yz_corr_filter/I')
out_tree.Branch('x_corr_filter', x_corr_filter, 'x_corr_filter/I')
out_tree.Branch('Cx', Cx, 'Cx/I')
out_tree.Branch('Cyz', Cyz, 'Cyz/I')

#yz_tree = RT.TTree('yz_tree', '')
yz_trees = [RT.TTree('yz_tree_%i'%i, '') for i in range(3)]
for yz_tree in yz_trees:
  yz_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
  yz_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
  yz_tree.Branch('hit_x', hit_x, 'hit_x/D')
  yz_tree.Branch('hit_y', hit_y, 'hit_y/D')
  yz_tree.Branch('hit_z', hit_z, 'hit_z/D')
  yz_tree.Branch('efield', efield_out, 'efield/D')
  yz_tree.Branch('resrange', resrange_out, 'resrange/D')
  yz_tree.Branch('range_bin', range_bin, 'range_bin/I')
  yz_tree.Branch('ke', ke, 'ke/D')
#x_tree = RT.TTree('x_tree', '')
x_trees = [RT.TTree('x_tree_%i'%i, '') for i in range(3)]
for x_tree in x_trees:
  x_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
  x_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
  x_tree.Branch('hit_x', hit_x, 'hit_x/D')
  x_tree.Branch('hit_y', hit_y, 'hit_y/D')
  x_tree.Branch('hit_z', hit_z, 'hit_z/D')
  x_tree.Branch('efield', efield_out, 'efield/D')
  x_tree.Branch('resrange', resrange_out, 'resrange/D')
  x_tree.Branch('range_bin', range_bin, 'range_bin/I')
  x_tree.Branch('ke', ke, 'ke/D')

dedx_tree = RT.TTree('dedx_tree', '')
dedx_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
dedx_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
dedx_tree.Branch('hit_x', hit_x, 'hit_x/D')
dedx_tree.Branch('hit_y', hit_y, 'hit_y/D')
dedx_tree.Branch('hit_z', hit_z, 'hit_z/D')
dedx_tree.Branch('efield', efield_out, 'efield/D')
dedx_tree.Branch('resrange', resrange_out, 'resrange/D')
dedx_tree.Branch('range_bin', range_bin, 'range_bin/I')
dedx_tree.Branch('ke', ke, 'ke/D')

nevent = 0
for e in tree:
  if nevent < args.nskip:
    nevent += 1
    print('skipping %i'%nevent)
    continue
  if args.n > 0 and nevent > args.n: break
  print(nevent, end='\r')
  nevent += 1
  for i in range(0, e.cross_trks):
    track_index = i*9000
    dedx_event_filter = filter_event_for_dedx(e, i)
    yz_corr_filter[0] = filter_event_for_yz_corr(e, i)
    x_cross_filter = e.trkstartx[i]*e.trkendx[i] < 0

    if (not dedx_event_filter) and (not yz_corr_filter[0]): continue

    testneg = (e.trkstartx[i]<-350 or e.trkendx[i]<-350)
    testpos = (e.trkstartx[i]>350 or e.trkendx[i]>350)

    theta_xz_deg = 180./RT.TMath.Pi()*e.trackthetaxz[i]
    theta_yz_deg = 180./RT.TMath.Pi()*e.trackthetayz[i]
    for hit_plane in range(0, 3):
      hit_index = hit_plane*3000
      hit_plane_out[0] = hit_plane
      if (hit_plane == 2 and
          ((abs(theta_xz_deg)>60 and abs(theta_xz_deg)<120) or
           (abs(theta_yz_deg)>80 and abs(theta_yz_deg)<100))): continue

      res = []
      dq = []
      nhits = e.ntrkhits[i*3 + hit_plane]
      for k in range(track_index + hit_index, track_index + hit_index + nhits):
        res.append(e.trkresrange[k])
        dq.append(e.trkdqdx[k])
      dedx_plane_filter = filter_plane_for_dedx(e, i, hit_plane, track_index, hit_index, res, dq)
      hit_limit = min(nhits, 3000)
      for j in range(0, hit_limit):
        hit_x[0] = e.trkhitx[track_index + hit_index + j]
        hit_y[0] = e.trkhity[track_index + hit_index + j]
        hit_z[0] = e.trkhitz[track_index + hit_index + j]

        if (hit_y[0] < 0 or hit_y[0] > 600 or hit_z[0] < 0 or hit_z[0] > 695):
          continue

        ##negative X
        if (hit_x[0] > -360 and hit_x[0] < 0.):
          hit_pos = False
          if hit_plane == 1 and abs(theta_xz_deg) < 140: continue
          elif hit_plane == 0 and abs(theta_xz_deg) > 40: continue
        elif (hit_x[0] > 0. and hit_x[0] < 360.):
          hit_pos = True
          if hit_plane == 1 and abs(theta_xz_deg) > 40: continue
          elif hit_plane == 0 and abs(theta_xz_deg) < 140: continue
        else: continue

        ke[0] = ke_spline.Eval(e.trkresrange[track_index + hit_index + j]);
        resrange_out[0] = e.trkresrange[track_index + hit_index + j]
        range_bin[0] = int(resrange_out[0]/5.) 
        hit_trkpitch = e.trkpitch[track_index + hit_index + j]

        ke_pitch_filter = (ke[0] > 250. and ke[0] < 450. and
                           hit_trkpitch > 0.5 and hit_trkpitch < 0.8)

        dedx_filter[0] = dedx_event_filter and dedx_plane_filter and ke_pitch_filter

        if (not dedx_filter[0]) and (not yz_corr_filter[0]): continue

        dq_dx[0] = e.trkdqdx[track_index + hit_index + j]
        if dedx_filter[0]: efield_out[0] = get_efield(hit_x[0], hit_y[0], hit_z[0])
        else: efield_out[0] = -999.
        
        if (j < nhits-1) and (j > 0):
          if testpos and hit_pos:
            good_x_hit = (e.trkhitx[track_index + hit_index + j+1] > 0 and e.trkhitx[track_index + hit_index + j-1] > 0)
          elif testneg and not hit_pos:
            good_x_hit = (e.trkhitx[track_index + hit_index + j+1] < 0 and e.trkhitx[track_index + hit_index + j-1] < 0)
          else: good_x_hit = False
        else: good_x_hit = False
        x_corr_filter[0] = good_x_hit and x_cross_filter

        if (not dedx_filter[0]) and (not yz_corr_filter[0]) and (not x_corr_filter[0]): continue

        if dedx_filter[0]: dedx_tree.Fill()
        if yz_corr_filter[0]:
          #yz_tree.Fill()
          yz_trees[hit_plane].Fill()
          y_bin = int(hit_y[0]/args.ywidth)
          z_bin = int(hit_z[0]/args.zwidth)
          #if y_bin < nbins_y and z_bin < nbins_z:
          #  if hit_x[0] < 0.:
          #    dqdx_neg[hit_plane][z_bin][y_bin].push_back(dq_dx[0])
          #    all_dqdx_neg[hit_plane].append(dq_dx[0])
          #  else:
          #    dqdx_pos[hit_plane][z_bin][y_bin].push_back(dq_dx[0])
          #    all_dqdx_pos[hit_plane].append(dq_dx[0])
        #if x_corr_filter[0]: x_tree.Fill()
        if x_corr_filter[0]: x_trees[hit_plane].Fill()
        out_tree.Fill()

fOut.cd()
out_tree.Write()
for t in yz_trees: t.Write()
for t in x_trees: t.Write()
#yz_tree.Write()
#x_tree.Write()
dedx_tree.Write()

'''
for i in range(len(dqdx_pos)):
  for j in range(len(dqdx_pos[i])):
    for k in range(len(dqdx_pos[i][j])):
      fOut.WriteObject(dqdx_pos[i][j][k], 'dqdx_pos_%i_%i_%i'%(i, j, k))
      fOut.WriteObject(dqdx_neg[i][j][k], 'dqdx_neg_%i_%i_%i'%(i, j, k))

median_dqdx_pos = RT.TVectorD(3*nbins_y*nbins_z)
median_dqdx_neg = RT.TVectorD(3*nbins_y*nbins_z)

for i in range(3):
  global_median_dqdx_pos = RT.TVectorD(1)
  global_median_dqdx_neg = RT.TVectorD(1)
  global_median_dqdx_pos[0] = median(all_dqdx_pos[i])
  global_median_dqdx_neg[0] = median(all_dqdx_neg[i])

  global_median_dqdx_pos.Write('global_median_dqdx_pos_%i'%i)
  global_median_dqdx_neg.Write('global_median_dqdx_neg_%i'%i)
  for j in range(len(dqdx_pos[i])):
    for k in range(len(dqdx_pos[i][j])):
      median_dqdx_pos[i*nbins_y*nbins_z + j*nbins_y + k] = median(dqdx_pos[i][j][k]) if len(dqdx_pos) >= 5 else -999.
median_dqdx_pos.Write('median_dqdx_pos')
median_dqdx_neg.Write('median_dqdx_neg')
'''

fOut.Close()
