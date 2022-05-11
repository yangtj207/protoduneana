import ROOT as RT
from array import array
from argparse import ArgumentParser as ap
from math import sqrt
parser = ap()

parser.add_argument('-o', type=str, required=True)
parser.add_argument('-i', type=str, required=True)
parser.add_argument('-x', type=str, required=True)
parser.add_argument('-e', type=str, required=True)
parser.add_argument('-yz', type=str, required=True)
parser.add_argument('--norms', nargs=3, default=[1., 1., 1.])
#parser.add_argument('-e', type=str, required=True)
args = parser.parse_args()


x_corr_file = RT.TFile(args.x, 'open')
yz_corr_file = RT.TFile(args.yz, 'open')
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

x_corr_hists = [x_corr_file.Get('dqdx_X_correction_hist_%i'%i) for i in range(0, 3)]
yz_corr_hists_neg = [yz_corr_file.Get('correction_dqdx_ZvsY_negativeX_hist_%i'%i) for i in range(0, 3)]
yz_corr_hists_pos = [yz_corr_file.Get('correction_dqdx_ZvsY_positiveX_hist_%i'%i) for i in range(0, 3)]

print(x_corr_hists)
print(yz_corr_hists_neg)
print(yz_corr_hists_pos)

with open(args.i, 'r') as f:
  input_files = [l.strip() for l in f.readlines()][0:10]
print(input_files)

tree = RT.TChain('michelremoving2/Event')
for i in input_files:
  tree.AddFile(i)
print(tree.GetEntries())

spline_range = [0.70437, 1.27937, 2.37894, 4.72636, 7.5788, 22.0917, 30.4441, 48.2235, 76.1461, 123.567, 170.845, 353.438, 441.476]
spline_ke = [10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000]
ke_spline = RT.TSpline3("Cubic Spline", array('d', spline_range), array('d', spline_ke), 13, "b2e2", 0, 0)

fOut = RT.TFile(args.o, 'recreate')
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
out_tree.Branch('Cx', Cx, 'Cx/I')
out_tree.Branch('Cyz', Cyz, 'Cyz/I')


nevent = 0
for e in tree:
  print(nevent, end='\r')
  nevent += 1
  for i in range(0, e.cross_trks):
    #print('\tTrack %i'%i)
    track_index = i*9000
    #if e.trkstartx[i]*e.trkendx[i] > 0: continue
    #if (e.peakT_min[i] < 100 or e.peakT_max[i] > 5900 or e.trklen[i] < 100 or
    #    e.trklen[i] > 700 or (e.trkendz[i] > 226 and e.trkendz[i] < 236) or
    #    (e.trkstartz[i] > 226 and e.trkstartz[i] < 236) or
    #    (e.trkendz[i] > 456 and e.trkendz[i] < 472) or
    #    (e.trkstartz[i] > 456 and e.trkstartz[i] < 472)): continue #filter for plane 2
    #if e.adjacent_hits[i] != 0 or e.dist_min[i] > 5: continue;
    if not filter_event_for_dedx(e, i):
      dedx_filter[0] = 0
      continue
    if not filter_event_for_yz_corr(e, i):
      yz_corr_filter[0] = 0

    if (not dedx_filter[0]) and (not yz_corr_filter[0]): continue


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
      if not filter_plane_for_dedx(e, i, hit_plane, track_index, hit_index, res, dq):
        dedx_filter[0] = 0
        continue
      '''
      if (hit_plane == 2 and (e.lastwire[i] <= 5 or e.lastwire[i] >= 475)):
        continue
      if (hit_plane == 3 and (e.lastwire[i] <= 5 or e.lastwire[i] >= 795)):
        continue

      res = []
      dq = []
      nhits = e.ntrkhits[i*3 + hit_plane]
      for k in range(track_index + hit_index, track_index + hit_index + nhits):
        res.append(e.trkresrange[k])
        dq.append(e.trkdqdx[k])

      if len(res) == 0: continue

      ##Skipping flipped tracks because they have no effect anyways
      if ((e.trkhity[track_index + hit_index + 1] < e.trkhity[track_index + hit_index] and
           e.trkresrange[track_index + hit_index + 1] > e.trkresrange[track_index + hit_index]) or
          (e.trkhity[track_index + hit_index] < e.trkhity[track_index + hit_index + nhits - 1] and
           e.trkresrange[track_index + hit_index ] > e.trkresrange[track_index + hit_index + nhits - 1])):
        continue

      max_res = max(res)
      first_5_dq = [q for r, q in zip(res, dq) if r < 5]
      last_5_dq = [q for r, q in zip(res, dq) if r > max_res - 5]
      if len(first_5_dq) < 5: continue

      med1 = RT.TMath.Median(len(first_5_dq), array('d', first_5_dq))
      med2 = RT.TMath.Median(len(last_5_dq), array('d', last_5_dq))
      if med1/med2 < 1.4: continue'''      
      for j in range(0, min(nhits, 3000)):
        hit_x[0] = e.trkhitx[track_index + hit_index + j]
        hit_y[0] = e.trkhity[track_index + hit_index + j]
        hit_z[0] = e.trkhitz[track_index + hit_index + j]

        if (hit_y[0] < 0 or hit_y[0] > 600 or hit_z[0] < 0 or hit_z[0] > 695):
          continue

        ##negative x
        if (hit_x[0] > -360 and hit_x[0] < 0.):
          if hit_plane == 1 and abs(theta_xz_deg) < 140: continue
          elif hit_plane == 0 and abs(theta_xz_deg) > 40: continue
          #yz_bin = yz_corr_hists_neg[hit_plane].FindBin(e.trkhitz[track_index + hit_index + j], e.trkhity[track_index + hit_index + j]) 
          #print(z_bin, int(hit_z[0]/5.)+1)
          #Cyz = yz_corr_hists_neg[hit_plane].GetBinContent(yz_bin)
        elif (hit_x[0] > 0. and hit_x[0] < 360.):
          if hit_plane == 1 and abs(theta_xz_deg) > 40: continue
          elif hit_plane == 0 and abs(theta_xz_deg) < 140: continue
          #yz_bin = yz_corr_hists_pos[hit_plane].FindBin(e.trkhitz[track_index + hit_index + j], e.trkhity[track_index + hit_index + j])
          #print(z_bin, int(hit_z[0]/5.)+1)
          #Cyz = yz_corr_hists_pos[hit_plane].GetBinContent(yz_bin)
        else: continue

        ke[0] = ke_spline.Eval(e.trkresrange[track_index + hit_index + j]);
        resrange_out[0] = e.trkresrange[track_index + hit_index + j]
        range_bin[0] = int(resrange_out[0]/5.) 
        hit_trkpitch = e.trkpitch[track_index + hit_index + j]

        ###de/dx filter
        if ((ke[0] < 250. or ke[0] > 450.) or
            (hit_trkpitch < 0.5 or hit_trkpitch > 0.8)):
          continue
          dedx_filter[0] = 0

        x_bin = x_corr_hists[hit_plane].FindBin(e.trkhitx[track_index + hit_index + j])
        if (x_bin != int((hit_x[0] + 360.)/5.)+1):
          print("ERROR")
          break 
        #Cx = x_corr_hists[hit_plane].GetBinContent(x_bin)
        #corrected_dq_dx[0] = e.trkdqdx[track_index + hit_index + j]*Cx*args.norms[hit_plane]*Cyz
        dq_dx[0] = e.trkdqdx[track_index + hit_index + j]
        efield_out[0] = get_efield(hit_x[0], hit_y[0], hit_z[0])
        out_tree.Fill()
fOut.cd()
out_tree.Write()
fOut.Close()
