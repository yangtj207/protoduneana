import ROOT as RT
from argparse import ArgumentParser as ap
parser = ap()

parser.add_argument("-i", type=str, help='Input file', default = "")
parser.add_argument("-p", type=str, help='Plus file', default = "")
parser.add_argument("-m", type=str, help='Minus file', default = "")
parser.add_argument("-o", type=str, help='Output file', default='beam_shift_results.root')
args = parser.parse_args()


names = ['grXSecThrowAbsUnderflow', 'grXSecThrowCexUnderflow',
         'grXSecThrowOtherInelUnderflow']

fN = RT.TFile(args.i, 'open')
gNs = [fN.Get('Throws/%s'%n) for n in names]

fP = RT.TFile(args.p, 'open')
gPs = [fP.Get('Throws/%s'%n) for n in names]

fM = RT.TFile(args.m, 'open')
gMs = [fM.Get('Throws/%s'%n) for n in names]

fOut = RT.TFile(args.o, 'recreate')

leg = RT.TLegend()
for i in range(0, len(gNs)):
  gN = gNs[i]
  gP = gPs[i]
  gM = gMs[i]

  c = RT.TCanvas('c%i'%i, '')
  c.SetTicks()

  gP.SetLineColor(RT.kRed)
  gP.SetMarkerColor(RT.kRed)
  gM.SetLineColor(RT.kBlue)
  gM.SetMarkerColor(RT.kBlue)

  gN.Draw('AP')
  gP.Draw('P same')
  gM.Draw('P same')
  gN.SetTitle(';Cross Section Bin;#sigma (mb)');

  if i == 0:
    leg.AddEntry(gN, 'Nominal Fit')
    leg.AddEntry(gP, '+1#sigma Fit')
    leg.AddEntry(gM, '-1#sigma Fit')
    leg.Write('leg')
  c.Write()

nom_vals = [gN.GetY()[i] for gN in gNs for i in range(0, gN.GetN())]
plus_vals = [gP.GetY()[i] for gP in gPs for i in range(0, gP.GetN())]
minus_vals = [gM.GetY()[i] for gM in gMs for i in range(0, gM.GetN())]
plus_diffs = [n - p for n, p in zip(nom_vals, plus_vals)]
minus_diffs = [n - p for n, p in zip(nom_vals, minus_vals)]
plus_cov_vals = [[di*dj for di in plus_diffs] for dj in plus_diffs]
minus_cov_vals = [[di*dj for di in minus_diffs] for dj in minus_diffs]

plus_cov = RT.TH2D('plus_cov', '', len(nom_vals), 0, len(nom_vals), len(nom_vals), 0, len(nom_vals))
minus_cov = RT.TH2D('minus_cov', '', len(nom_vals), 0, len(nom_vals), len(nom_vals), 0, len(nom_vals))
ave_cov = RT.TH2D('ave_cov', '', len(nom_vals), 0, len(nom_vals), len(nom_vals), 0, len(nom_vals))

for i in range(0, len(nom_vals)):
  for j in range(0, len(nom_vals)):
    plus_cov.SetBinContent(i+1, j+1, plus_cov_vals[i][j])
    minus_cov.SetBinContent(i+1, j+1, minus_cov_vals[i][j])
    ave_cov.SetBinContent(i+1, j+1, (plus_cov_vals[i][j] + minus_cov_vals[i][j])/2.)

plus_cov.SetTitle(';Cross Section Bin;Cross Section Bin')
plus_cov.Write()
minus_cov.SetTitle(';Cross Section Bin;Cross Section Bin')
minus_cov.Write()
ave_cov.SetTitle(';Cross Section Bin;Cross Section Bin')
ave_cov.Write()

fOut.Close()
