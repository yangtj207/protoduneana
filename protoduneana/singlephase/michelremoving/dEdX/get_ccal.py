import ROOT as RT
from math import sqrt
from argparse import ArgumentParser as ap

parser = ap()
parser.add_argument('--i0', type=str, required=True)
parser.add_argument('--i1', type=str, required=True)
parser.add_argument('--i2', type=str, required=True)

args = parser.parse_args()
f0 = RT.TFile(args.i0)
f1 = RT.TFile(args.i1)
f2 = RT.TFile(args.i2)

c0 = f0.Get('Scans/C_cal')
c0.Fit('pol2')
f = c0.GetFunction('pol2')
a = f.GetParameter(2)
b = f.GetParameter(1)
c = f.GetParameter(0)

x0 = f.GetMinimumX()
chi2 = f.GetMinimum()
x1 = (-b + sqrt(b*b - 4*a*(c - chi2 - 1)))/(2*a)
print(x1)
print(x0)
x2 = (-b - sqrt(b*b - 4*a*(c - chi2 - 1)))/(2*a)
print(x2)
print(x1 - x0)

with open('calconst_0.txt', 'w') as fout:
  fout.write('%f %f'%(x0, x1-x0))

c1 = f1.Get('Scans/C_cal')
c1.Fit('pol2')
f = c1.GetFunction('pol2')
a = f.GetParameter(2)
b = f.GetParameter(1)
c = f.GetParameter(0)

x0 = f.GetMinimumX()
chi2 = f.GetMinimum()
x1 = (-b + sqrt(b*b - 4*a*(c - chi2 - 1)))/(2*a)
print(x1)
print(x0)
x2 = (-b - sqrt(b*b - 4*a*(c - chi2 - 1)))/(2*a)
print(x2)
print(x1 - x0)

with open('calconst_1.txt', 'w') as fout:
  fout.write('%f %f'%(x0, x1-x0))

c2 = f2.Get('Scans/C_cal')
c2.Fit('pol2')
f = c2.GetFunction('pol2')
a = f.GetParameter(2)
b = f.GetParameter(1)
c = f.GetParameter(0)

x0 = f.GetMinimumX()
chi2 = f.GetMinimum()
x1 = (-b + sqrt(b*b - 4*a*(c - chi2 - 1)))/(2*a)
print(x1)
print(x0)
x2 = (-b - sqrt(b*b - 4*a*(c - chi2 - 1)))/(2*a)
print(x2)
print(x1 - x0)

with open('calconst_2.txt', 'w') as fout:
  fout.write('%f %f'%(x0, x1-x0))

