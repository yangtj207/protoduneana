#!/bin/env python3

import subprocess
import argparse
import os 
from glob import glob as ls

parser = argparse.ArgumentParser(description = 'Submission script for fits')
parser.add_argument('--config', type=str, help='Which config',
                    default='/dune/app/users/calcuttj/larsoft-protoduneana/srcs/protoduneana/protoduneana/Utilities/FitUtils/fcl_cfg_files/asimov.cfg')
parser.add_argument('--type', type=str, help='Which type', default='asimov')

parser.add_argument('--output_dir', type=str, help='Output top dir', default=None)
parser.add_argument('--config_dump', action='store_true', help='Tell fife_launch to do a config_dump')
parser.add_argument('--dry_run', action='store_true', help='Tell fife_launch to do a dry_run')
#parser.add_argument('--fcl', type=str, required=True)
parser.add_argument('--output_name', type=str, default='hadd_wrapper_test.root')
parser.add_argument('--lifetime', type=str, default='3h')
parser.add_argument('--input_file', type=str, required=True)
parser.add_argument('--copy-input', action='store_true')
parser.add_argument('--cov', type=str, default='')
parser.add_argument('--data_input', type=str, default='')
parser.add_argument('-N', type=int, default=1000)
parser.add_argument('--sites', type=str, nargs='+')
parser.add_argument('--memory', type=str, default=None)
parser.add_argument('--blacklist', type=str, nargs='+')

parser.add_argument('--pduneana_tar', type=str, default='',
                    help='Optional Protoduneana tarball to be set up before NTupleProd')

parser.add_argument('--list-types', action='store_true')
args = parser.parse_args()
##Begin fife_launch cmd 
cmd = ['fife_launch', '-c', args.config]


'''
fcl_names = {
  'asimov':'grid_asimov.fcl',
  'split':'grid_split.fcl',
  'toy':'grid_toy.fcl',
  'toy_prot_syst':'grid_toy_prot_syst.fcl',
  'g4rw_prot':'grid_g4rw_prot.fcl',
  'g4rw_prot_high':'grid_g4rw_prot_high.fcl',
  'lowp':'grid_lowp.fcl',
  'alt_sce':'grid_alt_sce.fcl',
  'data':'grid_data.fcl',
  'data_thresh':'data_fit_thresh.fcl',
  'data_thresh_10pct_endz':'data_fit_thresh_10pct_endz.fcl',
  'data_thresh_match':'data_fit_thresh_match.fcl',
  'data_thresh_match_no_beam_cut':'data_fit_thresh_match_no_beam_cut.fcl',
  'data_thresh_match_no_beam_cut_binning_test':'data_fit_thresh_match_no_beam_cut_binning_test.fcl',
  'data_no_thresh':'data_fit_no_thresh.fcl',
  'data_thresh_no_beam_cut':'data_fit_thresh_no_beam_cut.fcl',
  'toy_thresh_match_no_beam_cut':'toy_thresh_match_no_beam_cut.fcl',
  'data_inclusive': 'data_inclusive.fcl',
  'toy_no_box_beam_no_beam_shift_norm_P': 'toy_no_box_beam_no_beam_shift_norm_P.fcl',
  'systs_test':'dEdX_test.fcl',
  'ediv_test':'ediv_test.fcl',
  'no_beam_match':'no_beam_match.fcl',
  'no_upstream':'no_upstream.fcl',
  'no_endz':'no_endz.fcl',
  'fluc':'grid_fluc.fcl',
  'set_toy':'set_toy.fcl',
  'set_toy_fluc': 'set_toy_fluc.fcl',
  'asimov_fine': 'asimov_fine.fcl',
  'fluc_fine': 'fluc_fine.fcl',
  'set_toy_fine': 'set_toy_fine.fcl',
  'set_toy_fine_fluc': 'set_toy_fine_fluc.fcl',
  'data_fine_plus_beam': 'data_fine_plus_beam.fcl',
  'data_fine_plus_2pct_beam': 'data_fine_plus_2pct_beam.fcl',
  'data_fine_minus_beam': 'data_fine_minus_beam.fcl',
  'data_fine_plus_cal': 'data_fine_plus_cal.fcl',
  'data_fine_minus_cal': 'data_fine_minus_cal.fcl',
}

if args.list_types:
  for n in fcl_names.keys():
    print(n)
  exit(0)

if args.type not in fcl_names.keys():
  print('Error: type %s invalid. Choose one of:'%args.type)
  print(fcl_names.keys())
  exit(1)
  '''

output_names = {
  'asimov':'%(process)s_asimov.root',
  'split':'%(process)s_split.root',
  'toy':'%(process)s_toy.root',
  'toy_prot_syst':'%(process)s_toy_prot_syst.root',
  'g4rw_prot':'%(process)s_g4rw_prot.root',
  'g4rw_prot_high':'%(process)s_g4rw_prot_high.root',
  'lowp':'%(process)s_lowp.root',
  'alt_sce':'%(process)s_alt_sce.root',
  'data':'%(process)s_data.root',
  'data_thresh':'%(process)s_data_thresh.root',
  'data_thresh_10pct_endz':'%(process)s_data_thresh_10pct_endz.root',
  'data_thresh_match':'%(process)s_data_thresh_match.root',
  'data_thresh_match_no_beam_cut':'%(process)s_data_thresh_match_no_beam_cut.root',
  'data_thresh_match_no_beam_cut_binning_test':'%(process)s_data_thresh_match_no_beam_cut_binning_test.root',
  'data_thresh_no_beam_cut':'%(process)s_data_thresh_no_beam_cut.root',
  'data_no_thresh':'%(process)s_data_no_thresh.root',
  'toy_thresh_match_no_beam_cut':'%(process)s_toy_thresh_match_no_beam_cut.root',
  'data_inclusive': '%(process)s_data_inclusive.root',
  'toy_no_box_beam_no_beam_shift_norm_P': '%(process)s_toy_no_box_beam_no_beam_shift_norm_P.root',
  'systs_test':'%(process)s_systs_test.root'
}

output_str = '*%s.root'%args.type
print('output_str:', output_str)
cmd += ['-Ojob_output.addoutput=%s'%output_str]

##Choose output dir
if args.output_dir:
  #cmd += ['-Oenv_pass.OUTPUT_DIR=%s'%args.output_dir]
  cmd += ['-Oglobal.output_dir=%s'%args.output_dir]

##Tell Fife_launch to do a dry run
if args.dry_run:
  cmd += ['--dry_run']

##Tell Fife_launch to do a dry run
if args.config_dump:
  cmd += ['--config_dump=after']

##Choose ntupleprod version
cmd += ['-Oglobal.protoduneana_version=%s'%os.getenv('PROTODUNEANA_VERSION')]

if args.cov != '':
  cmd += ['-Oglobal.cov_file=%s'%args.cov]
  cmd += ['-Oglobal.cov_name=%s'%args.cov.split('/')[-1]]

##Miscellanea
#cmd += ['-Oglobal.output_name=%s'%output_names[args.type]]
cmd += ['-Oglobal.output_name=%(process)s_' + '%s.root'%args.type]
#cmd += ['-Oglobal.fcl_name=%s'%fcl_names[args.type]]
cmd += ['-Oglobal.fcl_name=%s.fcl'%args.type]
cmd += ['-Osubmit.expected-lifetime=%s'%args.lifetime]
cmd += ['-Oglobal.input_file=%s'%args.input_file]
cmd += ['-Osubmit.f_2=dropbox://%(fcl)s']

if 'toy' in args.type or 'alt_sce' in args.type:
  cmd += ['-Osubmit.N=%i'%args.N]

if args.data_input != '':
  cmd += ['-Oexecutable.arg_7=-d',
          '-Oexecutable.arg_8=%s'%args.data_input]

##Special commands for overriding some setup stuff
if not args.pduneana_tar == '':
  cmd += ['-Ojob_setup.setup_local=True',
          '-Osubmit.tar_file_name=%s'%args.pduneana_tar, 
          #'-Ojob_setup.setup=delete',
         ]
else:
  cmd += ['-Ojob_setup.setup=protoduneana %(protoduneana_version)s -q e20:prof']
cmd += ['-Ojob_setup.prescript_4=ups active',]

if args.sites and len(args.sites) > 0:
  print("Sites:", args.sites)
  cmd += ['-Osubmit.site=%s'%','.join(args.sites)]
if args.blacklist and len(args.blacklist) > 0:
  print("Blacklist:", args.blacklist)
  cmd += ['-Osubmit.blacklist=%s'%','.join(args.blacklist)]
if args.memory:
  cmd += ['-Osubmit.memory=%s'%args.memory]
if args.copy_input:
  cmd += ['-Osubmit.f_4=%(input_file)s',
          '-Ojob_setup.prescript_3=ln -s ${CONDOR_DIR_INPUT}/eventSelection_mc_all.root eventSelection_mc_all.root']
print(cmd)

subprocess.run(cmd)
