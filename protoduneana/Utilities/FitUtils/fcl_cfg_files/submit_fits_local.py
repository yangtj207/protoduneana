#!/bin/env python3

import subprocess
import argparse
import os 
from glob import glob as ls

parser = argparse.ArgumentParser(description = 'Submission script for fits')
parser.add_argument('--config', type=str, help='Which config',
                    #default='/dune/app/users/calcuttj/larsoft-protoduneana/srcs/protoduneana/protoduneana/Utilities/FitUtils/fcl_cfg_files/asimov.cfg')
                    default='%s/cfg_files/asimov.cfg'%(os.environ['PROTODUNEANA_DIR']))
parser.add_argument('--type', type=str, help='Which type', default='asimov')

parser.add_argument('--output_dir', type=str, help='Output top dir', default=None)
parser.add_argument('--config_dump', action='store_true', help='Tell fife_launch to do a config_dump')
parser.add_argument('--dry_run', action='store_true', help='Tell fife_launch to do a dry_run')
#parser.add_argument('--fcl', type=str, required=True)
parser.add_argument('--output_name', type=str, default='hadd_wrapper_test.root')
parser.add_argument('--lifetime', type=str, default='6h')
parser.add_argument('--input_file', type=str, required=True)
parser.add_argument('--copy-input', action='store_true')
#parser.add_argument('--cov', type=str, default='')
parser.add_argument('--data_input', type=str, default='')
parser.add_argument('-N', type=int, default=1000)
parser.add_argument('--sites', type=str, nargs='+')
parser.add_argument('--memory', type=str, default=None)
parser.add_argument('--blacklist', type=str, nargs='+')
parser.add_argument('--multiple', action='store_true')
parser.add_argument('--extra_life', type=str, default=None)

parser.add_argument('--pduneana_tar', type=str, default='',
                    help='Optional Protoduneana tarball to be set up before NTupleProd')
parser.add_argument('--dropbox', action='store_true',
                    help='Optional, use dropbox uri for tar')

parser.add_argument('--list-types', action='store_true')
args = parser.parse_args()
##Begin fife_launch cmd 
cmd = ['fife_launch', '-c', args.config]


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

#if args.cov != '':
#  cmd += ['-Oglobal.cov_file=%s'%args.cov]
#  cmd += ['-Oglobal.cov_name=%s'%args.cov.split('/')[-1]]

##Miscellanea
cmd += ['-Oglobal.output_name=%(process)s_' + '%s.root'%args.type]
cmd += ['-Oglobal.fcl_name=%s.fcl'%args.type]
cmd += ['-Oexecutable.arg_2=%(fcl_name)s']
cmd += ['-Osubmit.expected-lifetime=%s'%args.lifetime]
cmd += ['-Oglobal.input_file=%s'%args.input_file]


if 'toy' in args.type or 'alt_sce' in args.type or args.multiple:
  cmd += ['-Osubmit.N=%i'%args.N]

if args.data_input != '':
  cmd += ['-Oexecutable.arg_7=-d',
          '-Oexecutable.arg_8=%s'%args.data_input]

##Special commands for overriding some setup stuff
if not args.pduneana_tar == '':
  extra = 'dropbox://' if args.dropbox else ''
  cmd += ['-Ojob_setup.setup_local=True',
          f'-Osubmit.tar_file_name={extra}{args.pduneana_tar}',
          #'-Ojob_setup.setup=delete',
         ]
else:
  cmd += ['-Ojob_setup.setup=protoduneana %(protoduneana_version)s -q e20:prof']
cmd += ['-Ojob_setup.prescript_4=ups active',]

if args.sites and len(args.sites) > 0:
  print("Sites:", args.sites)
  cmd += ['-Osubmit.site=%s'%','.join(args.sites)]

blacklist = ['SURFsara']
if args.blacklist:
  blacklist += args.blacklist

print("Blacklist:", blacklist)
cmd += ['-Osubmit.blacklist=%s'%','.join(blacklist)]

'''
if args.blacklist and len(args.blacklist) > 0:
  blacklist=['SURFsara'] + args.blacklist
  print("Blacklist:", blacklist)
  cmd += ['-Osubmit.blacklist=%s'%','.join(blacklist)]
'''
if args.memory:
  cmd += ['-Osubmit.memory=%s'%args.memory]
if args.copy_input:
  cmd += ['-Osubmit.f_4=%(input_file)s',
          '-Ojob_setup.prescript_3=ln -s ${CONDOR_DIR_INPUT}/eventSelection_mc_all.root eventSelection_mc_all.root']
if args.extra_life:
  if 'h' in args.extra_life:
    extra_life = 3600*int(args.extra_life.replace('h', ''))
  elif 'd' in args.extra_life:
    extra_life = 3600*24*int(args.extra_life.replace('d', ''))
  cmd +=[f'-Oglobal.extra_life={extra_life}']

print(cmd)

subprocess.run(cmd)
