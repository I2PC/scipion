#!/usr/bin/env python

# basic setup, import all environment and custom tools
import os
from os.path import join, abspath, dirname
import sys
import platform

import SCons.Script


env = Environment(ENV=os.environ, toolpath=['scons/tools'],
                  tools=['default', 
                         'disttar', 
                         "URLDownload", 
                         "Unpack", 
                         'Make', 
                         'AutoConfig'
                         ]
                  )

PYTHON_URL = "http://scipionwiki.cnb.csic.es/files/scipion/software/python/Python-2.7.6.tgz"
PYTHON_DIR = 'Python-2.7.6'

SOFT_DIR = join(os.environ['SCIPION_HOME'], 'software')
TMP_DIR = join(SOFT_DIR, 'tmp')

env['URLDOWNLOAD_DIRECTORY'] = TMP_DIR
env['UNPACK']['EXTRACTDIR'] = TMP_DIR
env['CROSS_BUILD'] = False

download = env.URLDownload('python-download', PYTHON_URL)
# Unpack the python module
env['UNPACK']['EXTRACTDIR'] = TMP_DIR
p = join('software', 'tmp', PYTHON_DIR)

unpack = env.Unpack('python-unpack', download)

dp = Dir(p)
configured = env.AutoConfig(dp,
			    AutoConfigSource='Makefile.pre.in',
			    AutoConfigTarget='Makefile',
                AutoConfigParams=['--prefix=%s' % SOFT_DIR])

maked = env.Make(source=None, 
                 target=[join(p, 'python'), join(p, 'libpython2.7.a')],
                 MakePath=dp)

#Depends(configured, unpack)

Depends(maked, configured)


last = maked


env.Default(last)