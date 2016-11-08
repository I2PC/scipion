#!/usr/bin/env python

import os
import sys

cleanAll = 'all' in sys.argv


def cmd(command):
    print command
    os.system(command)

if cleanAll:
    cmd('rm -rf software/lib/* ')
    cmd('rm -rf software/bin/*  ')
    cmd('rm -rf software/include/* ')
    cmd('rm -rf software/man/* ')
    cmd('rm -rf software/share/* ')
    cmd('rm -rf software/tmp/*')
    cmd('rm -rf software/log/*')

for ext in ['so', 'os', 'o']:
    cmd('find software/em/xmipp -name "*.%s" -exec rm -rf {} \;' % ext)
