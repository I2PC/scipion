#!/usr/bin/env python

import os

RELEASE_ORI = '1.0.1'
RELEASE_END = '1.1'

EMPACK_DIRS_CONVERT = [
    {RELEASE_END:'ctffind-3.5', RELEASE_ORI: 'ctffind_V3.5'},
    {RELEASE_END:'ctffind4-4.0.15', RELEASE_ORI: 'ctffind_V4.0.15'},
    {RELEASE_END:'eman-2.11', RELEASE_ORI: 'eman2.11.linux64'},
    {RELEASE_END:'eman-2.12', RELEASE_ORI: 'eman2.12.linux64'},
    {RELEASE_END:'frealign-9.07', RELEASE_ORI: 'frealign_v9.07'},
    {RELEASE_END:'motioncorr-2.1', RELEASE_ORI: 'motioncorr_v2.1'},
    {RELEASE_END:'relion-1.4f', RELEASE_ORI: 'relion-1.4_float'},
    {RELEASE_END:'resmap-1.1.5s2', RELEASE_ORI: 'resmap-1.1.5-s2'},
    {RELEASE_END:'simple-2.1', RELEASE_ORI: 'simple2'},
    {RELEASE_END:'spider-21.13', RELEASE_ORI: 'spider-web-21.13'},
    {RELEASE_END:'summovie-1.0.2', RELEASE_ORI: 'summovie_1.0.2'},
    {RELEASE_END:'unblur-1.0.15', RELEASE_ORI: 'unblur_1.0_150529'}
]

print "This script is going to upgrade EM packages names from v1.0.1 to v1.1.0"

# Script folder
scd = os.path.dirname(__file__)
# Software EM folder
semd = os.path.join(scd, '..', 'software', 'em')

print "Script folder %s" %  scd
print "SEM folder %s" % semd

# Move to Software EM folder
os.chdir(semd)

for emp in EMPACK_DIRS_CONVERT:
    print "Check if %s exists" % emp[RELEASE_ORI]
    if os.path.exists(emp[RELEASE_ORI]):
        print "Create symlink between %s and %s" % (emp[RELEASE_END], emp[RELEASE_ORI])
        os.symlink(emp[RELEASE_ORI], emp[RELEASE_END])

# Move to back Scipion home
os.chdir("../..")

print "EM packages installed after the upgrade"
os.system('./scipion install --help')