#!/usr/bin/python

# Watch out: not "/usr/bin/env python", because we will overwrite the
# python in SCIPION_HOME/software/bin and we cannot use itself to run it.

# **************************************************************************
# *
# * Authors:     J. Burguet Castell     (jburguet@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

#
# To use after compiling with something like:
#   LDFLAGS=-Wl,-rpath,REPLACE_ME_WITH_FUNNY_ ./scipion install 
#

"""
Replaces in-place a text in binary files.
"""

import sys
import os
from tempfile import TemporaryFile


def main():
    if len(sys.argv) < 2:
        sys.exit('usage: %s <dir1> [<dir2> [...]]' % sys.argv[0])

    replacements = 0
    for subdir in sys.argv[1:]:
        print 'Replacing RPATH in elf files under %s ...' % subdir
        for dirpath, dirnames, filenames in os.walk(subdir):
            for fname in filenames:
                fpath = os.path.join(dirpath, fname)
                if isElf(fpath):
                    replacements += replace(fpath)

    print 'Total replacements: %d' % replacements
    if replacements == 0:
        print 'Warning: You probably did not really compile. Clean and repeat.'


def isElf(fname):
    "Is fname a file in elf format?"
    return os.path.isfile(fname) and open(fname).read(4) == '\x7fELF'


def replace(fname,
            txt_from='REPLACE_ME_WITH_FUNNY_',
            txt_to='$ORIGIN:$ORIGIN/../lib'):
    "Change txt_from to txt_to in binary file fname. Return 1 if changed."

    # First add write permission to the file if it does not have it.
    if os.stat(fname).st_mode & 0o200 == 0:  # mode has no write permission
        print 'Adding write permission to %s ...' % fname
        os.chmod(fname,  os.stat(fname).st_mode + 0o200)

    changed = 0
    with TemporaryFile('wb+') as fout:
        with open(fname, 'rb') as fin:
            for line in fin:
                if txt_from in line:
                    print 'Replacing in %s' % fname
                    line = line.replace(txt_from, txt_to)
                    changed += 1
                fout.write(line)
        fout.seek(0)
        open(fname, 'wb').writelines(fout)  # overwrite fname

    return changed



if __name__ == '__main__':
    main()
