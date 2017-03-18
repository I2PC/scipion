#!/usr/bin/env python

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

"""
Show all the system packages on which scipion depends.
"""

import sys
import os
from subprocess import check_output, CalledProcessError
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add = parser.add_argument  # short notation
    add('dirs', metavar='DIR', nargs='+', help="directories with binaries")
    add('--packages', action='store_true', help="show system packages")
    add('--pkg-command', default='dpkg -S',
        help="command used to find which package a file belongs to")
    args = parser.parse_args()

    # Get all the library dependencies.
    neededLibs = set()
    for subdir in args.dirs:
        for dirpath, dirnames, filenames in os.walk(subdir):
            for fname in filenames:
                fpath = os.path.join(dirpath, fname)
                if isElf(fpath):
                    neededLibs.update(libDeps(fpath))

    # Show where they are in the system.
    existingLibs = {path: os.listdir(path)
                    for path in searchPaths() if os.path.exists(path)}

    shown = set()
    for lib in neededLibs:
        for libPath in existingLibs:
            if lib not in existingLibs[libPath]:
                continue
            if not args.packages:
                libFull = '%s/%s' % (libPath, lib)
                if libFull not in shown:
                    print libFull
                    shown.add(libFull)
            elif 'software/lib' not in libPath:
                command = '%s %s/%s' % (args.pkg_command, libPath, lib)
                try:
                    pkg = check_output(command.split()).split(':')[0].strip()
                    if pkg not in shown:
                        print pkg
                        shown.add(pkg)
                except CalledProcessError:
                    print '%s/%s -> no package found' % (libPath, lib)
            break
        else:
            print '%s -> not found in any path' % lib


def libDeps(fpath):
    "Return set of libraries that the file at fpath links to"
    libs = set()
    for line in check_output(['objdump', '-p', fpath]).splitlines():
        fields = line.split()
        if fields and fields[0] == 'NEEDED':
            libs.add(fields[-1])
    return libs

    
def isElf(fname):
    "Is fname a file in elf format?"
    return os.path.isfile(fname) and open(fname).read(4) == '\x7fELF'


def searchPaths():
    "Return the list of paths were the libraries are searched in the system"
    paths = []

    if 'LD_LIBRARY_PATH' in os.environ:
        paths += os.environ.get('LD_LIBRARY_PATH', '').split(':')
    
    for fname in os.listdir('/etc/ld.so.conf.d'):
        if fname.endswith('.conf'):
            for line in open('/etc/ld.so.conf.d/%s' % fname):
                if not line.startswith('#'):
                    paths.append(line.strip())
    paths = [x for x in paths if 'i386' not in x]  # FIXME: something nicer
    return paths



if __name__ == '__main__':
    main()
