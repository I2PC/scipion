# **************************************************************************
# *
# * Authors: J. Burguet Castell (jburguet@cnb.csic.es)
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
# **************************************************************************
"""
Check the local configuration files, and/or create them if requested
or if they do not exist.
"""

import sys
import os
from os.path import join, exists, dirname, basename
import time
import argparse
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser  # Python 3


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add = parser.add_argument  # shortcut
    add('--overwrite', action='store_true',
        help=('Rewrite the configuration files using the original templates.'))
    args = parser.parse_args()

    try:
        settingsDir = join(os.environ['SCIPION_HOME'], 'settings')
        for fpath, tmplt in [
                (os.environ['SCIPION_CONFIG'], 'scipion.conf'),
                (os.environ['SCIPION_PROTOCOLS'], 'protocols.conf'),
                (os.environ['SCIPION_HOSTS'], 'hosts.conf')]:
            if not exists(fpath) or args.overwrite:
                createConf(fpath, join(settingsDir, tmplt))
            else:
                checkConf(fpath, join(settingsDir, tmplt))
        # After all, check some extra things are fine in scipion.conf
        checkPaths(os.environ['SCIPION_CONFIG'])
        # TODO: say that we are checking scipion.conf, and be more
        # nice to the user.
    except Exception:
        # This way of catching exceptions works with Python 2 & 3
        sys.stderr.write('Error: %s\n' % sys.exc_info()[1])
        sys.exit(1)


def createConf(fpath, ftemplate):
    "Create config file in fpath following the template in ftemplate"
    dname = dirname(fpath)
    print('')
    if not exists(dname):
        os.makedirs(dname)
    elif exists(fpath):
        if not exists(join(dname, 'backups')):
            os.makedirs(join(dname, 'backups'))
        backup = join(dname, 'backups',
                      '%s.%d' % (basename(fpath), int(time.time())))
        print('* Creating backup: %s' % backup)
        os.rename(fpath, backup)
    print('* Creating local configuration file: %s' % fpath)
    open(fpath, 'w').write(open(ftemplate).read())  # cp ftemplate fpath
    print('Edit it to reflect the configuration of your system.')


def checkPaths(conf):
    "Check that some paths in the config file actually make sense"
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(conf) != [], 'Missing file %s' % conf
    for var in ['MPI_LIBDIR', 'MPI_INCLUDE', 'MPI_BINDIR',
                'JAVA_HOME', 'JAVA_BINDIR']:
        path = cf.get('BUILD', var)
        if not os.path.isdir(path):
            print('Path to %s (%s) should exist but it does not.' % (var, path))
    # TODO: also check that some libraries and header files are actually there.


def checkConf(fpath, ftemplate):
    "Check that all the variables in the template are in the config file too"
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(fpath) != [], 'Missing file %s' % fpath
    ct = ConfigParser()
    ct.optionxform = str
    assert ct.read(ftemplate) != [], 'Missing file %s' % ftemplate
    df = dict([(s, set(cf.options(s))) for s in cf.sections()])
    dt = dict([(s, set(ct.options(s))) for s in ct.sections()])
    if df == dt:
        print('The configuration file %s looks fine.' % fpath)
    else:
        print('Found differences between the configuration file\n  %s\n'
              'and the template file\n  %s' % (fpath, ftemplate))
        sf = set(df.keys())
        st = set(dt.keys())
        for s in sf - st:
            print('Section %s exists in the configuration file but '
                  'not in the template.' % s)
        for s in st - sf:
            print('Section %s exists in the template but '
                  'not in the configuration file.' % s)
        for s in st & sf:
            for o in df[s] - dt[s]:
                print('In section %s, option %s exists in the configuration '
                      'file but not in the template.' % (s, o))
            for o in dt[s] - df[s]:
                print('In section %s, option %s exists in the template '
                      'but not in the configuration file.' % (s, o))


if __name__ == '__main__':
    main()
