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
import optparse
# We use optparse instead of argparse because we want this script to
# be compatible with python >= 2.3

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser  # Python 3


def main():
    parser = optparse.OptionParser(description=__doc__)
    add = parser.add_option  # shortcut
    add('--overwrite', action='store_true',
        help=('Rewrite the configuration files using the original templates.'))
    options, args = parser.parse_args()

    globalIsLocal = (os.environ['SCIPION_CONFIG'] ==
                     os.environ['SCIPION_LOCAL_CONFIG'])  # if we used --config
    if globalIsLocal:
        localSections = []
    else:
        localSections = ['DIRS_LOCAL']

    try:
        templatesDir = join(os.environ['SCIPION_HOME'], 'config', 'templates')
        # Global installation configuration files.
        for fpath, tmplt in [
                (os.environ['SCIPION_CONFIG'], 'scipion'),
                (os.environ['SCIPION_PROTOCOLS'], 'protocols'),
                (os.environ['SCIPION_HOSTS'], 'hosts')]:
            if not exists(fpath) or options.overwrite:
                createConf(fpath, join(templatesDir, tmplt + '.template'),
                           remove=localSections)
            else:
                checkConf(fpath, join(templatesDir, tmplt + '.template'),
                          remove=localSections)

        if not globalIsLocal:  # which is normally the case
            # Local user configuration files (well, only "scipion.conf").
            if not exists(os.environ['SCIPION_LOCAL_CONFIG']):
                #  It might make sense to add   "or options.overwrite" ...
                createConf(os.environ['SCIPION_LOCAL_CONFIG'],
                           join(templatesDir, 'scipion.template'),
                           keep=localSections)
            else:
                checkConf(os.environ['SCIPION_LOCAL_CONFIG'],
                          join(templatesDir, 'scipion.template'),
                          keep=localSections)

        # After all, check some extra things are fine in scipion.conf
        print('Checking paths in %s ...' % os.environ['SCIPION_CONFIG'])
        checkPaths(os.environ['SCIPION_CONFIG'])
    except Exception:
        # This way of catching exceptions works with Python 2 & 3
        sys.stderr.write('Error: %s\n' % sys.exc_info()[1])
        sys.exit(1)


def createConf(fpath, ftemplate, remove=[], keep=[]):
    "Create config file in fpath following the template in ftemplate"
    # Remove from the template the sections in "remove", and if "keep"
    # is used only keep those sections.
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
    print('* Creating configuration file: %s' % fpath)
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(ftemplate) != [], 'Missing file: %s' % ftemplate
    for section in set(remove) - set(keep):
        cf.remove_section(section)
    if keep:
        for section in set(cf.sections()) - set(keep):
            cf.remove_section(section)
    cf.write(open(fpath, 'w'))
    print('Edit it to reflect the configuration of your system.\n')


def checkPaths(conf):
    "Check that some paths in the config file actually make sense"
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(conf) != [], 'Missing file: %s' % conf
    for var in ['MPI_LIBDIR', 'MPI_INCLUDE', 'MPI_BINDIR',
                'JAVA_HOME', 'JAVA_BINDIR']:
        path = cf.get('BUILD', var)
        if not os.path.isdir(path):
            print('Path to %s (%s) should exist but it does not.' % (var, path))
    # TODO: also check that some libraries and header files are actually there.


def checkConf(fpath, ftemplate, remove=[], keep=[]):
    "Check that all the variables in the template are in the config file too"
    # Remove from the checks the sections in "remove", and if "keep"
    # is used only check those sections.
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(fpath) != [], 'Missing file %s' % fpath
    ct = ConfigParser()
    ct.optionxform = str
    assert ct.read(ftemplate) != [], 'Missing file %s' % ftemplate

    # Keep only the sections we want to compare from the files.
    for section in set(remove) - set(keep):
        ct.remove_section(section)
        cf.remove_section(section)
    if keep:
        for section in set(ct.sections()) - set(keep):
            ct.remove_section(section)
            cf.remove_section(section)

    df = dict([(s, set(cf.options(s))) for s in cf.sections()])
    dt = dict([(s, set(ct.options(s))) for s in ct.sections()])
    # That funny syntax to create the dictionaries works with old pythons.

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
