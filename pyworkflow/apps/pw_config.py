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
from os.path import join, exists, basename
import time
import optparse
# We use optparse instead of argparse because we want this script to
# be compatible with python >= 2.3

try:
    from ConfigParser import ConfigParser, Error
except ImportError:
    from configparser import ConfigParser, Error  # Python 3

def ansi(n):
    "Return function that escapes text with ANSI color n."
    return lambda txt: '\x1b[%dm%s\x1b[0m' % (n, txt)

black, red, green, yellow, blue, magenta, cyan, white = map(ansi, range(30, 38))
# We don't take them from pyworkflow.utils because this has to run
# with all python versions (and so it is simplified).


def main():
    parser = optparse.OptionParser(description=__doc__)
    add = parser.add_option  # shortcut
    add('--overwrite', action='store_true',
        help=("Rewrite the configuration files using the original templates."))
    options, args = parser.parse_args()

    if args:  # no args which aren't options
        sys.exit(parser.format_help())

    globalIsLocal = (os.environ['SCIPION_CONFIG'] ==
                     os.environ['SCIPION_LOCAL_CONFIG'])  # if we used --config
    if globalIsLocal:
        localSections = []
    else:
        localSections = ['DIRS_LOCAL', 'PACKAGES', 'VARIABLES']

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
        checkPaths(os.environ['SCIPION_CONFIG'])
        if (not globalIsLocal and '[BUILD]' in
            [x.strip() for x in open(os.environ['SCIPION_LOCAL_CONFIG'])]):
            print(red("Found a BUILD section in the local configuration file %s"
                      "\nthat would override %s -- Not checking it." %
                      (os.environ['SCIPION_LOCAL_CONFIG'],
                       os.environ['SCIPION_CONFIG'])))
    except Exception:
        # This way of catching exceptions works with Python 2 & 3
        sys.stderr.write('Error: %s\n' % sys.exc_info()[1])
        sys.exit(1)


def createConf(fpath, ftemplate, remove=[], keep=[]):
    "Create config file in fpath following the template in ftemplate"
    # Remove from the template the sections in "remove", and if "keep"
    # is used only keep those sections.

    # Create directory and backup if necessary.
    dname = os.path.dirname(fpath)
    if not exists(dname):
        os.makedirs(dname)
    elif exists(fpath):
        if not exists(join(dname, 'backups')):
            os.makedirs(join(dname, 'backups'))
        backup = join(dname, 'backups',
                      '%s.%d' % (basename(fpath), int(time.time())))
        print(yellow("* Creating backup: %s" % backup))
        os.rename(fpath, backup)

    # Read the template configuration file.
    print(yellow("* Creating configuration file: %s" % fpath))
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(ftemplate) != [], 'Missing file: %s' % ftemplate
    for section in set(remove) - set(keep):
        cf.remove_section(section)
    if keep:
        for section in set(cf.sections()) - set(keep):
            cf.remove_section(section)

    # Update with our guesses.
    if 'BUILD' in cf.sections():
        for options in [guessJava(), guessMPI()]:
            for key in options:
                if key in cf.options('BUILD'):
                    cf.set('BUILD', key, options[key])

    # Create the actual configuration file.
    cf.write(open(fpath, 'w'))
    print("Please edit it to reflect the configuration of your system.\n")


def checkPaths(conf):
    "Check that some paths in the config file actually make sense"

    print("Checking paths in %s ..." % conf)
    cf = ConfigParser()
    cf.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
    assert cf.read(conf) != [], 'Missing file: %s' % conf
    def get(var):
        try:
            return cf.get('BUILD', var)
        except Exception:
            _, e = sys.exc_info()[:2]
            print(red("While getting '%s' in section BUILD: %s" % (var, e)))
            return '/'
    allOk = True
    for var in ['MPI_LIBDIR', 'MPI_INCLUDE', 'MPI_BINDIR',
                'JAVA_HOME', 'JAVA_BINDIR']:
        if not os.path.isdir(get(var)):
            print("  Path to %s (%s) should exist but it doesn't." %
                  (var, red(get(var))))
            allOk = False
    for fname in [join(get('JAVA_BINDIR'), 'java'),
                  get('JAVAC'), get('JAR'),
                  join(get('MPI_BINDIR'), get('MPI_CC')),
                  join(get('MPI_BINDIR'), get('MPI_CXX')),
                  join(get('MPI_BINDIR'), get('MPI_LINKERFORPROGRAMS')),
                  join(get('MPI_INCLUDE'), 'mpi.h')]:
        if not exists(fname):
            print("  Cannot find file: %s" % red(fname))
            allOk = False
    if allOk:
        print(green("All seems fine with %s" % conf))
    else:
        print(red("Errors found."))
        print("Please edit %s and check again." % conf)
        print("To regenerate the config files trying to guess the paths, you "
              "can run: scipion config --overwrite")


def checkConf(fpath, ftemplate, remove=[], keep=[]):
    "Check that all the variables in the template are in the config file too"
    # Remove from the checks the sections in "remove", and if "keep"
    # is used only check those sections.

    # Read the config file fpath and the template ftemplate
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
        print(green("All the expected sections and options found in " + fpath))
    else:
        print("Found differences between the configuration file\n  %s\n"
              "and the template file\n  %s" % (fpath, ftemplate))
        sf = set(df.keys())
        st = set(dt.keys())
        for s in sf - st:
            print("Section %s exists in the configuration file but "
                  "not in the template." % red(s))
        for s in st - sf:
            print("Section %s exists in the template but "
                  "not in the configuration file." % red(s))
        for s in st & sf:
            for o in df[s] - dt[s]:
                print("In section %s, option %s exists in the configuration "
                      "file but not in the template." % (red(s), red(o)))
            for o in dt[s] - df[s]:
                print("In section %s, option %s exists in the template "
                      "but not in the configuration file." % (red(s), red(o)))


def guessJava():
    "Guess the system's Java installation, return a dict with the Java keys"

    options = {}
    candidates = []

    # First check if the system has a favorite one.
    if 'JAVA_HOME' in os.environ:
        candidates.append(os.environ['JAVA_HOME'])

    # Add also all the ones related to a "javac" program.
    for d in os.environ.get('PATH', '').split(':'):
        if not os.path.isdir(d) or 'javac' not in os.listdir(d):
            continue
        javaBin = os.path.realpath(join(d, 'javac'))
        if javaBin.endswith('/bin/javac'):
            javaHome = javaBin[:-len('/bin/javac')]
            candidates.append(javaHome)
            if javaHome.endswith('/jre'):
                candidates.append(javaHome[:-len('/jre')])

    # Check in order if for any of our candidates, all related
    # directories and files exist. If they do, that'd be our best guess.
    for javaHome in candidates:
        allExist = True
        for path in ['include', join('bin', 'javac'), join('bin', 'jar')]:
            if not exists(join(javaHome, path)):
                allExist = False
        if allExist:
            options['JAVA_HOME'] = javaHome
            break
        # We could instead check individually for JAVA_BINDIR, JAVAC
        # and so on, as we do with MPI options, but we go for an
        # easier and consistent case instead: everything must be under
        # JAVA_HOME, which is the most common case for Java.

    if not options:
        print(red("Warning: could not detect a suitable JAVA_HOME."))
        if candidates:
            print(red("Our candidates were:\n  %s" % '\n  '.join(candidates)))

    return options


def guessMPI():
    "Guess the system's MPI installation, return a dict with MPI keys"
    # Returns MPI_LIBDIR, MPI_INCLUDE and MPI_BINDIR as a dictionary.

    options = {}
    candidates = []

    # First check if the system has a favorite one.
    for prefix in ['MPI_', 'MPI', 'OPENMPI_', 'OPENMPI']:
        if '%sHOME' % prefix in os.environ:
            candidates.append(os.environ['%sHOME' % prefix])

    # Add also all the ones related to a "mpicc" program.
    for d in os.environ.get('PATH', '').split(':'):
        if not os.path.isdir(d) or 'mpicc' not in os.listdir(d):
            continue
        mpiBin = os.path.realpath(join(d, 'mpicc'))
        if 'MPI_BINDIR' not in options:
            options['MPI_BINDIR'] = os.path.dirname(mpiBin)
        if mpiBin.endswith('/bin/mpicc'):
            mpiHome = mpiBin[:-len('/bin/mpicc')]
            candidates.append(mpiHome)

    # Add some extra directories that are commonly around.
    candidates += ['/usr/lib/openmpi', '/usr/lib64/mpi/gcc/openmpi']

    # Check in order if for any of our candidates, all related
    # directories and files exist. If they do, that'd be our best guess.
    for mpiHome in candidates:
        if (exists(join(mpiHome, 'include', 'mpi.h')) and
            'MPI_INCLUDE' not in options):
            options['MPI_INCLUDE'] = join(mpiHome, 'include')
        if (exists(join(mpiHome, 'lib', 'libmpi.so')) and
            'MPI_LIBDIR' not in options):
            options['MPI_LIBDIR'] = join(mpiHome, 'lib')
        if (exists(join(mpiHome, 'bin', 'mpicc')) and
            'MPI_BINDIR' not in options):
            options['MPI_BINDIR'] = join(mpiHome, 'bin')

    return options



if __name__ == '__main__':
    main()
