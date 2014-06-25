#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *              J. Burguet Castell (jburguet@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************

"""
Scipion data synchronization.
"""

from __future__ import division

import sys
import os
from os.path import join, isdir, exists, relpath
import subprocess
import argparse
import hashlib
from urllib2 import urlopen
import getpass

from pyworkflow.utils.path import makePath, makeFilePath, copyFile, cleanPath


scipion_logo = """
QQQQQQQQQT!^'::""?$QQQQQQ  S   S   S
QQQQQQQY`          ]4QQQQ  C   C   C
QQQQQD'              "$QQ  I   I   I
QQQQP                 "4Q  P   P   P
QQQP        :.,        -$  I   I   I
QQD       awQQQQwp      )  O   O   O
QQ'     qmQQQWQQQQg,   jm  N   N   N
Qf     QQQD^   -?$QQp jQQ ################################################
Q`    qQQ!        4WQmQQQ # Integrating image processing packages for EM #
F     QQ[          ~)WQQQ ################################################
[    ]QP             4WQQ
f    dQ(             -$QQ Data synchronization script
'    QQ              qQQQ
.   )QW            _jQQQQ
-   =QQ           jmQQQQQ
/   -QQ           QQQQQQQ
f    4Qr    jQk   )WQQQQQ
[    ]Qm    ]QW    "QQQQQ
h     $Qc   jQQ     ]$QQQ
Q,  :aQQf qyQQQL    _yQQQ
QL jmQQQgmQQQQQQmaaaQWQQQ
"""



def downloadFile(datasetName, fname, workingCopy=None, askMsg="download it?",
                 url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
                 verbose=False):
    """ Download data set from a url, return 1 if it worked, 0 if not """

    datasetFolder = join(workingCopy if workingCopy else os.environ['SCIPION_TESTS'],
                         datasetName)

    if verbose:
        while True:  # so verbose also interrupts the flow... uhm...
            answer = raw_input("\t%s ([y]/n): " % askMsg)
            if answer in ['N', 'n']:
                return 0
            elif answer in ['Y', 'y', '']:
                break

    if verbose:
        print "\tDownloading file..."

    makeFilePath(join(datasetFolder, fname))
    try:
        data = urlopen('%s/%s/%s' % (url, datasetName, fname)).read()
        open(join(datasetFolder, fname), 'w').write(data)
        if verbose:
            print "\t...done."
        return 1
    except Exception as e:
        print "\tError downloading %s (%s)" % (fname, e)
        print "URL: %s/%s/%s" % (url, datasetName, fname)
        print "destination: %s" % join(datasetFolder, fname)
        return 0


def downloadDataset(datasetName, destination=None,
                    url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
                    verbose=False, onlyManifest=False):
    """ Download all the data files mentioned in url/datasetName/MANIFEST """

    destination = destination if destination else os.environ['SCIPION_TESTS']

    # First make sure that we ask for a known dataset.
    datasetsAvailable = [x.strip('./') for x in
                         urlopen('%s/MANIFEST' % url).read().splitlines()]
    if datasetName not in datasetsAvailable:
        print 'Unknown dataset "%s".' % datasetName
        print 'Use --list to see the available datasets.'
        return

    # Retrieve the dataset's MANIFEST file.
    # It contains a list of "file md5sum" of all files included in the dataset.
    datasetFolder = join(destination, datasetName)
    makePath(datasetFolder)
    manifest = join(destination, datasetName, 'MANIFEST')
    try:
        if verbose:
            print "Retrieving MANIFEST file"
        data = urlopen('%s/%s/MANIFEST' % (url, datasetName)).read()
        open(manifest, 'w').write(data)
    except Exception as e:
        print "ERROR reading %s/%s/MANIFEST (%s)" % (url, datasetName, e)
        return

    if onlyManifest:
        print "\t...done\n"
        return

    # Now retrieve all of the files mentioned in MANIFEST, and check their md5.
    print "Fetching dataset %s files..." % datasetName
    manifestLines = open(manifest).readlines()
    totalNumber = len(manifestLines)
    prevPercent = 0
    for number, line in enumerate(manifestLines):
        fname, md5InLine = line.strip().split()
        makeFilePath(join(datasetFolder, fname))
        try:
            # Create file with downloaded content
            data = urlopen('%s/%s/%s' % (url, datasetName, fname)).read()
            open(join(datasetFolder, fname), 'w').write(data)

            percent = ((number+1)/(totalNumber*1.0))*100
            if verbose:
                print "\t%s ...OK (%02d %%) " % (fname, percent)
            else:
                progress = percent - prevPercent
                sys.stdout.write("#" * int(1 + progress))
                sys.stdout.flush()
            prevPercent = percent

            md5 = md5sum(join(datasetFolder, fname))
            assert md5 == md5InLine, \
                'Bad md5. Expected: %s Computed: %s' % (md5InLine, md5)
        except Exception as e:
            print "\tError in %s (%s)" % (fname, e)
            print "URL: %s/%s/%s" % (url, datasetName, fname)
            print "Destination: %s" % join(datasetFolder, fname)
            while True:
                answer = raw_input("Continue downloading? (y/[n]): ")
                if not answer or answer in ['N', 'n']:
                    break
                if answer in ['Y', 'y']:
                    continue

        if verbose:
            print "md5 verification...OK\n"


def md5sum(fname):
    """ Return the md5 hash of file fname """

    mhash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(128 * mhash.block_size), ''):
            mhash.update(chunk)
    return mhash.hexdigest()


def checkForUpdates(datasetName, workingCopy=None,
                    url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
                    verbose=False):
    """ Compare md5 of files in url/datasetName/MANIFEST with their local counterpart """

    # Get default values for variables if we got none
    workingCopy = workingCopy or os.environ['SCIPION_TESTS']

    # Read contents of *remote* MANIFEST file, and create a dict {fname: md5}
    manifestRaw = urlopen('%s/%s/MANIFEST' % (url, datasetName)).read()
    md5sRemote = dict(x.split() for x in manifestRaw.splitlines())

    # Read contents of *local* MANIFEST file, and create a dict {fname: md5}
    datasetFolder = join(workingCopy, datasetName)
    md5sLocal = dict(x.split() for x in open(join(datasetFolder, 'MANIFEST')))

    # Check that all the files mentioned in MANIFEST are up-to-date
    print "Verifying MD5s..."

    vlog = lambda txt: sys.stdout.write(txt) if verbose else None  # verbose log
    filesUpdated = 0  # number of files that have been updated
    taintedMANIFEST = False  # can MANIFEST be out of sync? (if so, say it)

    for fname in md5sRemote:
        vlog('\t%s' % fname)

        if not exists(join(datasetFolder, fname)):
            vlog("\n\tFile %s doesn't exist.\n" % fname)
            filesUpdated += downloadFile(datasetName, fname, workingCopy,
                                         askMsg='download it?', url=url, verbose=verbose)
        else:
            if md5sLocal[fname] == md5sRemote[fname]:
                vlog('\r\tOK  %s\n' % fname)
                pass  # just to emphasize that we do nothing in this case
            else:
                vlog('\r\tBAD %s  \t==> checksum differs\n' % fname)
                taintedMANIFEST = True  # if we don't update, it can be wrong
                filesUpdated += downloadFile(datasetName, fname, workingCopy,
                                             askMsg='update it?', url=url, verbose=verbose)

    print '\t...done. Updated files: %d' % filesUpdated

    # Save the new MANIFEST file in the folder of the downloaded dataset.
    open(join(datasetFolder, 'MANIFEST'), 'w').write(manifestRaw)

    if taintedMANIFEST:
        print 'If you stopped any files from being updated, please run '
        print 'with the --format option to regenerate the MANIFEST file.'


def Cmd(command):
    print ">>>>>", command
    subprocess.call(command, shell=True)


def main():
    # Arguments parsing
    parser = argparse.ArgumentParser(description=__doc__)
    g = parser.add_mutually_exclusive_group()
    g.add_argument('--download', action='store_true', help="Download dataset.")
    g.add_argument('--upload', action='store_true',
        help=("Upload local dataset to the scipion server. The dataset name must"
              "be the name of its folder relative to the SCIPION_TESTS folder."))
    g.add_argument('--list', action='store_true',
        help=('List local datasets (from $SCIPION_TESTS) and remote ones '
              '(remote url can be specified with --url).'))
    g.add_argument('--format', action='store_true',
        help='Create a MANIFEST file with checksums in the datasets folders.')
    add = parser.add_argument  # shortcut
    add('datasets', metavar='DATASET', nargs='*', help='Name of a dataset.')
    add('--delete', action='store_true',
        help=('Delete remote files not present in local. It leaves the remote '
              'scipion data directory as it is in the local one. Extremely '
              'dangerous! Only valid with -r.'))
    add('-u', '--url',
        default='http://scipionwiki.cnb.csic.es/files/scipion/data/tests',
        help='URL where remote datasets will be looked for.')
    add('-q', '--query-for-modifications', action='store_true',
        help=("Look for last_m.txt file. There, there will be (1) the IP of "
              "the last modifier and date or (2) nothing. If the script finds "
              "(1) it moves it to modifications.log and returns 0. If it "
              "finds (2) it checks whether the file was modified in the last "
              "30 minutes. If not, it returns 1, if yes it returns 2."))
    add('--last-mod-file', default='last_m.txt',
        help=('File which lists the computers that have last changed the '
              'Scipion tests data (using this script).'))
    add('--mod-log-file', default='modifications.log',
        help=('File with the whole modifications log to keep track of what '
              'has been done in the Scipion tests data. The path must be '
              'relative to the Scipion folder.'))
    add('-v', '--verbose', action='store_true', help='Print more details.')
    args = parser.parse_args()

    #print scipion_logo

    # Dispatch the easy cases first (list and query), and then take care of
    # the more complex ones.
    if args.list:
        listDatasets(args.url)
        sys.exit(0)

    if args.query_for_modifications:
        last_mfile = join(os.environ['HOME'], 'Scipion', args.last_mod_file)
        log_mfile = join(os.environ['HOME'], 'Scipion', args.mod_log_file)
        queryModifications(last_mfile, log_mfile)
        sys.exit(0)

    if not args.datasets:
        sys.exit("At least --list, -q or datasets needed. See --help for more info.")

    print "Selected datasets: %s" % " ".join(args.datasets)

    if args.format:
        for dataset in args.datasets:
            print "Formatting %s" % dataset
            if not exists(join(os.environ['SCIPION_TESTS'], dataset)):
                sys.exit("ERROR: %s folder doesn't exist in datasets folder %s." %
                         (dataset, os.environ['SCIPION_TESTS']))
            createMANIFEST(join(os.environ['SCIPION_TESTS'], dataset))
        sys.exit(0)

    if args.download:
        # Download datasets.
        for dataset in args.datasets:
            if exists(join(os.environ['SCIPION_TESTS'], dataset)):
                print 'Local copy of dataset %s detected.' % dataset
                print 'Checking for updates...'
                checkForUpdates(dataset, url=args.url, verbose=args.verbose)
            else:
                print "Dataset %s not in local machine. Downloading..." % dataset
                downloadDataset(dataset, url=args.url, verbose=args.verbose)
        sys.exit(0)

    if args.upload:
        # Upload datasets. This is the tricky part.
        for dataset in args.datasets:
            deleteFlag = '--delete' if args.delete else ''
            localFolder = join(os.environ['SCIPION_TESTS'], dataset)
            remoteUser = 'scipion'
            remoteServer = 'ramanujan.cnb.csic.es'
            remoteFolder = join('/home', 'twiki', 'public_html', 'files', 'scipion', 'data', 'tests')

            if not exists(localFolder):
                sys.exit("ERROR: local folder %s doesn't exist." % localFolder)

            print 'Warning: Uploading, please BE CAREFUL! This is potentially dangerous.'
            print ('You are going to be connected to server "%s" as user "%s" to write '
                   'in folder "%s" the dataset "%s".' %
                   (remoteServer, remoteUser, remoteFolder, dataset))
            while True:
                ans = raw_input('Continue? (y/n): ')
                if ans in ['N', 'n']:
                    sys.exit(0)
                elif ans in ['Y', 'y']:
                    break

            # First make sure we have our MANIFEST file up-to-date
            createMANIFEST(join(os.environ['SCIPION_TESTS'], dataset))

            # Upload the dataset files (with rsync)
            Cmd('rsync -av %s %s %s@%s:%s' % (deleteFlag, localFolder, remoteUser, remoteServer, remoteFolder))

            # Regenerate remote MANIFEST (which contains a list of datasets)
            print "Regenerating remote MANIFEST file..."
            Cmd('ssh %s@%s "cd %s && find -type d -mindepth 1 -maxdepth 1 > MANIFEST"' %
                (remoteUser, remoteServer, remoteFolder))
            # This is a file that just contains the name of the directories
            # in remoteFolder. Nothing to do with the MANIFEST files in
            # the datasets, which contain file names and md5s.

            # Now do stuff related to last_m.txt - useful for buildbot
            lastmFile = join("Scipion", 'last_m.txt')
            localHostname = " ".join(os.uname())
            localUser = getpass.getuser()

            print "Registering modification attempt in last_m.txt file"
            Cmd("ssh " + remoteUser + '@' + remoteServer + " \"echo '++++' >> " +
                lastmFile + " && echo 'Modification to " + dataset +
                " dataset made at' >> " + lastmFile + " && date >> " +
                lastmFile + " && echo 'by " + localUser + " at " + localHostname +
                "' >> " + lastmFile + " && echo '----' >> " + lastmFile + "\"")
            print "...done."
            # Leave last_m.txt file indicating modifications have been done, to let buildbot trigger its automatic testing system
            #TODO: this step
        sys.exit(0)

    parser.print_usage()


def createMANIFEST(path):
    """ Create a MANIFEST file in path with the md5 of all files below """

    with open(join(path, 'MANIFEST'), 'w') as manifest:
        for root, dirs, files in os.walk(path):
            for filename in files:
                if filename != 'MANIFEST':  # do not include ourselves
                    fn = join(root, filename)  # file to check
                    manifest.write('%s %s\n' % (relpath(fn, path), md5sum(fn)))


def listDatasets(url):
    """ Print a list of local and remote datasets """

    tdir = os.environ['SCIPION_TESTS']
    print "Local datasets in %s" % tdir
    for folder in sorted(os.listdir(tdir)):
        if isdir(join(tdir, folder)):
            if exists(join(tdir, folder, 'MANIFEST')):
                print "  * " + folder
            else:
                print "  * " + folder + ' (not in dataset format)'

    try:
        print "\nRemote datasets in %s" % url
        data = urlopen(url + '/MANIFEST').read()
        for line in sorted(data.splitlines()):
            print "  * " + line.strip('./')
    except Exception as e:
        print 'Error reading %s (%s)' % (url, e)


def queryModifications(last_mfile, log_mfile):
    """ Show modifications stored in last_mfile, and update log_mfile if there are """

    print "Querying the modifications log file (%s)..." % last_mfile
    if exists(last_mfile):
        print "File %s exists. Checking its content..." % last_mfile
        if os.stat(last_mfile).st_size != 0: #File contains data
            print "File is not empty. Copying the content to log file " + log_mfile
            modif_file = open(log_mfile, 'a')
            file_content = open(last_mfile).read()
            modif_file.write(file_content)
            print "Last modifications file shows following content:"
            print file_content
            # TODO: In case we want to add the responsible of the modifications
            # to the blame list, this is the place to do it
            modif_file.close()
            open(last_mfile, 'w').close()  # wipe out contents
        else: #File's empty
            print "File is empty, so no modifications since last check."
            sys.exit(1)
            # FIXME: document this strange exit code - something to do with buildbot??
    else:
        print "File %s doesn't exist. Creating it..." % last_mfile
        try:
            os.makedirs(os.path.dirname(last_mfile))
            open(last_mfile, 'w').close()  # "touch"
            sys.exit(2)
            # We return with 2, to let Buildbot know that no modification was made
            # (when failure there was modification)
        except OSError as e:
            print 'Error creating %s: %s' % (last_mfile, e)



if __name__ == "__main__":
    main()
