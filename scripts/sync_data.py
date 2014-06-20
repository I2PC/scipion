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
from os.path import join, isdir, exists
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


        
def initDownloadBar(len):
    sys.stdout.write("[%s]" % (" " * len))
    sys.stdout.flush()
    sys.stdout.write("\b" * (len+1))

    
def backDownloadBar(number):
    sys.stdout.write("\b" * (number))
    sys.stdout.flush()


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
    destination = destination if destination else os.environ['SCIPION_TESTS']

    datasetFolder = join(destination, datasetName)
    makePath(datasetFolder)
    manifest = join(destination, datasetName, 'MANIFEST')
    try:
        if verbose:
            print "retrieving MANIFEST file"
        data = urlopen('%s/%s/MANIFEST' % (url, datasetName)).read()
        open(manifest, 'w').write(data)
    except Exception as e:
        if verbose:
            print "URL: " + url+'/'+datasetName+'/MANIFEST'
            print "destination:" + manifest
        raise Exception("URL %s could not be retrieved (%s)" % (url, e))
    if not onlyManifest:
        manifestFile = open(manifest, 'r+')
        manifestLines = manifestFile.readlines()
        print "Fetching dataset %s files..." % datasetName
        totalNumber = len(manifestLines)
        percent = prevPercent = 0
        downloadBarWidth = 100
        if not verbose:
            initDownloadBar(downloadBarWidth)
        for number, lineExt in enumerate(manifestLines):
            line = os.path.normpath(lineExt.replace("\n","").split(" ")[0])
            md5InLine = lineExt.replace("\n","").split(" ")[1] 
            makeFilePath(join(datasetFolder, line))
            try:
                data = urlopen('%s/%s/%s' % (url, datasetName, line)).read()
                open(join(datasetFolder, line), 'w').write(data)
                percent = ((number+1)/(totalNumber*1.0))*100
                if verbose:
                    print "\t%s ...OK (%02d %%) " % (line, percent)
                else:
                    progress = percent - prevPercent
                    remaining = downloadBarWidth-progress
                    sys.stdout.write("#" * int(1 + progress))
                    sys.stdout.flush()
                prevPercent = percent
            except:
                print "\t "+ line+" ...ERROR"
                print "URL: "+url+'/'+datasetName+line
                print "destination: "+join(datasetFolder, line)
                while True:
                    answer = raw_input("continue downloading? (y/[n]): ")
                    if not answer or answer.lower() == 'n':
                        sys.exit(1)
                    if answer in 'Yy':
                        break

            md5sum = md5Sum(join(datasetFolder, os.path.normpath(line.replace("\n","").split(" ")[0])))
                #data = fileToCheck.read()
                #md5sum = hashlib.md5(data).hexdigest()
            if verbose:
                print "\t\tmd5 verification...%s %s" % (md5sum, md5InLine)
            if md5sum != md5InLine:
                print "ERROR in md5 verification for file %s" % line
                print "md5 sum calculated for downloaded file is %s" % md5sum
                print "md5 sum that file should have is %s" % md5InLine
                print
                while True:
                    answer = raw_input("continue downloading? (y/[n]): ")
                    if not answer or answer.lower() == "n":
                        sys.exit(1)
                    if answer in 'Yy':
                        break
            elif verbose:
                print "md5 verification...OK"
                print
        if not verbose:
            print
    print "\t...done"
    print


def md5Sum(file):
    md5sum = 0
    md5 = hashlib.md5()
    with open(file,'r+') as fileToCheck:
        for chunk in iter(lambda: fileToCheck.read(128*md5.block_size), b''):
            md5.update(chunk)
    md5sum = md5.hexdigest()
    return md5sum


def checkForUpdates(datasetName, workingCopy=None,
                    url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
                    verbose=False):
    # Get default values for variables if we got none
    workingCopy = workingCopy or os.environ['SCIPION_TESTS']

    # We need to download the remote manifest file
    datasetFolderTmp = join(os.environ['SCIPION_TMP'], datasetName)
    manifest = join(os.environ['SCIPION_TMP'], datasetName, 'MANIFEST')
    manifestFileTmp = open(manifest, 'r+')
    manifestLinesTmp = manifestFileTmp.readlines()

    # and check it with the local copy
    datasetFolder = join(workingCopy, datasetName)
    manifestFile = open(manifest)
    manifestLines = manifestFile.readlines()
    
    filesUpdated = 0
    print "Verifying MD5..."
    for number, line in enumerate(manifestLinesTmp):
        fname = os.path.normpath(line.replace("\n","").split(" ")[0])
        if verbose:
            print '\t%s' % fname,
        if exists(join(datasetFolder, fname)):
            lineInManifest = findLineInFile(line, manifest)
            if lineInManifest == -1:
                if verbose:
                    print "\n\tERROR: file %s not found in MANIFEST" % fname
                sys.exit(1)
            md5fcalc = (manifestLines[lineInManifest]).split(" ")[1]
            md5Tmpfcalc = line.split(" ")[1]
            if md5fcalc == md5Tmpfcalc:
                if verbose:
                    print "\r\tOK  %s" % fname
            else:
                if verbose:
                    print "\r\tBAD %s  \t==> checksum differs" % fname
                filesUpdated += downloadFile(datasetName, fname, workingCopy, askMsg="update it?", url=url, verbose=verbose)
        else: #file does not exist, show option for downloading it
            if verbose:
                print "\n\t file %s doesn't exist." % fname
            filesUpdated += downloadFile(datasetName, fname, workingCopy, askMsg="download it?", url=url, verbose=verbose)
    copyFile(join(datasetFolderTmp, 'MANIFEST'), join(datasetFolder, 'MANIFEST'))
    if filesUpdated == 0:
        print "\t ...done. Nothing changed."
    elif filesUpdated == 1:
        print "\t ...done. 1 file was updated."
    else:
        print "\t ...done. %d files were updated." % filesUpdated
    print


def findLineInFile(text, filePath):
    "Return the line number where text first appears in a file, or -1 if nowhere"
    for i, line in enumerate(open(filePath)):
        if text in line:
            return i  # yeah, found it!
    return -1  # oh, we didn't find it :(


def Cmd(command):
    print ">>>>>", command
    subprocess.call(command, shell=True)


def main():
    # Arguments parsing
    parser = argparse.ArgumentParser(description=__doc__)
    add = parser.add_argument  # shortcut
    add('--datasets', nargs='+', help='Datasets to use.')
    add('-l', '--list', action='store_true',
        help=('Look for local datasets in $SCIPION_TESTS and for remote '
              'datasets in http://scipionwiki.cnb.csic.es/files/tests.'))
    add('-s', '--syncfile', default='sync_data', help='Sync Bash Script.')
    add('--last-mod-file', default='last_m.txt', 
        help=('File which lists the computers that have last changed the '
              'Scipion tests data (using this script).'))
    add('--mod-log-file', default='modifications.log',
        help=('File with the whole modifications log to keep track of what '
              'has been done in the Scipion tests data. The path must be '
              'relative to the Scipion folder.'))
    add('-u', '--url',
        default='http://scipionwiki.cnb.csic.es/files/scipion/data/tests',
        help='URL where remote datasets will be looked for.')
    add('--delete', action='store_true',
        help=('Delete remote files not present in local. It leaves the remote '
              'scipion data directory as it is in the local one. Extremely '
              'dangerous! Only valid with -r.'))
    add('-q', '--query-for-modifications', action='store_true',
        help=("Look for last_m.txt file. There, there will be (1) the IP of "
              "the last modifier and date or (2) nothing. If the script finds "
              "(1) it moves it to modifications.log and returns 0. If it "
              "finds (2) it checks whether the file was modified in the last "
              "30 minutes. If not, it returns 1, if yes it returns 2."))
    add('-r', '--reverse-sync', action='store_true',
        help=("Synchronize from the local data to the scipion machine. When "
              "the wildcard 'all' is given, it will synchronize all the local "
              "folder with the remote one. When a set of locations is given, "
              "they're synchronized one by one against the remote scipion "
              "server. File path must be given from the SCIPION_TESTS folder."))
    add('-v', '--verbose', action='store_true', help='Print more details.')
    args = parser.parse_args()

    #print scipion_logo

    # Dispatch the easy cases first (list and query), and then take care of
    # the more complex ones.
    if args.list:
        listDatasets(args.url)
        sys.exit(0)
    elif args.query_for_modifications:
        last_mfile = join(os.environ['HOME'], 'Scipion', args.last_mod_file)
        log_mfile = join(os.environ['HOME'], 'Scipion', args.mod_log_file)
        queryModifications(last_mfile, log_mfile)
        sys.exit(0)

    if not args.datasets:
        sys.exit("At least -l, -q or --datasets needed. See --help for more info.")

    print "Selected datasets: %s" % " ".join(args.datasets)

    if args.datasets == ["all"] and args.reverse_sync:
        print "You've selected to synchronize all your data with the remote data."
        print "You will squash all remote data."
        while True:
            ans = raw_input("Are you sure you want to continue? (y/n): ")
            if ans in ['N', 'n']:
                sys.exit(1)  # FIXME: why return with 1? comment please
            elif ans in ['Y', 'y']:
                break
        # FIXME: and then what? what happens next if we actually selected "all"?

    if not args.reverse_sync:
        # Download datasets.
        for dataset in args.datasets:
            if exists(join(os.environ['SCIPION_TESTS'], dataset)):
                print "Dataset %s working copy detected. Checking checksum for updates..." % dataset
                downloadDataset(dataset, destination=os.environ['SCIPION_TMP'], url=args.url, verbose=args.verbose, onlyManifest=True)
                checkForUpdates(dataset, url=args.url, verbose=args.verbose)
                cleanPath(join(os.environ['SCIPION_TMP'], dataset))
            else:
                print "Dataset %s not in local machine, trying to download..." % dataset
                downloadDataset(dataset, url=args.url, verbose=args.verbose)
    else:
        # Upload datasets. This is the tricky part.
        for dataset in args.datasets:
            deleteFlag = '--delete' if args.delete else ''
            localFolder = join(os.environ['SCIPION_TESTS'], dataset)
            remoteUser = 'scipion'
            remoteServer = 'ramanujan.cnb.csic.es'
            remoteFolder = join('/home', 'twiki', 'public_html', 'files', 'scipion', 'data', 'tests')
            lastmFile = join("Scipion", 'last_m.txt')
            localHostname = " ".join(os.uname())
            localUser = getpass.getuser()

            if not exists(localFolder):
                sys.exit("ERROR: local folder %s doesn't exist." % localFolder)

            print "Reverse synchronizing, BE CAREFUL!!! OPERATION EXTREMELY DANGEROUS!!!"
            print ("You're going to be connected to %s server as %s user to write "
                   "in %s folder for %s dataset. Only authorized users shall pass." %
                   (remoteServer, remoteUser, remoteFolder, dataset))
            while True:
                ans = raw_input('Continue? (y/n): ')
                if ans in ['N', 'n']:
                    sys.exit(1)  # FIXME: why exit with 1? why is it an error? please explain
                elif ans in ['Y', 'y']:
                    break
            # If the folder is not in the proper format, create the format and then upload
            Cmd('%s scripts/generate_md5.py %s %s' % (os.environ['SCIPION_PYTHON'], dataset, os.environ['SCIPION_TESTS']))
            # Synchronize files
            Cmd('rsync -av %s %s %s@%s:%s' % (deleteFlag, localFolder, remoteUser, remoteServer, remoteFolder))
            # Regenerate remote MANIFEST
            print "Regenerating remote MANIFEST file..."
            Cmd('ssh %s@%s "cd %s && find -maxdepth 1 -type d -type d ! -iname \'.\' > MANIFEST"' % (remoteUser, remoteServer, remoteFolder))

            print "Registering modification attempt in last_m.txt file"
            Cmd("ssh " + remoteUser + '@' + remoteServer + " \"echo '++++' >> " + lastmFile + " && echo 'Modification to " + dataset + " dataset made at' >> " + lastmFile + " && date >> " + lastmFile + " && echo 'by " + localUser + " at " + localHostname +"' >> " + lastmFile + " && echo '----' >> " + lastmFile + "\"")
            print "...done."
            # Leave last_m.txt file indicating modifications have been done, to let buildbot trigger its automatic testing system
            #TODO: this step


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
        data = urlopen(url + '/MANIFEST').read()
        print "\nRemote datasets in %s" % url
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
