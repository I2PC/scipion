#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
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

from __future__ import division

import os
import subprocess
import getopt
import sys
import argparse
import hashlib
import types
from urllib2 import urlopen

from pyworkflow.utils.path import moveFile, makePath, makeFilePath, copyFile, cleanPath


def scipion_logo():
    return """

QQQQQQQQQT!^'::\"\"?$QQQQQQ  S   S   S
QQQQQQQY`          ]4QQQQ  C   C   C
QQQQQD'              \"$QQ  I   I   I
QQQQP                 \"4Q  P   P   P
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
[    ]Qm    ]QW    \"QQQQQ
h     $Qc   jQQ     ]$QQQ
Q,  :aQQf qyQQQL    _yQQQ
QL jmQQQgmQQQQQQmaaaQWQQQ
"""


def backupFile(basePath, fileName='MANIFEST'):
    filePath = os.path.join(basePath, fileName)
    backPath = os.path.join(basePath, "."+fileName+".backup")
    if os.path.exists(filePath):
        moveFile(filePath, backPath)

        
def restoreFile(basePath, fileName='MANIFEST'):
    filePath = os.path.join(basePath, fileName)
    backPath = os.path.join(basePath, "."+fileName+".backup")
    if os.path.exists(backPath):
        if os.path.exists(filePath):
            os.remove(filePath)
        moveFile(backPath, filePath)

            
def delBackFile(basePath, fileName='MANIFEST'):
    filePath = os.path.join(basePath, fileName)
    backPath = os.path.join(basePath, "."+fileName+".backup")
    if os.path.exists(backPath):
        os.remove(backPath)

        
def initDownloadBar(len):
    sys.stdout.write("[%s]" % (" " * len))
    sys.stdout.flush()
    sys.stdout.write("\b" * (len+1))

    
def backDownloadBar(number):
    sys.stdout.write("\b" * (number))
    sys.stdout.flush()


def downloadFile(datasetName, fname,
                 workingCopy=os.environ['SCIPION_TESTS'],
                 tmpMd5copy=os.environ['SCIPION_TMP'],
                 askMsg="download it?",
                 url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
                 verbose=False):
    fileDownloaded = 0
    datasetFolder = os.path.join(workingCopy, datasetName)
    datasetFolderTmp = os.path.join(tmpMd5copy, datasetName)
    answer = "-"
    if verbose:
        while answer != "y" and answer != "n" and answer != "":
            answer = raw_input("\t "+askMsg + " ([y]/n): ")
    if answer == "y" or answer == "" or not verbose:
        if verbose:
            print "\t Downloading file..."
        makeFilePath(os.path.join(datasetFolder, fname))
        try:
            data = urlopen('%s/%s/%s' % (url, datasetName, fname)).read()
            open(os.path.join(datasetFolder, fname), 'w').write(data)
            if verbose:
                print "\t ...done."
        except:
            print "\t "+ fname+" ...ERROR"
            print "URL: "+url+'/'+datasetName+fname
            print "destination: "+os.path.join(datasetFolder, fname)
        fileDownloaded = 1
    return fileDownloaded


def downloadDataset(datasetName, destination=os.environ['SCIPION_TESTS'], url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests", verbose=False, onlyManifest=False):
    datasetFolder = os.path.join(destination, "%(dataset)s" % ({'dataset': datasetName}))
    makePath(datasetFolder)
    manifest = os.path.join(destination, datasetName, 'MANIFEST')
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
            makeFilePath(os.path.join(datasetFolder, line))
            try:
                data = urlopen('%s/%s/%s' % (url, datasetName, line)).read()
                open(os.path.join(datasetFolder, line), 'w').write(data)
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
                print "destination: "+os.path.join(datasetFolder, line)
                while True:
                    answer = raw_input("continue downloading? (y/[n]): ")
                    if not answer or answer.lower() == 'n':
                        sys.exit(1)
                    if answer in 'Yy':
                        break

            md5sum = md5Sum(os.path.join(datasetFolder, os.path.normpath(line.replace("\n","").split(" ")[0])))
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


def checkForUpdates(datasetName, workingCopy=None, tmpMd5copy=None,
                    url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
                    verbose=False):
    # Get default values for variables if we got none
    workingCopy = workingCopy or os.environ['SCIPION_TESTS']
    tmpMd5copy = tmpMd5copy or os.environ['SCIPION_TMP']

    # We need to download the remote manifest file
    datasetFolderTmp = os.path.join(tmpMd5copy, datasetName)
    manifest = os.path.join(os.environ['SCIPION_TMP'], datasetName, 'MANIFEST')
    manifestFileTmp = open(manifest, 'r+')
    manifestLinesTmp = manifestFileTmp.readlines()

    # and check it with the local copy
    datasetFolder = os.path.join(workingCopy, datasetName)
    manifestFile = open(manifest)
    manifestLines = manifestFile.readlines()
    
    filesUpdated = 0
    print "Verifying MD5..."
    for number, line in enumerate(manifestLinesTmp):
        fname = os.path.normpath(line.replace("\n","").split(" ")[0])
        if verbose:
            print '\t%s' % fname,
        if os.path.exists(os.path.join(datasetFolder, fname)):
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
                filesUpdated += downloadFile(datasetName, fname, workingCopy, tmpMd5copy, askMsg="update it?", url=url, verbose=verbose)
        else: #file does not exist, show option for downloading it
            if verbose:
                print "\n\t file %s doesn't exist." % fname
            filesUpdated += downloadFile(datasetName, fname, workingCopy, tmpMd5copy, askMsg="download it?", url=url, verbose=verbose)
    copyFile(os.path.join(datasetFolderTmp, 'MANIFEST'), os.path.join(datasetFolder, 'MANIFEST'))
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
    parser = argparse.ArgumentParser(description="Sync Data Python Bash Executor (WRAPPER)".center(5, '*'))
    exclusive = parser.add_mutually_exclusive_group()

    parser.add_argument('-s', '--syncfile', 
        default="sync_data", 
        help="Sync Bash Script")

    parser.add_argument('-t', '--last-mod-file', 
        default="last_m.txt", 
        help="File that contain the IP/s of the computer/s that acceeded last time to change the Scipion tests data.")

    parser.add_argument('-m', '--mod-log-file', 
        default="modifications.log",
        help="File that contain the whole modifications log to keep a tracking of what has been done in the Scipion tests data. The path given must be relative to $SCIPION_HOME")
    
    parser.add_argument('-u', '--url',
        default="http://scipionwiki.cnb.csic.es/files/scipion/data/tests",
        help="String storing the url where remote datasets will be looked for")

    exclusive.add_argument('-q', '--query-for-modifications', 
        action="store_true",
        help="Look for last_m.txt file. There, there will be or (1) the IP of the last modifier and date or (2) nothing. If the script finds (1) it moves it to modifications.log file and returns 0. If it finds (2) it check whether the file was modified in the last 30 minutes. If not, it returns 1, if yes it returns 2.")

    parser.add_argument('-d', '--delete', 
        action="store_true",
        help="Only valid when -r option enabled. This deletes remote files not existing in local. In the nutshell, it leaves the remote scipion data directory as it is the local one. Extremely dangerous option!.")

    exclusive.add_argument('-r', '--reverse-sync', 
        action='store_true',
        help="Synchronize from the local data to scipion machine. When wildcard 'all' is given, it will synchronize all the local folder with the remote one. When a set of locations is given, they're synchronized one by one against the remote scipion server. File path must be given from $SCIPION_HOME/tests folder")
    
    parser.add_argument('-f', '--dataset',
        nargs='+',
        help="Determine the dataset to use. The selected operation will be applied to this dataset.")
    
    parser.add_argument('-l', '--list',
        action="store_true",
        help="Look for local datasets in $SCIPION_TESTS folder and for remote datasets in http://scipionwiki.cnb.csic.es/files/tests.")
    
    parser.add_argument('-v', '--verbose',
        action="store_true",
        help="Verbose mode. This will print more detailed messages")

    args = parser.parse_args()


    # Depending on the arguments selected, doing one thing or another

    #print scipion_logo()
    if args.list:
        print "List of local datasets in %(datasetsFolder)s" % ({'datasetsFolder': os.environ['SCIPION_TESTS']})
        folders = os.listdir(os.environ['SCIPION_TESTS'])
        for folder in folders:
            if os.path.isdir(os.path.join(os.environ['SCIPION_TESTS'], folder)):
                if os.path.exists(os.path.join(os.environ['SCIPION_TESTS'], folder, 'MANIFEST')):
                    print "\t * " + folder
                else:
                    print "\t * " + folder + ' (not in dataset format)'
        print
        print "updating remote info..."
        manifest = os.path.join(os.environ['SCIPION_TMP'], 'MANIFEST')
        backupFile(manifest)
        try:
            data = urlopen(args.url + '/MANIFEST').read()
            open(manifest, 'w').write(data)
        except:
            print "MANIFEST"+" ...ERROR"
            print "URL: "+args.url+'/MANIFEST'
            print "destination: "+manifest
            restoreFile(manifest)
        if os.path.exists(manifest):
            print "List of remote datasets in %(urlAddress)s" % ({'urlAddress': args.url})
            for line in open(manifest):
                print "\t * "+os.path.split(line.replace("\n",""))[1]
        delBackFile(manifest)

    elif args.query_for_modifications:
        print "Querying the modifications log file..."
        if os.path.exists(os.path.join(os.environ['SCIPION_HOME'],args.last_mod_file)):
            print "File " + args.last_mod_file + " exists. Checking its content..."
            if os.stat(os.path.join(os.environ['SCIPION_HOME'],args.last_mod_file)).st_size != 0: #File contains data
                print "File's not empty. Copying the content to log file " + args.mod_log_file
                last_file = open(os.path.join(os.environ['SCIPION_HOME'],args.last_mod_file), 'r+')
                modif_file = open(os.path.join(os.environ['SCIPION_HOME'],args.mod_log_file), 'a')
                file_content = last_file.read()
                modif_file.write(file_content)
                print "Last modifications file shows following content:"
                print file_content
                #TODO: In case we want to add the responsible of the modifications to the blame list, this is the place to do it
                last_file.close()
                modif_file.close()
                last_file = open(os.path.join(os.environ['SCIPION_HOME'],args.last_mod_file), 'w').close()
            else: #File's empty
                print "File's empty, so no modification since last check."
                sys.exit(1)
        else:
            print "File " + args.last_mod_file + " doesn't exist. Creating it..."
            open(os.path.join(os.environ['SCIPION_HOME'],args.last_mod_file), 'w').close()
            sys.exit(2) #We return with 2, to let Buildbot know that no modification was made (when failure there was modification)
    elif args.dataset:
        if len(args.dataset) != 1:
            print "Selected datasets: %(datasets)s" % ({'datasets': " ".join(args.dataset)})
        else:
            print "Selected datasets: %(datasets)s" % ({'datasets': args.dataset[0]})
        ans = ""
        if len(args.dataset) == 1 and " ".join(args.dataset) == "all" and args.reverse_sync:
            while ans != "y" and ans != "Y" and ans != "n" and ans != "N":
                ans = raw_input("You've selected to synchronize all your data with the remote data. You will squash all remote data. Are you sure you want to continue? (y/n): ")
            if ans != ("y" or "Y"):
                sys.exit(1)
        for dataset in args.dataset:
            if args.reverse_sync:
                deleteFlag=""
                if args.delete:
                    deleteFlag="--delete "
                localFolder = os.path.join(os.environ['SCIPION_TESTS'], dataset)
                remoteUser = 'scipion'
                remoteServer = 'ramanujan.cnb.csic.es'
                remoteFolder = os.path.join('/home', 'twiki', 'public_html', 'files', 'scipion', 'data', 'tests')
                lastmFile = os.path.join("Scipion", 'last_m.txt')
                localUser = None
                localHostname = " ".join(os.uname())
                try:
                    import pwd
                except ImportError:
                    import getpass
                    pwd = None
                if pwd:
                    localUser = pwd.getpwuid(os.geteuid()).pw_name
                else:
                    localUser = getpass.getuser()
                
                if not os.path.exists(localFolder):
                    print "ERROR: local folder %(localFolder)s doesn't exist." % ({'localFolder': localFolder})
                    sys.exit(1)
                print "Reverse synchronizing, BE CAREFUL!!! OPERATION EXTREMELY DANGEROUS!!!"
                ans = ''
                while ans not in 'YyNn':
                    ans = raw_input("You're going to be connected to %(remoteServer)s server as %(remoteUser)s user to write in %(remoteFolder)s folder for %(dataset)s dataset. Only authorized users shall pass. Continue? (y/n): " % ({'remoteServer': remoteServer, 'remoteUser': remoteUser, 'remoteFolder': remoteFolder, 'dataset': dataset}))
                if ans not in 'Yy':
                    sys.exit(1)
                
                # If the folder is not in the proper format, create the format and then upload
                command = '%(scipionPython)s scripts/generate_md5.py %(dataset)s %(datasetPath)s' % ({'scipionPython': os.environ['SCIPION_PYTHON'], 'dataset': dataset, 'datasetPath': os.environ['SCIPION_TESTS']})
                Cmd(command)
                # Synchronize files
                command = 'rsync -av '+ deleteFlag + localFolder + ' ' + remoteUser + '@' + remoteServer + ':' + remoteFolder
                Cmd(command)
                # Regenerate remote MANIFEST
                print "Regenerating remote MANIFEST file..."
                command = "ssh " + remoteUser + '@' + remoteServer + " \"cd " + remoteFolder + ' && ' + "find -maxdepth 1 -type d -type d ! -iname '.' > MANIFEST" + "\""
                Cmd(command)
                
                print "Registering modification attempt in last_m.txt file"
                command = "ssh " + remoteUser + '@' + remoteServer + " \"echo '++++' >> " + lastmFile + " && echo 'Modification to " + dataset + " dataset made at' >> " + lastmFile + " && date >> " + lastmFile + " && echo 'by " + localUser + " at " + localHostname +"' >> " + lastmFile + " && echo '----' >> " + lastmFile + "\""
                Cmd(command)
                print "...done."
                # Leave last_m.txt file indicating modifications have been done, to let buildbot trigger its automatic testing system
                #TODO: this step
            else:
                datasetFolder = os.path.join(os.environ['SCIPION_TESTS'],"%(dataset)s" % ({'dataset': dataset}))
                if os.path.exists(datasetFolder):
                    print "Dataset %s working copy detected. Checking checksum for updates..." % dataset
                    downloadDataset(dataset, destination=os.environ['SCIPION_TMP'], url=args.url, verbose=args.verbose, onlyManifest=True)
                    checkForUpdates(dataset, url=args.url, verbose=args.verbose)
                    cleanPath(os.path.join(os.environ['SCIPION_TMP'], "%(dataset)s" % ({'dataset': dataset})))
                else:
                    print "dataset %(dataset)s not in local machine, trying to download..." % ({'dataset': dataset})
                    downloadDataset(dataset, url=args.url, verbose=args.verbose)
                
        

if __name__ == "__main__":
    main()
