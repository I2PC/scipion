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
import urllib
import hashlib
from pyworkflow.utils.path import moveFile, makePath, makeFilePath, copyFile, cleanPath

def scipion_logo():
    print ""
    print "QQQQQQQQQT!^'::\"\"?$QQQQQQ" + "  S   S   S"
    print "QQQQQQQY`          ]4QQQQ"   + "  C   C   C"
    print "QQQQQD'              \"$QQ"  + "  I   I   I"
    print "QQQQP                 \"4Q"  + "  P   P   P"
    print "QQQP        :.,        -$"   + "  I   I   I"
    print "QQD       awQQQQwp      )"   + "  O   O   O"
    print "QQ'     qmQQQWQQQQg,   jm"   + "  N   N   N"
    print "Qf     QQQD^   -?$QQp jQQ"   + " ################################################"
    print "Q`    qQQ!        4WQmQQQ"   + " # Integrating image processing packages for EM #"
    print "F     QQ[          ~)WQQQ"   + " ################################################"
    print "[    ]QP             4WQQ"   + ""
    print "f    dQ(             -$QQ"   + " Data synchronization script"
    print "'    QQ              qQQQ"
    print ".   )QW            _jQQQQ"
    print "-   =QQ           jmQQQQQ"
    print "/   -QQ           QQQQQQQ"
    print "f    4Qr    jQk   )WQQQQQ"
    print "[    ]Qm    ]QW    \"QQQQQ"
    print "h     $Qc   jQQ     ]$QQQ"
    print "Q,  :aQQf qyQQQL    _yQQQ"
    print "QL jmQQQgmQQQQQQmaaaQWQQQ"
    print ""

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
    
def forwardDownloadbar(percent, progress):
    for x in range(0, (int)(progress//1)):
        sys.stdout.write("#")
        sys.stdout.flush()
        progress-=1
    if progress != 0:
        if ((percent-progress)//1) != ((percent)//1):
            sys.stdout.write("#")
            sys.stdout.flush()
    
def downloadFile(datasetName, file, fileMd5, workingCopy=os.environ['SCIPION_TESTS'], tmpMd5copy=os.environ['SCIPION_TMP'], askMsg="download it?", url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests/NACHOTESTS_PLEASEDONOTREMOVE/tests", verbose=False):
    fileDownloaded = 0
    datasetFolder = os.path.join(workingCopy, "dataset_%(dataset)s" % ({'dataset': datasetName}))
    datasetFolderTmp = os.path.join(tmpMd5copy, "dataset_%(dataset)s" % ({'dataset': datasetName}))
    answer = "-"
    if verbose:
        while answer != "y" and answer != "n" and answer != "":
            answer = raw_input("\t "+askMsg + " ([y]/n): ")
    if answer == "y" or answer == "" or not verbose:
        if verbose:
            print "\t Downloading file..."
        makeFilePath(os.path.join(datasetFolder, file))
        try:
            urllib.urlretrieve(url+'/dataset_'+datasetName+'/'+file, os.path.join(datasetFolder, file))
            if verbose:
                print "\t ...done."
        except:
            print "\t "+ file+" ...ERROR"
            print "URL: "+url+'/dataset_'+datasetName+file
            print "destination: "+os.path.join(datasetFolder, file)
        if verbose:
            print "\t Copying md5..."
        makeFilePath(os.path.join(datasetFolder, fileMd5))
        copyFile(os.path.join(datasetFolderTmp, fileMd5), os.path.join(datasetFolder, fileMd5))
        if verbose:
            print "\t ...done."
        fileDownloaded = 1
    return fileDownloaded

def downloadDataset(datasetName, destination=os.environ['SCIPION_TESTS'], url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests/NACHOTESTS_PLEASEDONOTREMOVE/tests", verbose=False, onlyMd5=False):
    datasetFolder = os.path.join(destination, "dataset_%(dataset)s" % ({'dataset': datasetName}))
    makePath(datasetFolder)
    try:
        if verbose:
            print "retreiving MANIFEST file"
        urllib.urlretrieve(url+'/dataset_'+datasetName+'/MANIFEST', getManifestPath(basePath=destination, dataset=datasetName))
    except:
        print "URL could not be retreived"
        if verbose:
            print "URL: " + url+'/dataset_'+datasetName+'/MANIFEST'
            print "destination:" + getManifestPath(basePath=destination, dataset=datasetName)
    manifestFile = open(getManifestPath(basePath=destination, dataset=datasetName), 'r+')
    manifestLines = manifestFile.readlines()
    print "Fetching dataset %(dataset)s files..." % ({'dataset': datasetName})
    totalNumber = len(manifestLines)*(2-1.0*onlyMd5)
    percent = prevPercent = 0
    downloadBarWidth = 100
    if not verbose:
        initDownloadBar(downloadBarWidth)
    for number, line in enumerate(manifestLines):
        fileOptions = []
        if not onlyMd5:
            fileOptions = [ os.path.normpath(line.replace("\n","")), os.path.normpath(line.replace("\n","") + ".md5") ]
        else:
            fileOptions = [ os.path.normpath(line.replace("\n","") + ".md5") ]
        for indx, fileOption in enumerate(fileOptions):
            file = fileOption
            makeFilePath(os.path.join(datasetFolder, file))
            try:
                urllib.urlretrieve(url+'/dataset_'+datasetName+'/'+file, os.path.join(datasetFolder, file))
                percent = ((number*(2-1.0*onlyMd5)+indx+1)/(totalNumber*1.0))*100
                if verbose:
                    print "\t "+file+" ...OK ( %(percent)02d " % ({'percent': percent}) + ' %)'
                else:
                    progress = percent-prevPercent
                    remaining = downloadBarWidth-progress
                    forwardDownloadbar(percent, progress)
                prevPercent = percent
            except:
                print "\t "+ file+" ...ERROR"
                print "URL: "+url+'/dataset_'+datasetName+file
                print "destination: "+os.path.join(datasetFolder, file)
                answer = "-"
                while answer != "y" and answer != "n" and answer != "":
                    answer = raw_input("continue downloading? (y/[n]): ")
                    if answer == "n" or answer == "":
                        sys.exit(1)
        if not onlyMd5:
            md5sum = 0
            md5 = hashlib.md5()
            with open(os.path.join(datasetFolder, os.path.normpath(line.replace("\n",""))),'r+') as fileToCheck:
                for chunk in iter(lambda: fileToCheck.read(128*md5.block_size), b''):
                    md5.update(chunk)
            md5sum = md5.hexdigest()
                #data = fileToCheck.read()
                #md5sum = hashlib.md5(data).hexdigest()
            md5file = open(os.path.join(datasetFolder, os.path.normpath(line.replace("\n","")+".md5")),'r+')
            md5calc = md5file.readlines()[0].split(" ")[0]
            if verbose:
                print "\t \t md5 verification...%(md5sum)s %(md5calc)s" % ({'md5sum': md5sum, 'md5calc': md5calc})
            if md5sum != md5calc:
                print "ERROR in md5 verification"
                answer = "-"
                while answer != "y" and answer != "n" and answer != "":
                    answer = raw_input("continue downloading? (y/[n]): ")
                    if answer == "n" or answer == "":
                        sys.exit(1)
            elif verbose:
                print "md5 verification...OK"
                print ""
    if not verbose:
        sys.stdout.write("\n")
    print "\t ...done"
    print ""

def checkForUpdates(datasetName, workingCopy=os.environ['SCIPION_TESTS'], tmpMd5copy=os.environ['SCIPION_TMP'], url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests/NACHOTESTS_PLEASEDONOTREMOVE/tests", verbose=False):
    
    # Check with local copy
    datasetFolder = os.path.join(workingCopy, "dataset_%(dataset)s" % ({'dataset': datasetName}))
    datasetFolderTmp = os.path.join(tmpMd5copy, "dataset_%(dataset)s" % ({'dataset': datasetName}))
    manifestFileTmp = open(getManifestPath(basePath=os.environ['SCIPION_TMP'], dataset=datasetName), 'r+')
    manifestLinesTmp = manifestFileTmp.readlines()
    filesUpdated = 0
    print "Verifying MD5..."
    for number, line in enumerate(manifestLinesTmp):
        file = os.path.normpath(line.replace("\n",""))
        fileMd5 = os.path.normpath(line.replace("\n","")+".md5")
        if verbose:
            print "\t checking file %(file)s..." % ({'file': file})
        if os.path.exists(os.path.join(datasetFolder, file)):
            if os.path.join(os.path.join(datasetFolder, fileMd5)):
                md5f = open(os.path.normpath(os.path.join(datasetFolder, fileMd5)))
                md5Tmpf = open(os.path.normpath(os.path.join(datasetFolderTmp, fileMd5)))
                md5fcalc = md5f.readlines()[0].split(" ")[0]
                md5Tmpfcalc = md5Tmpf.readlines()[0].split(" ")[0]
                if md5fcalc == md5Tmpfcalc:
                    if verbose:
                        print "\t md5 checked and OK!"
                else:
                    if verbose:
                        print ""
                        print "\t file %(file)s checksum differs." % ({'file': file})
                    filesUpdated += downloadFile(datasetName, file, fileMd5, workingCopy, tmpMd5copy, askMsg="update it?", url=url, verbose=verbose)
            else: #file exists, but its md5 file doesn't
                print "ERROR: file %(fileMd5)s doesn't exists. Dataset corrupt, please delete it and download it again." % ({'fileMd5': fileMd5})
                sys.exit(1)
        else: #file does not exist, show option for downloading it
            if verbose:
                print ""
                print "\t file %(file)s doesn't exist." % ({'file': file})
            filesUpdated += downloadFile(datasetName, file, fileMd5, workingCopy, tmpMd5copy, askMsg="download it?", url=url, verbose=verbose)
    copyFile(os.path.join(datasetFolderTmp, 'MANIFEST'), os.path.join(datasetFolder, 'MANIFEST'))
    if filesUpdated == 0:
        print "\t ...done. Nothing changed."
    elif filesUpdated == 1:
        print "\t ...done. 1 file was updated."
    else:
        print "\t ...done. %(filesUpdated)d files were updated." % ({'filesUpdated': filesUpdated})

def getManifestPath(basePath=os.environ['SCIPION_TESTS'], dataset=""):
    datasetName = "dataset_"+ dataset
    if dataset == "":
        datasetName = ""
    return os.path.join(basePath, datasetName, 'MANIFEST')

def main(argv):
    # Arguments parsing
    parser = argparse.ArgumentParser(description="*"*5 + "Sync Data Python Bash Executor (WRAPPER)" + "*"*5)
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
        default="http://scipionwiki.cnb.csic.es/files/scipion/data/tests/NACHOTESTS_PLEASEDONOTREMOVE/tests",
        help="String storing the url where remote datasets will be looked for")

    exclusive.add_argument('-q', '--query-for-modifications', 
        action="store_true",
        help="Look for last_m.txt file. There, there will be or (1) the IP of the last modifier and date or (2) nothing. If the script finds (1) it moves it to modifications.log file and returns 0. If it finds (2) it check whether the file was modified in the last 30 minutes. If not, it returns 1, if yes it returns 2.")

    parser.add_argument('-d', '--delete', 
        action="store_true",
        help="Only valid when -r option enabled. This deletes remote files not existing in local. In the nutshell, it leaves the remote scipion data directory as it is the local one. Extremely dangerous option!.")

    exclusive.add_argument('-r', '--reverse-sync', 
        nargs='+',
        help="Synchronize from the local data to scipion machine. When wildcard 'all' is given, it will synchronize all the local folder with the remote one. When a set of locations is given, they're synchronized one by one against the remote scipion server. File path must be given from $SCIPION_HOME/tests folder")
    
    parser.add_argument('-f', '--dataset',
        help="Determine the dataset to use. The selected operation will be applied to this dataset.")
    
    parser.add_argument('-l', '--list',
        action="store_true",
        help="Look for local datasets in $SCIPION_TESTS folder and for remote datasets in http://scipionwiki.cnb.csic.es/files/tests.")
    
    parser.add_argument('-v', '--verbose',
        action="store_true",
        help="Verbose mode. This will print more detailed messages")

    args = parser.parse_args()


    # Depending on the arguments selected, doing one thing or another

    deleteFlag=""
    if args.delete:
        deleteFlag=" -f"
    if args.list:
        print "List of local datasets in %(datasetsFolder)s" % ({'datasetsFolder': os.environ['SCIPION_TESTS']})
        folders = os.listdir(os.environ['SCIPION_TESTS'])
        for folder in folders:
            if os.path.isdir(os.path.join(os.environ['SCIPION_TESTS'], folder)):
                print "\t * " + folder.replace("dataset_", "")
        print "updating remote info..."
        backupFile(getManifestPath(basePath=os.environ['SCIPION_TMP']))
        try:
            urllib.urlretrieve(args.url+'/MANIFEST', getManifestPath(basePath=os.environ['SCIPION_TMP']))
            print "done."
        except:
            print "MANIFEST"+" ...ERROR"
            print "URL: "+args.url+'/MANIFEST'
            print "destination: "+getManifestPath(basePath=os.environ['SCIPION_TMP'])
            restoreFile(getManifestPath(basePath=os.environ['SCIPION_TMP']))
        if os.path.exists(getManifestPath(basePath=os.environ['SCIPION_TMP'])):
            print "List of remote datasets in %(urlAddress)s" % ({'urlAddress': args.url})
            manifestFile = open(getManifestPath(basePath=os.environ['SCIPION_TMP']), 'r+')
            manifestLines = manifestFile.readlines()
            for line in manifestLines:
                print "\t * "+os.path.split(line.replace("\n","").replace("dataset_",""))[1]
        delBackFile(getManifestPath(basePath=os.environ['SCIPION_TMP']))
            
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
        print "Selected dataset: %(dataset)s" % ({'dataset': args.dataset})
        if args.reverse_sync:
            print "Reverse synchronizing, BE CAREFUL!!! OPERATION EXTREMELY DANGEROUS!!!"
            ans = ""
            if "".join(args.reverse_sync) == "all":
                while ans != "y" and ans != "Y" and ans != "n" and ans != "N":
                    ans = raw_input("You're about to synchronize all your data with the remote data. Are you sure you want to continue? (y/n): ")
            else:
                print "Following files/folders are going to be upload to the remote scipion repository:"
                for sfile in args.reverse_sync:
                    print sfile
                while ans != "y" and ans != "Y" and ans != "n" and ans != "N":
                    ans = raw_input("You're about to synchronize all of them. Are you sure you want to smash 'em all? (y/n): ")
            if ans == ("y" or "Y"):
                print "You've chosen to proceed. Executing bash script " + args.syncfile + " -r ..."
                print "bash " + args.syncfile + deleteFlag + " -r " + " ".join(args.reverse_sync)
                #subprocess.call("bash " + args.syncfile + " -l " + args.last_mod_file + deleteFlag + " -r " + " ".join(args.reverse_sync), shell=True)
            else:
                print "You've chosen to abort. Goodbye!."
                sys.exit(3)
        else:
            datasetFolder = os.path.join(os.environ['SCIPION_TESTS'],"dataset_%(dataset)s" % ({'dataset': args.dataset}))
            if(os.path.exists(datasetFolder)):
                scipion_logo()
                print "Dataset working copy detected. Checking checksum for updates..."
                downloadDataset(args.dataset, destination=os.environ['SCIPION_TMP'], url=args.url, verbose=args.verbose, onlyMd5=True)
                checkForUpdates(args.dataset, url=args.url, verbose=args.verbose)
                cleanPath(os.path.join(os.environ['SCIPION_TMP'], "dataset_%(dataset)s" % ({'dataset': args.dataset})))
            else:
                scipion_logo()
                print "dataset not in local machine, trying to download..."
                downloadDataset(args.dataset, url=args.url, verbose=args.verbose)
                
        

if __name__ == "__main__":
    main(sys.argv[1:])
