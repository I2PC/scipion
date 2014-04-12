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
from pyworkflow.utils.path import moveFile, makePath, makeFilePath

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
        
def downloadDataset(datasetName, destination=os.path.join(os.environ['SCIPION_USER_DATA'],'data','tests'), url="http://scipionwiki.cnb.csic.es/files/scipion/data/tests/NACHOTESTS_PLEASEDONOTREMOVE/tests", verbose=False):
    datasetFolder = os.path.join(destination, "dataset_%(dataset)s" % ({'dataset': datasetName}))
    makePath(datasetFolder)
    try:
        if verbose:
            print "retreiving MANIFEST file"
        urllib.urlretrieve(url+'/dataset_'+datasetName+'/MANIFEST', getManifestPath(dataset=datasetName))
    except:
        print "URL could not be retreived"
        if verbose:
            print "URL: " + url+'/dataset_'+datasetName+'/MANIFEST'
            print "destination:" + getManifestPath(dataset=datasetName)
    manifestFile = open(getManifestPath(dataset=datasetName), 'r+')
    manifestLines = manifestFile.readlines()
    print "Fetching dataset %(dataset)s files..." % ({'dataset': datasetName})
    totalNumber = len(manifestLines)*2.0
    percent = prevPercent = 0
    downloadBarWidth = 100
    if not verbose:
        sys.stdout.write("[%s]" % (" " * downloadBarWidth))
        sys.stdout.flush()
        sys.stdout.write("\b" * (downloadBarWidth+1))
    for number, line in enumerate(manifestLines):
        for indx, fileOption in enumerate([ os.path.normpath(line.replace("\n","")), os.path.normpath(os.path.join('md5', os.path.relpath(line.replace("\n",""), 'data') + ".md5")) ]):
            file = fileOption
            md5 = ''
            if indx == 1:
                md5 = 'md5'
            makeFilePath(os.path.join(datasetFolder, file))
            try:
                urllib.urlretrieve(url+'/dataset_'+datasetName+'/'+file, os.path.join(datasetFolder, file))
                percent = ((number*2+indx+1)/(totalNumber*1.0))*100
                if verbose:
                    print "\t "+file+" ...OK ( %(percent)02d " % ({'percent': percent}) + ' %)'
                else:
                    progress = percent-prevPercent
                    remaining = downloadBarWidth-progress
                    for x in range(0, (int)(progress//1)):
                        sys.stdout.write("#")
                        sys.stdout.flush()
                        progress-=1
                    if progress != 0:
                        if ((percent-progress)//1) != ((percent)//1):
                            sys.stdout.write("#")
                            sys.stdout.flush()
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
        md5sum = 0
        with open(os.path.join(datasetFolder, os.path.normpath(line.replace("\n",""))),'r+') as fileToCheck:
            data = fileToCheck.read()
            md5sum = hashlib.md5(data).hexdigest()
        md5file = open(os.path.normpath(os.path.join(datasetFolder, 'md5', os.path.relpath(line.replace("\n",""), 'data') + ".md5")),'r+')
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
    print "done"
    print ""
    #print "Executing bash script " + args.syncfile + "..."
    #print "bash " + args.syncfile + " -l " + args.last_mod_file + deleteFlag
    #subprocess.call("bash " + args.syncfile, shell=True)

def getManifestPath(basePath=os.path.join(os.environ['SCIPION_USER_DATA'],'data','tests'), dataset=""):
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
    
    parser.add_argument('-l', '--list-local-datasets',
        action="store_true",
        help="Look for local datasets in $SCIPION_HOME/data/tests folder.")
    
    parser.add_argument('-x', '--list-remote-datasets',
        action="store_true",
        help="Look for remote datasets in http://scipionwiki.cnb.csic.es/files/tests folder.")
    
    parser.add_argument('-v', '--verbose',
        action="store_true",
        help="Verbose mode. This will print more detailed messages")

    args = parser.parse_args()


    # Depending on the arguments selected, doing one thing or another

    deleteFlag=""
    if args.delete:
        deleteFlag=" -f"
    if args.list_local_datasets:
        if os.path.exists(getManifestPath()):
            print "List of local datasets in %(datasetsFolder)s" % ({'datasetsFolder': os.path.join(os.environ['SCIPION_USER_DATA'], 'data', 'tests')})
            manifestFile = open(getManifestPath(), 'r+')
            manifestLines = manifestFile.readlines()
            for line in manifestLines:
                print "\t * "+os.path.split(line.replace("\n","").replace("dataset_",""))[1]
            print "-"*40
    elif args.list_remote_datasets:
        print "updating remote info..."
        backupFile(getManifestPath(basePath=os.environ['SCIPION_TMP']))
        try:
            urllib.urlretrieve(args.url+'/MANIFEST', getManifestPath(basePath=os.environ['SCIPION_TMP']))
            print "done."
        except:
            traceback.print_exc()
            restoreFile(getManifestPath(basePath=os.environ['SCIPION_TMP']))
        if os.path.exists(getManifestPath(basePath=os.environ['SCIPION_TMP'])):
            print "List of remote datasets in %(urlAddress)s" % ({'urlAddress': args.url})
            manifestFile = open(getManifestPath(basePath=os.environ['SCIPION_TMP']), 'r+')
            manifestLines = manifestFile.readlines()
            for line in manifestLines:
                print "\t * "+os.path.split(line.replace("\n","").replace("dataset_",""))[1]
            print "-"*40
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
            datasetFolder = os.path.join(os.environ['SCIPION_USER_DATA'],'data','tests',"dataset_%(dataset)s" % ({'dataset': args.dataset}))
            if(os.path.exists(datasetFolder)):
                print "checking checksum for updates..."
                
                print "local checksum: " 
            else:
                scipion_logo()
                print "dataset not in local machine, trying to download..."
                downloadDataset(args.dataset, url=args.url, verbose=args.verbose)
                
        

if __name__ == "__main__":
    main(sys.argv[1:])
