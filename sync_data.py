#!/usr/bin/env python
import os
import subprocess
import getopt
import sys
import argparse

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

    exclusive.add_argument('-q', '--query-for-modifications', 
        action="store_true",
        help="Look for last_m.txt file. There, there will be or (1) the IP of the last modifier and date or (2) nothing. If the script finds (1) it moves it to modifications.log file and returns 0. If it finds (2) it check whether the file was modified in the last 30 minutes. If not, it returns 1, if yes it returns 2.")

    parser.add_argument('-d', '--delete', 
        action="store_true",
        help="Only valid when -r option enabled. This deletes remote files not existing in local. In the nutshell, it leaves the remote scipion data directory as it is the local one. Extremely dangerous option!.")

    exclusive.add_argument('-r', '--reverse-sync', 
        nargs='+',
        help="Synchronize from the local data to scipion machine. When wildcard 'all' is given, it will synchronize all the local folder with the remote one. When a set of locations is given, they're synchronized one by one against the remote scipion server. File path must be given from $SCIPION_HOME/tests folder")

    args = parser.parse_args()


    # Depending on the arguments selected, doing one thing or another

    deleteFlag=""
    if args.delete:
        deleteFlag=" -f"

    if args.query_for_modifications:
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

    elif args.reverse_sync:
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
            subprocess.call("bash " + args.syncfile + " -l " + args.last_mod_file + deleteFlag + " -r " + " ".join(args.reverse_sync), shell=True)
        else:
            print "You've chosen to abort. Goodbye!."
            sys.exit(3)

    else:
        scipion_logo()
        print "Executing bash script " + args.syncfile + "..."
        print "bash " + args.syncfile + " -l " + args.last_mod_file + deleteFlag
        subprocess.call("bash " + args.syncfile, shell=True)
        

if __name__ == "__main__":
    main(sys.argv[1:])
