#!/usr/bin/env python
import os
import subprocess
import getopt
import sys

def main(argv):
    syncFile = "sync_data"
    try:
        opts, args = getopt.getopt(argv, "hs:", ["help", "syncfile="])
    except getopt.GetoptError:
        print "ERROR: Bad syntax"
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-s", "--syncfile"):
            syncFile = arg

    print "Executing bash script " + syncFile + "..."
    subprocess.call("bash "+syncFile, shell=True)

def usage():
    print "#" * 100
    print "#  Sync Data Python Bash Executor (WRAPPER)"
    print "#" * 100

if __name__ == "__main__":
    main(sys.argv[1:])
