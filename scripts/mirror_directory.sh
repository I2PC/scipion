#!/bin/sh

#run for three days
TIMEOUT=259200

#Script Usage
help()
{
    printf "\n\nmirror_directory.sh:  synchronize origin directory with target directory.\n (target  may be a remote machine)"
    printf "\n\nUsage: timeout runTime $0 sourceDir targetDir"    
    printf   "\nExample1: timeout $TIMEOUT $0 MyProjecDir MyTargetDir"    
    printf   "\nExample2: timeout $TIMEOUT $0 MyProjecDir MyUserName@MyHost:\n\n"    
}

#call help if script is called without parameters
if [ $# -eq 0 ]
 then
   help
   exit
fi

#set source and target directories,
SOURCE=$1
TARGET=$2

#retry after SEC seconds
SEC=300

#main loop
while true; do
   echo rsync -h -r -P -t $SOURCE $TARGET
   rsync -h -r -P -t $SOURCE $TARGET
   date
   printf "sleeping $SEC seconds\n"
   sleep $SEC
done

#-h: human readable numbers
#-r: recurse into directories
#-P: --partial (keep partially transferred files) +
#        --progress (show progress during transfer)
#-t: preserve modification times
