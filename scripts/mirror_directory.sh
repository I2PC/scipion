#!/bin/bash
# See some rcomendations at the end of this script
#run for TIMEOUT seconds, 3 days= 259200, 7 days=604800
TIMEOUT=259200
#retry after SEC seconds
SEC=300

# contains(string, substring)
#
# Returns 0 if the specified string contains the specified substring,
# otherwise returns 1.
contains() {
    string="$1"
    substring="$2"
    if test "${string#*$substring}" != "$string"
    then
        result=0    # $substring is in $string
    else
        result=1    # $substring is not in $string
    fi
}

#Script Usage
#return string with help
#
help()
{
    printf "\n\nmirror_directory.sh:  synchronize origin directory with target directory.\n (target  may be a remote machine)"
    printf "\n\nUsage: $0 dataDir projectName targetDir"    
    printf   "\nExample1: $0 MyDataDir projectName MyTargetDir"    
    printf   "\nExample2: $0 MyDataDir projectName MyUserName@MyHost:\n\n"    
}

#call help and exit if script is called without parameters
if [ $# -eq 0 ]
 then
   help
   exit
fi

#set source and target directories,
SOURCE=$1
PROJECT=$2
TARGET=$3

#if target a remote computer?,
# result of contains function is stored in global variable $result
contains $TARGET "@"

#create origin and target directories (full relative path)
FROMDIRECTORY=$SOURCE/$PROJECT
TODIRECTORY=$TARGET/$PROJECT

if [ $result -eq 1 ]; then #if TARGET is remote machine skip this check
                           # and do not create directory
  if [ -d "$TODIRECTORY" ]; then
    echo "Target directory $TODIRECTORY exists. I should not overwrite it. I quit"
    exit 1
  fi
  mkdir -p $TODIRECTORY
fi

#main loop
#SECONDS has the number of seconds the script has been running
#it is a variable defined by the system
END=`expr $SECONDS + $TIMEOUT`
while [ $SECONDS -lt $END ]; do
   if [ $result -eq 1 ]; then
      if [ ! -d "$TODIRECTORY" ]; then
         echo "Target directory $TODIRECTORY does not exist. Did you remove the disk? I abort"
         break
      fi
   fi
   echo rsync -r -P -t -l $FROMDIRECTORY $TARGET
   rsync -r -P -t -l  $FROMDIRECTORY $TARGET
   date
   printf "sleeping $SEC sec.\n"
   ETA=`expr $END - $SECONDS`
   printf "Remaining time= $ETA sec.\n"
   sleep $SEC
done

#rsync flags
#-r: recurse into directories
#-P: --partial (keep partially transferred files) +
#        --progress (show progress during transfer)
#-t: preserve modification times
#-l copy symlinks as symlinks

# RECOMENDATIONS***********
# 1) Do not use external disk formated with vfat since it cannot
# handle symbolic links. Use some linux format as ext4 o xsf
#
# 2) If you format the disk as root remember to change the ownership
#    to scipion user or make it writtable by anybody
#
# 3) In first form remember to double click the directory
#
# 4) Many checks are made in local directories to
# but nothing is checked in remote ones
#
# 5) When done and disk is moved to another computer remember to link project
# in ScionUserData (alternatively there is a import project option)
#
# 6) fix links with command ~/Scipion/scipion_box/scipion python ~/Scipion/scipion_box/scripts/fix_links.py 20170227_uu_nn 20170227_uu_nn
