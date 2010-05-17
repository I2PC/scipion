#!/usr/bin/env python
"""This is a script to isolate each qsub implementation from
the user. It needs to be modified for each cluster:

std_file is a template pbs file
The command line options and replace are machine dependent
"""
# modification Roberto Marabini 9-oct-2009 
#create temp name for script files so several can be sent to a queue

import os, string
import sys
import optparse
#local QSUB
QSUB = 'qsub'
#check command line
            
def main():
    
    std_file="""
#!/bin/bash
### Inherit all current environment variables
#PBS -V
### Job name
#PBS -N XXXjobIDXXX
### Queue name
#PBS -q high
### Standard output and standard error messages
#PBS -k eo
### Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=XXXnodesXXX:ppn=XXXppnXXX
### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=XXXhoursXXX:00:00
#################################
### Switch to the working directory;
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR
# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`
# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`
### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes
echo PBS_NODEFILE:
cat $PBS_NODEFILE
#################################

python XXXscriptXXX
#delete temporal file
rm XXXscriptXXX
    """.strip()#note strip removes first newline
    std_file += "\n"    
    #check command line

    ####Crunchy####
    parser = optparse.OptionParser("usage: %prog [options] protocols_file.py")
    parser.add_option("-t", "--time", dest="wallClockTime",
                      default="72", type="int",
                      help="Maximum execution time in wallclock hours (default=24)")
    parser.add_option("-m", "--memory", dest="maximumMemory", default='16g',
                      type="string", help="Memory per mpi process (default=16g)")
    parser.add_option("-i", "--id", dest="jobID", default='-1',
                      type="string", help="Job ID (default=-1)")

    #options are the options ;-)
    #arg the values that do not requiere '-flags', that is, the python script
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        return
    
    wallClockTime = options.wallClockTime
    maximumMemory = options.maximumMemory
    jobID         = options.jobID
    
    inFileName       = args[1]

    #get NumberOfMpiProcesses and NumberOfThreads from protocol script itself

    pyFile = open(inFileName,"r")
    doParallel = False
    numberOfThreads = 1
    numberOfNodes = 1
    workingDir='xmipp'
    while 1:
        line = (pyFile.readline()).replace('\n','')
        if (line.find("WorkingDir") == 0):
            myline=line.split("=")
            workingDir=myline[1].replace('/','_')
        elif (line.find("NumberOfThreads") == 0):
            myline=line.split("=")
            numberOfThreads = int(myline[1])
        elif (line.find("DoParallel=True") == 0):
            doParallel = True
        elif (line.find("NumberOfMpiProcesses=") == 0):
            myline=line.split("=")
            numberOfNodes = int(myline[1])
        elif (line.find("end-of-header") > 0):
            break

    if (not doParallel):
        numberOfNodes = 1

    if (jobID == '-1'):
        jobID=workingDir.replace("'","")
        jobID=jobID.replace('"','')
     
    #ask for a unique name
    tmp_inFileName=uniquefile(inFileName)
    #copy standard name file to tempporal file
    import shutil
    shutil.copy(inFileName,tmp_inFileName)
    #Replace old words with new ones
    outFileName=tmp_inFileName.replace(".py",".pbs");
    o = open(outFileName,"w")
    std_file = std_file.replace("XXXjobIDXXX",jobID)
    std_file = std_file.replace("XXXnodesXXX",str(numberOfNodes))
    std_file = std_file.replace("XXXppnXXX",str(numberOfThreads))
    std_file = std_file.replace("XXXhoursXXX",str(wallClockTime))
    std_file = std_file.replace("XXXscriptXXX",tmp_inFileName)
    std_file = std_file.replace("XXXmemXXX",maximumMemory)
    

    #create command and add it (as a comment) to the pbs file
    command = QSUB + " " 
    command += outFileName
    o.write(std_file)
    o.close()
    
    # Give the user some feedback
    print "---------------------------------------"
    print "qsub.py made job with:"
    print " #PBS -N "+jobID
    print " #PBS -l nodes="+str(numberOfNodes)+":ppn="+str(numberOfThreads)
    print " #PBS -l walltime="+str(wallClockTime)+":00:00"
    print " #PBS -l mem="+maximumMemory
    print "and submitted this script using command:"
    print command
    print "---------------------------------------"

    # Really submit the job to the queue
    os.system(command)
    
##################
#create unique file names
##################
def uniquefile(filename=""):
    counter = 0
    filename = filename.replace(".py","_")
    while 1:
        counter = counter + 1
        file = "%s%05d.py" % (filename,counter)
        if not os.path.exists(file):
            return file


if __name__ == '__main__':
    main()
                
                
                
