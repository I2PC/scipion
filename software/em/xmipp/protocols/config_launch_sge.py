#--------------------------------------------------------------------------------
# This is a sample config file for launching jobs to a queue system
# 
# There should be some variables in this file:
#
# 1.- FileTemplate: string holding the template for launch script
#      Following vars are availables to template:
#      - %(jobName)s     : job name 
#      - %(nodes)d       : number of mpi nodes
#      - %(threads)d     : number of threads
#      - %(hours)d       : limit of running hours
#      - %(memory)d      : limit of memory used
#      - %(command)s     : command to be executed
#
#
# 2.- ArgsTemplate: string with the template arguments passed to launch Program
#      Following vars are availables to template:
#      - %(file)s     : file to be launched
#
# 3.- MpiArgsTemplate: string with the template arguments to run MPI program
#      Following vars are availables to template:
#      - %(nodes)d       : number of mpi nodes
#      - %(command)s     : command to be executed in parallel
#
#
#--------------------------------------------------------------------------------

#Program to launch jobs
Program = "qsub"
#Arguments template to launch
ArgsTemplate = "%(file)s"
# Command to stop a job
StopCommand = "canceljob"
StopArgsTemplate = "%(jobid)d"
# Command to query about job status
QueryCommand = "qstat"
QueryArgsTemplate = "%(jobid)d"

# Mpi executable
MpiProgram = "mpirun"
# Mpi run arguments template
MpiArgsTemplate = "-np %(nodes)d %(command)s"

# Other templates for other environments: 
# Using nodefile setted as environment variable or at home
# hostfile = os.environ.get('HOSFILE')
# hostfile = os.environ.get('HOME') + "/machinefile.txt" 
# MpiArgsTemplate = "-n %(nodes)d -hostfile " + hostfile

# Queue submition file template
FileTemplate = """
#!/bin/bash
### Inherit all current environment variables
#$ -V
### Job name
#$ -N %(jobName)s
### Specify the number of nodes and threads for your job.
#PBS -l nodes=%(nodes)d:ppn=%(threads)d
### Tell the expected run-time for your job (HH:MM:SS)
#$ -l h_rt=%(hours)d:00:00
### Number of processors to use
#$ -pe orte* %(nodes)d
# Use as working dir the path where qsub was launched
WORKDIR=$SGE_O_WORKDIR
#################################
### Set environment varible to know running mode is non interactive
export XMIPP_IN_QUEUE=1
### Switch to the working directory;
cd $WORKDIR
#################################
echo Running on host `hostname`
echo Time is `date`
echo Working directory is `pwd`
#################################

%(command)s
"""

