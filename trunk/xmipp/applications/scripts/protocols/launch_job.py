# The job should be launched from the working directory!
def launch_job(programname,
               params,
               log,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour,
	       onlyBuildCommand=False):

    import os,sys

    if not DoParallel:
        command = programname + ' ' + params

    else:
        mpiprogramname=programname.replace('xmipp','xmipp_mpi')

        if (SystemFlavour=='SLURM-MPICH'): # like BSCs MareNostrum, LaPalma etc
            mpicommand = 'srun '

        elif (SystemFlavour=='TORQUE-OPENMPI'): # like our crunchy
            if (int(NumberOfThreads) > 1):
                mpicommand = 'mpirun --bynode -np ' + str(NumberOfMpiProcesses) 
            else:
                mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) 

        elif (SystemFlavour=='SGE-OPENMPI'): # like cluster at imp.ac.at (no variable nr_cpus yet...)
            mpicommand = 'mpiexec -n ' + str(NumberOfMpiProcesses) 

        elif (SystemFlavour=='PBS'): # like in Vermeer and FinisTerrae
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -hostfile ' + os.environ.get('PBS_NODEFILE')

        elif (SystemFlavour=='XMIPP_MACHINEFILE'): # environment variable $XMIPP_MACHINEFILE points to machinefile
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -machinefile ' + os.environ.get('XMIPP_MACHINEFILE')

        elif (SystemFlavour=='HOME_MACHINEFILE'): # machinefile is called $HOME/machines.dat
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -machinefile ' + os.environ.get('HOME') + '/machinefile.dat'

        elif (SystemFlavour==''):
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses)

        else:
            message= "Error: unrecognized SystemFlavour: ", SystemFlavour
            print '* ',message
            print '*********************************************************************'
            log.info(message)
            sys.exit(1)

        command = mpicommand + ' `which '+ str(mpiprogramname) +'` ' + params

    print '* ',command,'\n'
    log.info(command)
    if (not onlyBuildCommand): os.system(command)
    return(command);
