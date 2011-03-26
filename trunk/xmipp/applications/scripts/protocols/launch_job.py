# The job should be launched from the working directory!
def launch_job(programname,
               params,
               log,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour,
               onlyBuildCommand=False):
    import os, sys
    command = buildCommand(programname,
               params,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour)
    if command == '':
            message = "Error: unrecognized SystemFlavour: ", SystemFlavour
            print '* ', message
            print '*********************************************************************'
            log.info(message)
            sys.exit(1)
    print '* ', command, '\n'
    log.info(command)
    from subprocess import call

    retcode = 0
    if (not onlyBuildCommand):
        try:
            retcode = call(command, shell=True)
            if retcode > 0:
                print >> sys.stderr, "Warning:", command, ' was terminated by signal', retcode
        except OSError, e:
            print >> sys.stderr, "Execution failed:", e
            exit(1)

    return(command, retcode);

def buildCommand(programname,
               params,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour):
    import os, sys

    if not DoParallel:
        command = programname + ' ' + params

    else:
        mpiprogramname = programname.replace('xmipp', 'xmipp_mpi')

        if (SystemFlavour == 'SLURM-MPICH'): # like BSCs MareNostrum, LaPalma etc
            mpicommand = 'srun '

        elif (SystemFlavour == 'TORQUE-OPENMPI'): # like our crunchy
            if (int(NumberOfThreads) > 1):
                mpicommand = 'mpirun -mca mpi_yield_when_idle 1 --bynode -np ' + str(NumberOfMpiProcesses)
            else:
                mpicommand = 'mpirun -mca mpi_yield_when_idle 1 -np ' + str(NumberOfMpiProcesses)

        elif (SystemFlavour == 'SGE-OPENMPI'): # like cluster at imp.ac.at (no variable nr_cpus yet...)
            mpicommand = 'mpiexec -n ' + str(NumberOfMpiProcesses)

        elif (SystemFlavour == 'PBS'): # like in Vermeer and FinisTerrae
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -hostfile ' + os.environ.get('PBS_NODEFILE')

        elif (SystemFlavour == 'XMIPP_MACHINEFILE'): # environment variable $XMIPP_MACHINEFILE points to machinefile
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -machinefile ' + os.environ.get('XMIPP_MACHINEFILE')

        elif (SystemFlavour == 'HOME_MACHINEFILE'): # machinefile is called $HOME/machines.dat
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -machinefile ' + os.environ.get('HOME') + '/machinefile.dat'

        elif (SystemFlavour == ''):
            mpicommand = 'mpirun -mca mpi_yield_when_idle 1 -np ' + str(NumberOfMpiProcesses)

        else:
            return ''
        command = mpicommand + ' `which ' + str(mpiprogramname) + '` ' + params
    return(command);
