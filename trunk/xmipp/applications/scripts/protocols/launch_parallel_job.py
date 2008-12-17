# The job should be launched from the working directory!
def launch_job(DoParallel,
               programname,
               mpiprogramname,
               params,
               log,
               MyNumberOfCPUs,
               MyMachineFile,
               RunInBackground):

    if DoParallel and MyNumberOfCPUs>1:
        launch_parallel_job(mpiprogramname,
                            params,
                            log,
                            MyNumberOfCPUs,
                            MyMachineFile,
                            RunInBackground)
    else:
        launch_sequential_job(programname,
                              params,
                              log,
                              RunInBackground)
        
def launch_sequential_job(programname,
                          params,
                          log,
                          RunInBackground):

    import os

    command=programname+' '+params
    if RunInBackground==True:
        command += ' &'
    print '* ',command,'\n'
    log.info(command)
    os.system(command)
   

def launch_parallel_job(mpiprogramname,
                        params,
                        log,
                        MyNumberOfCPUs,
                        MyMachineFile,
                        RunInBackground):

    import os
    if (MyMachineFile.find('USE_SLURM') > -1):
        command = 'srun `which '+ str(mpiprogramname) +'` ' + params
    else:
        if (len(MyMachineFile)==0):
            machineParam=""
            nr_cpus=MyNumberOfCPUs
        else:
            if (MyMachineFile[0]=="$"):
                machinefile= os.environ.get(MyMachineFile[1:])
            else:
                machinefile= MyMachineFile
            machinefile_flag = str(os.environ.get('XMIPP_MACHINEFILE_FLAG','machinefile'))
            machineParam=' -' + machinefile_flag + ' ' +  machinefile
            #machineParam=' -machinefile ' +  machinefile

            if (MyNumberOfCPUs<0):
                # use the number of entries in the machinefile
                fh = open(machinefile)
                lines = fh.readlines()
                fh.close()
                nr_cpus=len(lines)
            else:
                nr_cpus=MyNumberOfCPUs

        command='mpirun -np ' + str(nr_cpus)+machineParam
        command+=' `which '+ str(mpiprogramname) +'` ' + params

    if (RunInBackground==True):
        command = command + ' &'

    print '* ',command,'\n'
    log.info(command)
    os.system(command)
