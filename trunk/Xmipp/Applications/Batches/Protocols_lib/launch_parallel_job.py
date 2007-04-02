def launch_parallel_job(mpiprogramname,
                        params,
                        ParallelScript,
                        LaunchJobCommand,
                        log,
                        MyNumberOfCPUs,
                        MyMachineFile,
                        WorkingDir):
    import os
    fh=open(ParallelScript,'r')
    lines=fh.readlines()
    newlines=[]
    for line in lines:
        line=line.replace('MyNumberOfCPUs',str(MyNumberOfCPUs))
        line=line.replace('MyMachineFile' ,str(MyMachineFile))
        newlines+=line 
    line="`which "+ str(mpiprogramname)+"` "+params
    newlines+=line
    scriptname=str(WorkingDir)+'.script'
    fh=open(scriptname,'w')
    fh.writelines(newlines)
    os.chmod(scriptname,0777)
    if (LaunchJobCommand==""):
        command=scriptname+' & \n' 
    else:
        command=LaunchJobCommand + ' ' + scriptname + ' & \n' 
    print '* ',command
    log.info(command)
    os.system(command)
