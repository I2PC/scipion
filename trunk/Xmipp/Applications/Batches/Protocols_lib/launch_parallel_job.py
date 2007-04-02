def launch_parallel_job(mpiprogramname,
                        params,
                        ParallelScript,
                        LaunchJobCommand,
                        log,
                        MyNumberOfCPUs,
                        MyMachineFile,
                        WorkingDir,
                        RunInBackground=True):
    import os

    fh=open(ParallelScript,'r')
    lines=fh.readlines()
    newlines=""
    for line in lines:
        line=line.replace('MyNumberOfCPUs',str(MyNumberOfCPUs))
        line=line.replace('MyMachineFile' ,str(MyMachineFile))
        if len(line)>0:
           newlines+=line 
    print newlines
    line="`which "+ str(mpiprogramname)+"` "+params
    print line
    newlines+=line
    scriptname=str(WorkingDir)+'/'+mpiprogramname+'.script'
    fh.close()
    fh=open(scriptname,'w')
    fh.writelines(newlines)
    fh.close()
    os.chmod(scriptname,0777)
    if (LaunchJobCommand==""):
        command=scriptname
        if RunInBackground==False:
           command = command + ' &'
        command = command + '\n'    
    else:
        command=LaunchJobCommand + ' ' + scriptname
        if RunInBackground==False:
           command = command + ' &'
        command = command + '\n'    
    print '* ',command
    log.info(command)
    os.system(command)
