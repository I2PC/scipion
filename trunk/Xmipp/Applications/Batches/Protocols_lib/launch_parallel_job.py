# The job should be launched from the working directory!
def launch_job(DoParallel,
               programname,
               mpiprogramname,
               params,
               ParallelScript,
               LaunchJobCommand,
               log,
               MyNumberOfCPUs,
               MyMachineFile,
               ScriptName,
               RunInBackground=True):
    import os

    if not DoParallel:
        command=programname+' '+params
        if RunInBackground==True:
            command = command + ' &'
        command = command + '\n'    

        print '* ',command
        self.log.info(command)
        os.system(command)
    else:
        launch_parallel_job(mpiprogramname,
                            params,
                            ParallelScript,
                            LaunchJobCommand,
                            log,
                            MyNumberOfCPUs,
                            MyMachineFile,
                            ScriptName,
                            RunInBackground)
        

def launch_parallel_job(mpiprogramname,
                        params,
                        ParallelScript,
                        LaunchJobCommand,
                        log,
                        MyNumberOfCPUs,
                        MyMachineFile,
                        ScriptName,
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
    # ParallelScript should end either in \ or without newline
    line="`which "+ str(mpiprogramname)+"` "+params
    newlines+=line
    scriptname=ScriptName+'.script'
    fh.close()
    fh=open(scriptname,'w')
    fh.writelines(newlines)
    fh.close()
    os.chmod(scriptname,0777)
    if (LaunchJobCommand==""):
        command=scriptname
        if RunInBackground==True:
           command = command + ' &'
        command = command + '\n'    
    else:
        command=LaunchJobCommand + ' ' + scriptname
        if RunInBackground==True:
           command = command + ' &'
        command = command + '\n'    
    print '* ',command
    log.info(command)
    os.system(command)
