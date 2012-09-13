
def expandCommentRun(allowContinue=False):
    list = "Resume, Restart"
    linesStr = '''
#------------------------------------------------------------------------------------------
# {section}{has_question} Comment
#------------------------------------------------------------------------------------------
# Display comment
DisplayComment = False

# {text} Write a comment:
Comment = """Describe your run here..."""
#-----------------------------------------------------------------------------
# {section} Run 
#-----------------------------------------------------------------------------
# RUN name:
""" 
This will identify your protocol run. It need to be unique for each protocol. 
You could have <run1>, <run2> for protocol X, but not two run1 for same protocol. 
This name together with the protocol output folder will determine the working
directory for this run.
"""
RunName = "run_001"

# {list}(%(list)s) Run behavior
""" 
Resume from the last step, restart the whole process 
or continue at a given step or iteration.
"""
Behavior = "Resume"
'''
    if allowContinue:
        list += ", Continue"
        linesStr += '''
# {condition}(Behavior=="Continue") Continue at step:
""" Set to a positive number N to continue the protocol run at step N. """
ContinueAtStep = 1
'''
    return linesStr % locals()

def expandParallel(threads=1, mpi=8, condition="", hours=72, jobsize=0):
    conditionValue = 'True'
    if len(condition) > 0:
        conditionValue = condition
        condition = "{condition}(%s)" % condition
    linesStr = ''' 
#------------------------------------------------------------------------------------------
# {section} %(condition)s Parallelization
#------------------------------------------------------------------------------------------
# {hidden} Parallel Condition
ParallelCondition = "%(conditionValue)s"
'''
    if threads > 0:
        linesStr += '''
# Number of threads
""" 
This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than <Xmipp>
"""
NumberOfThreads = %(threads)d
'''
    if mpi > 0:
        linesStr += '''
# Number of MPI processes
""" 
This option provides the number of independent processes spawned 
in parallel by <mpirun> command in a cluster, usually throught
a queue system. This will require that you have compile <Xmipp>
with <mpi> support.
"""
NumberOfMpi = %(mpi)d
'''
    if jobsize > 0:
        linesStr += '''
        #MPI job size 
"""
Minimum size of jobs in mpi processes. 
Set to 1 for large images (e.g. 500x500)
and to 10 for small images (e.g. 100x100)
"""
MpiJobSize ='%(jobsize)d'
'''
        
    linesStr += '''
# Submit to queue ? 
"""Submit to queue
"""
SubmitToQueue = True

# {expert}{condition}(SubmitToQueue) Queue name
"""Name of the queue to submit the job
"""
QueueName = "default"

# {condition}(SubmitToQueue) Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the
queue system
"""
QueueHours = %(hours)d
''' 
    return linesStr % locals()

def expandExpert():
    return '''    
# {hidden} Show expert options
"""If True, expert options will be displayed """
ShowExpertOptions = False
'''
    
def expandJavaMemory():
    return '''
# {expert} Memory to use (In Gb)
"""
Amount of memory passed to the JVM for the application.
If you have very large micrographs (more than 10k pixels width)
is recommended to use 2Gb or more
"""
Memory = 2
'''
