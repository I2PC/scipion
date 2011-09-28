#------------------------------------------------------------------------------------------
# {section} Parallelization 
#------------------------------------------------------------------------------------------
# Number of threads
""" 
This option provides shared-memory parallelization on multi-core machines.
It does not require any additional software, other than <Xmipp>
"""
NumberOfThreads = 1

# Number of MPI processes
""" 
This option provides the number of independent processes spawned 
in parallel by <mpirun> command in a cluster, usually throught
a queue system. This will require that you have compile <Xmipp>
with <mpi> support.
"""
NumberOfMpi = 3

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
QueueHours = 72
