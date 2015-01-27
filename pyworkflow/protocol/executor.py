# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module have the classes for execution of protocol steps.
The basic one will run steps, one by one, after completion.
There is one based on threads to execute steps in parallel
using different threads and the last one with MPI processes.
"""

import time
import datetime
import traceback
import threading

import pyworkflow.utils.process as process
import pyworkflow.protocol.constants as cts


class StepExecutor():
    """ Run a list of Protocol steps. """
    def __init__(self, hostConfig):
        self.hostConfig = hostConfig 
    
    def runJob(self, log, programName, params,           
           numberOfMpi=1, numberOfThreads=1, 
           env=None, cwd=None):
        """ This function is a wrapper around runJob, 
        providing the host configuration. 
        """
        process.runJob(log, programName, params,
                       numberOfMpi, numberOfThreads, 
                       self.hostConfig,
                       env=env, cwd=cwd)
    
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback):
        """ Simply iterate over the steps and run each one. """
        for s in steps:
            if not s.isFinished():
                s.setRunning()
                stepStartedCallback(s)
                s.run()
                doContinue = stepFinishedCallback(s)
                if not doContinue:
                    break


class StepThread(threading.Thread):
    """ Thread to run Steps in parallel. """
    def __init__(self, thId, step):
        threading.Thread.__init__(self)
        self.thId = thId
        self.step = step

    def run(self):
        try:
            self.step._run()  # not self.step.run() , to avoid race conditions
            self.step.setStatus(cts.STATUS_FINISHED)
        except Exception as e:
            self.step.setFailed(str(e))
            traceback.print_exc()
        finally:
            self.step.endTime.set(datetime.datetime.now())


class ThreadStepExecutor(StepExecutor):
    """ Run steps in parallel using threads. """
    def __init__(self, hostConfig, nThreads):
        StepExecutor.__init__(self, hostConfig)
        self.numberOfProcs = nThreads
        
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback):
        """ Create threads and synchronize the steps execution.
        n: the number of threads.
        """
        # Keep a list of busy slots ids (used to number threads and mpi)
        freeNodes = range(self.numberOfProcs)
        runningSteps = {}

        def getRunnable():
            """ Return the first step that is 'new' and all its
            dependencies have been finished.
            """
            for s in steps:
                if (s.getStatus() == cts.STATUS_NEW and
                        all(steps[i-1].isFinished() for i in s._prerequisites)):
                    return s
            return None

        while True:
            notRunning = [ns for ns in runningSteps.iteritems() if not ns[1].isRunning()]
            for nodeId, step in notRunning:
                runningSteps.pop(nodeId)
                freeNodes.append(nodeId)
                stepFinishedCallback(step)

            if freeNodes:
                step = getRunnable()
                if step is not None:
                    # We found a step to work in, so let's start a new
                    # thread to do the job and book it.
                    step.setRunning()
                    stepStartedCallback(step)
                    nodeId = freeNodes.pop()
                    runningSteps[nodeId] = step
                    StepThread(nodeId, step).start()
                elif all(not s.isRunning() for s in steps):
                    break  # yeah, we are done, either failed or finished :)
            time.sleep(0.1)
        # wait for all now
        for t in threading.enumerate():
            if t is not threading.current_thread():
                t.join()


class MPIStepExecutor(ThreadStepExecutor):
    """ Run steps in parallel using threads.
    But execution the runJob statement through MPI workers
    """
    def __init__(self, hostConfig, nMPI, comm):
        ThreadStepExecutor.__init__(self, hostConfig, nMPI)
        self.comm = comm
    
    def runJob(self, log, programName, params,           
           numberOfMpi=1, numberOfThreads=1, 
           env=None, cwd=None):
        from pyworkflow.utils.mpi import runJobMPI
        node = threading.current_thread().thId + 1
        print "==================calling runJobMPI=============================="
        print " to node: ", node
        runJobMPI(programName, params, self.comm, node,
                  numberOfMpi, hostConfig=self.hostConfig,
                  env=env, cwd=cwd)

    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback):
        ThreadStepExecutor.runSteps(self, steps, stepStartedCallback, stepFinishedCallback)

        # We import mpi here just to avoid failing in case MPI4py
        # was not properly compiled, we still can run in parallel with threads
        from pyworkflow.utils.mpi import TAG_RUN_JOB
        # Send special command 'None' to MPI slaves to notify them
        # that there is not more jobs to do and they can finish.
        for i in range(1, self.numberOfProcs+1):
            self.comm.send('None', dest=i, tag=(TAG_RUN_JOB+i))
