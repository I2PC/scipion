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
import constants as cts


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
        
    def _getRunnable(self, steps):
        """ Return the first step that is 'new' and all its
        dependencies have been finished, or None if none ready.
        """
        for s in steps:
            if (s.getStatus() == cts.STATUS_NEW and
                    all(steps[i-1].isFinished() for i in s._prerequisites)):
                return s
        return None
    
    def _arePending(self, steps):
        """ Return True if there are pending steps (either running or waiting)
        that can be done and thus enable other steps to be executed.
        """
        return any(s.isRunning() or s.isWaiting() for s in steps)
    
    def runSteps(self, steps, 
                 stepStartedCallback, 
                 stepFinishedCallback,
                 stepsCheckCallback):
        # Even if this will run the steps in a single thread
        # let's follow a similar approach than the parallel one
        # In this way we can take into account the steps graph
        # dependency and also the case when using streamming

        while True:
            # Get an step to run, if there is one
            step = self._getRunnable(steps)
            if step is not None:
                # We found a step to work in, so let's start a new
                # thread to do the job and book it.
                step.setRunning()
                stepStartedCallback(step)
                step.run()
                doContinue = stepFinishedCallback(step)
            
                if not doContinue:
                    break
                
            elif not self._arePending(steps):  # nothing running
                break  # yeah, we are done, either failed or finished :)

            time.sleep(3)  # be gentle on the main thread
            stepsCheckCallback()

        stepsCheckCallback() # one last check to finalize stuff


class StepThread(threading.Thread):
    """ Thread to run Steps in parallel. """
    def __init__(self, thId, step, lock):
        threading.Thread.__init__(self)
        self.thId = thId
        self.step = step
        self.lock = lock

    def run(self):
        error = None
        try:
            self.step._run()  # not self.step.run() , to avoid race conditions
        except Exception as e:
            error = str(e)
            traceback.print_exc()
        finally:
            with self.lock:
                if error is None:
                    self.step.setStatus(cts.STATUS_FINISHED)
                else:
                    self.step.setFailed(error)
                self.step.endTime.set(datetime.datetime.now())


class ThreadStepExecutor(StepExecutor):
    """ Run steps in parallel using threads. """
    def __init__(self, hostConfig, nThreads):
        StepExecutor.__init__(self, hostConfig)
        self.numberOfProcs = nThreads
        
    def runSteps(self, steps, 
                 stepStartedCallback, 
                 stepFinishedCallback,
                 stepsCheckCallback):
        """ Create threads and synchronize the steps execution.
        n: the number of threads.
        """
        sharedLock = threading.Lock()

        runningSteps = {}  # currently running step in each node ({node: step})
        freeNodes = range(self.numberOfProcs)  # available nodes to send mpi jobs

        while True:
            # See which of the runningSteps are not really running anymore.
            # Update them and freeNodes, and call final callback for step.
            with sharedLock:
                nodesFinished = [node for node, step in runningSteps.iteritems()
                                 if not step.isRunning()]
            doContinue = True
            for node in nodesFinished:
                step = runningSteps.pop(node)  # remove entry from runningSteps
                freeNodes.append(node)  # the node is available now
                doContinue = stepFinishedCallback(step)  # and do final work on the finished step
                if not doContinue:
                    break
            if not doContinue:
                break

            # If there are available nodes, send next runnable step.
            with sharedLock:
                if freeNodes:
                    step = self._getRunnable(steps)
                    if step is not None:
                        # We found a step to work in, so let's start a new
                        # thread to do the job and book it.
                        step.setRunning()
                        stepStartedCallback(step)
                        node = freeNodes.pop()  # take an available node
                        runningSteps[node] = step
                        t = StepThread(node, step, sharedLock)
                        t.daemon = True  # won't keep process up if main thread ends
                        t.start()
                    elif not self._arePending(steps):  # nothing running
                        break  # yeah, we are done, either failed or finished :)

            time.sleep(3)  # be gentle on the main thread
            stepsCheckCallback()

        stepsCheckCallback()

        # Wait for all threads now.
        for t in threading.enumerate():
            if t is not threading.current_thread():
                t.join()


class MPIStepExecutor(ThreadStepExecutor):
    """ Run steps in parallel using threads.
    But call runJob through MPI workers.
    """
    def __init__(self, hostConfig, nMPI, comm):
        ThreadStepExecutor.__init__(self, hostConfig, nMPI)
        self.comm = comm
    
    def runJob(self, log, programName, params,
               numberOfMpi=1, numberOfThreads=1, env=None, cwd=None):
        # Import mpi here so if MPI4py was not properly compiled
        # we can still run in parallel with threads.
        from pyworkflow.utils.mpi import runJobMPI
        node = threading.current_thread().thId + 1
        runJobMPI(programName, params, self.comm, node,
                  numberOfMpi, hostConfig=self.hostConfig, env=env, cwd=cwd)

    def runSteps(self, steps, 
                 stepStartedCallback, 
                 stepFinishedCallback,
                 checkStepsCallback):
        ThreadStepExecutor.runSteps(self, steps, 
                                    stepStartedCallback, 
                                    stepFinishedCallback,
                                    checkStepsCallback)

        # Import mpi here so if MPI4py was not properly compiled
        # we can still run in parallel with threads.
        from pyworkflow.utils.mpi import TAG_RUN_JOB

        # Send special command 'None' to MPI slaves to notify them
        # that there are no more jobs to do and they can finish.
        for node in range(1, self.numberOfProcs+1):
            self.comm.send('None', dest=node, tag=(TAG_RUN_JOB+node))
