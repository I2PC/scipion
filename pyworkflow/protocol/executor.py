# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from threading import Thread, Condition, Event, current_thread
import pyworkflow.utils.process as process

STATUS_FINISHED = "finished"  # successfully finished
STATUS_READY = "ready" # The step is ready for execution, i.e. all requirements are done
STATUS_WAITING_OTHERS = "waiting_others" # There are some prerequisites steps that are not done yet

NO_READY_STEPS = -1 # no ready steps at this moment, should wait for it
NO_MORE_STEPS = -2  # all steps were done and nothing else to do.


class StepExecutor():
    """ Run a list of Protocol steps. """
    def __init__(self, hostConfig):
        self.hostConfig = hostConfig 
    
    def runJob(self, log, programName, params,           
           numberOfMpi=1, numberOfThreads=1, 
           runInBackground=False):
        """ This function is a wrapper around runJob, 
        providing the host configuration. 
        """
        process.runJob(log, programName, params,
                       numberOfMpi, numberOfThreads, 
                       runInBackground, self.hostConfig)
    
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback):
        """ Simply iterate over the steps and run each one. """
        for s in steps:
            s.setRunning()
            stepStartedCallback(s)
            s.run()
            doContinue = stepFinishedCallback(s)
            if not doContinue:
                break


class StepThread(Thread):
    """ Thread to run Steps in parallel. 
    If there is not work to do in this moment, the thread
    should be waiting in his Event variable.
    When the event is set to True, can happens two thing:
    1. The step variable is None and the thread should exit
       since there is not more work to do, or
    2. The step will be run and reported back after completion.
    """
    def __init__(self, thId, condition):
        Thread.__init__(self)
        self.thId = thId
        self.condition = condition
        self.event = Event() # Wait for work or exit
        self.step = None
        
    def isReady(self):
        ready = not self.event.is_set()
        return ready 
    
    def setStep(self, step):
        self.step = step
        self.event.set() # Work to do!!!
    
    def run(self):
        while True:
            # Wait for work
            self.event.wait()
            if self.step is None:
                break
            self.step.run()
            # Notify finished step
            self.condition.acquire()
            self.event.clear()
            self.condition.notify()
            self.condition.release()
            
        
class ThreadStepExecutor(StepExecutor):
    """ Run steps in parallel using threads. 
    """
    def __init__(self, hostConfig, nThreads):
        StepExecutor.__init__(self, hostConfig)
        self.numberOfThreads = nThreads
        
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback):
        """ Create threads and syncronize the steps execution. 
        n: the number of threads.
        """
        self.stepStartedCallback = stepStartedCallback
        self.stepFinishedCallback = stepFinishedCallback
        self.steps = steps
        self.stepsLeft = len(steps)
        self.condition = Condition() # Condition over global state
        
        for s in steps:
            s.status.set(STATUS_WAITING_OTHERS)
        
        self.thList = []
        
        #print "main: creating %d threads. " % self.numberOfThreads
        for i in range(self.numberOfThreads):
            th = StepThread(i, self.condition)
            self.thList.append(th)
            th.start()
            
        self.condition.acquire()
        
        while self.stepsLeft:
            self._launchThreads() # Check ready steps and launch threads
            self.condition.wait() # Wait to some steps completed
            self._updateThreads() # Check stepsLeft and clean finished threads status

        self.condition.release() # Not needed
        
        #print "main: Waiting for threads..."
        # Wait until all threads finish
        for th in self.thList:
            self.finishThread(th) # Let thread exit
            
        #print "main: Exit..."
         
    def finishThread(self, th):
        th.setStep(None)
        th.join()
        
    def _getReadyThread(self):
        """ Get the first thread waiting for work. """
        for th in self.thList:
            #print "main: checking thread ready for th: ", th.thId
            if th.isReady():
                return th
        return None
    
    def _updateThreads(self):
        """ Check which threads are done with theirs job. """
        for th in self.thList:
            if th.isReady(): # Waiting for work
                if th.step is not None: # Step finished
                    self.stepFinishedCallback(th.step)
                    self.stepsLeft -= 1
                    th.step = None # clean the thread step

    def _isStepReady(self, step):
        """ Check if a step has all prerequisites done. """
        if step.status != STATUS_WAITING_OTHERS:
            return False
        for i in step._prerequisites:
            if self.steps[i-1].status != STATUS_FINISHED:
                #print "main: prerequisite ", i, " is not finished!!!"
                return False
        return True
                    
    def _launchThreads(self):
        """ Check ready steps and awake threads to work. """
        for s in self.steps:
            #print "main, step ", i
            if self._isStepReady(s):
                #print " step ready."
                th = self._getReadyThread()
                if th is None:
                    #print "main: no thread available"
                    break # exit if no available threads to work
                s.setRunning()
                self.stepStartedCallback(s)
                #print "main: ", "awaking thread: ", th.thId
                th.setStep(s) # Awake thread to work 


class MPIStepExecutor(ThreadStepExecutor):
    """ Run steps in parallel using threads.
    But execution the runJob statement through MPI workers
    """
    def __init__(self, hostConfig, nMPI, comm):
        ThreadStepExecutor.__init__(self, hostConfig, nMPI)
        #self.runJob = runJob
        self.comm = comm
    
    def runJob(self, log, programName, params,           
           numberOfMpi=1, numberOfThreads=1, 
           runInBackground=False):
        from pyworkflow.utils.mpi import runJobMPI
        node = current_thread().thId + 1
        print "==================calling runJobMPI=============================="
        print " to node: ", node
        runJobMPI(log, programName, params, self.comm, node,
                  numberOfMpi, numberOfThreads, 
                  runInBackground, self.hostConfig)
        
    def finishThread(self, th):
        from pyworkflow.utils.mpi import TAG_RUN_JOB
        th.setStep(None)
        self.comm.send('None', dest=th.thId+1, tag=TAG_RUN_JOB)
        th.join()
        
        
