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
from protocol import STATUS_READY, STATUS_WAITING_OTHERS, STATUS_FINISHED

class StepExecutor():
    """ Run a list of Protocol steps. """
    
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback):
        """ Simply iterate over the steps and run each one. """
        for s in steps:
            s.setRunning()
            stepStartedCallback(s)
            s.run()
            doContinue = stepFinishedCallback(s)
            if not doContinue:
                break


NO_READY_STEPS = -1 # no ready steps at this moment, should wait for it
NO_MORE_STEPS = -2  # all steps were done and nothing else to do.

class ThreadStepExecutor(StepExecutor):
    """ Run steps in parallel using threads. """
    
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback, n):
        """ Create threads and syncronize the steps execution. 
        n: the number of threads.
        """
        from threading import Thread, Condition
        self.stepStartedCallback = stepStartedCallback
        self.stepFinishedCallback = stepFinishedCallback
        self.steps = steps
        self.cond = Condition()
        
        self._updateReadySteps()
        
        thList = []
        for i in range(n):
            th = Thread(target=self._threadMain, args=(i))
            thList.append(th)
            th.start()
            
        # Wait until all threads finish
        for th in thList:
            th.join()
            
    def _isReady(self, step):
        """ Check if a step has all prerequisites done. """
        for i in step._prerequisites:
            if self.steps[i].status != STATUS_FINISHED:
                return False
        return True
            
    def _updateReadySteps(self):
        """ Check which steps become READY after some 
        of theirs prerequisites have finished.
        Return the number of new ready steps 
        """
        
        newReady = 0
        for s in self.steps:
            if s.status == STATUS_WAITING_OTHERS and self._isReady(s):
                s.status.set(STATUS_READY)
                newReady += 1
        return newReady
    
    def _getNextStep(self):
        """ Get next step to execute. 
        This function should protect access to step list
        since will be called concurrently.
        The index of the step in the list is returned, or:
        NO_READY_STEPS: no ready steps at this moment, should wait for it
        NO_MORE_STEPS: all steps were done and nothing else to do.
        """
        finished = True
        
        for i, s in enumerate(self.steps):
            if s.status == STATUS_READY:
                return i
            if s.status == STATUS_WAITING_OTHERS:
                finished = False
        
        if finished:
            return NO_MORE_STEPS
        else:
            return NO_READY_STEPS
        
    def _threadMain(self, thNumber):
        """ Function that will be executed in parallel by each thread. """
        while True:
            self.cond.acquire()
            
            next = self._getNextStep()
            while next == NO_READY_STEPS:
                self.cond.wait()
                next = self._getNextStep()
            
            self.cond.release()
            
            if next == NO_MORE_STEPS: # No more work to do
                break
            
            step = self.steps[next]
            step.run()
            self.stepFinishedCallback(step)
    
    def _reportStepFinished(self, step):
        """ Updating step status should also need to be syncronized. """
        self.cond.acquire()
        
        self.stepFinishedCallback(step)
        newReady = self._updateReadySteps()
        self.cond.notify(newReady)
                
        self.cond.release()
            
            