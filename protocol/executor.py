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
    
    def runSteps(self, steps, stepStartedCallback, stepFinishedCallback, n=12):
        """ Create threads and syncronize the steps execution. 
        n: the number of threads.
        """
        from threading import Thread, Condition
        self.stepStartedCallback = stepStartedCallback
        self.stepFinishedCallback = stepFinishedCallback
        self.steps = steps
        self.stepsDone = []
        stepsLeft = len(self.steps)
        self.cond = Condition() # Condition over to-do steps
        self.cond2 = Condition() # Condition over done steps
        
        for s in steps:
            s.status.set(STATUS_WAITING_OTHERS)
            
        self._updateReadySteps()
        
        thList = []
        for i in range(n):
            th = Thread(target=self._threadMain, args=(i,))
            thList.append(th)
            th.start()
            
        # On steps completion, notify calling the callback
        while stepsLeft:
            #print"main: cond2 acquiring...", stepsLeft
            self.cond2.acquire()
            while stepsLeft and not len(self.stepsDone):
                self.cond2.wait()
            finished = self.stepsDone
            self.stepsDone = []
            self.cond2.release()
            #print"main: cond2 released...", stepsLeft
            n = len(finished)
            #print"main:  len(finished)", n
            if n > 0:
                self._reportStepsFinished(finished)
                
                stepsLeft -= n
                #print"main: Steps left: ", stepsLeft
            
        self.cond.acquire()
        self.cond.notify_all()
        self.cond.release()
        
        #print"main: Waiting for threads..."
        # Wait until all threads finish
        for th in thList:
            th.join()
            
        #print"main: Exit..."
         
         
    def _reportStepsFinished(self, finished):
        """ Updating step status should also need to be syncronized. """
        
        for s in finished:
            self.stepFinishedCallback(s)
        
        self.cond.acquire()
        newReady = self._updateReadySteps()
        #print"main: reporting finished, new: ", newReady
        self.cond.notify(newReady)
                
        self.cond.release()
           
    def _isReady(self, step):
        """ Check if a step has all prerequisites done. """
        for i in step._prerequisites:
            if self.steps[i-1].status != STATUS_FINISHED:
                #print"main: prerequisite ", i, " is not finished!!!"
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
            #printthNumber, ": cond acquiring..."
            self.cond.acquire()
            
            next = self._getNextStep()
            #printthNumber, ":     next: ", next
            while next == NO_READY_STEPS:
                self.cond.wait()
                next = self._getNextStep()
            
            if next != NO_MORE_STEPS:
                step = self.steps[next]
                step.setRunning()
                
            self.cond.release()
            #printthNumber, ": cond released..."
            
            if next == NO_MORE_STEPS: # No more work to do
                break
            
            step.run()
            
            # Get lock on done jobs
            #printthNumber, ": cond2 acquiring..."
            self.cond2.acquire()
            self.stepsDone.append(step)
            #printthNumber, "        step done."
            self.cond2.notify()
            self.cond2.release()
            #printthNumber, ": cond2 released..."
            

            
            