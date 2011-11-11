#!/usr/bin/env xmipp_python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

from protlib_base import XmippProject, getProtocolFromModule
from protlib_utils import printLog
from protlib_xmipp import XmippScript
from mpi4py import MPI
from protlib_sql import SqliteDb

JOB_REQUEST = 0
JOB_RESPONSE = 1
NO_JOBS = -1
NO_AVAIL_JOB = -2
WAIT_TIME = 30 # 5 min


class ScriptParallelStepRunner(XmippScript):
    def __init__(self):
        XmippScript.__init__(self, True)
        
    def defineParams(self):
        self.addUsageLine('Run steps from protocol database in parallel')
        ## params
        self.addParamsLine('--script <protocol_script>     : Protocol run script')
        self.addParamsLine('--range <from_step> <to_step>  : Range of steps to run')
            
    def run(self):
        try:
            script = self.getParam('--script')
            fromStep = self.getIntParam('--range', 0)
            toStep = self.getIntParam('--range', 1)
            project = XmippProject()
            project.load()
            # Create protocol instance from script and setup database
            protocol = getProtocolFromModule(script, project)
            protocol.runSetup(isMainLoop=False)
            self.db = protocol.Db
            # All nodes retrieve parallel steps to work on 
            self.steps = self.db.getStepsRange(fromStep, toStep)
            

            comm = MPI.COMM_WORLD
            self.rank = comm.Get_rank()
            self.size = comm.Get_size()
            
            if self.rank:
                self.runWorker(comm)
            else:
                self.runMaster(comm)
            
            #comm.Disconnect()
        
        except Exception, e:
            printLog("Stopping MPI process because of error %s" % e, self.db.Log, out=True, err=True, isError=True)
            self.db.updateRunState(SqliteDb.RUN_FAILED)
            exit(1)
            
    def createStepsDict(self):
        stepsDict = {}
        for step in self.steps:
            stepId = step['step_id']
            stepsDict[stepId] = StepInfo(id=stepId, step=step, free=True, finish=False)
        return stepsDict
            
    def runMaster(self, comm):
        ''' Master will distribute steps to run and write to db'''
        # Build a dictionary with steps dependencies
        stepsDict = self.createStepsDict()
        
        #f = open('nodo%d.log' % self.rank, 'w')
        workingNodes = self.size - 1
        remainingJobs = len(self.steps)
        
        while workingNodes:
            #print >> f, "workingNodes %d, remainingJobs %d " % (workingNodes, remainingJobs)
            # Wait for a step request
            status = MPI.Status()
            jobId = NO_JOBS
            ##print >> f, "Waiting for request"
            jobId = comm.recv(None, source=MPI.ANY_SOURCE, tag=JOB_REQUEST, status=status)
            #print >> f, "Request received from nodo %d" % status.source
            if jobId != NO_JOBS:
                #Update job completion
                si = stepsDict[jobId]
                si.finish = True
                try:
                    self.db.endSingleStep(si.step, si.info)
                except Exception, e:
                    #print >> f, "ERROR: ", e
                    pass
                #break
                
            # Try to find next step to execute
            try:
                if remainingJobs:
                    jobId = NO_AVAIL_JOB
                    #print >> f, "Searching for job"
                    for step in self.steps:
                        step_id = step['step_id']
                        parent_id = step['parent_step_id']
                        si = stepsDict[step_id]
                        #print >> f, "id: ", step_id, si.free, si.finish, 
                        #print >> f, "parent: ", parent_id
                        if si.free and (not stepsDict.has_key(parent_id) or stepsDict[parent_id].finish):
                            #Select this step to send to worker
                            si.free = False
                            jobId = step_id
                            remainingJobs -= 1
                            try:
                                si.info = self.db.beginSingleStep(si.step)
                            except Exception, e:
                                #print >> f, "ERROR: ", e
                                pass
                            break
                else:
                    jobId = NO_JOBS
                    workingNodes -= 1 # Decrement number of working nodes            
                # Send response to worker
                #print >> f, "job returned %d to node %d" % (jobId, status.source)
                comm.send(jobId, dest=status.source, tag=JOB_RESPONSE)
            except Exception, e:
                #print >> f, "ERROR: ", e
                comm.send(NO_JOBS, dest=status.source, tag=JOB_RESPONSE)
                pass
        #print >> f, "MASTER FINISH: No more jobs to do"
    
    def runWorker(self, comm):
        '''Worker will ask for steps to run'''
        moreJobs = True
        jobId = NO_JOBS
        stepsDict = self.createStepsDict()
        #f = open('nodo%d.log' % self.rank, 'w')
        
        while moreJobs:
            #print >> f, "Asking for job, sending REQUEST"
            comm.send(jobId, dest=0, tag=JOB_REQUEST)
            #print >> f, "Waiting for job, expecting RESPONSE"
            jobId = comm.recv(None, source=0, tag=JOB_RESPONSE)
            #print >> f, "Job obtained: %d" % jobId
            if jobId == NO_AVAIL_JOB:
                # Sleep a time before requesting new jobs
                from time import sleep
                #print >> f, "No jobs at this moment, sleeping..."
                sleep(WAIT_TIME)
                jobId = NO_JOBS # Set this as expected by master when not job done
            elif jobId == NO_JOBS:
                moreJobs = False
            else: # Job obtained
                #print >> f, "Node %d: running step: %d" % (self.rank, jobId)
                try:
                    #pass
                    self.db.execSingleStep(stepsDict[jobId].step)
                except Exception, e:
                    pass
                    #print >> f, "ERROR: ", e
                #print >> f, "Node %d: step: %d run finished" % (self.rank, jobId)
    
class StepInfo():
    def __init__(self, **args):
        for k, v in args.iteritems():
            setattr(self, k, v)

if __name__ == '__main__':
    ScriptParallelStepRunner().tryRun()
    
