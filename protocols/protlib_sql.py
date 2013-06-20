from sqlite3 import dbapi2 as sqlite
from sqlite3 import IntegrityError
import pickle
import os, sys
from os.path import exists
from config_protocols import projectDefaults
from protlib_utils import reportError, getScriptPrefix, printLog
from protlib_xmipp import blueStr, headerStr, greenStr
#The following imports are not directly used, but are common operations
#that will be performed by running steps on database
from protlib_utils import runJob
from protlib_filesystem import *

runColumns = ['run_id',
              'run_name',
              'run_state',
              'script',
              'init',
              'last_modfied',
              'protocol_name',
              'comment',
              'group_name', 'pid', 'jobid']

def getRunDict(run):
    ''' Convert the sqlite3 row to a python dict '''
    zrun = dict(zip(runColumns, run))
    zrun['source'] = zrun['script']        
    return zrun   

NO_MORE_GAPS = 0 #no more gaps to work on
NO_AVAIL_GAP = 1 #no available gaps now, retry later
STEP_GAP = 2 #step gap to work on
DB_TIMEOUT = 1000 # ms for timeout waiting

def existsDB(dbName):
    """check if database has been created by checking if the table tableruns exist"""
    result = False
    try:
        connection = sqlite.Connection(dbName, timeout=DB_TIMEOUT)
        connection.row_factory = sqlite.Row
        cur = connection.cursor()
        _sqlCommand = "SELECT count(*) from sqlite_master where tbl_name = '%(TableRuns)s';" % projectDefaults
        cur.execute(_sqlCommand)
        result = True
    except sqlite.Error, e:
        print "Could not connect to database %s\nERROR: %S" % (dbName, str(e))
    return result

class SqliteDb:
    #Some constants of run state
    RUN_SAVED = 0
    RUN_LAUNCHED = 1
    RUN_STARTED = 2
    RUN_FINISHED = 3
    RUN_FAILED = 4
    RUN_ABORTED = 5
    
    EXEC_PARALLEL = 2 # Any value greater than 1 will be parallel
    EXEC_MAINLOOP = 1
    EXEC_ALWAYS = 0
    
    NO_JOBID = -1
    UNKNOWN_JOBID = 0
    
    StateNames = ['Saved', 'Launched', 'Running', 'Finished', 'Failed', 'Aborted']
    
    def execSqlCommand(self, sqlCmd, errMsg):
        """Helper function to execute sqlite commands"""
        try:
            self.cur.executescript(sqlCmd)
        except sqlite.Error, e:
            print errMsg, e.args[0]
            if(e.args[0].find('database is locked') != -1):
                print 'consider deleting the database (%s)' % self.dbName
            sys.exit(1)  
        self.connection.commit()
    
    def getRunId(self, protName, runName):
        self.sqlDict['protocol_name'] = protName
        self.sqlDict['run_name'] = runName        
        _sqlCommand = """ SELECT run_id 
                          FROM %(TableRuns)s 
                          WHERE run_name = '%(run_name)s'
                            AND protocol_name = '%(protocol_name)s' """ % self.sqlDict
                            
        self.cur.execute(_sqlCommand)
        result = self.cur.fetchone()
        if result:
            result = result['run_id']
        return result
    
    def updateRunState(self, runState, runId=None, cursor=None, connection=None):
        if runId:
            self.sqlDict['run_id'] = runId
        self.sqlDict['run_state'] = runState
        _sqlCommand = """UPDATE %(TableRuns)s SET
                            run_state = %(run_state)d
                        WHERE run_id = %(run_id)d""" % self.sqlDict
        if cursor is None:
            cursor = self.cur
            connection = self.connection
        cursor.execute(_sqlCommand)
        connection.commit()
        
    def updateRunPid(self, run):
        self.sqlDict.update(run)
        _sqlCommand = """UPDATE %(TableRuns)s SET
                            pid = %(pid)s
                        WHERE run_id = %(run_id)d""" % self.sqlDict
        self.cur.execute(_sqlCommand)
        self.connection.commit()
        
    def updateRunJobid(self, run):
        self.sqlDict.update(run)
        _sqlCommand = """UPDATE %(TableRuns)s SET
                            jobid = %(jobid)s
                        WHERE run_id = %(run_id)d""" % self.sqlDict
        self.cur.execute(_sqlCommand)
        self.connection.commit()   
        
    def getRunJobid(self, runId):
        self.sqlDict['run_id'] = runId
        _sqlCommand = """ SELECT jobid 
                          FROM %(TableRuns)s 
                          WHERE run_id = %(run_id)d""" % self.sqlDict                            
        self.cur.execute(_sqlCommand)
        result = self.cur.fetchone()
        if result:
            result = result['jobid']
        return result    
        
class XmippProjectDb(SqliteDb):
    LAST_STEP = -1
    FIRST_STEP = 1
    FIRST_ITER = 1
    BIGGEST_STEP = 99999
            
    def createProtocolTables(self):
        """Create protocol related tables:
           protocols
           groups
           protocol_group
        """
        self.execSqlCommand('pragma foreign_keys=ON',"Foreing key activation failed")
        
        _sqlCommand = """CREATE TABLE IF NOT EXISTS %(TableGroups)s
             (group_name TEXT PRIMARY KEY);""" % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableGroups)s' table: " % self.sqlDict)

        _sqlCommand = """CREATE TABLE IF NOT EXISTS %(TableProtocols)s
                     (protocol_name TEXT PRIMARY KEY);""" % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableProtocols)s' table: " % self.sqlDict)
        
        _sqlCommand = """CREATE TABLE IF NOT EXISTS %(TableProtocolsGroups)s
                     (protocol_name TEXT,
                      group_name TEXT,
                      PRIMARY KEY(protocol_name, group_name));""" % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableProtocolsGroups)s' table: " % self.sqlDict)        
              
        
    def __init__(self, dbName):
        try:
            self.dbName = dbName
            self.connection = sqlite.Connection(dbName, timeout=DB_TIMEOUT)
            self.connection.row_factory = sqlite.Row
            self.cur = self.connection.cursor()
            self.sqlDict = projectDefaults
            self.sqlDict['execution_parallel'] = SqliteDb.EXEC_PARALLEL
            self.sqlDict['execution_mainloop'] = SqliteDb.EXEC_MAINLOOP
            self.sqlDict['execution_always'] = SqliteDb.EXEC_ALWAYS
            #enable foreign keys must be executed BEFORE table creation
            self.execSqlCommand('pragma foreign_keys=ON',"Foreing key activation failed")
            self.createProtocolTables()
                            
            _sqlCommand = """CREATE TABLE IF NOT EXISTS %(TableRuns)s
                         (run_id INTEGER PRIMARY KEY AUTOINCREMENT,
                          run_name TEXT,  -- label 
                          run_state INT  DEFAULT 0,  -- state of the run, possible values are:
                                          -- 0 - Saved (Never executed)
                                          -- 1 - Launched (Submited to queue)
                                          -- 2 - Running (Directly or from queue)
                                          -- 3 - Finished (Run finish correctly)
                                          -- 4 - Failed (Run produced an error)
                                          -- 5 - Aborted
                          script TEXT,    -- strip full name
                          init DATE,      -- run started at
                          last_modified DATE, --last modification (edition usually)
                          protocol_name TEXT REFERENCES %(TableProtocols)s(protocol_name), -- protocol name
                          comment TEXT,             -- user defined comment
                          pid  INTEGER DEFAULT -1,  -- process id
                          jobid INTEGER DEFAULT -1, -- this will be different of -1 of queue launched jobs
                          CONSTRAINT unique_workingdir UNIQUE(run_name, protocol_name));""" % self.sqlDict
            self.execSqlCommand(_sqlCommand, "Error creating '%(TableRuns)s' table: " % self.sqlDict)
            _sqlCommand = """ CREATE TABLE IF NOT EXISTS %(TableSteps)s 
                         (step_id INTEGER DEFAULT 0, -- primary key (weak entity)
                         command TEXT,               -- command
                         parameters TEXT,            -- wrapper parameters
                         init DATE,                  -- process started at
                         finish DATE,                -- process finished at
                         verifyFiles TEXT,           -- list with files to modify
                         iter INTEGER DEFAULT 1,     -- for iterative scripts, iteration number
                                                     -- useful to resume at iteration n
                         execution_mode INTEGER,     -- Possible values are: 0 - DoAlways, 1 - Mainloop, >1 - Parallel
                                                     -- an external program that will run it
                         passDb BOOL,                -- Should the script pass the database handler
                         run_id INTEGER REFERENCES %(TableRuns)s(run_id)  ON DELETE CASCADE,
                                                     -- key that unify all processes belonging to a run 
                         parent_step_id INTEGER       -- do not set this as a foreign key
                                                      -- because can not reffer to step_id in the same row
                         ,PRIMARY KEY(step_id, run_id)
                         )""" % self.sqlDict
            self.execSqlCommand(_sqlCommand, "Error creating %(TableSteps)s table: " % self.sqlDict)

            _sqlCommand = """CREATE TABLE IF NOT EXISTS %(TableParams)s 
                            (parameters TEXT,
                            run_id INTEGER REFERENCES %(TableRuns)s(run_id) ON DELETE CASCADE);""" % self.sqlDict
            self.execSqlCommand(_sqlCommand, "Error creating '%(TableParams)s' table: " % self.sqlDict)
                    
            _sqlCommand = """CREATE TRIGGER IF NOT EXISTS increment_step_id 
                             AFTER INSERT ON %(TableSteps)s FOR EACH ROW  
                             BEGIN 
                                UPDATE steps SET step_id = (SELECT MAX(step_id) + 1 
                                                           FROM %(TableSteps)s 
                                                           WHERE run_id = NEW.run_id)
                                WHERE step_id = 0 AND run_id = NEW.run_id; 
                             END """ % self.sqlDict
            self.execSqlCommand(_sqlCommand, "Error creating trigger: increment_step_id ")
        except sqlite.Error, e:
            reportError("database initialization failed: " + e.args[0])
            
    def insertGroup(self, groupName):
        self.sqlDict['group'] = groupName
        _sqlCommand = "INSERT OR IGNORE INTO %(TableGroups)s VALUES('%(group)s');" % self.sqlDict
        self.cur.execute(_sqlCommand)
        
    def insertProtocol(self, groupName, protName):
        self.sqlDict['group'] = groupName
        self.sqlDict['protocol_name'] = protName
        #check if protocol exists
        _sqlCommand = "SELECT COUNT(*) FROM %(TableProtocols)s WHERE protocol_name = '%(protocol_name)s'" % self.sqlDict
        self.cur.execute(_sqlCommand)
        if self.cur.fetchone()[0] == 0:
            _sqlCommand = "INSERT OR IGNORE INTO %(TableProtocols)s VALUES('%(protocol_name)s')" % self.sqlDict
            self.cur.execute(_sqlCommand)
        _sqlCommand = "INSERT OR IGNORE INTO %(TableProtocolsGroups)s VALUES('%(protocol_name)s', '%(group)s')" % self.sqlDict
        self.cur.execute(_sqlCommand)
          
    def insertRun(self, run):#run_name, script, comment=''):
        self.sqlDict.update(run)
        _sqlCommand = """INSERT INTO %(TableRuns)s values(
                            NULL, 
                            '%(run_name)s', 
                            0,
                            '%(script)s', 
                            datetime('now'), 
                            datetime('now'), 
                            '%(protocol_name)s',
                            '%(comment)s', -1, -1); """ % self.sqlDict
        self.cur.execute(_sqlCommand)
        run['run_id'] = self.cur.lastrowid
        run['run_state'] = SqliteDb.RUN_SAVED
        self.connection.commit()
        
    def suggestRunName(self, protName, prefix=None):
        if not prefix:
            prefix = self.sqlDict['RunsPrefix']
        self.sqlDict['prefix'] = prefix
        self.sqlDict['protocol_name'] = protName
        _sqlCommand = """SELECT COALESCE(MAX(run_name), '%(prefix)s') AS run_name 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocols)s
                         WHERE run_name LIKE '%(prefix)s%%' 
                           AND protocol_name = '%(protocol_name)s'""" % self.sqlDict
        self.cur.execute(_sqlCommand)         
        lastRunName = self.cur.fetchone()[0]    
        prefix, suffix = getScriptPrefix(lastRunName)
        n = 1
        if suffix:
            n = int(suffix) + 1
        runName = "%s_%03d" % (prefix, n)
        return runName     
        
    def updateRun(self, run):
        self.sqlDict.update(run)
        _sqlCommand = """UPDATE %(TableRuns)s SET
                            run_state = %(run_state)d,
                            last_modified = datetime('now'),
                            comment = '%(comment)s'
                        WHERE run_id = %(run_id)d""" % self.sqlDict
                         
        self.execSqlCommand(_sqlCommand, "Error updating run: %(run_name)s" % run)  
        self.connection.commit()
        
    def deleteRun(self, run):
        self.sqlDict.update(run)
        _sqlCommand = "DELETE FROM %(TableRuns)s WHERE run_id = %(run_id)d " % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error deleting run: %(run_name)s" % run)
        #this should be working automatically with the ON DELETE CASCADE
        #but isn't working, so...
        _sqlCommand = "DELETE FROM %(TableSteps)s WHERE run_id = %(run_id)d " % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error deleting steps of run: %(run_name)s" % run)        
        
    def selectRunsCommand(self):
        sqlCommand = """SELECT run_id, run_name, run_state, script, 
                               datetime(init, 'localtime') as init, 
                               datetime(last_modified, 'localtime') as last_modified,
                               protocol_name, comment,
                               group_name, 
                               pid, jobid 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocolsGroups)s """ % self.sqlDict
        return sqlCommand
                         
    def selectRuns(self, groupName='All', order='DESC'):
        self.connection.commit()
        self.sqlDict['group'] = groupName
        sqlCommand = self.selectRunsCommand()
        if groupName != 'All':
            sqlCommand += "WHERE group_name = '%(group)s'"
        sqlCommand += " ORDER BY last_modified " + order
        self.cur.execute(sqlCommand % self.sqlDict) 
        return self.cur.fetchall()
    
    def selectRunByName(self, protocol_name, runName):
        self.sqlDict['run_name'] = runName
        self.sqlDict['protocol_name'] = protocol_name
        sqlCommand = self.selectRunsCommand() + """WHERE protocol_name = '%(protocol_name)s'
                                                     AND run_name = '%(run_name)s'
                                                   ORDER BY last_modified DESC """ % self.sqlDict
        self.cur.execute(sqlCommand) 
        return self.cur.fetchone()
    
    def selectRunsByProtocol(self, protocol_name, state=None):
        ''' Select runs of a give protocol, if the state arguments is passed
        will restrict select to thouse in that state.
        Example of use: select all FINISH runs from ML2D
        '''
        self.sqlDict['protocol_name'] = protocol_name
        sqlCommand = self.selectRunsCommand() + "WHERE protocol_name = '%(protocol_name)s'"
        if state:
            self.sqlDict['run_state'] = state
            sqlCommand += "AND run_state == %(run_state)d "
        sqlCommand += "ORDER BY last_modified DESC "
        self.cur.execute(sqlCommand % self.sqlDict) 
        return self.cur.fetchall()
    
    def getRunProgress(self, run):
        self.sqlDict['run_id'] = run['run_id']
        sqlCommand = """ SELECT COUNT(step_id) FROM %(TableSteps)s WHERE run_id=%(run_id)d""" % self.sqlDict 
        steps_total = self.cur.execute(sqlCommand).fetchone()[0]
        steps_done = self.cur.execute(sqlCommand + ' AND finish IS NOT NULL').fetchone()[0]
        if steps_done > 0 and steps_done == steps_total:
            self.updateRunState(SqliteDb.RUN_FINISHED, run['run_id'])
            #run['state'] = SqliteDb.RUN_FINISHED
        return (steps_done, steps_total)

    def getRunStateByName(self, protocol_name, runName):
        run = self.selectRunByName(protocol_name, runName)
        self.sqlDict['run_id'] = run['run_id']
        sqlCommand = "SELECT run_state FROM %(TableRuns)s WHERE run_id = %(run_id)d" % self.sqlDict
        self.cur.execute(sqlCommand)
        return self.cur.fetchone()[0]
    
    def getRunSteps(self, run):
        self.sqlDict['run_id'] = run['run_id']
        sqlCommand = """ SELECT step_id, iter, init, finish, command, 
                            parameters, verifyFiles, parent_step_id 
                         FROM %(TableSteps)s 
                         WHERE (run_id = %(run_id)d)
                         ORDER BY step_id """ % self.sqlDict 
        self.cur.execute(sqlCommand % self.sqlDict) 
        return self.cur.fetchall()       
     
class XmippProtocolDb(SqliteDb):
    def __init__(self, protocol, script, isMainLoop=True):
        self.ContinueAtStep = getattr(protocol, 'ContinueAtStep', 0) 
        self.runBehavior = getattr(protocol, 'Behavior', 'Resume')
        self.dbName = protocol.project.dbName
        self.Import = protocol.Import  
        self.Log = protocol.Log  
        self.NumberOfMpi = getattr(protocol, 'NumberOfMpi', 1) 
        self.sqlDict = projectDefaults
        self.connection = sqlite.Connection(self.dbName, timeout=DB_TIMEOUT)
        self.connection.row_factory = sqlite.Row
        self.cur = self.connection.cursor()
        self.cur_aux = self.connection.cursor()
        self.lastStepId = XmippProjectDb.FIRST_STEP
        self.iter = XmippProjectDb.FIRST_ITER
        self.ProjDir = "."
        self.execSqlCommand('pragma foreign_keys=ON',"Foreing key activation failed")
        self.protocolScript = script
        #get run_id
        run_id = self.getRunId(protocol.Name, protocol.RunName)
        if not run_id:
            reportError("Protocol run '%(run_name)s' has not been registered in project database" % self.sqlDict)
        self.sqlDict['run_id'] = run_id

        if isMainLoop:
            # Restart or resume, only meaningless for execution on main protocol loop
            if self.runBehavior == "Restart":
                self.insertStatus = True
                _sqlCommand = 'DELETE FROM %(TableSteps)s WHERE run_id = %(run_id)d' % self.sqlDict
                self.execSqlCommand(_sqlCommand, "Error cleaning table: %(TableSteps)s" % self.sqlDict)
            else:
                #This will select steps for comparision in resume mode when insertStep is invoked
                sqlCommand = """ SELECT step_id, iter, command, parameters, verifyFiles 
                                 FROM %(TableSteps)s 
                                 WHERE (run_id = %(run_id)d)
                                 ORDER BY step_id """ % self.sqlDict
                self.cur.execute(sqlCommand)
                self.insertStatus = False

    def setIteration(self, iteration):
        self.iter = iteration

    def getRunState(self, cursor):
        sqlCommand = "SELECT run_state FROM %(TableRuns)s WHERE run_id = %(run_id)d" % self.sqlDict
        cursor.execute(sqlCommand)
        return cursor.fetchone()[0]
    
    def getRunIter(self):
        sqlCommand = """ SELECT COALESCE(MAX(iter),1) FROM %(TableSteps)s WHERE run_id=%(run_id)d AND finish IS NOT NULL""" % self.sqlDict 
        return self.cur.execute(sqlCommand).fetchone()[0]

    def checkRunOk(self, cursor):
        return self.getRunState(cursor) == SqliteDb.RUN_STARTED 
    
    def runAborted(self):
        ''' Return true if the run has been aborted '''
        state = self.getRunState(self.cur)
        print "CHECKING RUN ABORTED, STATE: %d", state
        return state == SqliteDb.RUN_ABORTED 
        
    def differentParams(self, Parameters, RowParameters):
        coreParameters=dict(Parameters)
        exclusion = ['NumberOfMpi', 'NumberOfThreads']
        for p in exclusion:
            if p in coreParameters: del coreParameters[p] 
            if p in RowParameters: del RowParameters[p]
        return coreParameters != RowParameters
    
    def insertStep(self, command,
                           verifyfiles=[],
                           parent_step_id=None,
                           execution_mode=SqliteDb.EXEC_MAINLOOP,
                           passDb=False,
                           **_Parameters):
        
        if not parent_step_id:
            parent_step_id = self.lastStepId

        parameters = pickle.dumps(_Parameters, 0)
        #Eventually verifyfiles will be dropped and replaced by verifyfilesDictionary

        verifyfilesString = pickle.dumps(verifyfiles, 0)

        if not self.insertStatus:
            #This will use previous select query in constructor
            row = self.cur.fetchone()
            if row is None:
                self.insertStatus = True
            else:
                if self.runBehavior == "Continue" and row['step_id'] >= self.ContinueAtStep:
                    self.insertStatus = True
                elif self.differentParams(_Parameters, pickle.loads(str(row['parameters']))) or row['verifyFiles'] != verifyfilesString:
                    self.insertStatus = True
                else:
                    for f in verifyfiles:
                        if not exists(f):
                            self.insertStatus = True
                            break
                self.lastStepId = row['step_id']
                if self.insertStatus:
                    self.sqlDict['step_id'] = row['step_id']
                    sqlCommand = """DELETE FROM %(TableSteps)s 
                                           WHERE run_id = %(run_id)d 
                                             AND step_id>=%(step_id)d""" % self.sqlDict
                    self.cur.execute(sqlCommand)
                    self.connection.commit()
        if self.insertStatus:
            try:

#                if parent_step_id ==-1:
#                    self.execSqlCommand('pragma foreign_keys=OFF',"Foreing key deactivation failed")
#                    parent_step_id = self.lastStepId
                    
                self.cur_aux.execute("""INSERT INTO 
                                    %(TableSteps)s(command,
                                                   parameters,
                                                   verifyfiles,
                                                   iter,
                                                   execution_mode,
                                                   passDb,
                                                   run_id,
                                                   parent_step_id)
                                     VALUES (?,?,?,?,?,?,?,?)""" % self.sqlDict,
                                     [command,
                                      parameters,
                                      verifyfilesString,
                                      self.iter,
                                      execution_mode,
                                      passDb,
                                      self.sqlDict['run_id'],
                                      parent_step_id])

#                if parent_step_id ==-1:
#                     self.execSqlCommand('pragma foreign_keys=ON',"Foreing key activation failed")
                    
                #select the last step_id inserted for this run
                self.cur_aux.execute("""SELECT MAX(step_id) 
                                        FROM %(TableSteps)s 
                                        WHERE run_id = %(run_id)d""" % self.sqlDict)

                #self.connection.commit()
                self.lastStepId = self.cur_aux.fetchone()[0]
                #fill table with verify files aliases, since they are linke to step_id in cascade I
                #do not need to worry about deleting them
            except sqlite.Error, e:
                reportError("Cannot insert command " + e.args[0])
        return self.lastStepId
    
    def getStepsRange(self, fromStep, toStep):
        self.sqlDict['from_step'] = fromStep
        self.sqlDict['to_step'] = toStep
        sqlCommand = """ SELECT step_id, iter, passDb, command, parameters, verifyFiles, parent_step_id
                         FROM %(TableSteps)s 
                         WHERE run_id = %(run_id)d 
                         AND step_id >= %(from_step)d
                         AND step_id < %(to_step)d
                         ORDER BY step_id """ % self.sqlDict
        self.cur.execute(sqlCommand)
        return self.cur.fetchall()
          
    def runSteps(self):
        #Update run state to STARTED
        self.updateRunState(SqliteDb.RUN_STARTED)
        #Clean init and finish for unfinished and doAlways steps
        sqlCommand = """ UPDATE %(TableSteps)s  SET init=NULL, finish=NULL
                         WHERE run_id = %(run_id)d 
                         AND (finish IS NULL OR execution_mode = %(execution_always)d)""" % self.sqlDict
        self.cur.execute(sqlCommand)
        self.connection.commit()
        #Select steps to run, include parallel ones
        sqlCommand = """ SELECT step_id, iter, passDb, command, parameters, verifyFiles, execution_mode
                         FROM %(TableSteps)s 
                         WHERE run_id = %(run_id)d 
                         AND finish IS NULL
                         ORDER BY step_id """ % self.sqlDict
        self.cur.execute(sqlCommand)
        steps = self.cur.fetchall()
        n = len(steps)
        msg = '***************************** Protocol STARTED mode: %s' % self.runBehavior
        printLog(msg, self.Log, out=True, err=True)
        
        i = 0
        final_status = (SqliteDb.RUN_FINISHED, 'FINISHED')
        
        while i < n:
            #Execute the step
            try:
                # Detect if there are parallel steps and execute them in parallel
                first = i
                while i < n and steps[i]['execution_mode'] > SqliteDb.EXEC_MAINLOOP:
                    i += 1
                if first < i: # There are parallel steps
                    mpiForParallelSteps = self.NumberOfMpi
                    fromStep = steps[first]['step_id']
                    if i < n: 
                        toStep = steps[i]['step_id']
                    else:
                        toStep = 99999
                    self.runParallelSteps(fromStep, toStep, mpiForParallelSteps)
                if i < n: # To be sure there are normal steps after parallels
                    self.runSingleStep(self.connection, self.cur, steps[i])
                i += 1
            except Exception as e:
                msg = "Stopping batch execution"
                if self.runAborted():
                    msg += " ABORTED by user request"
                    final_status = (SqliteDb.RUN_ABORTED, 'ABORTED')
                else:
                    msg += " since one of the steps could not be performed: %s" % str(e)
                    final_status = (SqliteDb.RUN_FAILED, 'FAILED')                    
                printLog(msg, self.Log, out=True, err=True, isError=True)
                break
        
        self.updateRunState(final_status[0])
        msg = '***************************** Protocol %s' % final_status[1]
        printLog(msg, self.Log, out=True, err=True)
        
    def runParallelSteps(self, fromStep, toStep, mpiForParallelSteps):
        '''This should run in parallel steps from fromStep and to toStep-1'''
        numberOfMpi = min(mpiForParallelSteps, toStep - fromStep+1)
        numberOfMpi = max(numberOfMpi,2)
        script = self.protocolScript
        retcode = runJob(self.Log, "xmipp_steps_runner", 
                         "--script %(script)s --range %(fromStep)d %(toStep)d" % locals(), numberOfMpi)
        if retcode != 0:
            raise Exception('xmipp_mpi_steps_runner execution failed')

    def verifyStepFiles(self, fileList):
        missingFilesStr = ''
        for f in fileList:
            if not xmippExists(f):
                missingFilesStr += ' ' + f

        if len(missingFilesStr) > 0:
            raise Exception("Missing result files: " + missingFilesStr)
        
    def runSingleStep(self, _connection, _cursor, stepRow):
        info = self._beginSingleStep(_connection, _cursor, stepRow)
        try:
            self._execSingleStep(stepRow, info)
        except Exception as e:
            err = "         Step finished with error: %s: %s" % (info.stepStr, e)
            printLog(err, self.Log, out=True, err=True, isError=True)
            raise
        self._endSingleStep(_connection, _cursor, stepRow, info)
         
    def _getStepRowInfo(self, stepRow):
        myDict = self.sqlDict.copy()
        myDict['step_id'] = stepRow['step_id']   
        myDict['iter'] = stepRow['iter']  
        info = StepRowInfo()   
        info.command = stepRow['command']
        info.args = pickle.loads(str(stepRow["parameters"]))
        info.stepStr = "%(step_id)d (iter=%(iter)d)" % myDict
        info.dict = myDict
        return info
        
    def beginSingleStep(self, stepRow):
        return self._beginSingleStep(self.connection, self.cur, stepRow)
        
    def _beginSingleStep(self, _connection, _cursor, stepRow):
        info = self._getStepRowInfo(stepRow)
        # Print
        from pprint import pformat        
        msg = blueStr("-------- Step start:  %s" % info.stepStr)
        printLog(msg, self.Log, out=True, err=True)
        printLog(headerStr((info.command.split())[-1]),self.Log, out=True)
        printLog(pformat(info.args, indent=4, width=20),self.Log, out=False)
        # pprint.PrettyPrinter(indent=4, width=20,stream=self.log).pprint(info.args)
        # Set init time
        sqlCommand = """UPDATE %(TableSteps)s SET init = CURRENT_TIMESTAMP 
                        WHERE step_id=%(step_id)d
                          AND run_id=%(run_id)d""" % info.dict
        _cursor.execute(sqlCommand)
        _connection.commit()
        return info
        
    def execSingleStep(self, stepRow):
        info = self._getStepRowInfo(stepRow)
        self._execSingleStep(stepRow, info)
        
    def _execSingleStep(self, stepRow, info):
        # Execute Python function
        exec(self.Import)
        if stepRow['passDb']:
            exec ( info.command + '(self, **info.args)')
        else:
            exec ( info.command + '(self.Log, **info.args)')
        # Check that expected result files were produced
        self.verifyStepFiles(pickle.loads(str(stepRow["verifyFiles"])))
        
    def endSingleStep(self, stepRow, info):
        self._endSingleStep(self.connection, self.cur, stepRow, info)
        
    def _endSingleStep(self, _connection, _cursor, stepRow, info):  
        # Report step finish and update database
        msg = greenStr("         Step finished: %s" % info.stepStr)
        printLog(msg, self.Log, out=True, err=True)
        sqlCommand = """UPDATE %(TableSteps)s SET finish = CURRENT_TIMESTAMP 
                        WHERE step_id=%(step_id)d
                          AND run_id=%(run_id)d""" % info.dict
        _cursor.execute(sqlCommand)
        _connection.commit()

class StepRowInfo():
    pass
#    def __init__(self):
#        self.args = None
#        self.command = None
#        self.dict = None

def escapeStr(str):
    return "'%s'" % str.replace("'", "''") 
          
class ProgramDb():
    ''' Class to handle Programs DB in the Xmipp installation folder '''
    def __init__(self, dbName=None):
        if dbName is None:
            from protlib_xmipp import getProgramsDbName
            dbName = getProgramsDbName()
        self.dbName = dbName
        self.connection = sqlite.Connection(dbName, timeout=DB_TIMEOUT)
        self.connection.row_factory = sqlite.Row
        self.cursor = self.connection.cursor()
        self.cursor.execute('pragma foreign_keys=ON')
            
    def create(self):
        self.createTables()
        
    def commit(self):
        self.connection.commit()
        
    def createTables(self):
            sqlCommand = """DROP TABLE IF EXISTS Category;
                            CREATE TABLE Category (
                               id INTEGER PRIMARY KEY ASC AUTOINCREMENT, 
                               name TEXT UNIQUE, 
                               desc TEXT, 
                               prefixes TEXT);
                            INSERT INTO Category VALUES(NULL, 'Classification', NULL, 'classify_ ml_ mlf_');
                            INSERT INTO Category VALUES(NULL, 'CTF', NULL, 'ctf_');
                            INSERT INTO Category VALUES(NULL, 'Images', NULL, 'image_ micrograph_');
                            INSERT INTO Category VALUES(NULL, 'Metadatas', NULL, 'metadata_');
                            INSERT INTO Category VALUES(NULL, 'Phantoms', NULL, 'phantom_ pdb_');
                            INSERT INTO Category VALUES(NULL, 'Angular assignment', NULL, 'angular_');
                            INSERT INTO Category VALUES(NULL, 'Tomography', NULL, 'tomo_ xray_');
                            INSERT INTO Category VALUES(NULL, 'Transformations', NULL, 'transform_');
                            INSERT INTO Category VALUES(NULL, 'Volumes', NULL, 'volume_ reconstruct_ resolution_');                             
                            
                            DROP TABLE IF EXISTS Program;
                            CREATE TABLE Program (
                               id INTEGER PRIMARY KEY ASC AUTOINCREMENT,
                               category_id INTEGER, 
                               name TEXT UNIQUE,
                               usage TEXT,
                               examples TEXT,
                               keywords TEXT);
                               
                            DROP TABLE IF EXISTS Label;
                            CREATE TABLE Label (
                               id INTEGER PRIMARY KEY ASC AUTOINCREMENT,
                               name TEXT UNIQUE,
                               type TEXT,
                               enum TEXT UNIQUE,
                               comment TEXT);
                         """
            self.cursor.executescript(sqlCommand)
            self.connection.commit()
            
    def insertProgram(self, program):
        program['desc'] = escapeStr(program['desc'])
        sqlCommand = """INSERT INTO Program VALUES (
                           NULL, %(category_id)d, %(name)s, %(desc)s, %(keywords)s);
                     """ % program
        self.cursor.execute(sqlCommand)
                     
    def selectPrograms(self, category=None):
        categoryWhere = ""
        if category:
            categoryWhere = "WHERE category_id=%d" % category['id']
        sqlCommand = "SELECT * FROM Program %s ORDER BY name;""" % categoryWhere
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def selectProgram(self, program_name):
        sqlCommand = "SELECT * FROM Program WHERE name='%s';""" % program_name
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchone()
        
    def selectCategories(self):
        sqlCommand = "SELECT * FROM Category;"
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def updateProgramCategory(self, program_name, category):
        sqlCommand = "UPDATE Program SET category_id = %d WHERE name='%s';" % (category['id'], program_name)
        self.cursor.execute(sqlCommand)
        self.connection.commit()
        
    def selectLabels(self):
        sqlCommand = "SELECT * FROM Label;"
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def insertLabel(self, labelData):
        labelData['comment'] = escapeStr(labelData['comment'])
        sqlCommand = """INSERT INTO Label VALUES (
                           NULL, '%(name)s', '%(type)s', '%(enum)s', %(comment)s);
                     """ % labelData
        self.cursor.execute(sqlCommand)