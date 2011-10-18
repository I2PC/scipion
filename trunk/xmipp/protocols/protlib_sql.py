from sqlite3 import dbapi2 as sqlite
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

NO_MORE_GAPS = 0 #no more gaps to work on
NO_AVAIL_GAP = 1 #no available gaps now, retry later
STEP_GAP   = 2 #step gap to work on

def existsDB(dbName):
    """check if database has been created by checking if the table tableruns exist"""
    result = False
    try:
        connection = sqlite.Connection(dbName)
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
    
    EXEC_GAP = 0
    EXEC_MAINLOOP = 1
    EXEC_ALWAYS = 2
    
    NO_JOBID = -1
    
    StateNames = ['Saved', 'Launched', 'Running', 'Finish', 'Failed', 'Aborted']
    
    def execSqlCommand(self, sqlCmd, errMsg):
        """Helper function to execute sqlite commands"""
        try:
            self.cur.executescript(sqlCmd)
        except sqlite.Error, e:
            print errMsg, e
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
                        WHERE run_id = %(run_id)d"""  % self.sqlDict
        if cursor is None:
            cursor = self.cur
            connection = self.connection
        cursor.execute(_sqlCommand)
        connection.commit()
        
    def updateRunPid(self, run):
        self.sqlDict.update(run)
        _sqlCommand = """UPDATE %(TableRuns)s SET
                            pid = %(pid)s
                        WHERE run_id = %(run_id)d"""  % self.sqlDict
        self.cur.execute(_sqlCommand)
        self.connection.commit()
        
    def updateRunJobid(self, run):
        self.sqlDict.update(run)
        _sqlCommand = """UPDATE %(TableRuns)s SET
                            jobid = %(jobid)s
                        WHERE run_id = %(run_id)d"""  % self.sqlDict
        self.cur.execute(_sqlCommand)
        self.connection.commit()   
        
class XmippProjectDb(SqliteDb):
    LAST_STEP  = -1
    FIRST_STEP = 1
    FIRST_ITER = 1
    BIGGEST_STEP = 99999
            
    def __init__(self, dbName):
        try:
            self.dbName = dbName
            self.connection = sqlite.Connection(dbName)
            self.connection.row_factory = sqlite.Row
            self.cur = self.connection.cursor()
            self.sqlDict = projectDefaults
            self.sqlDict['execution_gap'] = SqliteDb.EXEC_GAP
            self.sqlDict['execution_mainloop'] = SqliteDb.EXEC_MAINLOOP
            self.sqlDict['execution_always'] = SqliteDb.EXEC_ALWAYS
            
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
                          comment TEXT,       -- user defined comment
                          pid TEXT,           -- process id
                          jobid INTEGER DEFAULT -1, -- this will be different of -1 of queue launched jobs
                          CONSTRAINT unique_workingdir UNIQUE(run_name, protocol_name));""" % self.sqlDict
            self.execSqlCommand(_sqlCommand, "Error creating '%(TableRuns)s' table: " % self.sqlDict)
            
            _sqlCommand = """ CREATE TABLE IF NOT EXISTS %(TableSteps)s 
                         (step_id INTEGER DEFAULT 0, -- primary key (weak entity)
                         command TEXT,               -- comment (NOT USED, DROP?)
                         parameters TEXT,            -- wrapper parameters
                         init DATE,                  -- process started at
                         finish DATE,                -- process finished at
                         verifyFiles TEXT,           -- list with files to modify
                         iter INTEGER DEFAULT 1,     -- for iterative scripts, iteration number
                                                     -- useful to resume at iteration n
                         execution_mode INTEGER,     -- Possible values are: 0 - Gap, 1-Mainloop, 2-DoAlways
                                                     -- an external program that will run it
                         passDb BOOL,                -- Should the script pass the database handler
                         run_id INTEGER REFERENCES %(TableRuns)s(run_id)  ON DELETE CASCADE,
                                                     -- key that unify all processes belonging to a run 
                         parent_step_id INTEGER REFERENCES %(TableSteps)s(step_id) ON DELETE CASCADE,
                                                      -- parent_step_id step must be executed before 
                                                      -- step_id may be executed
                         PRIMARY KEY(step_id, run_id))""" % self.sqlDict
        
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
            reportError("database initialization failed: " + e)

    def insertGroup(self, groupName):
        self.sqlDict['group'] = groupName
        _sqlCommand = "INSERT INTO %(TableGroups)s VALUES('%(group)s');" % self.sqlDict
        self.cur.execute(_sqlCommand)
        
    def insertProtocol(self, groupName, protName):
        self.sqlDict['group'] = groupName
        self.sqlDict['protocol_name'] = protName
        #check if protocol exists
        _sqlCommand = "SELECT COUNT(*) FROM %(TableProtocols)s WHERE protocol_name = '%(protocol_name)s'" % self.sqlDict
        self.cur.execute(_sqlCommand)
        if self.cur.fetchone()[0] == 0:
            _sqlCommand = "INSERT INTO %(TableProtocols)s VALUES('%(protocol_name)s')" % self.sqlDict
            self.cur.execute(_sqlCommand)
        _sqlCommand = "INSERT INTO %(TableProtocolsGroups)s VALUES('%(protocol_name)s', '%(group)s')" % self.sqlDict
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
                            '%(comment)s',
                            -1, 0 --pid and type                            
                            ); """  % self.sqlDict
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
        prefix, suffix  = getScriptPrefix(lastRunName)
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
                        WHERE run_id = %(run_id)d"""  % self.sqlDict
                         
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
                         
    def selectRuns(self, groupName=None):
        self.connection.commit()
        self.sqlDict['group'] = groupName
        sqlCommand = self.selectRunsCommand()
        if groupName:
            sqlCommand += "WHERE group_name = '%(group)s'"
        sqlCommand += " ORDER BY last_modified DESC "
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
    
    def selectRunsByProtocol(self, protocol_name):
        self.sqlDict['protocol_name'] = protocol_name
        sqlCommand = self.selectRunsCommand() + """WHERE protocol_name = '%(protocol_name)s'
                                                   ORDER BY last_modified DESC """ % self.sqlDict
        self.cur.execute(sqlCommand) 
        return self.cur.fetchall()
    
    def getRunProgress(self, run):
        self.sqlDict['run_id'] = run['run_id']
        sqlCommand = """ SELECT COUNT(step_id) FROM %(TableSteps)s WHERE run_id=%(run_id)d""" % self.sqlDict 
        steps_total = self.cur.execute(sqlCommand).fetchone()[0]
        steps_done = self.cur.execute(sqlCommand + ' AND finish IS NOT NULL').fetchone()[0]
        return (steps_done, steps_total)

    def getRunStateByName(self,protocol_name, runName):
        run=self.selectRunByName(protocol_name, runName)
        self.sqlDict['run_id'] = run['run_id']
        sqlCommand = "SELECT run_state FROM %(TableRuns)s WHERE run_id = %(run_id)d" % self.sqlDict
        self.cur.execute(sqlCommand)
        return self.cur.fetchone()[0]
     
class XmippProtocolDb(SqliteDb):
    def __init__(self, protocol, isMainLoop=True):
        self.ContinueAtStep = getattr(protocol, 'ContinueAtStep', 0) 
        self.runBehavior = getattr(protocol, 'Behaviour', 'Resume')
        self.dbName = protocol.project.dbName
        self.Import = protocol.Import  
        self.Log = protocol.Log             
        self.sqlDict = projectDefaults
        self.connection = sqlite.Connection(self.dbName)
        self.connection.row_factory = sqlite.Row
        self.cur = self.connection.cursor()
        self.cur_aux = self.connection.cursor()
        self.lastStepId = XmippProjectDb.FIRST_STEP
        self.iter = XmippProjectDb.FIRST_ITER
        self.ProjDir = "."

        #get run_id
        run_id = self.getRunId(protocol.Name, protocol.RunName)
        if not run_id:
            reportError("Protocol run '%(run_name)s' has not been registered in project database" % self.sqlDict)
        self.sqlDict['run_id'] = run_id

        if isMainLoop:
            # Restart or resume, only meaningless for execution on main protocol loop
            if self.runBehavior=="Restart":
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

    def setIteration(self,iteration):
        self.iter=iteration

    def getRunState(self,cursor):
        sqlCommand = "SELECT run_state FROM %(TableRuns)s WHERE run_id = %(run_id)d" % self.sqlDict
        cursor.execute(sqlCommand)
        return cursor.fetchone()[0]

    def checkRunOk(self, cursor):
        return self.getRunState(cursor) == SqliteDb.RUN_STARTED 
        
    def insertStep(self, command,
                           verifyfiles=[],
                           parent_step_id=None, 
                           execution_mode = SqliteDb.EXEC_MAINLOOP,
                           passDb=False,
                           **_Parameters):
        if not parent_step_id:
            parent_step_id=self.lastStepId

        parameters = pickle.dumps(_Parameters, 0)
        verifyfilesString = pickle.dumps(verifyfiles, 0)
        
        if not self.insertStatus:
            #This will use previous select query in constructor
            row = self.cur.fetchone()
            if row is None:
                self.insertStatus = True
            else:
                if self.runBehavior=="Continue" and row['step_id']>=self.ContinueAtStep:
                    self.insertStatus = True
                if row['parameters'] != parameters or row['verifyFiles'] != verifyfilesString:
                    self.insertStatus = True
                else:
                    for f in verifyfiles:
                        if not exists(f):
                            self.insertStatus = True
                            break
                self.lastStepId=row['step_id']
                if self.insertStatus:
                    self.sqlDict['step_id'] = row['step_id']
                    sqlCommand = """DELETE FROM %(TableSteps)s 
                                           WHERE run_id = %(run_id)d 
                                             AND step_id>=%(step_id)d""" % self.sqlDict
                    self.cur.execute(sqlCommand)
                    self.connection.commit()
        if self.insertStatus:
            try:
                self.cur_aux.execute("""INSERT INTO 
                                    %(TableSteps)s(command,parameters,verifyFiles,iter,execution_mode,passDb,run_id,parent_step_id)
                                     VALUES (?,?,?,?,?,?,?,?)""" % self.sqlDict,
                                     [command, parameters, verifyfilesString, self.iter,execution_mode,passDb,self.sqlDict['run_id'], parent_step_id])
                #select the last step_id inserted for this run
                self.cur_aux.execute("""SELECT MAX(step_id) 
                                        FROM %(TableSteps)s 
                                        WHERE run_id = %(run_id)d""" % self.sqlDict)
                self.lastStepId = self.cur_aux.fetchone()[0]
            except sqlite.Error, e:
                reportError( "Cannot insert command: %s" % e.args[0])
        return self.lastStepId
    def runSteps(self):
        #Update run state to STARTED
        self.updateRunState(SqliteDb.RUN_STARTED)
        #Clean init and finish for unfinished and doAlways steps
        sqlCommand = """ UPDATE %(TableSteps)s  SET init=NULL, finish=NULL
                         WHERE run_id = %(run_id)d 
                         AND (finish IS NULL OR execution_mode = %(execution_always)d)""" % self.sqlDict
        self.cur.execute(sqlCommand)
        self.connection.commit()
        #Select steps to run
        sqlCommand = """ SELECT step_id, iter, passDb, command, parameters, verifyFiles 
                         FROM %(TableSteps)s 
                         WHERE run_id = %(run_id)d 
                         AND finish IS NULL AND execution_mode > %(execution_gap)d
                         ORDER BY step_id """ % self.sqlDict
        self.cur.execute(sqlCommand)
        commands = self.cur.fetchall()
        n = len(commands)
        msg='***************************** Protocol STARTED mode: %s'%self.runBehavior
        printLog(msg, self.Log, out=True, err=True)
        
        for i in range(n):
            #Execute the step
            try:
                self.runSingleStep(self.connection, self.cur, commands[i])
            except Exception as e:
                msg = "Stopping batch execution since one of the steps could not be performed: %s" % e
                printLog(msg, self.Log, out=True, err=True, isError=True)
                self.updateRunState(SqliteDb.RUN_FAILED)
                raise
        msg='***************************** Protocol FINISHED'
        printLog(msg, self.Log, out=True, err=True)
        self.updateRunState(SqliteDb.RUN_FINISHED)

    def verifyStepFiles(self, fileList):
        missingFilesStr = ''
        for f in fileList:
            if not exists(f):
                missingFilesStr += ' ' + f

        if len(missingFilesStr) > 0:
            raise Exception("Missing result files: " + missingFilesStr)
        
    def runSingleStep(self, _connection, _cursor, stepRow):
        exec(self.Import)
        step_id = stepRow['step_id']
        myDict = self.sqlDict.copy()
        myDict['step_id'] = step_id
        myDict['iter'] = stepRow['iter']
        command = stepRow['command']
        dict = pickle.loads(str(stepRow["parameters"]))
        
        # Print
        import pprint, time
        stepStr = "%d (iter=%d)" % (step_id, stepRow['iter'])
        msg = blueStr("-------- Step start:  %s" % stepStr)
        printLog(msg, self.Log, out=True, err=True)
        print headerStr((command.split())[-1])
        pprint.PrettyPrinter(indent=4,width=20).pprint(dict)

        # Set init time
        sqlCommand = """UPDATE %(TableSteps)s SET init = CURRENT_TIMESTAMP 
                        WHERE step_id=%(step_id)d
                          AND run_id=%(run_id)d""" % myDict
        _cursor.execute(sqlCommand)
        _connection.commit()
        
        # Execute Python function
        try:
            if stepRow['passDb']:
                exec ( command + '(self, **dict)')
            else:
                exec ( command + '(self.Log, **dict)')
            # Check that expected result files were produced
            self.verifyStepFiles(pickle.loads(str(stepRow["verifyFiles"])))
        except Exception as e:
            msg = "         Step finish with error: %(stepStr)s: %(e)s" % locals()
            printLog(msg, self.Log, out=True, err=True, isError=True)
            raise
        
        # Report step finish and update database
        msg = greenStr("         Step finish: %s" % stepStr)
        printLog(msg, self.Log, out=True, err=True)
        sqlCommand = """UPDATE %(TableSteps)s SET finish = CURRENT_TIMESTAMP 
                        WHERE step_id=%(step_id)d
                          AND run_id=%(run_id)d""" % myDict
        _cursor.execute(sqlCommand)
        _connection.commit()

    # Function to get the first avalaible gap to run 
    # it will return pair (state, stepRow)
    # if state is:
    # NO_MORE_GAPS, stepRow is None and there are not more gaps to work on
    # NO_AVAIL_GAP, stepRow is Nonew and not available gaps now, retry later
    # STEP_GAP, stepRow is a valid step row to work on
    def getStepGap(self, cursor):
        #The following query is just to start a transaction, pysqlite doesn't provide a better way
        cursor.execute("UPDATE %(TableSteps)s SET step_id = -1 WHERE step_id < 0" % self.sqlDict)
        if self.checkRunOk(cursor) and self.countStepGaps(cursor) > 0:
            sqlCommand = """ SELECT child.step_id, child.iter, child.passDb, child.command, child.parameters,child.verifyFiles 
                        FROM %(TableSteps)s parent, %(TableSteps)s child
                        WHERE (parent.step_id = child.parent_step_id) 
                          AND (child.step_id < %(next_step_id)d)
                          AND (child.init IS NULL)
                          AND (parent.finish IS NOT NULL)
                          AND (child.execution_mode = %(execution_gap)d)
                          AND (child.run_id=%(run_id)d)
                          AND (parent.run_id =%(run_id)d) 
                        LIMIT 1""" % self.sqlDict
            cursor.execute(sqlCommand)
            row = cursor.fetchone()
            if row is None:
                result = (NO_AVAIL_GAP, None)
            else:
                result = (STEP_GAP, row)
        else:
            result = (NO_MORE_GAPS, None)
        return result
    
    def countStepGaps(self, cursor):
        #Select first the next step id in main loop to limit the gaps
        sqlCommand = """SELECT COALESCE(MIN(step_id), 99999) 
                         FROM %(TableSteps)s 
                        WHERE run_id=%(run_id)d 
                          AND init IS NULL AND execution_mode > %(execution_gap)d""" % self.sqlDict
        cursor.execute(sqlCommand)
        # Count the gaps before the next step in main loop
        self.sqlDict['next_step_id'] = cursor.fetchone()[0]
        sqlCommand = """SELECT COUNT(*) FROM %(TableSteps)s 
                        WHERE (step_id < %(next_step_id)d)
                          AND (finish IS NULL)
                          AND (execution_mode = %(execution_gap)d)
                          AND (run_id=%(run_id)d) """ % self.sqlDict
        
        cursor.execute(sqlCommand)
        result =  cursor.fetchone()[0] 
        return result

# Function to fill gaps of step in database
#this will use mpi process
# this will be usefull for parallel processing, i.e., in threads or with MPI
# step_id is the step that will launch the process to fill previous gaps
def runStepGapsMpi(db, script, NumberOfMpi=1):
    if NumberOfMpi > 1:
        retcode = runJob(db.Log, "xmipp_steps_runner",  script, NumberOfMpi)
        if retcode != 0:
            raise Exception('xmipp_mpi_steps_runner execution failed')
    else:
        runThreadLoop(db, db.connection, db.cur) 
           
def runStepGaps(db, NumberOfThreads=1):
    # If run in separate threads, each one should create 
    # a new connection and cursor
    if NumberOfThreads > 1: 
        threads = []
        for i in range(NumberOfThreads):
            thr = ThreadStepGap(db)
            thr.start()
            threads.append(thr)
        #wait until all threads finish
        for thr in threads:
            thr.join() 
    else:
        runThreadLoop(db, db.connection, db.cur)
        
def runThreadLoop(db, connection, cursor):
    import time
    counter = 0
    while True:
        if counter > 100:
            reportError("runThreadLoop: more than 100 times here, this is probably a bug")
        state, stepRow = db.getStepGap(cursor)
        counter += 1
        if state == STEP_GAP: #database will be unlocked after commit on init timestamp
            db.runSingleStep(connection, cursor, stepRow)
        else:
            connection.rollback() #unlock database
            if state == NO_AVAIL_GAP:
                time.sleep(1)
            else:
                break
 
# RunJobThread
from threading import Thread
class ThreadStepGap(Thread):
    def __init__(self, db):
        Thread.__init__(self)
        self.db = db
    def run(self):
        try:
            conn = sqlite.Connection(self.db.dbName)
            conn.row_factory = sqlite.Row
            cur = conn.cursor()
            runThreadLoop(self.db, conn, cur)
        except Exception, e:
            printLog("Stopping threads because of error %s"%e,self.db.Log,out=True,err=True,isError=True)
            self.db.updateRunState(SqliteDb.RUN_FAILED, cursor=cur, connection=conn)
            
def escapeStr(str):
    return "'%s'" % str.replace("'", "''") 
          
class ProgramDb():
    def __init__(self, dbName=None):
        if dbName is None:
            from protlib_xmipp import getProgramsDbName
            dbName = getProgramsDbName()
        self.dbName = dbName
        self.connection = sqlite.Connection(dbName)
        self.connection.row_factory = sqlite.Row
        self.cursor = self.connection.cursor()
            
    def create(self):
        self.createTables()
        
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
        sqlCommand = "SELECT * FROM Category;"""
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def updateProgramCategory(self, program_name, category):
        sqlCommand = "UPDATE Program SET category_id = %d WHERE name='%s';" % (category['id'], program_name)
        self.cursor.execute(sqlCommand)
        self.connection.commit()        
