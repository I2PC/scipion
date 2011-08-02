from pysqlite2 import dbapi2 as sqlite
import pickle
import os, sys
from config_protocols import projectDefaults
from protlib_utils import reportError, getScriptPrefix, printLog, printLogError, runJob
from protlib_filesystem import createDir, deleteDir

runColumns = ['run_id',
              'run_name',
              'script',
              'init',
              'last_modfied',
              'protocol_name',
              'comment',
              'group_name']

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
                          script TEXT,    -- strip full name
                          init DATE,      -- run started at
                          last_modified DATE, --last modification (edition usually)
                          protocol_name TEXT REFERENCES %(TableProtocols)s(protocol_name), -- protocol name
                          comment TEXT, -- user defined comment
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
                         execute_mainloop BOOL,      -- Should the script execute this wrapper in main loop?
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
                            '%(script)s', 
                            datetime('now'), 
                            datetime('now'), 
                            '%(protocol_name)s',
                            '%(comment)s');"""  % self.sqlDict
        self.cur.execute(_sqlCommand)
        run['run_id'] = self.cur.lastrowid
        self.connection.commit()
        
    def suggestRunName(self, protName):
        self.sqlDict['protocol_name'] = protName
        _sqlCommand = """SELECT COALESCE(MAX(run_name), '%(RunsPrefix)s') AS run_name 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocols)s
                         WHERE protocol_name = '%(protocol_name)s'""" % self.sqlDict
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
                            last_modified = datetime('now'),
                            comment = '%(comment)s'
                        WHERE run_id = %(run_id)d"""  % self.sqlDict
                         
        self.execSqlCommand(_sqlCommand, "Error updating run: %(run_name)s" % run)  
        
    def deleteRun(self, run):
        self.sqlDict.update(run)
        _sqlCommand = "DELETE FROM %(TableRuns)s WHERE run_id = %(run_id)d " % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error deleting run: %(run_name)s" % run)
        #this should be working automatically with the ON DELETE CASCADE
        #but isn't working, so...
        _sqlCommand = "DELETE FROM %(TableSteps)s WHERE run_id = %(run_id)d " % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error deleting steps of run: %(run_name)s" % run)        
        
    def selectRunsCommand(self):
        sqlCommand = """SELECT run_id, run_name, script, 
                               datetime(init, 'localtime') as init, 
                               datetime(last_modified, 'localtime') as last_modified,
                               protocol_name, comment,
                               group_name 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocolsGroups)s """ % self.sqlDict
        return sqlCommand
                         
    def selectRuns(self, groupName):
        self.sqlDict['group'] = groupName
        sqlCommand = self.selectRunsCommand() + """WHERE group_name = '%(group)s'
                                                   ORDER BY last_modified DESC """ % self.sqlDict
        self.cur.execute(sqlCommand) 
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
     
class XmippProtocolDb(SqliteDb): 
    def __init__(self, protocol):
        self.runBehavior = protocol.Behavior
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

        # Restart or resume
        if self.runBehavior=="Restart":
            _sqlCommand = 'DELETE FROM %(TableSteps)s WHERE run_id = %(run_id)d' % self.sqlDict
            self.execSqlCommand(_sqlCommand, "Error cleaning table: %(TableSteps)s" % self.sqlDict)
            self.insertStatus=True
        else:
            sqlCommand = """ SELECT step_id, iter, command, parameters, verifyFiles 
                             FROM %(TableSteps)s 
                             WHERE (run_id = %(run_id)d)
                             ORDER BY step_id """ % self.sqlDict
            self.cur.execute(sqlCommand)            
            self.insertStatus=False

    def setIteration(self,iter):
        self.iter=iter
        
    def insertStep(self, command,
                           verifyfiles=[],
                           parent_step_id=None, 
                           execute_mainloop = True,
                           passDb=False,
                           **_Parameters):
        if not parent_step_id:
            parent_step_id=self.lastStepId

        parameters = pickle.dumps(_Parameters, 0)
        verifyfilesString = pickle.dumps(verifyfiles, 0)
        if self.insertStatus==False:
            row=self.cur.fetchone()
            if row is None:
                self.insertStatus=True
            else:
                if row['parameters']!=parameters or row['verifyFiles']!=verifyfilesString:
                    self.insertStatus=True
                else:
                    for file in verifyfiles:
                        if not os.path.exists(file):
                            self.insertStatus=True
                            break
                self.lastStepId=row['step_id']
                if self.insertStatus==True:
                    self.sqlDict['step_id'] = row['step_id']
                    self.cur.execute("""DELETE FROM %(TableSteps)s WHERE run_id = %(run_id)d AND step_id>=%(step_id)d""" % self.sqlDict)
                    self.connection.commit()
        if self.insertStatus==True:
            try:
                self.cur_aux.execute("""INSERT INTO 
                                    %(TableSteps)s(command,parameters,verifyFiles,iter,execute_mainloop,passDb,run_id,parent_step_id)
                                     VALUES (?,?,?,?,?,?,?,?)""" % self.sqlDict,
                                     [command, parameters, verifyfilesString, self.iter,execute_mainloop,passDb,self.sqlDict['run_id'], parent_step_id])
                #select the last step_id inserted for this run
                self.cur_aux.execute("""SELECT MAX(step_id) FROM %(TableSteps)s WHERE run_id = %(run_id)d""" % self.sqlDict)
                self.lastStepId = self.cur_aux.fetchone()[0]
            except sqlite.Error, e:
                reportError( "Cannot insert command: %s" % e.args[0])
        return self.lastStepId

    def runSteps(self):
        self.connection.commit()
        sqlCommand = """ SELECT step_id, iter, passDb, command, parameters, verifyFiles 
                         FROM %(TableSteps)s 
                         WHERE run_id = %(run_id)d 
                         AND finish IS NULL AND execute_mainloop = 1
                         ORDER BY step_id """ % self.sqlDict
        self.cur.execute(sqlCommand)
        commands = self.cur.fetchall()
        n = len(commands)
        printLog(self.Log,'***************************** Protocol STARTED mode: %s'%self.runBehavior,True)
        
        for i in range(n):
            if i < n - 1:
                self.nextStepId = commands[i+1]['step_id']
            else:
                self.nextStepId = XmippProjectDb.BIGGEST_STEP
            self.runSingleStep(self.connection, self.cur, commands[i])
        printLog(self.Log,'***************************** Protocol FINISHED',True)

    def runSingleStep(self, _connection, _cursor, stepRow):
        exec(self.Import)
        step_id = stepRow['step_id']
        myDict = self.sqlDict.copy()
        myDict['step_id'] = step_id
        myDict['iter'] = stepRow['iter']
        command = stepRow['command']
        dict = pickle.loads(str(stepRow["parameters"]))
        
        # Print
        import pprint
        from protlib_utils import blueStr, headerStr
        print blueStr("-------- Step: %d (iter=%d)" % (step_id,stepRow['iter']))
        print headerStr((command.split())[-1])
        pprint.PrettyPrinter(indent=4,width=20).pprint(dict)

        # Set init time
        sqlCommand = """UPDATE %(TableSteps)s SET init = CURRENT_TIMESTAMP 
                        WHERE step_id=%(step_id)d
                          AND run_id=%(run_id)d""" % myDict
        _cursor.execute(sqlCommand)
        _connection.commit()
        
        # Execute Python function
        if stepRow['passDb']:
            exec ( command + '(self, **dict)')
        else:
            exec ( command + '(self.Log, **dict)')
        
        # Check verify files
        fileList =pickle.loads(str(stepRow["verifyFiles"]))
        print "About to check", fileList
        missingFiles=False
        for file in fileList:
            print "Checking",file
            if not os.path.exists(file):
                myDict['file'] = file
                printLogError(self.Log, "ERROR at  step: %(step_id)d, file %(file)s has not been created." % myDict)
                missingFiles=True
        if missingFiles:
            exit(1)
     
        # Update finish
        sqlCommand = """UPDATE %(TableSteps)s SET finish = CURRENT_TIMESTAMP 
                        WHERE step_id=%(step_id)d
                          AND run_id=%(run_id)d""" % myDict
        _cursor.execute(sqlCommand)
        _connection.commit()
        
        print "Step %d finished\n" % step_id

    # Function to get the first avalaible gap to run 
    # it will return pair (state, stepRow)
    # if state is:
    # NO_MORE_GAPS, stepRow is None and there are not more gaps to work on
    # NO_AVAIL_GAP, stepRow is Nonew and not available gaps now, retry later
    # STEP_GAP, stepRow is a valid step row to work on
    def getStepGap(self, cursor):
        #The following query is just to start a transaction, pysqlite doesn't provide a better way
        cursor.execute("UPDATE %(TableSteps)s SET step_id = -1 WHERE step_id < 0" % self.sqlDict)
        if self.countStepGaps(cursor) > 0:
            sqlCommand = """ SELECT child.step_id, child.iter, child.passDb, child.command, child.parameters,child.verifyFiles 
                        FROM %(TableSteps)s parent, %(TableSteps)s child
                        WHERE (parent.step_id = child.parent_step_id) 
                          AND (child.step_id < %(step_id)d)
                          AND (child.init IS NULL)
                          AND (parent.finish IS NOT NULL)
                          AND (child.execute_mainloop =  0)
                          AND (child.run_id=%(run_id)d)
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
        self.sqlDict['step_id'] = self.nextStepId
        sqlCommand = """SELECT COUNT(*) FROM %(TableSteps)s 
                        WHERE (step_id < %(step_id)d) 
                          AND (finish IS NULL)
                          AND (execute_mainloop =  0)
                          AND (run_id=%(run_id)d) """ % self.sqlDict
        
        cursor.execute(sqlCommand)
        result =  cursor.fetchone()[0] 
        return result
# Function to fill gaps of step in database
# this will be usefull for parallel processing, i.e., in threads or with MPI
# step_id is the step that will launch the process to fill previous gaps
def runStepGaps(db, NumberOfThreads):
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
        conn = sqlite.Connection(self.db.dbName)
        conn.row_factory = sqlite.Row
        cur = conn.cursor()
        runThreadLoop(self.db, conn, cur)

