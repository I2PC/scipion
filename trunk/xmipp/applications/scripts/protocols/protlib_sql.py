from pysqlite2 import dbapi2 as sqlite
import os,sys
import pickle
from config import *
from protlib_utils import *
from protlib_filesystem import *

runColumns = ['run_id',
              'run_name',
              'script',
              'init',
              'last_modfied',
              'protocol_name',
              'comment']

class XmippProtocolDbStruct(object):
    doAlways = 99999
    lastStep  = -1
    firstStep =  0

def existsDB(dbName):
    """check if database has been created by checking if the table tableruns exist"""
    result = False
    try:
        connection = sqlite.Connection(dbName)
        connection.row_factory = sqlite.Row
        cur = connection.cursor()
        _sqlCommand = "SELECT count(*) from sqlite_master where tbl_name = '%(TableRuns)s';"%projectDefaults
        cur.execute(_sqlCommand)
        result=True
    except sqlite.Error, e:
        print "database %s is missing"% dbName
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
    
    def getRunId(self, protName,runName):
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
            
            def createStepTable(tableName):            
                _sqlCommand = (" CREATE TABLE IF NOT EXISTS " + tableName +  """
                             (step_id INTEGER DEFAULT 0, -- primary key (weak entity)
                             command TEXT,               -- comment (NOT USED, DROP?)
                             parameters TEXT,            -- wrapper parameters
                             init DATE,                  -- process started at
                             finish DATE,                -- process finished at
                             verified BOOL DEFAULT 'f',  -- auxiliary column to mark verified files (NOT USED, DROP?)
                             fileNameList TEXT,          -- list with files to modify
                             iter INTEGER DEFAULT 1,     -- for iterative scripts, iteration number
                                                         -- useful to restart at iteration n
                             execute BOOL,               -- Should the script execute this wrapper or is there
                                                         -- an external program taht will run it
                             run_id         INTEGER REFERENCES        %(TableRuns)s(run_id)  ON DELETE CASCADE,
                                                         -- key that unify all processes belonging to a run 
                             parent_step_id INTEGER REFERENCES """ + tableName +"""(step_id) ON DELETE CASCADE,
                                                          -- parent_step_id step must be executed before 
                                                          -- step_id may be executed
                             PRIMARY KEY(step_id, run_id))""") % self.sqlDict
            
                self.execSqlCommand(_sqlCommand, ("Error creating " + tableName + " table: ") % self.sqlDict)
                
            createStepTable('%(TableSteps)s')
            createStepTable('%(TableStepsRestart)s')
            
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
        

        
    def getLastRunName(self, protName):
        self.sqlDict['protocol_name'] = protName
        _sqlCommand = """SELECT COALESCE(MAX(run_name), '%(RunsPrefix)s') AS run_name 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocols)s
                         WHERE protocol_name = '%(protocol_name)s'""" % self.sqlDict
        self.cur.execute(_sqlCommand) 
        return self.cur.fetchone()[0]     
        
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
        
    def selectRuns(self, groupName):
        self.sqlDict['group'] = groupName
        _sqlCommand = """SELECT * 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocolsGroups)s
                         WHERE group_name = '%(group)s'""" % self.sqlDict
        self.cur.execute(_sqlCommand) 
        return self.cur.fetchall()
    
class XmippProtocolDb(SqliteDb): 
    
    def __init__(self, dbName, continueAt, isIter, protocol):
        """Constructor of the Sqlite database
        dbName    -- The filename of the database, the full path will be created if not exists
        continueAt -- at wich point to continue
        isIter     -- if True continueAt refers to iteration, otherwise refers to one step
        """
        self.sqlDict = projectDefaults
        self.ContinueAtIteration = continueAt        
        self.dbName = dbName
        self.connection = sqlite.Connection(dbName)
        self.connection.row_factory = sqlite.Row
        self.cur                = self.connection.cursor()
        self.cur_aux            = self.connection.cursor()
        self.lastid = XmippProtocolDbStruct.firstStep
        self.iter = self.lastid + 1
        self.parentCase=XmippProtocolDbStruct.lastStep
        # print wrapper name
        self.PrintWrapperCommand=True
        # print wrapper parameters
        self.PrintWrapperParameters=False
        #project dir
        self.ProjDir="."
        #verify output files
        self.verify=True
        #show file veification
        self.viewVerifyedFiles=False
        #constant
        self.SystemFlavour = "None"
        #get run_id
        run_id = self.getRunId(protocol.Name, protocol.runName)
        if not run_id:
            reportError("Protocol run '%(run_name)s' has not been registered in project database" % self.sqlDict)
        self.sqlDict['run_id'] = run_id
                        
        #check if protocol has ben run previosluy (that is, there are finished steps)
        _sqlCommand = """ SELECT COUNT(*) 
                          FROM %(TableRuns)s NATURAL JOIN %(TableSteps)s
                          WHERE finish IS NOT NULL 
                            AND run_id = %(run_id)d""" % self.sqlDict
        
        self.cur.execute(_sqlCommand)
        stepNo = self.cur.fetchone()[0]
            
        self.createRestartTable = stepNo>0 and self.ContinueAtIteration != 1

        if self.createRestartTable:
            self.sqlDict['TableStepsCurrent'] = self.sqlDict['TableStepsRestart']
        else:
            self.sqlDict['TableStepsCurrent'] = self.sqlDict['TableSteps']
            
        _sqlCommand = 'DELETE FROM %(TableStepsCurrent)s WHERE run_id = %(run_id)d' % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error cleaning table: %(TableStepsCurrent)s" % self.sqlDict)
        #Auxiliary string to insert/UPDATE data
        self.sqlInsertcommand = """ INSERT INTO 
                                    %(TableStepsCurrent)s(command,parameters,iter,execute,run_id,parent_step_id)
                                     VALUES (?,?,?,?,?,?)""" % self.sqlDict
        
        self.sqlInsertVerify  = " UPDATE %(TableStepsCurrent)s SET fileNameList= ? WHERE step_id=?"% self.sqlDict
        #Calculate the step at which should starts
        self.setStartingStep(isIter)
        #set to null time in original table for step >= self.StartAtStepN
        self.sqlDict['step_id'] = self.StartAtStepN
        _sqlCommand = """UPDATE %(TableSteps)s SET finish=NULL 
                         WHERE step_id >= %(step_id)d
                           AND run_id = %(run_id)d""" % self.sqlDict
        self.execSqlCommand(_sqlCommand, "Error reseting finish date: ")

    # print wrapper name
    def setPrintWrapperCommand(self,value):
        self.PrintWrapperCommand=value
        
    # print wrapper parameters
    def setPrintWrapperParameters(self,value):
        self.PrintWrapperParameters=value

    def setVerify(self,value,viewVerifyedFiles):
        self.verify=value
        self.viewVerifyedFiles=viewVerifyedFiles
        
    def setStartingStep(self,isIter):
        #if(self.ContinueAtIteration==-1 and self.StartAtStepN>0):
        if(self.ContinueAtIteration>0 and not isIter):
            self.StartAtStepN = self.ContinueAtIteration
        elif (self.ContinueAtIteration > 0 and isIter):
            self.sqlDict['iter'] = self.ContinueAtIteration
            _sqlCommand = """SELECT MIN(step_id) 
                             FROM %(TableStepsCurrent)s
                             WHERE iter = %(iter)d
                               AND run_id = %(run_id)d"""  % self.sqlDict
            self.cur_aux.execute(_sqlCommand)
            self.StartAtStepN = self.cur_aux.fetchone()[0]
        elif (self.ContinueAtIteration < 0):
            self.StartAtStepN =self.getStartingStepVerify(isIter)
        else:
            raise Exception("self.ContinueAtIteration must be !=0")

    def getStartingStepVerify(self,isIter):
        _sqlCommand = """SELECT step_id, iter, command, fileNameList 
                        FROM %(TableSteps)s
                        WHERE finish IS NULL
                             AND fileNameList IS NOT NULL
                             AND run_id = %(run_id)d
                        ORDER BY id """ % self.sqlDict
        #print "getstart", sqlCommand
        self.cur_aux.execute(_sqlCommand)
                    
        def getMinId(row):
            if isIter:
                self.sqlDict['iter'] = row['iter']
                _sqlCommand= """ SELECT MIN(step_id)
                                 FROM %(TableSteps)s
                                 WHERE iter = %(iter)d
                                   AND run_id = %(run_id)d""" % self.sqlDict
                self.cur_aux.execute(_sqlCommand)
                return(self.cur_aux.fetchone()[0])
            else:
                return (row['step_id'])
            
        for row in self.cur_aux:
            _list = pickle.loads(str(row["fileNameList"]))
            for i in _list:
                return getMinId(row)
        return getMinId(row)
    
    def saveParameters(self, _log, SystemFlavour):
        """save a dictionary to an auxiliary table"""
        if self.SystemFlavour == SystemFlavour:
            return
        cur_aux = self.connection.cursor()
        sqlCommand = """DELETE FROM %(TableParams)s
                               WHERE run_id = %(run_id)d""" % self.sqlDict
        cur_aux.execute(sqlCommand)
        sqlCommand = """INSERT INTO %(TableParams)s(parameters, run_id) VALUES(?, %(run_id)d)"""% self.sqlDict
        self.SystemFlavour = SystemFlavour
        dict = { 
          'SystemFlavour':self.SystemFlavour
        }
        cur_aux.execute(sqlCommand, [pickle.dumps(dict, 0)])
        self.connection.commit()
        
    def loadParameters(self, _log):
        """load a dictionary from an auxiliary table"""
        sqlCommand = """ SELECT parameters FROM %(TableParams)s WHERE run_id = %(run_id)d """ % self.sqlDict
        try:
            self.cur_aux.execute(sqlCommand)
        except sqlite.Error, e:
            print "loadParameters: Can not access to parameters computed in previous iteration:", e.args[0]
            print "you may need to set ContinueAtIteration=1"
            exit(1)
        dict = pickle.loads(str(self.cur_aux.fetchone()[0]))
        print dict
        self.SystemFlavour=dict['SystemFlavour']

    def compareParameters (self):
        """return 0 if new execution of script (tableName2) is a subset of and old execution(tableName1)
        for those steps in with finish is not null. This is interesting for continue at iteration N"""
        
        _sqlCommand = """ SELECT count(*) FROM
                              (SELECT command,parameters 
                               FROM %(TableSteps)s 
                               WHERE finish IS NOT NULL 
                                 AND command <> 'self.saveParameters'
                                 AND run_id = %(run_id)d
                               
                               except
                               
                               SELECT command,parameters 
                               FROM %(TableStepsRestart)s 
                               WHERE run_id = %(run_id)d
                              )""" % self.sqlDict
        #cur = self.connection.cursor()
        self.cur_aux.execute(_sqlCommand)
        result = self.cur_aux.fetchone()[0]
        #if original table is not a subset of restart then return error, i.e. result !=0
        # else overwrite original with restart for all those values that
        # has finish set to null
        if(not result):#original table is a subset of restart
            _sqlCommand = """ DELETE FROM %(TableSteps)s WHERE finish IS NULL AND run_id = %(run_id)d;
                              INSERT INTO %(TableSteps)s 
                                 SELECT * 
                                 FROM %(TableStepsRestart)s 
                                 WHERE step_id > (SELECT MAX(step_id) FROM %(TableSteps)s WHERE run_id = %(run_id)d)
                                   AND run_id = %(run_id)d""" % self.sqlDict
        else:
            #do the query again and print result
            _sqlCommand =   """SELECT command,parameters 
                               FROM %(TableSteps)s 
                               WHERE finish IS NOT NULL and command <> 'self.saveParameters'
                                 AND run_id = %(run_id)d
                                 
                               except
                               
                               SELECT command, parameters 
                               FROM %(TableStepsRestart)s 
                               WHERE run_id = %(run_id)d""" % self.sqlDict

            self.cur_aux.execute(_sqlCommand)
            for i in self.cur_aux:
                print i
            
        self.cur_aux.executescript(_sqlCommand)
        self.connection.commit()
        return result

    def setParentDefault(self, value):
        self.parentCase = value
        
    def setIteration(self,iter):
        self.iter=iter
        
    def insertAction(self, command,
                           verifyfiles=None,
                           parent_step_id=None, 
                           execute=True,
                           **_Parameters):
        if not parent_step_id:
            if self.parentCase == XmippProtocolDbStruct.lastStep:
                parent_step_id=self.lastid
            elif self.parentCase == XmippProtocolDbStruct.firstStep:
                parent_step_id=XmippProtocolDbStruct.firstStep
        if not execute:
            execute=True
        parameters = pickle.dumps(_Parameters, 0)#Dict
        self.cur_aux = self.connection.cursor()
        #print self.sqlInsertcommand, [command, parameters, iter]
        self.cur_aux.execute(self.sqlInsertcommand, [command, parameters, self.iter,execute,
                                                     self.sqlDict['run_id'],parent_step_id])
        self.lastid = self.cur_aux.lastrowid
        if (verifyfiles):
            verifyfiles = pickle.dumps(verifyfiles, 0)#Dict
            self.cur_aux.execute(self.sqlInsertVerify, [verifyfiles,self.lastid])
        return self.lastid

    def runActions(self, _log, _import):
       
        #print "kk", bcolors.OKBLUE,"kk"
        import pprint
        exec(_import)
        #check if tableName and tablename_aux are identical if not abort
        if self.createRestartTable:
            if self.compareParameters():
                ##########################Restore original table from backup
                print "ERROR: Can not continue from old execution, parameters do not match. Relaunch execution from begining"
                exit(1)
        self.sqlDict['step_id'] = self.StartAtStepN
        self.sqlDict['iter'] = XmippProtocolDbStruct.doAlways
        sqlCommand = """ SELECT iter, step_id, command, parameters,fileNameList 
                        FROM %(TableSteps)s 
                        WHERE (step_id >= %(step_id)d OR iter = %(iter)d)
                          AND (run_id = %(run_id)d)
                       ORDER BY step_id """ % self.sqlDict
        self.cur.execute(sqlCommand)

        kommands = self.cur.fetchall()
        for row in kommands:
            self.sqlDict['step_id'] = row['step_id']
            self.sqlDict['iter'] = row['iter']
            command = row['command']
            dict = pickle.loads(str(row["parameters"]))
            if(self.PrintWrapperCommand):
                if iter == 99999:
                    self.sqlDict['siter'] = 'N/A'
                else:
                    self.sqlDict['siter'] = str(row['iter'])
                print bcolors.OKBLUE,"--------\nExecution of wrapper: %(step_id)d (iter=%(siter)s)" % self.sqlDict
                print bcolors.HEADER,(command.split())[-1],bcolors.ENDC
            #print in column format rather tahn in raw, easier to read 
            if(self.PrintWrapperParameters):
                pp = pprint.PrettyPrinter(indent=4,width=20)
                pp.pprint(dict)

            sqlCommand = "UPDATE %(TableSteps)s SET init = CURRENT_TIMESTAMP WHERE step_id=%(step_id)d" % self.sqlDict
            self.connection.execute(sqlCommand)
            self.connection.commit()
            
            exec ( command + '(_log, **dict)')
            
            if self.verify and row["fileNameList"]:
                _list =pickle.loads(str(row["fileNameList"]))
                for i in _list:
                    if not os.path.exists(i):
                        self.sqlDict['file'] = i
                        print "ERROR at  step: %(step_id)d, file %(file)s has not been created." % self.sqlDict
                        exit(1)
                    elif self.viewVerifyedFiles:
                        print "Verified file:", i
            
            sqlCommand = "UPDATE %(TableSteps)s SET finish = CURRENT_TIMESTAMP WHERE step_id=%(step_id)d" %  self.sqlDict
            self.connection.execute(sqlCommand)
            self.connection.commit()
            
            if self.PrintWrapperCommand:
                print "Wrapper step: %(step_id)d finished\n" % self.sqlDict
        print '********************************************************'
        print ' Protocol FINISHED'
