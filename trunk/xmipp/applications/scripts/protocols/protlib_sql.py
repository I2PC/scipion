from pysqlite2 import dbapi2 as sqlite
import os,sys
import pickle
from  bcolors import *
from config import *
from protlib_utils import *
from protlib_filesystem import *

class XmippProtocolDbStruct(object):
    doAlways = 99999

        
class XmippProjectDb:
        
    def execSqlCommand(self, sqlCmd, errMsg):
        '''Helper function to execute sqlite commands'''
        try:
           self.cur.executescript(sqlCmd)
        except sqlite.Error, e:
            print errMsg, e
            if(e.args[0].find('database is locked') != -1):
                print 'consider deleting the database (%s)' % LogName
            sys.exit(1)  
        self.connection.commit()
        
    def __init__(self, dbName):
        self.dbName = dbName
        self.connection = sqlite.Connection(dbName)
        self.connection.row_factory = sqlite.Row
        self.cur = self.connection.cursor()
        self.defaults = projectDefaults
        
        _sqlCommand = '''CREATE TABLE IF NOT EXISTS %(TableGroups)s
             (group_name TEXT PRIMARY KEY);''' % self.defaults
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableGroups)s' table: " % self.defaults)

        _sqlCommand = '''CREATE TABLE IF NOT EXISTS %(TableProtocols)s
                     (protocol_name TEXT PRIMARY KEY);''' % self.defaults
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableProtocols)s' table: " % self.defaults)
        
        _sqlCommand = '''CREATE TABLE IF NOT EXISTS %(TableProtocolsGroups)s
                     (protocol_name TEXT,
                      group_name TEXT,
                      PRIMARY KEY(protocol_name, group_name));''' % self.defaults
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableProtocolsGroups)s' table: " % self.defaults)        
                
        _sqlCommand = '''CREATE TABLE IF NOT EXISTS %(TableRuns)s
                     (id INTEGER PRIMARY KEY AUTOINCREMENT,
                      run_name TEXT UNIQUE,
                      script TEXT,
                      init DATE, 
                      last_modified DATE,
                      protocol_name TEXT REFERENCES %(TableProtocols)s(protocol_name),
                      comment TEXT);''' % self.defaults
        self.execSqlCommand(_sqlCommand, "Error creating '%(TableRuns)s' table: " % self.defaults)
        
        _sqlCommand = '''CREATE TABLE IF NOT EXISTS %s
                     (step_id INTEGER DEFAULT 0,
                     command TEXT, 
                     parameters TEXT,
                     init DATE, 
                     finish DATE,
                     verified BOOL,
                     fileNameList TEXT,
                     iter INTEGER,
                     run_id INTEGER REFERENCES runs(id) ON DELETE CASCADE,
                     PRIMARY KEY(step_id, run_id))'''
        self.execSqlCommand((_sqlCommand % '%(TableSteps)s') % self.defaults, 
                            "Error creating '%(TableSteps)s' table: " % self.defaults)
        self.execSqlCommand((_sqlCommand % '%(TableStepsRestart)s') % self.defaults, 
                            "Error creating '%(TableStepsRestart)s' table: " % self.defaults)
        
        _sqlCommand = '''CREATE TABLE IF NOT EXISTS parameters 
                        (parameters TEXT,
                        run_id INTEGER REFERENCES runs(id) ON DELETE CASCADE);'''
        self.execSqlCommand(_sqlCommand, "Error creating 'parameters' table: ")
                
        _sqlCommand = '''CREATE TRIGGER IF NOT EXISTS increment_step_id 
                         AFTER INSERT ON steps FOR EACH ROW  
                         BEGIN 
                            UPDATE steps SET step_id = SELECT MAX(step_id) + 1 
                                                       FROM steps 
                                                       WHERE run_id = NEW.run_id
                            WHERE step_id = 0 AND run_id = NEW.run_id; 
                         END'''
        
    def insertGroup(self, groupName):
        self.defaults['group'] = groupName
        _sqlCommand = "INSERT INTO %(TableGroups)s VALUES('%(group)s');" % self.defaults
        self.cur.execute(_sqlCommand)
        
    def insertProtocol(self, groupName, protName):
        self.defaults['group'] = groupName
        self.defaults['protocol'] = protName
        #check if protocol exists
        _sqlCommand = "SELECT COUNT(*) FROM %(TableProtocols)s WHERE protocol_name = '%(protocol)s'" % self.defaults
        self.cur.execute(_sqlCommand)
        if self.cur.fetchone()[0] == 0:
            _sqlCommand = "INSERT INTO %(TableProtocols)s VALUES('%(protocol)s')" % self.defaults
            self.cur.execute(_sqlCommand)
        _sqlCommand = "INSERT INTO %(TableProtocolsGroups)s VALUES('%(protocol)s', '%(group)s')" % self.defaults
        self.cur.execute(_sqlCommand)
          
    def insertRun(self, run):#Name, script, comment=''):
        _sqlCommand = """INSERT INTO runs values(
                            NULL, 
                            '%(run_name)s', 
                            '%(script)s', 
                            datetime('now'), 
                            datetime('now'), 
                            '%(protocol_name)s',
                            '%(comment)s');"""  % run
        self.cur.execute(_sqlCommand)
        run['id'] = self.cur.lastrowid
        self.connection.commit()
        
    def getLastRunName(self, protName):
        self.defaults['protocol'] = protName
        _sqlCommand = """SELECT COALESCE(MAX(run_name), '%(RunsPrefix)s') AS run_name 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocols)s
                         WHERE protocol_name = '%(protocol)s'""" % self.defaults
        self.cur.execute(_sqlCommand) 
        return self.cur.fetchone()[0]     
        
    def updateRun(self, run):
        _sqlCommand = """UPDATE runs SET
                            last_modified = datetime('now'),
                            comment = %(comment)s
                        WHERE id = %(id)d"""  % run
                         
        self.execSqlCommand(_sqlCommand, "Error updating run: %s" % run.name)  
        
    def deleteRun(self, run):
        _sqlCommand = "DELETE FROM runs WHERE id = %d " % run.id 
        self.execSqlCommand(_sqlCommand, "Error deleting run: %s" % run.name)
        
    def selectRuns(self, groupName):
        self.defaults['group'] = groupName
        _sqlCommand = """SELECT * 
                         FROM %(TableRuns)s NATURAL JOIN %(TableProtocolsGroups)s
                         WHERE group_name = '%(group)s'""" % self.defaults
        self.cur.execute(_sqlCommand) 
        return self.cur.fetchall()
    
class XmippProtocolDb: 
    
    def __init__(self, dbName, tableName, continueAt, isIter, run_id):
        '''Constructor of the Sqlite database
        dbName    -- The filename of the database, the full path will be created if not exists
        tableName -- The name of the table to be used
        continueAt -- at wich point to continue
        isIter     -- if True continueAt refers to iteration, otherwise refers to one step
        '''
        self.tableInsertOriginal = 'steps'
        self.tableInsertRestart = 'steps_restart'
        self.ContinueAtIteration = continueAt        
        self.dbName = dbName
        self.run_id = run_id
        self.connection = sqlite.Connection(dbName)
        self.connection.row_factory = sqlite.Row

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
        
        #check if table already exists
        _sqlCommand = """ SELECT COUNT(*) 
                          FROM runs 
                          WHERE finish IS NOT NULL 
                            AND run_id = %d""" % run_id
                            
        self.cur                = self.connection.cursor()
        self.cur_aux            = self.connection.cursor()
        
        self.cur.execute(_sqlCommand, [self.tableInsertOriginal])
        self.createRestartTable = self.cur.fetchone()[0] == 1 and self.ContinueAtIteration != 1

        if self.createRestartTable:
            self.tableInsert = self.tableInsertRestart
        else:
            self.tableInsert = self.tableInsertOriginal
        _sqlCommand = 'DELETE FROM ' + self.tableInsert +\
                             'WHERE run_id = %d' % run_id
        self.execSqlCommand(_sqlCommand, "Error cleaning table: " + self.tableInsert)
        #Auxiliary string to insert/UPDATE data
        self.sqlInsertcommand = " INSERT INTO " + self.tableInsert +\
                                " (command,parameters,iter,run_id) VALUES (?,?,?,?)"
        self.sqlInsertVerify  = " UPDATE " + self.tableInsert + " set fileNameList= ? WHERE id=?"
        #Calculate the step at which should starts
        self.setStartingStep(isIter)
        #set to null time in original table for step >= self.StartAtStepN
        _sqlCommand = """UPDATE %s
                              SET finish=NULL WHERE id >= %d
                                                AND run_id = %d""" \
                                                % (self.tableInsertOriginal, self.StartAtStepN, run_id)
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
            _sqlCommand = """SELECT MIN(id) 
                             FROM %s
                             WHERE iter = %d
                               AND run_id = %d""" \
                          % (self.tableInsert, self.ContinueAtIteration, self.run_id)
            self.cur_aux.execute(_sqlCommand)
            self.StartAtStepN = self.cur_aux.fetchone()[0]
        elif (self.ContinueAtIteration < 0):
            self.StartAtStepN =self.getStartingStepVerify(isIter)
        else:
            raise Exception("self.ContinueAtIteration must be !=0")

    def getStartingStepVerify(self,isIter):
        _sqlCommand = '''SELECT id, iter, command, fileNameList 
                        FROM %s
                        WHERE finish IS NULL
                             AND fileNameList IS NOT NULL
                             AND run_id = %d
                       ORDER BY id ''' \
                       % (self.tableInsertOriginal, self.run_id)
        #print "getstart", sqlCommand
        self.cur_aux.execute(_sqlCommand)
                    
        def getMinId(row):
            if isIter:
                _sqlCommand= ''' SELECT MIN(id)
                                 FROM %s
                                 WHERE iter = %d
                                   AND run_id = %d'''\
                                 % (self.tableInsertOriginal, row['iter'], self.run_id)
                self.cur_aux.execute(_sqlCommand)
                return(self.cur_aux.fetchone()[0])
            else:
                return (row['id'])
            
        for row in self.cur_aux:
            _list = pickle.loads(str(row["fileNameList"]))
            for i in _list:
                return getMinId(row)
        return getMinId(row)
    
    def saveParameters(self, _log, SystemFlavour):
        '''save a dictionary to an auxiliary table'''
        if self.SystemFlavour == SystemFlavour:
            return
        cur_aux = self.connection.cursor()
        sqlCommand = '''DELETE FROM parameters
                               WHERE run_id = %d''' % self.run_id
        cur_aux.execute(sqlCommand)
        sqlCommand = '''INSERT INTO parameters(parameters, run_id) VALUES(?, ?)'''
        self.SystemFlavour = SystemFlavour
        dict = { 
          'SystemFlavour':self.SystemFlavour
        }
        cur_aux.execute(sqlCommand, [pickle.dumps(dict, 0)])
        self.connection.commit()
        
    def loadParameters(self, _log):
        '''load a dictionary from an auxiliary table'''
        sqlCommand = '''SELECT parameters FROM parameters WHERE run_id = %d''' % self.run_id
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
        '''return 0 if new execution of script (tableName2) is a subset of and old execution(tableName1)
        for those steps in with finish is not null. This is interesting for continue at iteration N'''
        
        sqlDict = {'id': self.run_id, 
                'tableOriginal': self.tableInsertOriginal, 
                'tableRestart': self.tableInsertRestart}
        
        _sqlCommand = '''SELECT count(*) FROM
                              (SELECT command,parameters 
                               FROM %(tableOriginal)s 
                               WHERE finish IS NOT NULL 
                                 AND command <> 'self.saveParameters'
                                 AND run_id = %(id)d
                               
                               except
                               
                               SELECT command,parameters 
                               FROM %(tableRestart)s 
                               WHERE run_id = %(id)d
                              )''' % sqlDict
        #cur = self.connection.cursor()
        self.cur_aux.execute(_sqlCommand)
        result = self.cur_aux.fetchone()[0]
        #if original table is not a subset of restart then return error, i.e. result !=0
        # else overwrite original with restart for all those values that
        # has finish set to null
        if(not result):#original table is a subset of restart
            _sqlCommand = ''' DELETE FROM %(tableOriginal)s WHERE finish IS NULL AND run_id = %(id)d;
                              INSERT INTO %(tableOriginal)s 
                                 SELECT * 
                                 FROM %(tableRestart)s 
                                 WHERE id > (SELECT max(id) FROM %(tableOriginal)s WHERE run_id = %(id)d
                                   AND run_id = %(id)d''' % sqlDict
        else:
            #do the query again and print result
            _sqlCommand =   '''SELECT command,parameters 
                               FROM %(tableOriginal)s 
                               WHERE finish IS NOT NULL and command <> 'self.saveParameters'
                                 AND run_id = %(id)d
                                 
                               except
                               
                               SELECT command, parameters 
                               FROM %(tableRestart)s 
                               WHERE run_id = %(id)d''' % sqlDict

            self.cur_aux.execute(_sqlCommand)
            for i in self.cur_aux:
                print i
            
        self.cur_aux.executescript(_sqlCommand)
        self.connection.commit()
        return result

    def insertAction(self, command, _Parameters,iter,verifyfiles=None):
        parameters = pickle.dumps(_Parameters, 0)#Dict
        self.cur_aux = self.connection.cursor()
        #print self.sqlInsertcommand, [command, parameters, iter]
        self.cur_aux.execute(self.sqlInsertcommand, [command, parameters, iter])
        lastid = self.cur_aux.lastrowid
        if (verifyfiles):
            verifyfiles = pickle.dumps(verifyfiles, 0)#Dict
            self.cur_aux.execute(self.sqlInsertVerify, [verifyfiles,lastid])

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
        sqlCommand = '''SELECT iter, id, command, parameters,fileNameList 
                        FROM %s 
                        WHERE (id >= %d OR iter = %d)
                          AND (run_id = %d)
                       ORDER BY id''' \
                       % (self.tableInsertOriginal, self.StartAtStepN, XmippProtocolDbStruct.doAlways, self.run_id)
        self.cur.execute(sqlCommand)

        kommands = self.cur.fetchall()
        for row in kommands:
            id=row['id']
            iter=row['iter']
            command=row['command']
            dict = pickle.loads(str(row["parameters"]))
            if(self.PrintWrapperCommand):
                if iter == 99999:
                    siter='N/A'
                else:
                    siter=str(iter)
                print bcolors.OKBLUE,"--------\nExecution of wrapper: %d (iter=%s)"%(id,siter),
                print bcolors.HEADER,(command.split())[-1],bcolors.ENDC
            #print in column format rather tahn in raw, easier to read 
            if(self.PrintWrapperParameters):
                pp = pprint.PrettyPrinter(indent=4,width=20)
                pp.pprint(dict)

#                for i in dict.items():
#                    print i

            sqlCommand = "UPDATE %s set init = CURRENT_TIMESTAMP WHERE id=%d" % (self.tableInsertOriginal, id)
            self.connection.execute(sqlCommand)
            print 'command: ', row["command"]
            exec (row["command"] + '(_log, **dict)')
            if(self.verify and row["fileNameList"]):
                _list =pickle.loads(str(row["fileNameList"]))
                for i in _list:
                    if not os.path.exists(i):
                        print "ERROR at  step: %d, file %s has not been created." % (id, i)
                        exit(1)
                    elif self.viewVerifyedFiles:
                        print "Verified file:", i
            sqlCommand = "UPDATE %s set finish = CURRENT_TIMESTAMP WHERE id=%d" % (self.tableInsertOriginal, id)
            self.connection.execute(sqlCommand)
            if(self.PrintWrapperCommand):
                print "Wrapper step: %d finished\n" % id
#        self.cur.execute(sqlCommand)
            self.connection.commit()
        print '********************************************************'
        print ' Protocol FINISHED'
