from pysqlite2 import dbapi2 as sqlite
import os,sys
import pickle
from  bcolors import *
#
#

from protlib_utils import *
from protlib_filesystem import *

class XmippProtocolDbStruct(object):
    ''' class to mimic and structure for verification'''
    doAlways = 99999
#    def __init__(self):
#        a = 0

class XmippProtocolDb:
 
    # print wrapper name
    PrintWrapperCommand=True
    # print wrapper parameters
    PrintWrapperParameters=False
    #project dir
    ProjDir="."
    #verify output files
    verify=True
    #show file veification
    viewVerifyedFiles=False
    #constant
    SystemFlavour = "None"
    
    def __init__(self, dbName, tableName, continueAt, isIter):
        '''Constructor of the Sqlite database
        dbName    -- The filename of the database, the full path will be created if not exists
        tableName -- The name of the table to be used
        continueAt -- at wich point to continue
        isIter     -- if True continueAt refers to iteration, otherwise refers to one step
        '''
        
        self.tableInsertOriginal = tableName
        self.tableInsertRestart = tableName + "Restart"
        self.ContinueAtIteration = continueAt        
        self.tableVerify = tableName + "Verify"
        self.dbName = dbName

        self.connection = sqlite.Connection(dbName)
        self.connection.row_factory = sqlite.Row
        #check if table already exists
        _sqlCommand = "SELECT count(*) from sqlite_master where tbl_name = ?;"
        self.cur = self.connection.cursor()
        self.cur_aux = self.connection.cursor()
        self.cur.execute(_sqlCommand, [self.tableInsertOriginal])
        self.createRestartTable = self.cur.fetchone()[0] == 1 and self.ContinueAtIteration != 1

        if self.createRestartTable:
            self.tableInsert = self.tableInsertRestart
        else:
            self.tableInsert = self.tableInsertOriginal
        _sqlCommand = '''create table if not exists ''' + self.tableInsert + '''
                     (id INTEGER PRIMARY KEY,
                     command text, 
                     parameters text,
                     init date, 
                     finish date,
                     verified bool,
                     fileNameList text,
                     iter int)
                     ;delete from '''  + self.tableInsert
        try:
            self.cur.executescript(_sqlCommand)
        except sqlite.Error, e:
            print "kk",e
            if(e.args[0].find('database is locked')!= -1):
                print 'consider deleting the database (',LogName,')'
            sys.exit(1)


        self.cur.executescript(_sqlCommand)
        self.sqlInsertcommand = " insert into " + self.tableInsert + " (command,parameters,iter)             VALUES (?,?,?)"
        self.sqlInsertVerify = "update " + self.tableInsert + " set fileNameList= ? where id=?"
        
        #Calculate the step at which should starts
        self.setStartingStep(isIter)

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
            sqlString="select min(id) from " + self.tableInsert + " where iter=" +str(self.ContinueAtIteration)
            self.cur_aux.execute(sqlString)
            self.StartAtStepN = self.cur_aux.fetchone()[0]
        elif (self.ContinueAtIteration<0):
            self.StartAtStepN =self.getStartingStepVerify(isIter)
        else:
            raise Exception("self.ContinueAtIteration must be !=0")

    def getStartingStepVerify(self,isIter):
        sqlCommand = '''select id, iter, command, fileNameList 
                       FROM ''' + self.tableInsertOriginal + \
                   ''' where finish IS NULL
                             AND fileNameList IS NOT NULL
                       ORDER BY id '''
        #print "getstart", sqlCommand
        self.cur_aux.execute(sqlCommand)
        for row in self.cur_aux:
            _list = pickle.loads(str(row["fileNameList"]))
            for i in _list:
#                print i, row['iter'],row['id']
                if not os.path.exists(i):
                    if(isIter):
                        sqlCommand=''' select min(id) from 
                        ''' + self.tableInsertOriginal + '''
                        WHERE iter=''' + str(row['iter'])
                        self.cur_aux.execute(sqlCommand)
                        return(self.cur_aux.fetchone()[0])
                    else:
#                        print "return",row['id']
                        return (row['id'])
        if(isIter):
            sqlCommand=''' select min(id) from 
            ''' + self.tableInsertOriginal + '''
            WHERE iter=''' + str(row['iter'])
        else:
            sqlCommand=''' select max(id) from 
            ''' + self.tableInsertOriginal 
        self.cur_aux.execute(sqlCommand)
        return(self.cur_aux.fetchone()[0])
    
    def saveParameters(self, _log, SystemFlavour):
        '''save a dictionary to an auxiliary table'''
        if self.SystemFlavour == SystemFlavour:
            return
        cur_aux = self.connection.cursor()
        sqlCommand = '''CREATE TABLE if not exists parameters (parameters text);
                        DELETE FROM parameters;'''
        cur_aux.executescript(sqlCommand,)
        sqlCommand = '''INSERT into parameters(parameters) VALUES(?)'''
        self.SystemFlavour = SystemFlavour
        dict = { 
          'SystemFlavour':self.SystemFlavour
        }
        cur_aux.execute(sqlCommand, [pickle.dumps(dict, 0)])
        self.commit()
        
    def loadParameters(self, _log):
        '''load a dictionary from an auxiliary table'''
        sqlCommand = '''SELECT parameters FROM parameters'''
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

        _sqlCommand = '''SELECT count(*) FROM
                              (SELECT command,parameters 
                               FROM ''' + self.tableInsertOriginal + ''' 
                               WHERE finish IS NOT NULL and command <> 'self.saveParameters'
                               
                               except
                               
                               SELECT command,parameters 
                               FROM ''' + self.tableInsertRestart + ''' 
                              )'''
        #cur = self.connection.cursor()
        self.cur_aux.execute(_sqlCommand)
        result = self.cur_aux.fetchone()[0]
        #if original table is not a subset of restart then return error, i.e. result !=0
        # else overwrite original with restart for all those values that
        # has finish set to null
        if(not result):#original table is a subset of restart
            _sqlCommand = ''' delete from ''' + self.tableInsertOriginal + ''' where finish IS NULL'''
            _sqlCommand += ''' ; INSERT into ''' + self.tableInsertOriginal + '''
                                 SELECT * 
                                 FROM ''' + self.tableInsertRestart + '''
                                 WHERE id > (SELECT max(id) FROM ''' + self.tableInsertOriginal + ')'
        else:
            #do the query again and print result
            _sqlCommand =   '''SELECT command,parameters 
                               FROM ''' + self.tableInsertOriginal + ''' 
                               WHERE finish IS NOT NULL and command <> 'self.saveParameters'
                               
                               except
                               
                               SELECT command,parameters 
                               FROM ''' + self.tableInsertRestart + ''' 
                             '''
            #cur = self.connection.cursor()
            self.cur_aux.execute(_sqlCommand)
            for i in self.cur_aux:
                print i
            
        self.cur_aux.executescript(_sqlCommand)
        self.commit()
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

    def commit(self):
        self.connection.commit()
        
    def runActions(self, _log, _import):
       
        #print "kk", bcolors.OKBLUE,"kk"
        import pprint
        exec(_import)
        #check if tableName and tablename_aux are identical if not abort
        if self.createRestartTable:
            if self.compareParameters():
                ##########################Restore original table from backup
                print "ERROR: Can not continue from old execution, parameters do not match. Relunch execution from begining"
                exit(1)
        sqlCommand = '''SELECT iter,id, command, parameters,fileNameList 
                        FROM %s 
                        WHERE id >= %d OR 
                             iter=%d
                       ORDER BY id''' % (self.tableInsertOriginal, self.StartAtStepN, XmippProtocolDbStruct.doAlways)
        self.cur.execute(sqlCommand)

        kommands=self.cur.fetchall()
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

            sqlCommand = "update %s set init = CURRENT_TIMESTAMP where id=%d" % (self.tableInsertOriginal, id)
            self.connection.execute(sqlCommand)
            print 'command: ', row["command"]
            exec (row["command"] + '(_log, **dict)')
            if(self.verify and row["fileNameList"]):
                _list =pickle.loads(str(row["fileNameList"]))
                for i in _list:
                    if not os.path.exists(i):
                        print "ERROR at  step: %d, file %d has not been created." % (id, i)
                        exit(1)
                    elif self.viewVerifyedFiles:
                        print "Verified file:", i
            sqlCommand = "update %s set finish = CURRENT_TIMESTAMP where id=%d" % (self.tableInsertOriginal, id)
            self.connection.execute(sqlCommand)
            if(self.PrintWrapperCommand):
                print "Wrapper step: %d finished\n" % id
#        self.cur.execute(sqlCommand)
            self.commit()
        print '********************************************************'
        print ' Protocol FINISHED'
