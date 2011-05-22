from pysqlite2 import dbapi2 as sqlite
import os,sys
import pickle
#
#

class dataBaseStruct(object):
    ''' class to mimic and structure for verification'''
    doAlways = 99999
#    def __init__(self):
#        a = 0

class dataBase:
 
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
    notInitRadius = -1000
    def __init__(self, projectdir, logdir, scriptname, workDirectory, _tableName, _ContinueAtIteration):
        self.tableInsertOriginal = _tableName
        self.tableInsertRestart = _tableName + "Restart"
        self.ContinueAtIteration = _ContinueAtIteration
        self.tableVerify = _tableName + "Verify"
        self.OuterRadius = self.notInitRadius
        self.NumberOfCtfGroups =  1

        if logdir[0] == '/':
            LogName = logdir
        else:
            LogName = projectdir + '/' + logdir
        if not LogName[-1] == '/':
            LogName += '/'
        if not os.path.exists(LogName):
            os.makedirs(LogName)
        scriptname = os.path.basename(scriptname)
        LogName += scriptname.replace('.py', '')
        if not (workDirectory == "."):
            LogName += '_'
            LogName += os.path.basename(workDirectory)
        LogName += '.sqlite'
        self.connection = sqlite.Connection(LogName)
        self.connection.row_factory = sqlite.Row
        #check if table already exists
        _sqlCommand = "SELECT count(*) from sqlite_master where tbl_name = ?;"
        self.cur = self.connection.cursor()
        self.cur_aux = self.connection.cursor()
        self.cur.execute(_sqlCommand, [self.tableInsertOriginal])
        self.createRestartTable = self.cur.fetchone()[0] == 1 and _ContinueAtIteration != 1

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
                print 'consider deleting the database (placed at Log directory and having extension sqlite'
            sys.exit(1)


        self.cur.executescript(_sqlCommand)
        self.sqlInsertcommand = " insert into " + self.tableInsert + " (command,parameters,iter)             VALUES (?,?,?)"
        self.sqlInsertVerify = "update " + self.tableInsert + " set fileNameList= ? where id=?"

    # print wrapper name
    def setPrintWrapperCommand(self,value):
        self.PrintWrapperCommand=value
        
    # print wrapper parameters
    def setPrintWrapperParameters(self,value):
        self.PrintWrapperParameters=value
        
    def setProjDir(self,value):
        self.ProjDir=value
        os.chdir(self.ProjDir)

    def setVerify(self,value,viewVerifyedFiles):
        self.verify=value
        self.viewVerifyedFiles=viewVerifyedFiles
    def getStartingStep(self,isIter):
        #if(self.ContinueAtIteration==-1 and self.StartAtStepN>0):
        if(self.ContinueAtIteration>0 and not isIter):
            self.StartAtStepN= self.ContinueAtIteration
        elif (self.ContinueAtIteration > 0 and isIter):
            sqlString="Select min(id) from " + self.tableInsert + " where iter=" +str(self.ContinueAtIteration)
            self.cur_aux.execute(sqlString)
            self.StartAtStepN= self.cur_aux.fetchone()[0]
        elif (self.ContinueAtIteration<0):
            self.StartAtStepN=self.getStartingStepVerify(isIter)
        else:
            print "self.ContinueAtIteration must be !=0"
            exit(1)
        return self.StartAtStepN

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
    
    def saveParameters(self, _log,dict):
        '''save a dictionary to an auxiliary table'''
        if self.notInitRadius == self.OuterRadius:
            return
        cur_aux = self.connection.cursor()
        sqlCommand = '''CREATE TABLE if not exists parameters (parameters text);
                        DELETE FROM parameters;'''
        cur_aux.executescript(sqlCommand,)
        sqlCommand = '''INSERT into parameters(parameters) VALUES(?)'''
        dict = {
          'OuterRadius':self.OuterRadius
        , 'NumberOfCtfGroups':self.NumberOfCtfGroups
        }
        cur_aux.execute(sqlCommand, [pickle.dumps(dict, 0)])
        self.commit()
        
    def loadParameters(self,_log,dict):
        '''load a dictionary from an auxiliary table'''
        if self.notInitRadius != self.OuterRadius:
            return
        sqlCommand = '''SELECT parameters FROM parameters'''
        try:
            self.cur_aux.execute(sqlCommand)
        except sqlite.Error, e:
            print "loadParameters: Can not access to parameters computed in previous iteration:", e.args[0]
            print "you may need to set ContinueAtIteration=1"
            exit(1)
        dict = pickle.loads(str(self.cur_aux.fetchone()[0]))
        print dict
        self.OuterRadius=dict['OuterRadius']
        self.NumberOfCtfGroups=dict['NumberOfCtfGroups']

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

    def insertCommand(self, command, _Parameters,iter,verifyfiles=None):
        parameters = pickle.dumps(_Parameters, 0)#Dict
        self.cur_aux = self.connection.cursor()
        self.cur_aux.execute(self.sqlInsertcommand, [command, parameters, iter])
        lastid = self.cur_aux.lastrowid
        if (verifyfiles):
            verifyfiles = pickle.dumps(verifyfiles, 0)#Dict
            self.cur_aux.execute(self.sqlInsertVerify, [verifyfiles,lastid])

    def commit(self):
        self.connection.commit()
    def mainLoop(self, _log, stepNumber, _import):
        import pprint
        exec(_import)
        #check if tableName and tablename_aux are identical if not abort
        if self.createRestartTable:
            if self.compareParameters():
                ##########################Restore original table from backup
                print "ERROR: Can not continue from old execution, parameters do not match. Relunch execution from begining"
                exit(1)
        sqlCommand = '''SELECT iter,id, command, parameters,fileNameList 
                        FROM ''' + self.tableInsertOriginal + '''
                        WHERE id >= ''' + str(stepNumber) +''' OR 
                             iter=''' + str(dataBaseStruct.doAlways)+'''
                       ORDER BY id'''
        self.cur.execute(sqlCommand)

        kommands=self.cur.fetchall()
        for row in kommands:
            id=row['id']
            iter=row['iter']
            command=row['command']
            dict = pickle.loads(str(row["parameters"]))
            if(self.PrintWrapperCommand):
                print "--------\nExecution of wrapper: %d (iter=%d)"%(id,iter), (command.split())[-1]
            #print in column format rather tahn in raw, easier to read 
            if(self.PrintWrapperParameters):
                pp = pprint.PrettyPrinter(indent=4,width=20)
                pp.pprint(dict)

#                for i in dict.items():
#                    print i

            sqlCommand = "update " + self.tableInsertOriginal + " set init   = CURRENT_TIMESTAMP where id=%d" % id
            self.connection.execute(sqlCommand)
            print 'row["command"]
            exec (row["command"] command','+ '(_log, dict)')
            if(self.verify and row["fileNameList"]):
                _list =pickle.loads(str(row["fileNameList"]))
                for i in _list:
                    if not os.path.exists(i):
                        print "ERROR at  step: ", id, ", file", i, " has not been created."
                        exit(1)
                    elif self.viewVerifyedFiles:
                        print "Verified file:", i
            sqlCommand = "update " + self.tableInsertOriginal + " set finish   = CURRENT_TIMESTAMP where id=%d" % id
            self.connection.execute(sqlCommand)
            if(self.PrintWrapperCommand):
                print "Wrapper step: %d"%id, "finished\n"
#        self.cur.execute(sqlCommand)
            self.commit()
        print '********************************************************'
        print ' Protocol FINISHED'
