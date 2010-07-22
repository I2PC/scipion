/***************************************************************************
 *
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
#include "metadata_sql.h"
//#define DEBUG

//This is needed for static memory allocation
int MDSql::table_counter = 0;
sqlite3 *MDSql::db;
MDSqlStaticInit MDSql::initialization;
char *MDSql::errmsg;
const char *MDSql::zLeftover;
int MDSql::rc;


int MDSql::getUniqueId()
{
    //    if (table_counter == 0)
    //        sqlBegin();
    return ++table_counter;
}

MDSql::MDSql(MetaData *md)
{
    tableId = getUniqueId();
    myMd = md;
    myCache = new MDCache();
}

MDSql::~MDSql()
{
    delete myCache;
}

bool MDSql::createMd()
{
    return createTable(&(myMd->activeLabels));
}

bool MDSql::clearMd()
{
    return dropTable();
}

long int MDSql::addRow()
{
    //Fixme: this can be done in the constructor of MDCache only once
    sqlite3_stmt * &stmt = myCache->addRowStmt;
    //sqlite3_stmt * stmt = NULL;

    if (stmt == NULL)
    {
        std::stringstream ss;
        ss << "INSERT INTO " << tableName(tableId) << " DEFAULT VALUES;";
        rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    }
    rc = sqlite3_reset(stmt);
    long int id = -1;
    if (execSingleStmt(stmt))
        id = sqlite3_last_insert_rowid(db);

    //sqlite3_finalize(stmt);
    return id;
}

bool MDSql::addColumn(MDLabel column)
{
    std::stringstream ss;
    ss << "ALTER TABLE " << tableName(tableId)
    << " ADD COLUMN " << MDL::label2SqlColumn(column) <<";";
    return execSingleStmt(ss);
}

bool MDSql::setObjectValue(const int objId, const MDValue &value)
{
    bool r = true;
    MDLabel column = value.label;
    std::stringstream ss;
    //Check cached statements for setObjectValue
    sqlite3_stmt * &stmt = myCache->setValueCache[column];
    //sqlite3_stmt * stmt = NULL;
    if (stmt == NULL)//if not exists create the stmt
    {
        std::string sep = (MDL::isString(column) || MDL::isVector(column)) ? "'" : "";
        ss << "UPDATE " << tableName(tableId)
        << " SET " << MDL::label2Str(column) << "=? WHERE objID=?;";
        rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    }
    rc = sqlite3_reset(stmt);
    bindValue(stmt, 1, value);
    rc = sqlite3_bind_int(stmt, 2, objId);
    rc = sqlite3_step(stmt);
    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::setObjectValue: " << std::endl
        << "   " << ss.str() << std::endl
        <<"    code: " << rc << " error: " << sqlite3_errmsg(db) << std::endl;
        r = false;
    }
    return r;
}

bool MDSql::getObjectValue(const int objId, MDValue  &value)
{
    std::stringstream ss;
    MDLabel column = value.label;
    sqlite3_stmt * &stmt = myCache->getValueCache[column];
    //sqlite3_stmt * stmt = NULL;

    if (stmt == NULL)//prepare stmt if not exists
    {
        //std::cerr << "Creating cache " << ++count <<std::endl;
        ss << "SELECT " << MDL::label2Str(column)
        << " FROM " << tableName(tableId)
        << " WHERE objID=?";// << objId << ";";
        rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    }
#ifdef DEBUG

    std::cerr << "getObjectValue: " << ss.str() <<std::endl;
#endif

    rc = sqlite3_reset(stmt);
    rc = sqlite3_bind_int(stmt, 1, objId);
    rc = sqlite3_step(stmt);

    if (rc == SQLITE_ROW)
        //|| rc== SQLITE_DONE)
    {
        extractValue(stmt, 0, value);
        rc = sqlite3_step(stmt);
    }
    else
    {
        return false;
    }
    //not finalize now because cached
    //sqlite3_finalize(stmt);
    return true;
}

void MDSql::selectObjects(std::vector<long int> &objectsOut, const MDQuery *queryPtr)
{
    std::stringstream ss;
    sqlite3_stmt *stmt;
    objectsOut.clear();
    long int id;

    ss << "SELECT objID FROM " << tableName(tableId);
    if (queryPtr != NULL)
    {
        ss << queryPtr->whereString();
        ss << queryPtr->orderByString();
        ss << queryPtr->limitString();
    }
    rc = sqlite3_prepare(db, ss.str().c_str(), -1, &stmt, &zLeftover);
#ifdef DEBUG

    std::cerr << "selectObjects: " << ss.str() <<std::endl;
#endif

    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
    {
        objectsOut.push_back(sqlite3_column_int(stmt, 0));
    }
    rc = sqlite3_finalize(stmt);
}

long int MDSql::deleteObjects(const MDQuery *queryPtr)
{
    std::stringstream ss;
    ss << "DELETE FROM " << tableName(tableId);
    if (queryPtr != NULL)
        ss << queryPtr->whereString();

    if (execSingleStmt(ss))
    {
        return sqlite3_changes(db);
    }
    return 0;

}

long int MDSql::copyObjects(MetaData *mdPtrOut, const MDQuery *queryPtr) const
{
    copyObjects(mdPtrOut->myMDSql, queryPtr);
}

long int MDSql::copyObjects(MDSql * sqlOut, const MDQuery *queryPtr) const
{
    //NOTE: Is assumed that the destiny table has
    // the same columns that the source table, if not
    // the INSERT will fail
    std::stringstream ss, ss2;
    ss << "INSERT INTO " << tableName(sqlOut->tableId);
    //Add columns names to the insert and also to select
    //* couldn't be used because maybe are duplicated objID's
    std::string sep = " ";
    int size = myMd->activeLabels.size();

    for (int i = 0; i < size; i++)
    {
        ss2 << sep << MDL::label2Str( myMd->activeLabels[i]);
        sep = ", ";
    }
    ss << "(" << ss2.str() << ") SELECT " << ss2.str();
    ss << " FROM " << tableName(tableId);
    if (queryPtr != NULL)
    {
        ss << queryPtr->whereString();
        ss << queryPtr->orderByString();
        ss << queryPtr->limitString();
    }
    if (sqlOut->execSingleStmt(ss))
    {
        return sqlite3_changes(db);
    }
    return 0;
}

void MDSql::aggregateMd(MetaData *mdPtrOut,
                        const std::vector<AggregateOperation> &operations,
                        MDLabel operateLabel)
{
    std::stringstream ss;
    std::stringstream ss2;
    std::string aggregateStr = MDL::label2Str(mdPtrOut->activeLabels[0]);
    ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId)
    << "(" << aggregateStr;
    ss2 << aggregateStr;
    //Start iterating on second label, first is the
    //aggregating one
    for (int i = 0; i < operations.size(); i++)
    {
        ss << ", " << MDL::label2Str(mdPtrOut->activeLabels[i+1]);
        ss2 << ", " ;
        switch (operations[i])
        {
        case AGGR_COUNT:
            ss2 << "COUNT";
            break;
        case AGGR_MAX:
            ss2 << "MAX";
            break;
        case AGGR_MIN:
            ss2 << "MIN";
            break;
        case AGGR_SUM:
            ss2 << "SUM";
            break;
        case AGGR_AVG:
            ss2 << "AVG";
            break;
        default:
            REPORT_ERROR(-66, "Invalid aggregate operation.");
        }
        ss2 << "(" << MDL::label2Str(operateLabel)
        << ") AS " << MDL::label2Str(mdPtrOut->activeLabels[i+1]);
    }
    ss << ") SELECT " << ss2.str();
    ss << " FROM " << tableName(tableId)
    << " GROUP BY " << aggregateStr
    << " ORDER BY " << aggregateStr << ";";

    execSingleStmt(ss);
}

void MDSql::indexModify(const MDLabel column, bool create)
{
    std::stringstream ss;
    if (create)
    {
        ss << "CREATE INDEX IF NOT EXISTS " << MDL::label2Str(column) << "_INDEX "
        << " ON " << tableName(tableId) << " (" << MDL::label2Str(column) << ")";
    }
    else
    {
        ss << "DROP INDEX IF EXISTS " << MDL::label2Str(column) << "_INDEX ";
    }
    execSingleStmt(ss);
}

long int MDSql::firstRow()
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MIN(objID), -1) AS MDSQL_FIRST_ID FROM "
    << tableName(tableId) << ";";
    return execSingleIntStmt(ss);
}

long int MDSql::lastRow()
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MAX(objID), -1) AS MDSQL_LAST_ID FROM "
    << tableName(tableId) << ";";
    execSingleIntStmt(ss);
}

long int MDSql::nextRow(long int currentRow)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MIN(objID), -1) AS MDSQL_NEXT_ID FROM "
    << tableName(tableId)
    << " WHERE objID>" << currentRow << ";";
    return execSingleIntStmt(ss);
}

long int MDSql::previousRow(long int currentRow)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MAX(objID), -1) AS MDSQL_PREV_ID FROM "
    << tableName(tableId)
    << " WHERE objID<" << currentRow << ";";
    return execSingleIntStmt(ss);
}

int MDSql::columnMaxLength(MDLabel column)
{
    std::stringstream ss;
    ss << "SELECT MAX(COALESCE(LENGTH("<< MDL::label2Str(column)
    <<"), -1)) AS MDSQL_STRING_LENGTH FROM "
    << tableName(tableId) << ";";
    return execSingleIntStmt(ss);
}

void MDSql::setOperate(MetaData *mdPtrOut, MDLabel column, SetOperation operation)
{
    std::stringstream ss, ss2;
    bool execStmt = true;
    int size;
    std::string sep = " ";
    switch (operation)
    {
    case UNION:
        copyObjects(mdPtrOut->myMDSql);
        execStmt = false;
        break;
    case UNION_DISTINCT: //unionDistinct
        //Create string with columns list
        size = myMd->activeLabels.size();
        for (int i = 0; i < size; i++)
        {
            ss2 << sep << MDL::label2Str( myMd->activeLabels[i]);
            sep = ", ";
        }
        ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId)
        << " (" << ss2.str() << ")"
        << " SELECT " << ss2.str()
        << " FROM " << tableName(tableId)
        << " WHERE "<< MDL::label2Str(column)
        << " NOT IN (SELECT " << MDL::label2Str(column)
        << " FROM " << tableName(mdPtrOut->myMDSql->tableId) << ");";
        break;
    case INTERSECTION:
    case SUBSTRACTION:
        ss << "DELETE FROM " << tableName(mdPtrOut->myMDSql->tableId)
        << " WHERE " << MDL::label2Str(column);
        if (operation == INTERSECTION)
            ss << " NOT";
        ss << " IN (SELECT " << MDL::label2Str(column)
        << " FROM " << tableName(tableId) << ");";
        break;
    }
    if (execStmt)
        execSingleStmt(ss);
}

void MDSql::setOperate(const MetaData *mdInLeft, const MetaData *mdInRight, MDLabel column, SetOperation operation)
{
    std::stringstream ss, ss2, ss3;
    int size;
    std::string join_type = "", sep = "";

    switch (operation)
    {
    case INNER_JOIN:
        join_type = " INNER ";
        break;
    case LEFT_JOIN:
        join_type = " LEFT OUTER ";
        break;
    case OUTER_JOIN:
        join_type = " OUTER ";
        break;
    }
    MDSql temp(NULL);
    temp.createTable(&(mdInRight->activeLabels), false);
    mdInRight->myMDSql->copyObjects(&temp);
    size = myMd->activeLabels.size();
    for (int i = 0; i < size; i++)
    {
        ss2 << sep << MDL::label2Str( myMd->activeLabels[i]);
        ss3 << sep;
        if (column == myMd->activeLabels[i])
            ss3 << tableName(mdInLeft->myMDSql->tableId) << ".";
        ss3 << MDL::label2Str( myMd->activeLabels[i]);
        sep = ", ";
    }
    ss << "INSERT INTO " << tableName(tableId)
    << " (" << ss2.str() << ")"
    << " SELECT " << ss3.str()
    << " FROM " << tableName(mdInLeft->myMDSql->tableId)
    << join_type << " JOIN " << tableName(temp.tableId)
    << " ON " << tableName(mdInLeft->myMDSql->tableId) << "." << MDL::label2Str(column)
    << "=" << tableName(temp.tableId) << "." << MDL::label2Str(column) << ";";

    std::cerr << "JOIN STMT: " << ss.str() <<std::endl;
    execSingleStmt(ss);
}

bool MDSql::operate(const std::string &expression)
{
    std::stringstream ss;
    ss << "UPDATE " << tableName(tableId) << " SET " << expression;

    return execSingleStmt(ss);
}

void MDSql::dumpToFile(const FileName &fileName)
{
    sqlite3 *pTo;
    sqlite3_backup *pBackup;

    sqlCommitTrans();
    rc = sqlite3_open(fileName.c_str(), &pTo);
    if( rc==SQLITE_OK )
    {
        pBackup = sqlite3_backup_init(pTo, "main", db, "main");
        if( pBackup )
        {
            sqlite3_backup_step(pBackup, -1);
            sqlite3_backup_finish(pBackup);
        }
        rc = sqlite3_errcode(pTo);
    }
    else
        REPORT_ERROR(-55, "dumpToFile: error opening db file");
    sqlite3_close(pTo);
    sqlBeginTrans();
}

bool MDSql::sqlBegin()
{
    if (table_counter > 0)
        return true;
    //std::cerr << "entering sqlBegin" <<std::endl;
    rc = sqlite3_open("", &db);

    sqlite3_exec(db, "PRAGMA temp_store=MEMORY",NULL, NULL, &errmsg);
    sqlite3_exec(db, "PRAGMA synchronous=OFF",NULL, NULL, &errmsg);
    sqlite3_exec(db, "PRAGMA count_changes=OFF",NULL, NULL, &errmsg);
    sqlite3_exec(db, "PRAGMA page_size=4092",NULL, NULL, &errmsg);

    return sqlBeginTrans();
}

bool MDSql::sqlEnd()
{
    sqlCommitTrans();
    sqlite3_close(db);
    //std::cerr << "Database sucessfully closed." <<std::endl;
}

bool MDSql::sqlBeginTrans()
{
    if (sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errmsg) != SQLITE_OK)
    {
        std::cerr << "Couldn't begin transaction:  " << errmsg << std::endl;
        return false;
    }
    return true;
}

bool MDSql::sqlCommitTrans()
{
    char *errmsg;

    if (sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errmsg) != SQLITE_OK)
    {
        std::cerr << "Couldn't commit transaction:  " << errmsg << std::endl;
        return false;
    }
    return true;
}

bool MDSql::dropTable()
{
    std::stringstream ss;
    ss << "DROP TABLE IF EXISTS " << tableName(tableId) << ";";
    return execSingleStmt(ss);
}

bool MDSql::createTable(const std::vector<MDLabel> * labelsVector, bool withObjID)
{
    std::stringstream ss;
    ss << "CREATE TABLE " << tableName(tableId) << "(";
    std::string sep = "";
    if (withObjID)
    {
        ss << "objID INTEGER PRIMARY KEY ASC AUTOINCREMENT";
        sep = ", ";
    }
    if (labelsVector != NULL)
    {
        for (int i = 0; i < labelsVector->size(); i++)
        {
            ss << sep << MDL::label2SqlColumn(labelsVector->at(i));
            sep = ", ";
        }
    }
    ss << ");";
    execSingleStmt(ss);
}

void MDSql::prepareStmt(const std::stringstream &ss, sqlite3_stmt *stmt)
{
    const char * zLeftover;
    const char * str = ss.str().c_str();
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
}

bool MDSql::execSingleStmt(const std::stringstream &ss)
{
#ifdef DEBUG
    //std::cerr << "execSingleStmt, stmt: '" << stmtStr << "'" <<std::endl;
#endif

    sqlite3_stmt * stmt;
    rc = sqlite3_prepare(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    bool r = execSingleStmt(stmt, &ss);
    rc = sqlite3_finalize(stmt);
    return r;
}

bool MDSql::execSingleStmt(sqlite3_stmt * &stmt, const std::stringstream *ss)
{
    bool r = true;
    rc = sqlite3_step(stmt);


    //if (ss != NULL)
    //    char * sqlStr = new char[1024];
    //    sqlStr = sqlite3_sql(stmt);
    //            std::cerr << "   " << sqlStr << std::endl;
    //    delete sqlStr[];
    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::execSingleStmt: " << std::endl;
        std::cerr <<"    code: " << rc << " error: " << sqlite3_errmsg(db) << std::endl;
        r = false;
        exit(1);
    }
    return r;
}

long int MDSql::execSingleIntStmt(const std::stringstream &ss)
{
    sqlite3_stmt * stmt;
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    rc = sqlite3_step(stmt);
    long int result = sqlite3_column_int(stmt, 0);

    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::execSingleIntStmt: error executing statemente, code " << rc <<std::endl;
        result = -1;
    }
    rc = sqlite3_finalize(stmt);
    return result;
}

std::string MDSql::tableName(const int tableId) const
{
    std::stringstream ss;
    ss <<  "MDTable_" << tableId;
    return ss.str();
}

int MDSql::bindValue(sqlite3_stmt *stmt, const int position, const MDValue &valueIn)
{
    //First reset the statement
    //rc  = sqlite3_reset(stmt);
    //std::cerr << "rc after reset: " << rc <<std::endl;
    switch (valueIn.type)
    {
    case LABEL_BOOL: //bools are int in sqlite3
        return sqlite3_bind_int(stmt, position, valueIn.boolValue ? 1 : 0);
    case LABEL_INT:
        return sqlite3_bind_int(stmt, position, valueIn.intValue);
    case LABEL_LONG:
        return sqlite3_bind_int(stmt, position, valueIn.longintValue);
    case LABEL_DOUBLE:
        return sqlite3_bind_double(stmt, position, valueIn.doubleValue);
    case LABEL_STRING:
        return sqlite3_bind_text(stmt, position, valueIn.stringValue.c_str(), -1, NULL);
    case LABEL_VECTOR:
        return sqlite3_bind_text(stmt, position, valueIn.toString().c_str(), -1, NULL);
    }
}

int MDSql::extractValue(sqlite3_stmt *stmt, const int position, MDValue &valueOut)
{
    std::stringstream ss;
    switch (valueOut.type)
    {
    case LABEL_BOOL: //bools are int in sqlite3
        valueOut.boolValue = sqlite3_column_int(stmt, position) == 1;
        break;
    case LABEL_INT:
        valueOut.intValue = sqlite3_column_int(stmt, position);
        break;
    case LABEL_LONG:
        valueOut.longintValue = sqlite3_column_int(stmt, position);
        break;
    case LABEL_DOUBLE:
        valueOut.doubleValue = sqlite3_column_double(stmt, position);
        break;
    case LABEL_STRING:
        ss << sqlite3_column_text(stmt, position);
        valueOut.stringValue = ss.str();
        break;
    case LABEL_VECTOR:
        //FIXME: Now are stored as string in DB
        ss << sqlite3_column_text(stmt, position);
        valueOut.fromStream(ss);
        break;
    }
}

MDCache::MDCache()
{
    this->addRowStmt = NULL;
    this->iterStmt = NULL;
}

MDCache::~MDCache()
{
    //Clear cached statements
    std::map<MDLabel, sqlite3_stmt*>::iterator it;
    //FIXME: This is a bit dirty here...should be moved to MDSQl
    for (it = setValueCache.begin(); it != setValueCache.end(); it++)
        sqlite3_finalize(it->second);
    setValueCache.clear();

    for (it = getValueCache.begin(); it != getValueCache.end(); it++)
        sqlite3_finalize(it->second);
    getValueCache.clear();

    if (iterStmt != NULL)
        sqlite3_finalize(iterStmt);

    if (addRowStmt != NULL)
        sqlite3_finalize(iterStmt);
}

