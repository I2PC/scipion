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
sqlite3_stmt *MDSql::stmt;

int MDSql::getMdUniqueId()
{
    if (table_counter == 0)
        sqlBegin();

    return ++table_counter;
}

bool MDSql::createMd(const MetaData *mdPtr)
{
    int mdId = mdPtr->tableId;
    if (mdId < 0 || mdId > table_counter)
    {
        REPORT_ERROR(-1, "MDSql: invalid 'mdId' for create MD");
        return false;
    }
    return createTable(mdId, &(mdPtr->activeLabels));
}

bool MDSql::clearMd(const MetaData *mdPtr)
{
    //For now just drop the table
    int mdId = mdPtr->tableId;
    dropTable(mdId);
}

long int MDSql::addRow(const MetaData *mdPtr)
{
    std::stringstream ss;
    ss << "INSERT INTO " << tableName(mdPtr->tableId) << " VALUES(NULL";
    int size = mdPtr->activeLabels.size();
    for (int i = 0; i < size; i++)
        ss << ", NULL";
    ss <<");";
    execSingleStmt(ss.str());
    long int id = sqlite3_last_insert_rowid(db);
    return id;
}

bool MDSql::addColumn(const MetaData *mdPtr, MDLabel column)
{
    std::stringstream ss;
    ss << "ALTER TABLE " << tableName(mdPtr->tableId)
    << " ADD COLUMN " << MDL::label2SqlColumn(column) <<";";
    execSingleStmt(ss.str());
}

bool MDSql::setObjectValue(const MetaData *mdPtr, const int objId, const MDValue &value)
{
    MDLabel column = value.label;
    std::string sep = (MDL::isString(column) || MDL::isVector(column)) ? "'" : "";

    std::stringstream ss;
    ss << "UPDATE " << tableName(mdPtr->tableId)
    << " SET " << MDL::label2Str(column) << "="
    << sep << value.toString() << sep << " WHERE objID=" << objId << ";";
    execSingleStmt(ss.str());
}

bool MDSql::getObjectValue(const MetaData *mdPtr, const int objId, MDValue  &value)
{
    std::stringstream ss;
    sqlite3_stmt *stmt;
    MDLabel column = value.label;

    //if (mdPtr->getValueCache[value.label] == NULL)
    {
        ss << "SELECT " << MDL::label2Str(column)
        << " FROM " << tableName(mdPtr->tableId)
        << " WHERE objID=?";// << objId << ";";
    }
#ifdef DEBUG

    std::cerr << "getObjectValue: " << ss.str() <<std::endl;
#endif

    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    rc = sqlite3_reset(stmt);
    rc = sqlite3_bind_int(stmt, 1, objId);
    rc = sqlite3_step(stmt);

    if (rc == SQLITE_ROW || rc== SQLITE_DONE)
    {
        extractValue(stmt, 0, value);
    }
    else
    {
        REPORT_ERROR(-1, sqlite3_errmsg(db));
    }
    return false;
}

void MDSql::selectObjects(const MetaData *mdPtr, std::vector<long int> &objectsOut, int limit, const MDQuery *queryPtr)
{
    std::stringstream ss;
    sqlite3_stmt *stmt;
    objectsOut.clear();
    long int id;

    ss << "SELECT objID FROM " << tableName(mdPtr->tableId);
    if (queryPtr != NULL)
    {
        ss << " WHERE " << queryPtr->queryString;
    }
    ss << " ORDER BY objID";
    ss << " LIMIT " << limit << ";";

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

long int MDSql::deleteObjects(const MetaData *mdPtr, const MDQuery *queryPtr)
{
    std::stringstream ss;
    ss << "DELETE FROM " << tableName(mdPtr->tableId);
    if (queryPtr != NULL)
        ss << " WHERE " << queryPtr->queryString << ";";
    execSingleStmt(ss.str());
    return sqlite3_changes(db);
}

long int MDSql::copyObjects(const MetaData *mdPtrIn, MetaData *mdPtrOut,
                            const MDQuery *queryPtr, const MDLabel sortLabel,
                            int limit, int offset)
{
    //NOTE: Is assumed that the destiny table has
    // the same columns that the source table, if not
    // the INSERT will fail
    std::stringstream ss, ss2;
    ss << "INSERT INTO " << tableName(mdPtrOut->tableId);
    //Add columns names to the insert and also to select
    //* couldn't be used because maybe are duplicated objID's
    std::string sep = " ";
    int size = mdPtrIn->activeLabels.size();

    for (int i = 0; i < size; i++)
    {
        ss2 << sep << MDL::label2Str( mdPtrIn->activeLabels[i]);
        sep = ", ";
    }

    ss << "(" << ss2.str() << ") SELECT " << ss2.str();
    ss << " FROM " << tableName(mdPtrIn->tableId);
    if (queryPtr != NULL)
        ss << " WHERE " << queryPtr->queryString;
    ss << " ORDER BY " << MDL::label2Str(sortLabel);
    ss << " LIMIT " << limit << " OFFSET " << offset << ";";
    //std::cerr << "copyObjects: " << ss.str() <<std::endl;
    execSingleStmt(ss.str());
    //std::cerr << "executed copyObjects:------- " <<std::endl;
    return sqlite3_changes(db);
}

void MDSql::aggregateMd(const MetaData *mdPtrIn, MetaData *mdPtrOut,
                        const std::vector<AggregateOperation> &operations,
                        MDLabel operateLabel)
{
    std::stringstream ss;
    std::stringstream ss2;

    std::string aggregateStr = MDL::label2Str(mdPtrOut->activeLabels[0]);

    ss << "INSERT INTO " << tableName(mdPtrOut->tableId)
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
    ss << " FROM " << tableName(mdPtrIn->tableId)
    << " GROUP BY " << aggregateStr
    << " ORDER BY " << aggregateStr << ";";

    execSingleStmt(ss.str());
}

void MDSql::indexModify(const MetaData *mdPtr, const MDLabel column, bool create)
{
    std::stringstream ss;
    if (create)
    {
        ss << "CREATE INDEX IF NOT EXISTS " << MDL::label2Str(column) << "_INDEX "
        << " ON " << tableName(mdPtr->tableId) << " (" << MDL::label2Str(column) << ")";
    }
    else
    {
        ss << "DROP INDEX IF EXISTS " << MDL::label2Str(column) << "_INDEX ";
    }
    execSingleStmt(ss.str());
}

long int MDSql::firstRow(const MetaData *mdPtr)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MIN(objID), -1) AS MDSQL_FIRST_ID FROM "
    << tableName(mdPtr->tableId) << ";";
    execSingleStmt(ss.str());
    return sqlite3_column_int(stmt, 0);
}

long int MDSql::lastRow(const MetaData *mdPtr)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MAX(objID), -1) AS MDSQL_LAST_ID FROM "
    << tableName(mdPtr->tableId) << ";";
    execSingleStmt(ss.str());
    return sqlite3_column_int(stmt, 0);
}

long int MDSql::nextRow(const MetaData *mdPtr, long int currentRow)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MIN(objID), -1) AS MDSQL_NEXT_ID FROM "
    << tableName(mdPtr->tableId)
    << " WHERE objID>" << currentRow << ";";
    execSingleStmt(ss.str());
    return sqlite3_column_int(stmt, 0);
}

long int MDSql::previousRow(const MetaData *mdPtr, long int currentRow)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MAX(objID), -1) AS MDSQL_PREV_ID FROM "
    << tableName(mdPtr->tableId)
    << " WHERE objID<" << currentRow << ";";
    execSingleStmt(ss.str());
    return sqlite3_column_int(stmt, 0);
}

int MDSql::columnMaxLength(const MetaData *mdPtr, MDLabel column)
{
    std::stringstream ss;
    ss << "SELECT MAX(COALESCE(LENGTH("<< MDL::label2Str(column)
    <<"), -1)) AS MDSQL_STRING_LENGTH FROM "
    << tableName(mdPtr->tableId) << ";";
#ifdef DEBUG

    std::cerr << ss.str() <<std::endl;
#endif

    execSingleStmt(ss.str());
    return sqlite3_column_int(stmt, 0);
}

void MDSql::setOperate(const MetaData *mdPtrIn, MetaData *mdPtrOut, MDLabel column, int operation)
{
    std::stringstream ss, ss2;

    if (operation == 1) //unionDistinct
    {
        //Create string with columns list
        std::string sep = " ";
        int size = mdPtrIn->activeLabels.size();
        for (int i = 0; i < size; i++)
        {
            ss2 << sep << MDL::label2Str( mdPtrIn->activeLabels[i]);
            sep = ", ";
        }
        ss << "INSERT INTO " << tableName(mdPtrOut->tableId)
        << " (" << ss2.str() << ")"
        << " SELECT " << ss2.str()
        << " FROM " << tableName(mdPtrIn->tableId)
        << " WHERE "<< MDL::label2Str(column)
        << " NOT IN (SELECT " << MDL::label2Str(column)
        << " FROM " << tableName(mdPtrOut->tableId) << ");";
    }
    else //difference or intersecction
    {

        ss << "DELETE FROM " << tableName(mdPtrOut->tableId)
        << " WHERE " << MDL::label2Str(column);
        if (operation == 3)
            ss << " NOT";
        ss << " IN (SELECT " << MDL::label2Str(column)
        << " FROM " << tableName(mdPtrIn->tableId) << ");";
    }
    execSingleStmt(ss.str());
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
        std::cerr << "errcode: " << rc <<std::endl;

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
    std::cerr << "entering sqlBegin" <<std::endl;
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
    std::cerr << "Database sucessfully closed." <<std::endl;
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

bool MDSql::dropTable(const int mdId)
{
    std::stringstream ss;
    ss << "DROP TABLE IF EXISTS " << tableName(mdId) << ";";
    execSingleStmt(ss.str());
}

bool MDSql::createTable(const int mdId, const std::vector<MDLabel> * labelsVector)
{
    std::stringstream ss;
    ss << "CREATE TABLE " << tableName(mdId) <<
    "(objID INTEGER PRIMARY KEY ASC AUTOINCREMENT";


    if (labelsVector != NULL)
    {
        for (int i = 0; i < labelsVector->size(); i++)
        {
            ss << ", " << MDL::label2SqlColumn(labelsVector->at(i));
        }
    }
    ss << ");";

    execSingleStmt(ss.str());
}

int MDSql::execSingleStmt(const std::string &stmtStr)
{
#ifdef DEBUG
    std::cerr << "execSingleStmt, stmt: '" << stmtStr << "'" <<std::endl;
#endif

    sqlite3_stmt * _stmt;
//    rc = sqlite3_reset(_stmt);
    rc = sqlite3_prepare(db, stmtStr.c_str(), -1, &_stmt, &zLeftover);
    bool r = false;
    rc = sqlite3_step(_stmt);
#ifdef DEBUG2

    std::cerr << "execSingleStmt, return code: " << rc <<std::endl;
    //#endif

    if (rc == SQLITE_MISUSE)
        std::cerr << "misuse: " << sqlite3_errmsg(db) << std::endl;
#endif
    rc=sqlite3_finalize(_stmt);
    return rc;
}

int MDSql::execSingleStmt(sqlite3_stmt *stmt)
{
    bool r = false;
    rc = sqlite3_step(stmt);
#ifdef DEBUG

    std::cerr << "execSingleStmt, return code: " << rc <<std::endl;
    //#endif

    if (rc == SQLITE_MISUSE)
    {
        std::cerr << "misuse: " << sqlite3_errmsg(db) << std::endl;
    }
#endif
    return rc;
}

std::string MDSql::tableName(const int tableId)
{
    std::stringstream ss;
    ss <<  "MDTable_" << tableId;
    return ss.str();
}

int MDSql::bindValue(sqlite3_stmt *stmt, const int position, const MDValue &valueIn)
{
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

