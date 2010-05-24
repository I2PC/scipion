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
#define DEBUG

//This is needed for static memory allocation
int MDSql::table_counter = 0;
sqlite3 *MDSql::db;
MDSqlStaticInit MDSql::initialization;
char *MDSql::errmsg;
const char *MDSql::zLeftover;
int MDSql::rc;

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

bool MDSql::copyMd(const MetaData *mdPtrIn, MetaData *mdPtrOut)
{}

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

    std::cerr << "setObjectValue: " << ss.str() <<std::endl;
    return id;
}

bool MDSql::addColumn(const MetaData *mdPtr, MDLabel column)
{
    std::stringstream ss;
    ss << "ALTER TABLE " << tableName(mdPtr->tableId)
    << " ADD COLUMN " << MDL::label2SqlColumn(column);
#ifdef DEBUG
    std::cerr << "addColumn: "<< ss.str() <<std::endl;
#endif
    execSingleStmt(ss.str());
}

bool MDSql::setObjectValue(const MetaData *mdPtr, const int objId, const MDLabel column, const MDValue &value)
{
    std::string sep = (MDL::isString(column)) ? "'" : "";

    std::stringstream ss;
    ss << "UPDATE " << tableName(mdPtr->tableId)
    << " SET " << MDL::label2Str(column) << "="
    << sep << MDL::value2Str(column, value)<< sep << " WHERE objID=" << objId;
#ifdef DEBUG
    std::cerr << "setObjectValue: " << ss.str() <<std::endl;
#endif
    execSingleStmt(ss.str());

    std::cerr << "setObjectValue: " << ss.str() <<std::endl;
}

bool MDSql::getObjectValue(const MetaData *mdPtr, const int objId, const MDLabel column, MDValue  &value)
{
    std::stringstream ss;
    sqlite3_stmt *stmt;

    ss << "SELECT " << MDL::label2Str(column)
    << " FROM " << tableName(mdPtr->tableId)
    << " WHERE objID=" << objId << ";";
#ifdef DEBUG
    std::cerr << "getObjectValue: " << ss.str() <<std::endl;
#endif
    rc = sqlite3_prepare(db, ss.str().c_str(), -1, &stmt, &zLeftover);

    rc = sqlite3_step(stmt);
    char * strPtr  = NULL;

    if (rc == SQLITE_ROW)
    {
        switch (MDL::labelType(column))
        {
        case LABEL_BOOL: //bools are int in sqlite3
            value.boolValue = sqlite3_column_int(stmt, 0) == 1;
            break;
        case LABEL_INT:
            value.intValue = sqlite3_column_int(stmt, 0);
            break;
        case LABEL_DOUBLE:
            value.doubleValue = sqlite3_column_double(stmt, 0);
            break;
        case LABEL_STRING:
            ss.str("");
            ss << sqlite3_column_text(stmt, 0);
            value.stringValue = ss.str();
            break;
        case LABEL_VECTOR:
            REPORT_ERROR(-55, "Metadata still not suport vectors");
            break;
        }
    }
    else
    {
        REPORT_ERROR(-1, sqlite3_errmsg(db));
    }
    rc = sqlite3_finalize(stmt);

    return false;
}

std::vector<long int> MDSql::selectObjects(const MetaData *mdPtr, int maxObjects, const MDQuery *queryPtr)
{
    std::stringstream ss;
    sqlite3_stmt *stmt;
    std::vector<long int> objects;
    long int id;

    ss << "SELECT objID from " << tableName(mdPtr->tableId);
    if (queryPtr != NULL)
    {
        //todo: add query string to select
    }
    ss << " ORDER BY objID";
    ss << " LIMIT " << maxObjects;

    rc = sqlite3_prepare(db, ss.str().c_str(), -1, &stmt, &zLeftover);

    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
    {
        objects.push_back(sqlite3_column_int(stmt, 0));
    }
    rc = sqlite3_finalize(stmt);

    //FIXME: NOW A COPY OF THE VECTOR IS DONE!!!!
    return objects;
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

    if (sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errmsg) != SQLITE_OK)
    {
        std::cerr << "Couldn't begin transaction:  " << errmsg << std::endl;
        return false;
    }

    return true;
}

bool MDSql::sqlEnd()
{
    sqlCommit();
    sqlite3_close(db);

    std::cerr << "Database sucessfully closed." <<std::endl;
}

bool MDSql::sqlCommit()
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
    ss << "DROP TABLE IF EXISTS " << tableName(mdId);
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
    std::cerr << ss.str() <<std::endl;
    execSingleStmt(ss.str());
}

bool MDSql::execSingleStmt(const std::string &stmtStr)
{
#ifdef DEBUG
    std::cerr << "execSingleStmt, stmt: '" << stmtStr << "'" <<std::endl;
#endif

    bool r = false;
    sqlite3_stmt *stmt;
    rc = sqlite3_prepare(db, stmtStr.c_str(), -1, &stmt, &zLeftover);
    rc = sqlite3_step(stmt);
#ifdef DEBUG

    std::cerr << "execSingleStmt, return code: " << rc <<std::endl;
#endif

    if (rc == SQLITE_MISUSE)
        std::cerr << "misuse: " << sqlite3_errmsg(db) << std::endl;
    rc = sqlite3_finalize(stmt);
    return r;
}

std::string MDSql::tableName(const int tableId)
{
    std::stringstream ss;
    ss <<  "MDTable_" << tableId;
    return ss.str();
}

