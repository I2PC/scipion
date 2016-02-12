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

#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include "metadata_sql.h"
#include "xmipp_threads.h"
#include <sys/time.h>
#include <regex.h>
//#define DEBUG

//This is needed for static memory allocation
int MDSql::table_counter = 0;
sqlite3 *MDSql::db;
MDSqlStaticInit MDSql::initialization;
char *MDSql::errmsg;
const char *MDSql::zLeftover;
int MDSql::rc;
Mutex sqlMutex; //Mutex to syncronize db access

std::stringstream MDSql::preparedStream;	// Stream.
sqlite3_stmt * MDSql::preparedStmt;

void sqlite_regexp(sqlite3_context* context, int argc, sqlite3_value** values) {
    int ret;
    regex_t regex;
    char* reg = (char*)sqlite3_value_text(values[0]);
    char* text = (char*)sqlite3_value_text(values[1]);

    if ( argc != 2 || reg == 0 || text == 0) {
        sqlite3_result_error(context, "SQL function regexp() called with invalid arguments.\n", -1);
        return;
    }

    ret = regcomp(&regex, reg, REG_EXTENDED | REG_NOSUB);
    if ( ret != 0 ) {
        sqlite3_result_error(context, "error compiling regular expression", -1);
        return;
    }

    ret = regexec(&regex, text , 0, NULL, 0);
    regfree(&regex);

    sqlite3_result_int(context, (ret != REG_NOMATCH));
}

int getBlocksInMetaDataFileDB(const FileName &inFile, StringVector& blockList)
{
    char **results;
    int rows;
    int columns;
    int rc;

    sqlite3 *db1;
    String sql = (String)"SELECT name FROM sqlite_master\
                 WHERE type='table';";
    if ((rc=sqlite3_open(inFile.c_str(), &db1)))
        REPORT_ERROR(ERR_MD_SQL,formatString("Error opening database code: %d message: %s",rc,sqlite3_errmsg(db1)));
    if ((rc=sqlite3_get_table (db1, sql.c_str(), &results, &rows, &columns, NULL)) != SQLITE_OK)
        REPORT_ERROR(ERR_MD_SQL,formatString("Error accessing table code: %d message: %s. SQL command %s",rc,sqlite3_errmsg(db1),sql.c_str()));
    //For tables, the type field will always be 'table' and the name field will
    //be the name of the table. So to get a list of all tables in the database,

    if (rows < 1)
    {
        std::cerr << "Empty Metadata" <<std::endl;
        return 0;
    }
    else
    {
        for (int i = 1; i <= rows; i++)
        {
            blockList.push_back((String)results[i]);
        }
    }
    sqlite3_free_table (results);
    sqlite3_close(db1);

    return rows;
}

int MDSql::getUniqueId()
{
    //    if (table_counter == 0)
    //        sqlBegin();
    return ++table_counter;
}

MDSql::MDSql(MetaData *md)
{
    sqlMutex.lock();
    tableId = getUniqueId();
    //std::cerr << ">>>> creating md with table id: " << tableId << std::endl;
    sqlMutex.unlock();
    myMd = md;
    myCache = new MDCache();

}

MDSql::~MDSql()
{
    delete myCache;
}

bool MDSql::createMd()
{
    sqlMutex.lock();
    //std::cerr << "creating md" <<std::endl;
    bool result = createTable(&(myMd->activeLabels));
    //std::cerr << "leave creating md" <<std::endl;
    sqlMutex.unlock();

    return result;
}

bool MDSql::clearMd()
{
    sqlMutex.lock();
    //std::cerr << "clearing md" <<std::endl;
    myCache->clear();
    bool result = dropTable();
    //std::cerr << "leave clearing md" <<std::endl;
    sqlMutex.unlock();

    return result;
}

size_t MDSql::addRow()
{
    //Fixme: this can be done in the constructor of MDCache only once
    sqlite3_stmt * &stmt = myCache->addRowStmt;
    //sqlite3_stmt * stmt = NULL;

    if (stmt == NULL)
    {
        std::stringstream ss;
        ss << "INSERT INTO " << tableName(tableId) << " DEFAULT VALUES;";
        rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
//#define DEBUG
#ifdef DEBUG
    std::cerr << "DEBUG_JM: addRow: " << ss.str() <<std::endl;
#endif
#undef DEBUG
    }
    rc = sqlite3_reset(stmt);
    size_t id = BAD_OBJID;
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

bool  MDSql::activateMathExtensions(void)
{
    const char* lib = "libXmippSqliteExt.so";
    sqlite3_enable_load_extension(db, 1);
    if( sqlite3_load_extension(db, lib, 0, 0)!= SQLITE_OK)
        REPORT_ERROR(ERR_MD_SQL,"Cannot activate sqlite extensions");
    else
        return true;
}

bool  MDSql::activateRegExtensions(void)
{
	if( sqlite3_create_function(db, "regexp", 2, SQLITE_ANY,0, &sqlite_regexp,0,0)!= SQLITE_OK)
        REPORT_ERROR(ERR_MD_SQL,"Cannot activate sqlite extensions");
    else
        return true;
}

bool  MDSql::deactivateThreadMuting(void)
{
    return ( sqlite3_config(SQLITE_CONFIG_MULTITHREAD) == SQLITE_OK);
}
bool  MDSql::activateThreadMuting(void)
{
    return ( sqlite3_config(SQLITE_CONFIG_SERIALIZED) == SQLITE_OK);
}


bool MDSql::renameColumn(const std::vector<MDLabel> oldLabel, const std::vector<MDLabel> newlabel)
{
    //1 Create an new table that matches your original table,
    // but with the changed columns.
    bool result;
    std::vector<MDLabel> v1(myMd->activeLabels);
    std::vector<MDLabel>::const_iterator itOld;
    std::vector<MDLabel>::const_iterator itNew;
    for( itOld = oldLabel.begin(), itNew = newlabel.begin();
         itOld < oldLabel.end();
         ++itOld, ++itNew )
        std::replace(v1.begin(), v1.end(), *itOld, *itNew);

    int oldTableId = tableId;
    sqlMutex.lock();
    tableId = getUniqueId();
    createTable(&v1);
    sqlMutex.unlock();
    //2 Now we can copy the original data to the new table:
    String oldLabelString=" objID";
    String newLabelString=" objID";
    for(std::vector<MDLabel>::const_iterator  it = (myMd->activeLabels)
            .begin();
        it != (myMd->activeLabels).end();
        ++it)
        oldLabelString += ", " + MDL::label2StrSql(*it);
    for(std::vector<MDLabel>
        ::const_iterator  it = v1.begin();
        it != v1.end();
        ++it)
        newLabelString += ", " + MDL::label2StrSql(*it);
    std::stringstream sqlCommand;
    sqlCommand << " INSERT INTO " + tableName(tableId)
    << " ("+ newLabelString +") "
    << " SELECT " + oldLabelString
    << " FROM " + tableName(oldTableId) ;
    execSingleStmt(sqlCommand);
    //drop old table
    sqlCommand.str(std::string());
    sqlCommand << "DROP TABLE " << tableName(oldTableId);
    execSingleStmt(sqlCommand);
    //rename new table
    sqlCommand.str(std::string());
    sqlCommand << " ALTER TABLE  " << tableName(tableId)
    << " RENAME TO " << tableName(oldTableId);
    result = execSingleStmt(sqlCommand);
    tableId=oldTableId;
    myMd->activeLabels=v1;
    return result;
}

size_t MDSql::size(void)
{
    std::stringstream ss;
    ss << "SELECT COUNT(*) FROM "<< tableName(tableId) << ";";
    return execSingleIntStmt(ss);
}

bool MDSql::setObjectValues( size_t id, const std::vector<MDObject*> columnValues, const std::vector<MDLabel> *desiredLabels)
{
    bool r = true;			// Return value.
    int i=0, j=0;			// Loop indexes.

    // Add values.
    if (desiredLabels==NULL)
    {
        bindValue( this->preparedStmt, 1, *(columnValues[0]));
        for (i=1; i<columnValues.size() ;i++)
        {
        	bindValue( this->preparedStmt, i+1, *(columnValues[i]));
        }
    }
    // Add only desired columns.
    else
    {
		for (i=0; i<desiredLabels->size() ;i++)
		{
			for (j=0; j<columnValues.size() ;j++)
			{
				if (columnValues[j]->label == (*desiredLabels)[i])
				{
					bindValue( this->preparedStmt, i+1, *(columnValues[j]));
					break;
				}
			}
		}
    }

    // id != -1 means there is a WHERE clause in the query sentence and id must be added.
    if (id != -1)
    {
    	sqlite3_bind_int( this->preparedStmt, i+1, id);
    }

    // Execute statement.
    rc = sqlite3_step( this->preparedStmt);
    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::setObjectValue(MDObject): " << std::endl
        << "   " << this->preparedStream.str() << std::endl
        <<"    code: " << rc << " error: " << sqlite3_errmsg(db) << std::endl;
        r = false;
    }

    // Reset statement and bindings.
    sqlite3_clear_bindings(this->preparedStmt);
    sqlite3_reset(this->preparedStmt);

    return r;
}

void MDSql::finalizePreparedStmt(void)
{
	if (this->preparedStmt != NULL)
	{
		sqlite3_finalize( this->preparedStmt);
		this->preparedStmt = NULL;
	}
}

//set column with a given value
bool MDSql::setObjectValue(const MDObject &value)
{
    bool r = true;
    MDLabel column = value.label;
    std::stringstream ss;
    sqlite3_stmt * stmt;
    ss << "UPDATE " << tableName(tableId)
    << " SET " << MDL::label2StrSql(column) << "=?;";
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    bindValue(stmt, 1, value);
    rc = sqlite3_step(stmt);
    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::setObjectValue(MDObject): " << std::endl
        << "   " << ss.str() << std::endl
        <<"    code: " << rc << " error: " << sqlite3_errmsg(db) << std::endl;
        r = false;
    }
    sqlite3_finalize(stmt);
    return r;
}

bool MDSql::setObjectValue(const int objId, const MDObject &value)
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
        << " SET " << MDL::label2StrSql(column) << "=? WHERE objID=?;";
//#define DEBUG
#ifdef DEBUG
        std::cerr << "DEBUG_JM: setObjectValue: " << ss.str() << std::endl;

#endif

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

bool MDSql::initializeSelect( bool addWhereObjId, std::vector<MDLabel> labels)
{
	int 	i=0;					// Loop counter.
	bool	createdOK=true;		// Return value.
	std::stringstream ss;		// Sentence string.

	// Initialize SELECT sentence.
	ss << "SELECT ";

	// Add columns names.
	if (labels.size() > 0)
	{
		// Add columns names.
		ss << MDL::label2StrSql(labels[0]);
		for (i=1; i<labels.size() ;i++)
		{
			if (labels[i] != MDL_STAR_COMMENT)
			{
				ss << "," << MDL::label2StrSql(labels[i]);
			}
		}
	}
	// Select all columns.
	else
	{
		ss << " * ";
	}

	// Add table.
	ss << " FROM " << tableName(tableId);

	// Check if add WHERE clause to select by id.
	if (addWhereObjId)
	{
		ss << " WHERE objID=?";
	}

	rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &this->preparedStmt, &zLeftover);
	if (rc != SQLITE_OK)
	{
		createdOK = false;
		printf( "could not prepare statement: %s\n", sqlite3_errmsg(db) );
		this->preparedStmt = NULL;
	}

	return(createdOK);
}

bool MDSql::initializeInsert(const std::vector<MDLabel> *labels, const std::vector<MDObject*> &values)
{
	int 	i=0;				// Loop counter.
	int		length=0;			// # labels.
	bool	createdOK=true;		// Return value.

	// Clear preparedStream.
	this->preparedStream.str(std::string());

	// Initialize SQL sentence.
	this->preparedStream << "INSERT INTO " << tableName(tableId);

	// Add columns.
	this->preparedStream << " (";

	// Execute branch if "labels" size is not zero.
	if (labels != NULL)
	{
		length = labels->size();
		this->preparedStream << MDL::label2StrSql((*labels)[0]);
		for (i=1; i<length ;i++)
		{
			this->preparedStream << "," << MDL::label2StrSql((*labels)[i]);
		}
	}
	// Execute branch if "values" size is not zero.
	else if ((length = values.size()) > 0)
	{
		this->preparedStream << MDL::label2StrSql((*values[0]).label);
		for (i=1; i<length ;i++)
		{
			this->preparedStream << "," << MDL::label2StrSql((*values[i]).label);
		}
	}

	// Add one '?' character for each label.
	this->preparedStream << ") VALUES (?";
	for (i=1; i<length ;i++)
	{
		this->preparedStream << ",?";
	}
	this->preparedStream << ");";

	// Prepare statement.
	//std::cout << this->preparedStream.str().c_str() << std::endl;
	rc = sqlite3_prepare_v2(db, this->preparedStream.str().c_str(), -1, &this->preparedStmt, &zLeftover);
	if (rc != SQLITE_OK)
	{
		printf( "initializeInsert: could not prepare statement: %s\n", sqlite3_errmsg(db) );
		this->preparedStmt = NULL;
		createdOK = false;
	}

	return(createdOK);
}

bool MDSql::initializeUpdate( std::vector<MDLabel> labels)
{
	int 	i=0;				// Loop counter.
	int		length=0;			// # labels.
	bool	createdOK=true;		// Return value.

	// Get # labels.
	length = labels.size();

	// Check there are labels.
	if (length > 0)
	{
		// Clear preparedStream.
		this->preparedStream.str(std::string());

		// Initialize SQL sentence.
		this->preparedStream << "UPDATE " << tableName(tableId) << " SET ";

		// Add labels.
		this->preparedStream << MDL::label2StrSql(labels[0]) << "=?";
		for (i=1; i<length ;i++)
		{
			this->preparedStream << ", " << MDL::label2StrSql(labels[i]) << "=?";
		}

		// Insert where clause.
		this->preparedStream << " WHERE objID=?;";

		// Prepare statement.
		rc = sqlite3_prepare_v2(db, this->preparedStream.str().c_str(), -1, &this->preparedStmt, &zLeftover);
		if (rc != SQLITE_OK)
		{
			printf( "initializeUpdate: could not prepare statement: %s\n", sqlite3_errmsg(db) );
			this->preparedStmt = NULL;
			createdOK = false;
		}
	}
	else
	{
		this->preparedStmt = NULL;
		createdOK = false;
	}

	return(createdOK);
}


bool MDSql::getObjectsValues( std::vector<MDLabel> labels, std::vector<MDObject> *values)
{
	bool ret=true;				// Return value.
	int i=0;					// Loop counter.

	// Execute statement.
	rc = sqlite3_step(this->preparedStmt);
	if (rc == SQLITE_ROW)
	{
		for (i=0; i<labels.size() ;i++)
		{
			if (labels[i] != MDL_STAR_COMMENT)
			{
				MDObject value(labels[i]);
				extractValue(this->preparedStmt, i, value);
				(*values).push_back(value);
			}
		}
	}
	// If no row retrieved then return false.
	else
	{
		ret = false;
	}

	return(ret);
}

bool MDSql::getObjectValue(const int objId, MDObject  &value)
{
    std::stringstream ss;
    MDLabel column = value.label;
    sqlite3_stmt * &stmt = myCache->getValueCache[column];

    if (stmt == NULL)//prepare stmt if not exists
    {
        //std::cerr << "Creating cache " << ++count <<std::endl;
        ss << "SELECT " << MDL::label2StrSql(column)
        << " FROM " << tableName(tableId)
        << " WHERE objID=?";// << objId << ";";
        rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    }

//#define DEBUG
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

void MDSql::selectObjects(std::vector<size_t> &objectsOut, const MDQuery *queryPtr)
{
    std::stringstream ss;
    sqlite3_stmt *stmt;
    objectsOut.clear();

    ss << "SELECT objID FROM " << tableName(tableId);
    if (queryPtr != NULL)
    {
        ss << queryPtr->whereString();
        ss << queryPtr->orderByString();
        ss << queryPtr->limitString();
    }
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
#ifdef DEBUG

    std::cerr << "selectObjects: " << ss.str() <<std::endl;
#endif

    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
    {
        objectsOut.push_back(sqlite3_column_int(stmt, 0));
    }
    rc = sqlite3_finalize(stmt);
}

size_t MDSql::deleteObjects(const MDQuery *queryPtr)
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

size_t MDSql::copyObjects(MetaData *mdPtrOut, const MDQuery *queryPtr) const
{
    return copyObjects(mdPtrOut->myMDSql, queryPtr);
}

size_t MDSql::copyObjects(MDSql * sqlOut, const MDQuery *queryPtr) const
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
        ss2 << sep << MDL::label2StrSql( myMd->activeLabels[i]);
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
                        const std::vector<MDLabel>            &operateLabel)
{
    std::stringstream ss;
    std::stringstream ss2;
    std::string aggregateStr = MDL::label2StrSql(mdPtrOut->activeLabels[0]);
    ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId)
    << "(" << aggregateStr;
    ss2 << aggregateStr;
    //Start iterating on second label, first is the
    //aggregating one
    for (size_t i = 0; i < operations.size(); i++)
    {
        ss << ", " << MDL::label2StrSql(mdPtrOut->activeLabels[i+1]);
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
            REPORT_ERROR(ERR_MD_SQL, "Invalid aggregate operation.");
        }
        ss2 << "(" << MDL::label2StrSql(operateLabel[i])
        << ") AS " << MDL::label2StrSql(mdPtrOut->activeLabels[i+1]);
    }
    ss << ") SELECT " << ss2.str();
    ss << " FROM " << tableName(tableId);
    ss << " GROUP BY " << aggregateStr;
    ss << " ORDER BY " << aggregateStr << ";";
    //std::cerr << "ss " << ss.str() <<std::endl;
    execSingleStmt(ss);
}


void MDSql::aggregateMdGroupBy(MetaData *mdPtrOut,
                               AggregateOperation operation,
                               const std::vector<MDLabel> &groupByLabels ,
                               MDLabel operateLabel,
                               MDLabel resultLabel)
{
    std::stringstream ss;
    std::stringstream ss2;
    std::stringstream groupByStr;

    groupByStr << MDL::label2StrSql(groupByLabels[0]);
    for (size_t i = 1; i < groupByLabels.size(); i++)
        groupByStr << ", " << MDL::label2StrSql(groupByLabels[i]);

    ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId) << "("
    << groupByStr.str() << ", " << MDL::label2StrSql(resultLabel) << ")";

    ss2 << groupByStr.str() << ", ";
    switch (operation)
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
        REPORT_ERROR(ERR_MD_SQL, "Invalid aggregate operation.");
    }
    ss2 << "(" << MDL::label2StrSql(operateLabel);
    ss2 << ") AS " << MDL::label2StrSql(resultLabel);

    ss << " SELECT " << ss2.str();
    ss << " FROM " << tableName(tableId);
    ss << " GROUP BY " << groupByStr.str();
    ss << " ORDER BY " << groupByStr.str() << ";";

    //std::cerr << "ss " << ss.str() <<std::endl;
    execSingleStmt(ss);
}


double MDSql::aggregateSingleDouble(const AggregateOperation operation,
                                    MDLabel operateLabel)
{
    std::stringstream ss;
    ss << "SELECT ";
    //Start iterating on second label, first is the
    //aggregating one
    switch (operation)
    {
    case AGGR_COUNT:
        ss << "COUNT";
        break;
    case AGGR_MAX:
        ss << "MAX";
        break;
    case AGGR_MIN:
        ss << "MIN";
        break;
    case AGGR_SUM:
        ss << "SUM";
        break;
    case AGGR_AVG:
        ss << "AVG";
        break;
    default:
        REPORT_ERROR(ERR_MD_SQL, "Invalid aggregate operation.");
    }
    ss << "(" << MDL::label2StrSql(operateLabel) << ")" ;
    ss << " FROM " << tableName(tableId);
    return (execSingleDoubleStmt(ss));
}

size_t MDSql::aggregateSingleSizeT(const AggregateOperation operation,
                                   MDLabel operateLabel)
{
    std::stringstream ss;
    ss << "SELECT ";
    //Start iterating on second label, first is the
    //aggregating one
    switch (operation)
    {
    case AGGR_COUNT:
        ss << "COUNT";
        break;
    case AGGR_MAX:
        ss << "MAX";
        break;
    case AGGR_MIN:
        ss << "MIN";
        break;
    case AGGR_SUM:
        ss << "SUM";
        break;
    case AGGR_AVG:
        ss << "AVG";
        break;
    default:
        REPORT_ERROR(ERR_MD_SQL, "Invalid aggregate operation.");
    }
    ss << "(" << MDL::label2StrSql(operateLabel) << ")" ;
    ss << " FROM " << tableName(tableId);
    return (execSingleIntStmt(ss));
}

void MDSql::indexModify(const std::vector<MDLabel> columns, bool create)
{
    std::stringstream ss,index_name,index_column;
    std::string sep1=" ";
    std::string sep2=" ";
    for (size_t i = 0; i < columns.size(); i++)
    {
        index_name << sep1 << tableName(tableId) << "_"
        << MDL::label2Str(columns.at(i));
        sep1 = "_";
        index_column << sep2 << MDL::label2StrSql(columns.at(i));
        sep2 = ", ";
    }

    if (create)
    {
        ss << "CREATE INDEX IF NOT EXISTS " << index_name.str() << "_INDEX "
        << " ON " << tableName(tableId) << " (" << index_column.str() << ")";
    }
    else
    {
        ss << "DROP INDEX IF EXISTS " << index_name.str() << "_INDEX ";
    }
    execSingleStmt(ss);
}

size_t MDSql::firstRow()
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MIN(objID), -1) AS MDSQL_FIRST_ID FROM "
    << tableName(tableId) << ";";
    return execSingleIntStmt(ss);
}

size_t MDSql::lastRow()
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MAX(objID), -1) AS MDSQL_LAST_ID FROM "
    << tableName(tableId) << ";";
    return execSingleIntStmt(ss);
}

size_t MDSql::nextRow(size_t currentRow)
{
    std::stringstream ss;
    ss << "SELECT COALESCE(MIN(objID), -1) AS MDSQL_NEXT_ID FROM "
    << tableName(tableId)
    << " WHERE objID>" << currentRow << ";";
    return execSingleIntStmt(ss);
}

size_t MDSql::previousRow(size_t currentRow)
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
    ss << "SELECT MAX(COALESCE(LENGTH("<< MDL::label2StrSql(column)
    <<"), -1)) AS MDSQL_STRING_LENGTH FROM "
    << tableName(tableId) << ";";
    return execSingleIntStmt(ss);
}

void MDSql::setOperate(MetaData *mdPtrOut, const std::vector<MDLabel> &columns, SetOperation operation)
{
    std::stringstream ss, ss2;
    bool execStmt = true;
    int size;
    std::string sep = " ";
    std::vector<MDLabel> * labelVector;

    switch (operation)
    {
    case UNION:
        copyObjects(mdPtrOut->myMDSql);
        execStmt = false;
        break;

    case UNION_DISTINCT: //unionDistinct
        //Create string with columns list
        size = mdPtrOut->activeLabels.size();
        //std::cerr << "LABEL" <<  MDL::label2StrSql(column) <<std::endl;
        for (int i = 0; i < size; i++)
        {
            ss2 << sep << MDL::label2StrSql( myMd->activeLabels[i]);
            sep = ", ";
        }
        ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId)
        << " (" << ss2.str() << ")"
        << " SELECT " << ss2.str()
        << " FROM " << tableName(tableId)
        << " WHERE ";
        for (size_t j=0; j<columns.size(); ++j)
        {
        	if (j>0)
        		ss << " AND ";
        	ss << MDL::label2StrSql(columns[j])
				<< " NOT IN (SELECT " << MDL::label2StrSql(columns[j])
				<< " FROM " << tableName(mdPtrOut->myMDSql->tableId) << ") ";
        }
        ss << ";";
        break;

    case DISTINCT:
    case REMOVE_DUPLICATE:
        //Create string with columns list
        size = mdPtrOut->activeLabels.size();
        sep = ' ';
        labelVector = &(mdPtrOut->activeLabels);
        for (int i = 0; i < size; i++)
        {
            ss2 << sep << MDL::label2StrSql( labelVector->at(i));
            sep = ", ";
        }
        if (operation == DISTINCT || columns[0] == MDL_UNDEFINED)
        {
            ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId)
            << " (" << ss2.str() << ")"
            << " SELECT DISTINCT " << ss2.str()
            << " FROM " << tableName(tableId) << ";";
        }
        else {
            // We need this special case for the REMOVE_DUPLICATE because when using a subset
            // of labels as key, the rest of values are taken randomly.
            // The following query ensures that the first ocurrence of the first group of rows
            // is used to take the remaining values
            ss << "INSERT INTO " << tableName(mdPtrOut->myMDSql->tableId)
            << " (ObjId," << ss2.str() << ")"
            << " SELECT M.* FROM (SELECT " << MDL::label2StrSql(columns[0]) << ", MIN(ObjId) AS first "
            << " FROM " << tableName(tableId) << " GROUP BY " << MDL::label2StrSql(columns[0])
            << " ) foo JOIN " << tableName(tableId) << " M ON foo.first = M.ObjId;";
        }
        break;

    case INTERSECTION:
    case SUBSTRACTION:
        ss << "DELETE FROM " << tableName(mdPtrOut->myMDSql->tableId)
        << " WHERE ";
        for (size_t j=0; j<columns.size(); ++j)
        {
        	if (j>0)
        		ss << " AND ";
        	ss << MDL::label2StrSql(columns[j]);
            if (operation == INTERSECTION)
                ss << " NOT";
			ss << " IN (SELECT " << MDL::label2StrSql(columns[j])
			   << " FROM " << tableName(tableId) << ") ";
        }
        ss << ";";
        break;
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Cannot use this operation for a set operation");
    }
    //std::cerr << "ss " << ss.str() <<std::endl;
    if (execStmt)
        execSingleStmt(ss);
}

bool MDSql::equals(const MDSql &op)
{
    std::vector<MDLabel> v1(myMd->activeLabels),v2(op.myMd->activeLabels);
    std::sort(v1.begin(),v1.end());
    std::sort(v2.begin(),v2.end());

    if(v1 != v2)
        return (false);
    int size  = myMd->activeLabels.size();
    std::stringstream sqlQuery,ss2,ss2Group;

    ss2 << MDL::label2StrSql(MDL_OBJID);
    ss2Group << MDL::label2StrSql(MDL_OBJID);
    int precision = myMd->precision;
    for (int i = 0; i < size; i++)
    {
        //when metadata  is double compare
        if(MDL::isDouble(myMd->activeLabels[i]))
        {
            ss2 << ", CAST (" << MDL::label2StrSql( myMd->activeLabels[i])
            << "*" << precision
            << " as INTEGER) as " << MDL::label2StrSql( myMd->activeLabels[i]);

        }
        else
        {
            ss2 << ", "       << MDL::label2StrSql( myMd->activeLabels[i]);
        }
        ss2Group<< ", "       << MDL::label2StrSql( myMd->activeLabels[i]);
    }
    sqlQuery
    << "SELECT count(*) FROM ("
    << "SELECT count(*) as result\
    FROM\
    (\
    SELECT " << ss2.str() << "\
    FROM " <<   tableName(tableId)
    <<      " UNION ALL \
    SELECT " << ss2.str() << "\
    FROM " << tableName(op.tableId)
    <<     ") tmp"
    << " GROUP BY " << ss2Group.str()
    << " HAVING COUNT(*) <> 2"
    << ") tmp1";
    return (execSingleIntStmt(sqlQuery)==0);
}

void MDSql::setOperate(const MetaData *mdInLeft,
                       const MetaData *mdInRight,
					   const std::vector<MDLabel> &columnsLeft,
					   const std::vector<MDLabel> &columnsRight,
                       SetOperation operation)
{
    std::stringstream ss, ss2, ss3;
    size_t size;
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
    case NATURAL_JOIN:
        /* We do not want natural join but natural join except for the obj-ID column */
        join_type = " INNER ";
        break;
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Cannot use this operation for a set operation");
    }
    if(operation==NATURAL_JOIN)
    {
        std::vector<MDLabel> intersectLabels;
        std::vector<MDLabel>::const_iterator left, right;
        for (right=(mdInRight->activeLabels).begin();
             right!=(mdInRight->activeLabels).end();
             ++right)
            for (left=(mdInLeft->activeLabels).begin();
                 left!=(mdInLeft->activeLabels).end();
                 ++left)
            {
                if (*left == *right)
                {
                    String labelStr = MDL::label2StrSql(*left);
                    intersectLabels.push_back(*left);
                }
            }
        mdInRight->addIndex(intersectLabels);
        mdInLeft->addIndex(intersectLabels);
    }
    else
    {
    	if (columnsRight.size()==1)
    	{
			mdInRight->addIndex(columnsRight[0]);
			mdInLeft->addIndex(columnsLeft[0]);
    	}
    }
    size = myMd->activeLabels.size();
    size_t sizeLeft = mdInLeft->activeLabels.size();

    for (size_t i = 0; i < size; i++)
    {
        ss2 << sep << MDL::label2StrSql( myMd->activeLabels[i]);
        ss3 << sep;
        if (i < sizeLeft && mdInLeft->activeLabels[i] == myMd->activeLabels[i])
            ss3 << tableName(mdInLeft->myMDSql->tableId) << ".";
        else
            ss3 << tableName(mdInRight->myMDSql->tableId) << ".";
        ss3 << MDL::label2StrSql( myMd->activeLabels[i]);
        sep = ", ";
    }
    ss << "INSERT INTO " << tableName(tableId)
    << " (" << ss2.str() << ")"
    << " SELECT " << ss3.str()
    << " FROM " << tableName(mdInLeft->myMDSql->tableId)
    << join_type << " JOIN " << tableName(mdInRight->myMDSql->tableId);

    if (operation != NATURAL_JOIN)
    {
        ss << " ON (";
        for (size_t j=0; j<columnsLeft.size(); ++j)
        {
        	if (j>0)
        		ss << " AND ";
        	ss << tableName(mdInLeft->myMDSql->tableId) << "." << MDL::label2StrSql(columnsLeft[j])
               << "=" << tableName(mdInRight->myMDSql->tableId) << "." << MDL::label2StrSql(columnsRight[j]);
        }
        ss << ") ";
    }
    else
    {
        sep = " ";
        ss << " WHERE ";
        for (size_t i = 0; i < mdInRight->activeLabels.size(); i++)
            for (size_t j = 0; j < sizeLeft; j++)
            {
                if(mdInRight->activeLabels[i] == mdInLeft->activeLabels[j])
                {
                    ss << sep
                    << tableName(mdInRight->myMDSql->tableId) << "."
                    << MDL::label2StrSql(mdInRight->activeLabels[i])
                    << " = "
                    << tableName(mdInLeft->myMDSql->tableId) << "."
                    << MDL::label2StrSql(mdInLeft->activeLabels[j]);
                    sep = " AND ";
                }
            }
    }
    ss << ";";
    //    std::cerr << "ss:" << ss.str() << std::endl;
    //    for (int j = 0; j < sizeLeft; j++)
    //     std::cerr << "mdInRight->activeLabels:" << mdInRight->activeLabels[0] << std::endl;
    //    for (int j = 0; j < sizeLeft; j++)
    //     std::cerr << "mdInLeft->activeLabels:"  << mdInLeft->activeLabels[1] << std::endl;
    execSingleStmt(ss);
    //std::cerr << "ss:" << ss.str() << std::endl;
    //dumpToFile("kk.sqlite");
    //exit(0);
}

bool MDSql::operate(const String &expression)
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
        REPORT_ERROR(ERR_MD_SQL, "dumpToFile: error opening db file");
    sqlite3_close(pTo);
    sqlBeginTrans();
}

void MDSql::copyTableFromFileDB(const FileName blockname,
                                const FileName filename,
                                const std::vector<MDLabel> *desiredLabels,
                                const size_t maxRows
                               )
{
    char **results;
    int rows;
    int columns;
    char *Labels;

    sqlite3 *db1;
    if (sqlite3_open(filename.c_str(), &db1))
        REPORT_ERROR(ERR_MD_SQL,formatString("Error opening database code: %d message: %s",rc,sqlite3_errmsg(db1)));
    String _blockname;
    String sql;
    if(blockname.empty())
    {
    	sql = (String)"SELECT name FROM sqlite_master\
    	                 WHERE type='table' LIMIT 1;";
        if ((rc=sqlite3_get_table (db1, sql.c_str(), &results, &rows, &columns, NULL)) != SQLITE_OK)
            REPORT_ERROR(ERR_MD_SQL,formatString("Error accessing table code: %d message: %s. SQL command %s",rc,sqlite3_errmsg(db1),sql.c_str()));
        else
        	_blockname=(String)results[1];
    }
    else
        _blockname=blockname;
    sql = (String)"PRAGMA table_info(" + _blockname +")";
    if (sqlite3_get_table (db1, sql.c_str(), &results, &rows, &columns, NULL) != SQLITE_OK)
        REPORT_ERROR(ERR_MD_SQL,formatString("Error accessing table code: %d message: %s. SQL command %s",rc,sqlite3_errmsg(db1),sql.c_str()));
    //This pragma returns one row for each column in the named table.
    //Columns in the result set include the column name,
    //data type, whether or not the column can be NULL, and the default value for the column.

    String activeLabel;
    MDLabel label;
    if (rows < 1)
        std::cerr << "Empty Metadata" <<std::endl;
    else if (desiredLabels != NULL)
    {

        myMd->activeLabels = *desiredLabels;
        for(std::vector<MDLabel>::const_iterator  it = desiredLabels->
                begin();
            it != desiredLabels->end();
            ++it)
        {
            activeLabel += *it + " ,";
        }
    }
    else
    {

        activeLabel="*";
        for (int i = 1; i <= rows; i++)
        {
            Labels = results[(i * columns) + 1];
            label = MDL::str2Label(Labels);

            if (label == MDL_UNDEFINED)
            {
                if(strcmp(Labels,"objID"))
                    std::cout << (String)"WARNING: Ignoring unknown column: " + Labels << std::endl;
            }
            else
            {
                myMd->activeLabels.push_back(label);
            }
        }
    }
    sqlite3_free_table (results);
    sqlite3_close(db1);

    //Copy table to memory
    //tableName(tableId);
    sqlCommitTrans();
    dropTable();
    createMd();

    String sqlCommand = formatString("ATTACH database '%s' AS load;", filename.c_str());

    if (sqlite3_exec(db, sqlCommand.c_str(), NULL, NULL, &errmsg) != SQLITE_OK)
    {
        std::cerr << "Couldn't attach or create table:  " << errmsg << std::endl;
        return;
    }
    String selectCmd = formatString("SELECT %s FROM load.%s", activeLabel.c_str(), _blockname.c_str());
    sqlCommand = formatString("INSERT INTO %s %s", tableName(tableId).c_str(), selectCmd.c_str());

    if (maxRows)
    {
        std::stringstream ss;
        ss << "SELECT COUNT(objId) FROM load." << _blockname;
        myMd->_parsedLines = execSingleIntStmt(ss);
        //std::cerr << ss.str() << " = " << myMd->_parsedLines << std::endl;

        sqlCommand += formatString(" LIMIT %lu", maxRows);
    }

    if (sqlite3_exec(db, sqlCommand.c_str(),NULL,NULL,&errmsg) != SQLITE_OK)
    {
        std::cerr << (String)"Couldn't write table: " << tableName(tableId)
        << " "                             << errmsg << std::endl
        << "sqlcommand " << sqlCommand << std::endl;
        return;
    }
    sqlite3_exec(db, "DETACH load",NULL,NULL,&errmsg);
    sqlBeginTrans();
}

void MDSql::copyTableToFileDB(const FileName blockname, const FileName &fileName)
{
    sqlCommitTrans();
    String _blockname;
    if(blockname.empty())
        _blockname=DEFAULT_BLOCK_NAME;
    else
        _blockname=blockname;
    String sqlCommand = (String)"ATTACH database '" + fileName+"' as save";
    sqlCommand += (String)";drop table  if exists save." + blockname;
    if (sqlite3_exec(db, sqlCommand.c_str(),NULL,NULL,&errmsg) != SQLITE_OK)
    {
        std::cerr << "Couldn't attach or create table:  " << errmsg << std::endl;
        return;
    }

    sqlCommand = (String)"create table save." +blockname
                 +" as select * from main."+tableName(tableId);
    //do not know how to create indexes in attached table
    //    sqlCommand += (String)";CREATE INDEX IF NOT EXISTS save.obj_Id_INDEX ON "
    //         + "save." +blockname +"(objID)";
    if (sqlite3_exec(db, sqlCommand.c_str(),NULL,NULL,&errmsg) != SQLITE_OK)
    {
        std::cerr << (String)"Couldn't write table: " << blockname
        << " "                             << errmsg << std::endl;
        return;
    }
    sqlite3_exec(db, "DETACH save",NULL,NULL,&errmsg);
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

void MDSql::sqlTimeOut(int miliseconds)
{
    if (sqlite3_busy_timeout(db, miliseconds) != SQLITE_OK)
    {
        std::cerr << "Couldn't not set timeOut:  " << std::endl;
        exit(0);

    }

}

void MDSql::sqlEnd()
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
        for (size_t i = 0; i < labelsVector->size(); i++)
        {
            ss << sep << MDL::label2SqlColumn(labelsVector->at(i));
            sep = ", ";
        }
    }
    ss << ");";
    return execSingleStmt(ss);
}

void MDSql::prepareStmt(const std::stringstream &ss, sqlite3_stmt *stmt)
{
    const char * zLeftover;
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
}

bool MDSql::execSingleStmt(const std::stringstream &ss)
{

    sqlite3_stmt * stmt;
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);

//#define DEBUG
#ifdef DEBUG

    std::cerr << "execSingleStmt, stmt: '" << ss.str() << "'" <<std::endl;
#endif
#undef DEBUG

    bool r = execSingleStmt(stmt, &ss);
    rc = sqlite3_finalize(stmt);
    return r;
}

bool MDSql::execSingleStmt(sqlite3_stmt * &stmt, const std::stringstream *ss)
{

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
        REPORT_ERROR(ERR_MD_SQL,formatString("Error code: %d message: %s\n  Sqlite query: %s",rc,sqlite3_errmsg(db), ss->str().c_str()));
    return true;
}

size_t MDSql::execSingleIntStmt(const std::stringstream &ss)
{
    sqlite3_stmt * stmt;
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    rc = sqlite3_step(stmt);
    size_t result = sqlite3_column_int(stmt, 0);

    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::execSingleIntStmt: error executing statement, code " << rc <<std::endl;
        result = -1;
    }
    rc = sqlite3_finalize(stmt);
    return result;
}

double MDSql::execSingleDoubleStmt(const std::stringstream &ss)
{
    sqlite3_stmt * stmt;
    rc = sqlite3_prepare_v2(db, ss.str().c_str(), -1, &stmt, &zLeftover);
    rc = sqlite3_step(stmt);
    double result = sqlite3_column_double(stmt, 0);

    if (rc != SQLITE_OK && rc != SQLITE_ROW && rc != SQLITE_DONE)
    {
        std::cerr << "MDSql::execSingleDoubleStmt: error executing statement, code " << rc <<std::endl;
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

bool MDSql::bindStatement( size_t id)
{
	bool success=true;		// Return value.

	// Clear current statement.
	sqlite3_clear_bindings(this->preparedStmt);
	sqlite3_reset(this->preparedStmt);

	// Bind object id.
	rc = sqlite3_bind_int(this->preparedStmt, 1, id);
	if (rc != SQLITE_OK)
	{
		success = false;
	}

	return(success);
}

int MDSql::bindValue(sqlite3_stmt *stmt, const int position, const MDObject &valueIn)
{
    //First reset the statement
    //rc  = sqlite3_reset(stmt);
    //std::cerr << "rc after reset: " << rc <<std::endl;
  if (valueIn.failed)
  {
    // If a value was wronly parsed, set NULL in their sqlite entry
    std::cerr << "WARNING!!! valueIn.failed = True, binding NULL" << std::endl;

    return sqlite3_bind_null(stmt, position);
  }
  else
  {
    switch (valueIn.type)
    {
    case LABEL_BOOL: //bools are int in sqlite3
        return sqlite3_bind_int(stmt, position, valueIn.data.boolValue ? 1 : 0);
    case LABEL_INT:
        return sqlite3_bind_int(stmt, position, valueIn.data.intValue);
    case LABEL_SIZET:
        return sqlite3_bind_int(stmt, position, valueIn.data.longintValue);
    case LABEL_DOUBLE:
        return sqlite3_bind_double(stmt, position, valueIn.data.doubleValue);
    case LABEL_STRING:
        return sqlite3_bind_text(stmt, position, valueIn.data.stringValue->c_str(), -1, SQLITE_TRANSIENT);
    case LABEL_VECTOR_DOUBLE:
    case LABEL_VECTOR_SIZET:
        return sqlite3_bind_text(stmt, position, valueIn.toString(false, true).c_str(), -1, SQLITE_TRANSIENT);
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Do not know how to handle this type");
    }
  }
}

void MDSql::extractValue(sqlite3_stmt *stmt, const int position, MDObject &valueOut)
{
    switch (valueOut.type)
    {
    case LABEL_BOOL: //bools are int in sqlite3
        valueOut.data.boolValue = sqlite3_column_int(stmt, position) == 1;
        break;
    case LABEL_INT:
        valueOut.data.intValue = sqlite3_column_int(stmt, position);
        break;
    case LABEL_SIZET:
        valueOut.data.longintValue = sqlite3_column_int(stmt, position);
        break;
    case LABEL_DOUBLE:
        valueOut.data.doubleValue = sqlite3_column_double(stmt, position);
        break;
    case LABEL_STRING:
    {
        std::stringstream ss;
        ss << sqlite3_column_text(stmt, position);
        valueOut.data.stringValue->assign(ss.str());

        break;
    }
    case LABEL_VECTOR_DOUBLE:
    case LABEL_VECTOR_SIZET:
    {
        std::stringstream ss;
        ss << sqlite3_column_text(stmt, position);
        valueOut.fromStream(ss);
        break;
    }
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Do not know how to extract a value from this type");
    }
}

MDCache::MDCache()
{
    this->addRowStmt = NULL;
    this->iterStmt = NULL;
}

MDCache::~MDCache()
{
    clear();
}

void MDCache::clear()
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
    {
        sqlite3_finalize(iterStmt);
        iterStmt = NULL;
    }

    if (addRowStmt != NULL)
    {
        sqlite3_finalize(addRowStmt);
        addRowStmt = NULL;
    }
}
