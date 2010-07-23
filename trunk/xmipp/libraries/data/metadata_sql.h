/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#ifndef METADATASQL_H
#define METADATASQL_H

#include <iostream>
#include <map>
#include "strings.h"
#include <external/sqlite-3.6.23/sqlite3.h>
#include "metadata_label.h"



class MDSqlStaticInit;
class MDQuery;
class MetaData;
class MDCache;

enum AggregateOperation
{
    AGGR_COUNT, AGGR_MAX, AGGR_MIN, AGGR_SUM, AGGR_AVG
};

enum SetOperation
{
    UNION, UNION_DISTINCT, INTERSECTION, SUBSTRACTION, INNER_JOIN, LEFT_JOIN, OUTER_JOIN
};

enum JoinType
{
    INNER=INNER_JOIN, LEFT=LEFT_JOIN, OUTER=OUTER_JOIN
};

#include "metadata.h"

/** MDSql this class will manage SQL database interactions
 * will be used to implement different MetaData functions
 * The class will be static(only static methods)
 * and only available for MetaData
 */
class MDSql
{
public:

    /**This will create the table to store
     * the metada objects, will return false
     * if the mdId table is already present
     */
    bool createMd();

    /**This function will drop the entire table
     * for use the metada again, a call
     * to 'createMd' should be done
     */
    bool clearMd();

    /**Add a new row and return the objId(rowId)
     *
     */
    long int addRow();

    /** Add a new column to a metadata
     */
    bool addColumn(MDLabel column);
    /**Set the value of an object in an specified column
     *
     */
    bool setObjectValue(const int objId, const MDValue &value);

    /**Get the value of an object
     *
     */
    bool getObjectValue(const int objId, MDValue &value);

    /** This function will do a select
     * the 'limit' is the maximum number of object
     * returned, if is -1, all will be returned
     * Also a query could be specified for selecting objects
     * if no query is provided by default all are returned
     */
    void selectObjects(std::vector<long int> &objectsOut, const MDQuery *queryPtr = NULL);

    /** This function will delete elements
     * that match the query
     * if not query is provided, all rows are deleted
     */
    long int deleteObjects(const MDQuery *queryPtr = NULL);

    /** Coppy the objects from a metada to other
     * return the number of objects copied
     * */
    long int copyObjects(MDSql * sqlOut,
                         const MDQuery *queryPtr = NULL) const;
    long int copyObjects(MetaData * mdPtrOut,
                         const MDQuery *queryPtr = NULL) const;

    /** This function is to perform aggregation operations
     *
     */
    void aggregateMd(MetaData *mdPtrOut,
                     const std::vector<AggregateOperation> &operations,
                     MDLabel operateLabel);

    /**
     *This function will be used to create o delete an index over a column
     *to improve searchs, but inserts become expensives
     */
    void indexModify(const MDLabel label, bool create=true);

    /** Some iteration methods
     *
     */
    long int firstRow();
    long int lastRow();
    long int nextRow(long int currentRow);
    long int previousRow(long int currentRow);

    int columnMaxLength(MDLabel column);

    /**Functions to implement set operations */
    void setOperate(MetaData *mdPtrOut, MDLabel column, SetOperation operation);
    void setOperate(const MetaData *mdInLeft, const MetaData *mdInRight, MDLabel column, SetOperation operation);
    /** Function to dump DB to file */
    bool operate(const std::string &expression);

    static void dumpToFile(const FileName &fileName);


    /** Constructor of MDSql
     * Now each MD should have an instance
     * of this class to interact with the DB
     */
    MDSql(MetaData *md);
    ~MDSql();


private:
    static int table_counter;
    static sqlite3 *db;

    static MDSqlStaticInit initialization; //Just for initialization

    ///Just call this function once, at static initialization
    static bool sqlBegin();
    static bool sqlEnd();
    static bool sqlBeginTrans();
    static bool sqlCommitTrans();
    /**Return an unique id for each metadata
        * this function should be called once for each
        * metada and the id will be used for operations
        */
    int getUniqueId();

    bool dropTable();
    bool createTable(const std::vector<MDLabel> * labelsVector = NULL, bool withObjID=true);
    bool insertValues(double a, double b);
    void prepareStmt(const std::stringstream &ss, sqlite3_stmt *stmt);
    bool execSingleStmt(const std::stringstream &ss);
    bool execSingleStmt(sqlite3_stmt *&stmt, const std::stringstream *ss = NULL);
    long int execSingleIntStmt(const std::stringstream &ss);
    std::string tableName(const int tableId) const;

    int bindValue(sqlite3_stmt *stmt, const int position, const MDValue &valueIn);
    int extractValue(sqlite3_stmt *stmt, const int position, MDValue &valueOut);

    static char *errmsg;
    static const char *zLeftover;
    static int rc;
    static sqlite3_stmt *stmt;

    ///Non-static attributes
    int tableId;
    MetaData *myMd;
    MDCache *myCache;

    friend class MDSqlStaticInit;
}
;//close class MDSql


class MDSqlStaticInit
{
private:
    MDSqlStaticInit()
    {
        MDSql::sqlBegin();
    }//close constructor

    ~MDSqlStaticInit()
    {
        MDSql::sqlEnd();
    }//close destructor

    friend class MDSql;
}
;//close class MDSqlStaticInit


///@defgroup MetaDataQuery represent queries to a metadata
//@ingroup DataLibrary

/** MDQuery this is the base class for queries, its abstract
 *@ingroup MetaDataQuery
 */
class MDQuery
{
public:
    int limit, offset;
    MDLabel orderLabel;

    MDQuery(int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID)
    {
        this->limit = limit;
        this->offset = offset;
        this->orderLabel = orderLabel;
    }

    std::string orderByString() const
    {
        return (std::string)" ORDER BY " + MDL::label2Str(orderLabel);
    }

    std::string limitString() const
    {
        std::stringstream ss;
        if (limit != -1)
            ss << " LIMIT " << limit << " ";
        if (offset > 0)
            ss << "OFFSET " << offset << " ";
        return ss.str();
    }

    virtual std::string queryStringFunc() const
    {
        return " ";
    }

    std::string whereString() const
    {
        std::string queryString = this->queryStringFunc();
        return (queryString == " ") ? " " : " WHERE " + queryString + " ";
    }
    /**This now is specific to the SQL implementation
     * and its requiered to all MDQuery subclasses
     * equal to 0 means that is ABSTRACT in this class
     * and only accesible for MetaData
     */

};

enum RelationalOp { EQ, NE, GT, LT, GE, LE };
/**MDValueRelational will compare a column with some value
 *@ingroup MetaDataQuery
 */
class MDValueRelational: public MDQuery
{
    MDValue *value;
    RelationalOp op;
public:


    MDValueRelational()
    {
        value = NULL;
    }

    template <class T>
    MDValueRelational(MDLabel label, const T &value, RelationalOp op, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        this->op = op;
        this->value = new MDValue(label, value);
    }

    MDValueRelational(const MDValue &value, RelationalOp op, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        this->op = op;
        this->value = new MDValue(value);
    }

    ~MDValueRelational()
    {
        if (value != NULL)
            delete value;
    }

    std::string opString() const
    {
        switch (op)
        {
        case EQ:
            return "=";
        case NE:
            return "!=";
        case GT:
            return ">";
        case LT:
            return "<";
        case GE:
            return ">=";
        case LE:
            return "<=";
        }
    }

    virtual std::string queryStringFunc() const
    {
        return (value == NULL) ? " " : MDL::label2Str(value->label) + opString() + value->toString(false, true);
    }

    template <class T>
    void setValue(T &value)
    {
        if (this->value != NULL)
            delete this->value;
        this->value = new MDValue(this->value->label, value);
    }
}
;//class MDValueEqual

/**MDValueEqual this will test if a label have a value
 *@ingroup MetaDataQuery
 */
class MDValueEQ: public MDValueRelational
{
public:
    MDValueEQ()
    {}
    template <class T>
    MDValueEQ(MDLabel label, const T &value, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, value, EQ, limit, offset, orderLabel)
    {}
}
;//class MDValueEQ

/**MDValueAbove this will test if a label have a value larger than a minimum
 *@ingroup MetaDataQuery
 */
class MDValueGE: public MDValueRelational
{
public:
    MDValueGE()
    {}
    template <class T>
    MDValueGE(MDLabel label, const T &valueMin, int limit = -1,int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, valueMin, GE, limit, offset, orderLabel)
    {}
}
;//class MDValueGE

/**MDValueAbove this will test if a label have a value smaller than a maximum
 *@ingroup MetaDataQuery
 */
class MDValueLE: public MDValueRelational
{
public:
    MDValueLE()
    {}
    template <class T>
    MDValueLE(MDLabel label, const T &valueMax, int limit = -1,int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, valueMax, LE, limit, offset, orderLabel)
    {}
}
;//class MDValueLE

/**MDMultiQuery this will combine many queries with AND and OR operations
 *@ingroup MetaDataQuery
 */
class MDMultiQuery: public MDQuery
{
private:
    std::vector<const MDQuery*> queries;
    std::vector<std::string> operations;

public:

    MDMultiQuery(int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        clear();
    }
    void addAndQuery(MDQuery &query)
    {
        queries.push_back(&query);
        operations.push_back("AND");
    }
    void addOrQuery(MDQuery &query)
    {
        queries.push_back(&query);
        operations.push_back("OR");
    }

    void clear()
    {
        queries.clear();
        operations.clear();
    }

    virtual std::string queryStringFunc() const
    {
        if (queries.size() > 0)
        {
            std::stringstream ss;
            ss << "(" << queries[0]->queryStringFunc() << ") ";
            for (int i = 1; i < queries.size(); i++)
                ss << operations[i] << " (" << queries[i]->queryStringFunc() << ") ";

            return ss.str();
        }
        return " ";
    }

}
;//class MDMultiQuery


/**MDValueRange this will test if a label have a value within a minimum and maximum
 *@ingroup MetaDataQuery
 */
class MDValueRange: public MDQuery
{
    MDValueRelational *query1, *query2;
public:
    MDValueRange()
    {}
    template <class T>
    MDValueRange(MDLabel label, const T &valueMin, const T &valueMax,
                 int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        query1 = new MDValueRelational(label, valueMin, GE);
        query2 = new MDValueRelational(label, valueMax, LE);
    }

    virtual std::string queryStringFunc() const
    {
        std::stringstream ss;
        ss << "(" << query1->queryStringFunc() << " AND " << query2->queryStringFunc() << ")";
        return ss.str();
    }

    ~MDValueRange()
    {
        delete query1;
        delete query2;
    }
}
;//class MDValueRange


/** Just a class to store some cached sql statements
 *
 */
class MDCache
{
public:
    sqlite3_stmt *iterStmt;
    std::map<MDLabel, sqlite3_stmt*> getValueCache;
    std::map<MDLabel, sqlite3_stmt*> setValueCache;
    sqlite3_stmt *addRowStmt;

    MDCache();
    ~MDCache();
};

template<class T>
MDValueEQ MDValueEqualSwig(MDLabel label, const T &value)
{
    return MDValueEQ(label, value);
}

template<class T>
MDValueRange MDValueRangeSwig(MDLabel label, const T &valueMin, const T &valueMax)
{
    return MDValueRange(label, valueMin, valueMax);
}

template<class T>
MDValueGE MDValueAboveSwig(MDLabel label, const T &valueMin)
{
    return MDValueGE(label, valueMin);
}

template<class T>
MDValueLE MDValueBelowSwig(MDLabel label, const T &valueMax)
{
    return MDValueLE(label, valueMax);
}
#endif
