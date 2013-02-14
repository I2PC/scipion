/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

#ifndef PROGRAM_SQL_H_
#define PROGRAM_SQL_H_

#include <map>
#include "xmipp_program.h"
#include "external/sqlite-3.6.23/sqlite3.h"

typedef std::map<const char*, String> DictDB;

String& escapeSqliteStr(String & str);

/** Class that will encapsulate the Xmipp objects representation
 * on Sqlite db. Program are a kind of this objects.
 * Programas db storing will be useful for handle meta-info.
 */
class ProgramDb: public Printer
{
private:
    sqlite3 * db;
    int rc;
    char * errmsg, *zLeftover;


public:
    /** Some initialization */
    void init(const FileName &dbName);
    /** Empty constructor */
    ProgramDb();
    /** Constructor, it will create the Sqlite db. */
    ProgramDb(const FileName &dbName);
    /** Destructor */
    virtual ~ProgramDb() {}
    /** Begin and end transaction */
    bool execStmt(const String &stmt, const String &error="");
    bool beginTrans();
    bool commitTrans();
    /** Create tables related with programs */
    bool createProgramTables();
    /** Insert a program into db, the id field will be filled */
    bool insertProgram(DictDB &program);
    /** Get from the db the comment for a label */
    String getLabelComment(MDLabel label);

    //Methods inherits from Program Printer
    virtual void printProgram(const ProgramDef &program, int v = 0);
    virtual void printSection(const SectionDef &section, int v = 0);
    virtual void printParam(const ParamDef &param, int v = 0);
    virtual void printArgument(const ArgumentDef & argument, int v = 0);
    virtual void printCommentList(const CommentList &comments, int v = 0);

}
;//end of class ProgramDB
//
///** Print wiki text */
//class SqlitePrinter: public Printer
//{
//protected:
//    ProgramDb *db;
//
//public:
//    /**Constructor */
//    SqlitePrinter(ProgramDb *db);
//};
#endif /* PROGRAM_SQL_H_ */
