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

#include "program_sql.h"

XmippDB::XmippDB(const FileName &dbName)
{
    rc = sqlite3_open(dbName.c_str(), &db);

    sqlite3_exec(db, "PRAGMA temp_store=MEMORY",NULL, NULL, &errmsg);
    sqlite3_exec(db, "PRAGMA synchronous=OFF",NULL, NULL, &errmsg);
    sqlite3_exec(db, "PRAGMA count_changes=OFF",NULL, NULL, &errmsg);
    sqlite3_exec(db, "PRAGMA page_size=4092",NULL, NULL, &errmsg);

}

bool XmippDB::execStmt(const std::string &stmt, const std::string &error)
{
    if (sqlite3_exec(db, stmt.c_str(), NULL, NULL, &errmsg) != SQLITE_OK)
    {
        std::cerr << error << ":  " << errmsg << std::endl;
        return false;
    }
    return true;
}

bool XmippDB::beginTrans()
{
    return execStmt("BEGIN TRANSACTION", "Couldn't begin transaction:  ");
}

bool XmippDB::commitTrans()
{
    return execStmt("COMMIT TRANSACTION", "Couldn't commit transaction:  ");
}

bool XmippDB::createCategoryTable()
{
    char * cmdStr =
        "DROP TABLE IF EXISTS Category;"
        "CREATE TABLE Category ("
        "id INTEGER PRIMARY KEY ASC AUTOINCREMENT,"
        "name TEXT UNIQUE, desc TEXT, prefixes TEXT);"
        "INSERT INTO Category VALUES(NULL, 'Micrograph', 'Programs to work with micrographs', 'micrograph_');"
        "INSERT INTO Category VALUES(NULL, 'Metadata', 'Selfiles, docfiles and metadatas', 'selfile_ docfile_ metadata_');"
        "INSERT INTO Category VALUES(NULL, 'Header', 'Header manipulation', 'header_');"
        "INSERT INTO Category VALUES(NULL, 'Classification', 'classification programs', 'classify_');"
        ;
    return execStmt(cmdStr, "Couldn't create Category table:  ");
}

/** Create tables related with programs */
bool XmippDB::createProgramTable()
{
    char * cmdStr =
        "DROP TABLE IF EXISTS Program;"
        "CREATE TABLE Program ("
        "id INTEGER PRIMARY KEY ASC AUTOINCREMENT,"
        "cat_id INTEGER, name TEXT UNIQUE, desc TEXT,"
        "keywords TEXT);";

    return execStmt(cmdStr, "Couldn't create Program table:  ");
}

/** Insert a program into db, the id field will be filled */
bool XmippDB::insertProgram(DbProgram * program)
{
    ///FIXME: remove single quote or other special characters to sqlite
    std::stringstream ss;
    ss << "INSERT INTO Program VALUES(NULL, NULL,'"
    << program->name << "','" << program->description
    << "', '" << program->keywords << "');";
    bool result = execStmt(ss.str(), "Couldn't insert program");
    program->id = sqlite3_last_insert_rowid(db);
    return result;
}

/** Update program data, id must be valid */
bool XmippDB::updateProgram(DbProgram * program)
{}

bool XmippDB::selectPrograms(std::vector<DbProgram*> &programs)
{

    sqlite3_stmt *stmt;
    char * cmdStr = "SELECT * FROM Program;";
    DbProgram * progData;

    rc = sqlite3_prepare_v2(db, cmdStr, -1, &stmt, NULL);
    programs.clear();

    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
    {
        progData = new DbProgram();
        progData->id = sqlite3_column_int(stmt, 0);
        progData->cat_id = sqlite3_column_int(stmt, 1);
        progData->name.assign((char*)sqlite3_column_text(stmt, 2));
        progData->description.assign((char*)sqlite3_column_text(stmt, 3));
        progData->keywords.assign((char*)sqlite3_column_text(stmt, 4));
        programs.push_back(progData);
    }
    rc = sqlite3_finalize(stmt);

    return rc == SQLITE_OK;
}
