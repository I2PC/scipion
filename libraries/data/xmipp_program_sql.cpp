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

#include "xmipp_program_sql.h"

void ProgramDb::init(const FileName &dbName)
{
  rc = sqlite3_open(dbName.c_str(), &db);
  sqlite3_exec(db, "PRAGMA temp_store=MEMORY",NULL, NULL, &errmsg);
  sqlite3_exec(db, "PRAGMA synchronous=OFF",NULL, NULL, &errmsg);
  sqlite3_exec(db, "PRAGMA count_changes=OFF",NULL, NULL, &errmsg);
  sqlite3_exec(db, "PRAGMA page_size=4092",NULL, NULL, &errmsg);
}

ProgramDb::ProgramDb(const FileName &dbName)
{
    init(dbName);
}

ProgramDb::ProgramDb()
{
  init(formatString("%s/.xmipp_programs.sqlite", getXmippPath()));
}

bool ProgramDb::execStmt(const String &stmt, const String &error)
{
    if (sqlite3_exec(db, stmt.c_str(), NULL, NULL, &errmsg) != SQLITE_OK)
    {
        std::cerr << error << ":  " << errmsg << std::endl;
        return false;
    }
    return true;
}

bool ProgramDb::beginTrans()
{
    return execStmt("BEGIN TRANSACTION", "Couldn't begin transaction:  ");
}

bool ProgramDb::commitTrans()
{
    return execStmt("COMMIT TRANSACTION", "Couldn't commit transaction:  ");
}

/** Create tables related with programs */
bool ProgramDb::createProgramTables()
{
    const char * cmdStr =
        "DROP TABLE IF EXISTS Category;"
        "CREATE TABLE Category ("
        "   id INTEGER PRIMARY KEY ASC AUTOINCREMENT, "
        "   name TEXT UNIQUE, "
        "   desc TEXT, "
        "   prefixes TEXT);"
        "INSERT INTO Category VALUES(NULL, 'Classification', NULL, 'classify_ ml_ mlf_');"
        "INSERT INTO Category VALUES(NULL, 'CTF', NULL, 'ctf_');"
        "INSERT INTO Category VALUES(NULL, 'Images', NULL, 'image_');"
        "INSERT INTO Category VALUES(NULL, 'Metadatas', NULL, 'metadata_');"
        "INSERT INTO Category VALUES(NULL, 'Phantoms', NULL, 'phantom_ pdb_');"
        "INSERT INTO Category VALUES(NULL, 'Angular assignment', NULL, 'angular_');"
        "INSERT INTO Category VALUES(NULL, 'Tomography', NULL, 'tomo_ xray_');"
        "INSERT INTO Category VALUES(NULL, 'Transformations', NULL, 'transform_');"
        "INSERT INTO Category VALUES(NULL, 'Volumes', NULL, 'volume_ reconstruct_ resolution_');"
        "DROP TABLE IF EXISTS Program;"
        "CREATE TABLE Program ("
        "id INTEGER PRIMARY KEY ASC AUTOINCREMENT,"
        "category_id INTEGER, name TEXT UNIQUE, usage TEXT, examples TEXT,"
        "keywords TEXT);";

    return execStmt(cmdStr, "Couldn't create Program table:  ");
}

String& escapeSqliteStr(String & str)
{
    size_t pos = 0;
    while ((pos = str.find_first_of("'", pos)) != String::npos)
    {
        str.replace(pos, 1, "''");
        pos += 2;
    }
    str = "'" + str + "'";
    return str;
}

/** Insert program into db, the id field will be filled */
bool ProgramDb::insertProgram(DictDB &program)
{
    ///Delete first the program if exist
    //deleteProgramByName(program["name"]);
    std::stringstream ss;
    String &progName = escapeSqliteStr(program["name"]);
    ss
    << "DELETE FROM Program WHERE name=" << progName
    << ";INSERT INTO Program VALUES(NULL, NULL,"
    << progName << ","
    << escapeSqliteStr(program["usage"]) << ", "
    << escapeSqliteStr(program["examples"]) << ", "
    << escapeSqliteStr(program["keywords"]) << ");";
    bool result = execStmt(ss.str(), "Couldn't insert program");
    //program["id"] = sqlite3_last_insert_rowid(db);
    return result;
}

/** Select programs from db **/
//bool ProgramDb::selectPrograms()



/** Update program data, id must be valid */
//bool ProgramDb::updateProgram(DbProgram * program)
//{}
//
//bool ProgramDb::selectPrograms(std::vector<DbProgram*> &programs)
//{
//
//    sqlite3_stmt *stmt;
//    char * cmdStr = "SELECT * FROM Program ORDER BY name;";
//    DbProgram * progData;
//
//    rc = sqlite3_prepare_v2(db, cmdStr, -1, &stmt, NULL);
//    programs.clear();
//
//    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
//    {
//        progData = new DbProgram();
//        progData->id = sqlite3_column_int(stmt, 0);
//        progData->cat_id = sqlite3_column_int(stmt, 1);
//        progData->name.assign((char*)sqlite3_column_text(stmt, 2));
//        progData->description.assign((char*)sqlite3_column_text(stmt, 3));
//        progData->keywords.assign((char*)sqlite3_column_text(stmt, 4));
//        programs.push_back(progData);
//    }
//    rc = sqlite3_finalize(stmt);
//
//    return rc == SQLITE_OK;
//}

String ProgramDb::getLabelComment(MDLabel label)
{
      sqlite3_stmt *stmt;
      String aux = MDL::label2Str(label);
      String cmd = formatString("SELECT * FROM Label WHERE name='%s';", aux.c_str());

      rc = sqlite3_prepare_v2(db, cmd.c_str(), -1, &stmt, NULL);

      if ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
      {
          aux.assign((char*)sqlite3_column_text(stmt, 4));
      }
      rc = sqlite3_finalize(stmt);

      return aux;
}

//--------- SQLITE  PRINTER -----------------------
void ProgramDb::printProgram(const ProgramDef &program, int v)
{
    //print program name and usage
    String usage, examples;
    DictDB dict;

    dict["name"] = program.name;
    dict["usage"] = "";
    dict["examples"] = "";
    dict["keywords"] = program.keywords;
    //print usage
    if (program.usageComments.size() > 0)
    {
        for (size_t i = 0; i < program.usageComments.size(); ++i)
          dict["usage"] += program.usageComments.comments[i] + '\n';
    }
    //print examples
    if (program.examples.size() > 0)
    {
        for (size_t i = 0; i < program.examples.size(); ++i)
          dict["examples"] += program.examples.comments[i] + '\n';
    }
    insertProgram(dict);
    //std::cerr << "DEBUG_JM: program.name: " << program.name << std::endl;

    //print sections and params
    if (program.sections.size() > 0)
    {
        for (size_t i = 0; i < program.sections.size(); ++i)
            printSection(*program.sections[i], v);
    }
}

void ProgramDb::printSection(const SectionDef &section, int v)
{
}

void ProgramDb::printParam(const ParamDef &param, int v)
{
}

void ProgramDb::printArgument(const ArgumentDef & argument, int v)
{
}

void ProgramDb::printCommentList(const CommentList &comments, int v)
{
}
