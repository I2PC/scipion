/***************************************************************************
 * Authors:     AUTHOR_NAME (josem@cnb.csic.es)
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

#include "xmipp_program.h"
#include <external/sqlite-3.6.23/sqlite3.h>

/** Class to represent data related to a program in a DB */
class DbProgram
{
public:
  int rank; //Used to compare and sorting
  int id;
  int cat_id;
  String name, description, keywords;

  bool operator < (const DbProgram & prog) { return rank < prog.rank; }

};//end of class DbProgram



/** Class to represent data of programs categories */
class DbCategory
{
public:
  int id;
  String name, description, prefixes;
};//end of class DbCategory

/** Class that will encapsulate the Xmipp objects representation
 * on Sqlite db. Program are a kind of this objects.
 * Programas db storing will be useful for handle meta-info.
 */
class XmippDB
{
private:
  sqlite3 * db;
  int rc;
  char * errmsg, *zLeftover;


public:
  /** Some initialization */
  void init(const FileName &dbName);
  /** Empty constructor */
  XmippDB();
  /** Constructor, it will create the Sqlite db. */
  XmippDB(const FileName &dbName);
  /** Begin and end transaction */
  bool execStmt(const String &stmt, const String &error="");
  bool beginTrans();
  bool commitTrans();
  /** Create program categories tables */
  bool createCategoryTable();
  /** Create tables related with programs */
  bool createProgramTable();
  /** Delete a program knowing the name */
  bool deleteProgramByName(const String &programName);
  /** Insert a program into db, the id field will be filled */
  bool insertProgram(DbProgram * program);
  /** Update program data, id must be valid */
  bool updateProgram(DbProgram * program);
  bool selectPrograms(std::vector<DbProgram*> &programs);

  /** Select all categories from DB */
  bool selectCategories(std::vector<DbCategory*> &categories);

};//end of class ProgramDB


#endif /* PROGRAM_SQL_H_ */
