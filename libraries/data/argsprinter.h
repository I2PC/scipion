/***************************************************************************
 * Authors:     J.M.de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your param) any later version.
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

#ifndef ARGSPRINTER_H_
#define ARGSPRINTER_H_

#include "argsparser.h"

#define XMIPP_MAJOR 3
#define XMIPP_MINOR 0

/**Define printers to show the arguments definitions.
 * This class is abstract and only define the basic
 * methods that a printer should have
 */
class Printer
{
public:
    virtual ~Printer() {}
    virtual void printProgram(const ProgramDef &program, int v = 0) = 0;
    virtual void printSection(const SectionDef &section, int v = 0) = 0;
    virtual void printParam(const ParamDef &param, int v = 0) = 0;
    virtual void printArgument(const ArgumentDef & argument, int v = 0) = 0;
    virtual void printCommentList(const CommentList &comments, int v = 0) = 0;
    virtual void printToken(ArgToken * token);
};


/** Just print to out stream */
class ConsolePrinter: public Printer
{
protected:
  std::ostream * pOut;
  void printRequiresList(StringVector requires);
public:
  bool color;
  /**Constructor */
  ConsolePrinter(std::ostream &out=std::cout, bool color = true);
    virtual void printProgram(const ProgramDef &program, int v = 0);
    virtual void printSection(const SectionDef &section, int v = 0);
    virtual void printParam(const ParamDef &param, int v = 0);
    virtual void printArgument(const ArgumentDef & argument, int v = 0);
    virtual void printCommentList(const CommentList &comments, int v = 0);
};

/** Print out to create Tk GUI */
class TkPrinter: public Printer
{
protected:
  FILE * output;
public:
  /** buffer to read the command line output */
  char readbuffer[1024];
  /** Constructor */
  TkPrinter();
  ~TkPrinter();
  virtual void printProgram(const ProgramDef &program, int v = 0);
  virtual void printSection(const SectionDef &section, int v = 0);
  virtual void printParam(const ParamDef &param, int v = 0);
  virtual void printArgument(const ArgumentDef & argument, int v = 0);
  virtual void printCommentList(const CommentList &comments, int v = 0){};
};

/** Print wiki text */
class WikiPrinter: public Printer
{
protected:
  std::ostream * pOut;
  void printRequiresList(StringVector requires);
public:
  /**Constructor */
  WikiPrinter(std::ostream &out=std::cout);
    virtual void printProgram(const ProgramDef &program, int v = 0);
    virtual void printSection(const SectionDef &section, int v = 0);
    virtual void printParam(const ParamDef &param, int v = 0);
    virtual void printArgument(const ArgumentDef & argument, int v = 0);
    virtual void printCommentList(const CommentList &comments, int v = 0);
};

/** Print out to create Protocol header script */
class ProtPrinter: public Printer
{
protected:
  FILE * output;
  String label, condition, parentName, exclusiveGroupName;
  StringVector stringBackup;
  size_t keyCounter;
  bool param_expert; //this will be used for arguments or expert parameters
  bool programGui; //this is for special case of use protocol gui for single programs outside project


public:
  /** buffer to read the command line output */
  char readbuffer[1024];
  /** Constructor */
  ProtPrinter(const char * scriptfile, bool programGui = false);
  virtual ~ProtPrinter();
  virtual void printProgram(const ProgramDef &program, int v = 0);
  virtual void printSection(const SectionDef &section, int v = 0);
  virtual void printParam(const ParamDef &param, int v = 0);
  virtual void printArgument(const ArgumentDef & argument, int v = 0);
  virtual void printCommentList(const CommentList &comments, int v = 0);
  void addCondition(const String &newcondition);

};

/** Print out to create Protocol header script */
class AutocompletePrinter: public Printer
{
protected:
  FILE * output;


public:
  /** Constructor */
  AutocompletePrinter(const char * scriptfile, bool programGui = false);
  virtual ~AutocompletePrinter();
  virtual void printProgram(const ProgramDef &program, int v = 0);
  virtual void printSection(const SectionDef &section, int v = 0);
  virtual void printParam(const ParamDef &param, int v = 0);
  virtual void printArgument(const ArgumentDef & argument, int v = 0);
  virtual void printCommentList(const CommentList &comments, int v = 0);
};

#endif /* ARGSPRINTER_H_ */
