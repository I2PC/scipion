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

#ifndef PROGRAM_H_
#define PROGRAM_H_

#include "argsparser.h"
#include "error.h"
#include "strings.h"

/** This class represent an Xmipp Program.
 * It have some of the basic functionalities of
 * the programs like argument parsing, cheking, usage printing.
 */
class XmippProgram
{
private:
  ArgLexer * progLexer;
  ProgramDef * progDef;

  /** Initialization function */
  void init();

protected:
  /** @name Functions to be implemented by subclasses.
   * @{
   */
  /** Function in which the param of each Program are defined. */
  virtual void defineParams();
  /** Function in which each program will read parameters that it need.
   * If some error occurs the usage will be printed out.
   * */
  virtual void readParams();
  /** This operator will be used in the function defineParams()
   *  to add lines to parser */
  void addParamsLine(const char * line);

  /** Get the argument of this param, first start at 0 position */
  const char * getParam(const char * param, int arg = 0);
  int getIntParam(const char * param, int arg = 0);
  double getDoubleParam(const char * param, int arg = 0);
  bool checkParam(const char * param);
public:
  /** @name Public common functions
   * The functions in this section are available for all programs
   * and should not be reimplemented in derived class, since they
   * are thought to have the same behaivor and it depends on the
   * definition of each program.
   * @{
   */
  /** Returns the name of the program
   * the name of the program is defined
   * by each subclass of this base class.
   */
  const char * name();
  /** Print the usage of the program, reduced version */
  void usage();
  /** Print detailed usage of the program */
  void extendedUsage();

  /** Return the version of the program */
  int version();
  /** Read the command line arguments
   * If an error occurs while reading arguments,
   * the error message will be printed and the
   * usage of the program showed. So the is no need
   * to do that in readParams();
   * */
  void read(int argc, char ** argv);

  /** @} */
  /** This function will be start running the program.
   * it also should be implemented by derived classes.
   */
  virtual void run() = 0;

  /** @name Constructors
   * @{
   */
  /** Constructor */
  XmippProgram();
  /** Constructor for read params */
  XmippProgram(int argc, char ** argv);
  /** Destructor */
  virtual ~XmippProgram();

};

#endif /* PROGRAM_H_ */
