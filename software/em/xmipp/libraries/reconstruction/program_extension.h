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

#ifndef PROGRAM_EXTENSION_H_
#define PROGRAM_EXTENSION_H_

#include <data/xmipp_program.h>

/** @defgroup Programs3 Some extensions related with Xmipp programs
 *  @ingroup DataLibrary
 *
 * Program extension routines.
 *
 * Here you will find some routines which might help you to use
 * programs from other programs. Also providing an accesible
 * interface for binding with Python and/or Java.
 *
 * @{
 */

/** Run a program either as a system call or as an internal call */
void runSystem(const String &program, const String &arguments, bool useSystem=true);

/** Run a program passing some arguments.
 * The program will be supplied as a pointer of XmippProgram
 * if destroy is true the program pointer will be freed */
int runProgram(XmippProgram *program, const String &arguments, bool destroy = true);
/** Run a program providing the program name.
 * Also the arguments should be supplied in a String.
 */
int runProgram(const String &programName, const String &arguments);
/** This function will create a pointer to XmippProgram taking the program name.
 * If the name is invalid, it will return NULL.
 * This function will be hardcoded since C++ has not easy way to do reflection,
 * so it need to be updated for new programs.
 */
XmippProgram * getProgramByName(const String &programName);
/** @} */

#endif /* PROGRAM_EXTENSION_H_ */
