/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosam@cnb.csic.es)
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

#include "program_extension.h"

//Needed includes for instanciate programs
#include "data/filters.h"
#include "angular_project_library.h"
#include "angular_discrete_assign.h"
#include "convert_pdb2vol.h"
#include "angular_continuous_assign.h"
#include "data/mask.h"


int runProgram(XmippProgram * program, const String &arguments, bool destroy)
{
  if (program == NULL)
          REPORT_ERROR(ERR_PARAM_INCORRECT, "Received a NULL as program pointer");
    program->read(arguments);
    int retCode = program->tryRun();
    if (destroy)
        delete program;
    return retCode;
}

int runProgram(const String &programName, const String &arguments)
{
    XmippProgram * program = getProgramByName(programName);
    return runProgram(program, arguments);
}

XmippProgram * getProgramByName(const String &programName)
{
//  if (programName == "xmipp_tranform_filter")
//    return new ProgFilter();

  if (programName == "xmipp_convert_pdb2vol")
    return new ProgPdbConverter();

  if (programName == "xmipp_angular_project_library")
    return new ProgAngularProjectLibrary();

  if (programName == "xmipp_mask")
    return new ProgMask();

  if (programName == "xmipp_angular_discrete_assign")
    return new ProgAngularDiscreteAssign();

  if (programName == "xmipp_angular_continuous_assign")
    return new ProgAngularContinuousAssign();



    return NULL;
}


