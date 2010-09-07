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

#include "program.h"

void XmippProgram::init()
{
    progLexer = new ArgLexer();
    progDef = new ProgramDef(progLexer);
    this->defineParams();
    progLexer->nextToken();
    progDef->parse();
}

XmippProgram::XmippProgram()
{
    progLexer = NULL;
    progDef = NULL;
}

XmippProgram::XmippProgram(int argc, char ** argv)
{
    init();
    read(argc, argv);
}

XmippProgram::~XmippProgram()
{
    delete progLexer;
    delete progDef;
}

void XmippProgram::defineParams()
{
    REPORT_ERROR(ERR_PROG_NOTDEF, "function 'defineParams'");
}

void XmippProgram::readParams()
{
    REPORT_ERROR(ERR_PROG_NOTDEF, "function 'readParams'");
}

void XmippProgram::read(int argc, char ** argv)
{
    if (progLexer == NULL || progDef == NULL)
        init();

    if (argc == 1)
        usage();

    try
    {
        //TODO: Check if the command line is correct
        progDef->read(argc, argv);
        //TODO: Check if only requested help message
        this->readParams();

    }
    catch (XmippError xe)
    {
        std::cerr << xe;
        usage();
    }

}

void XmippProgram::addParamsLine(const char * line)
{
    progLexer->addLine((std::string)line);
}

const char * XmippProgram::getParam(const char * param, int arg)
{
    return progDef->getParam(param, arg);
}

int XmippProgram::getIntParam(const char * param, int arg)
{
    return textToInteger(progDef->getParam(param, arg));
}

double XmippProgram::getDoubleParam(const char * param, int arg)
{
    return textToFloat(progDef->getParam(param, arg));
}

const char * XmippProgram::name()
{
    return progDef->name.c_str();
}

void XmippProgram::usage()
{
    ConsolePrinter cp;
    cp.printProgram(*progDef);
    exit(1);
}

void XmippProgram::extendedUsage()
{
    ConsolePrinter cp;
    cp.printProgram(*progDef, 1);
    exit(1);
}
