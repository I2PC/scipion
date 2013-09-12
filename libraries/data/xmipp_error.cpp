/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <stdlib.h>
#include "xmipp_error.h"
#include "xmipp_color.h"

/* Exception handling ------------------------------------------------------ */
void _Xmipp_error(const ErrorType nerr, const String &what,
                  const String &file, const long line)
{
    XmippError xe(nerr, what, file, line);
    std::cerr << xe << std::endl;
    //std::cout << nerr << ": " << what << std::endl
    //<< "File: " << file << " line: " << line << std::endl;
    exit(nerr);
}

// Object Constructor
XmippError::XmippError(const ErrorType nerr, const String &what,
                       const String &fileArg, const long lineArg)
{
    __errno = nerr;
    msg = colorString(what.c_str(), RED);
    file = fileArg;
    line = lineArg;
//    std::cerr << "DEBUG_JM: " << getMessage() << std::endl;

    //Store information about the stack calls
//#ifdef LINUX
//    void *array[10];
//    size = backtrace(array, 10);
//    strings = backtrace_symbols(array, size);
//#endif
}

//Object Destructor
XmippError::~XmippError()
{
    //free(strings);
}

// Show message
std::ostream& operator << (std::ostream& o, XmippError& xe)
{
    String error = formatString("XMIPP_ERROR %d: %s", xe.__errno, xe.getDefaultMessage().c_str());
    o << colorString(error.c_str(), RED) << std::endl;
    o << colorString(xe.msg.c_str(), RED) << std::endl;
    error = formatString("File: %s line: %ld", xe.file.c_str(), xe.line);
    o << colorString(error.c_str(), RED) << std::endl;
    return o;
}

String XmippError::getMessage() const
{
      String error = formatString("XMIPP_ERROR %d: %s\n   ", __errno, getDefaultMessage().c_str());
      error += msg;
      error += formatString("\n   File: %s line: %ld\n", file.c_str(), line);
      return error;
}

String XmippError::getDefaultMessage() const
{
    return getDefaultMessage(__errno);
}

String XmippError::getDefaultMessage(ErrorType e)
{
    switch (e)
    {
    case ERR_ARG_BADCMDLINE:
        return "Errors on command line parameters";
    case ERR_ARG_INCORRECT:
        return " Incorrect argument received";
    case ERR_ARG_MISSING:
        return " Argument missing";
    case ERR_ARG_DEPENDENCE:
        return "Error with some arguments dependecies";

    case ERR_PROG_NOTDEF:
        return "Requiered function not implemented in derived class";
    case ERR_DEBUG_TEST:
        return " Just an error for debugging purpose";
    case ERR_DEBUG_IMPOSIBLE:
        return " Just for debugging: situation that can't happens";

    case ERR_DOCFILE:
        return " Error in docfile format";

    case ERR_GRID:
        return " Grid general error";
    case ERR_GRID_SIZE:
        return " Incorrect number of GRID volumes or shapes.";

    case ERR_IMG_NOREAD:
        return " Image cannot be read from file.";
    case ERR_IMG_NOWRITE:
        return " Image cannot be written to file.";
    case ERR_IMG_UNKNOWN:
        return " Unknown image type";

    case ERR_INDEX_OUTOFBOUNDS:
        return " Index out of bounds";

    case ERR_IO:
        return " Input/output general error.";
    case ERR_IO_NOCLOSED:
        return " File cannot be closed.";
    case ERR_IO_NOPATH:
        return " Environment PATH cannot be read.";
    case ERR_IO_NOPERM:
        return " Insufficient permissions to perform operation.";
    case ERR_IO_NOREAD:
        return " Couldn't read from file";
    case ERR_IO_NOTDIR:
        return " It is not a directory";
    case ERR_IO_NOTEXIST:
        return " File or directory does not exists";
    case ERR_IO_NOTFILE:
        return " It is not a file";
    case ERR_IO_NOTOPEN:
        return "File cannot be open.";
    case ERR_IO_NOWRITE:
        return " Couldn't write to file";
    case ERR_IO_SIZE:
        return " Incorrect file size.";

    case ERR_MATRIX:
        return " Matrix error";
    case ERR_MATRIX_DIM:
        return " Incorrect matrix dimensions.";
    case ERR_MATRIX_EMPTY:
        return " The matrix is empty";
    case ERR_MATRIX_SIZE:
        return " Incorrect matrix size";

    case ERR_MD:
        return " MetaData error.";
    case ERR_MD_BADLABEL:
        return " Unexpected label.";
    case ERR_MD_MISSINGLABEL:
        return " Missing expected label";
    case ERR_MD_BADTYPE:
        return " Bad label type.";
    case ERR_MD_NOACTIVE:
        return " No active object in MetaData.";
    case ERR_MD_NOOBJ:
        return " No exist requested object.";
    case ERR_MD_SQL:
        return " Error in SQL of MetaData operations";
    case ERR_MD_OBJECTNUMBER:
        return " Bad number of objects in MetaData";
    case ERR_MD_UNDEFINED:
        return " Undefined label.";
    case ERR_MD_BADBLOCK:
        return " Block not existing.";

    case ERR_MEM_BADREQUEST:
        return " Bad amount of memory requested.";
    case ERR_MEM_NOTENOUGH:
        return " There is not enough memory for allocation.";
    case ERR_MEM_NOTDEALLOC:
        return " Memory has not been deallocated.";
    case ERR_MEM_NULLPOINTER:
        return " Null pointer passed as parameter";

    case ERR_MMAP:
        return " Global mmap error.";
    case ERR_MMAP_NOTADDR:
        return " Map addressing of file has failed.";

    case ERR_MULTIDIM_DIM:
        return " Incorrect MultidimArray dimensions.";
    case ERR_MULTIDIM_EMPTY:
        return " MultidimArray is empty.";
    case ERR_MULTIDIM_SIZE:
        return " Incorrect MultidimArray size.";

    case ERR_NOT_IMPLEMENTED:
        return " Algorithm not implemented yet.";

    case ERR_NUMERICAL:
        return " Error related to numerical calculation.";

    case ERR_PARAM_INCORRECT:
        return " Parameter incorrect.";
    case ERR_PARAM_MISSING:
        return " Parameter missing.";

    case ERR_PLANS_NOCREATE:
        return " FFT plan cannot be created.";

    case ERR_SELFILE:
        return " Error in selfile format";

    case ERR_THREADS_NOTINIT:
        return " Threads cannot be initiated.";

    case ERR_TYPE_INCORRECT:
        return " Incorrect type received";

    case ERR_UNCLASSIFIED:
        return " Not classified error.";

    case ERR_VALUE_EMPTY:
        return " Empty value.";
    case ERR_VALUE_INCORRECT:
        return " Incorrect received value.";
    case ERR_VALUE_NOTSET:
        return " Value has not been set.";

    default:
        return "Unrecognized error code";
    }
}

void reportWarning(const String& what)
{
    String error = formatString("=== XMIPP_WARNING ===\n%s", what.c_str());
    std::cerr << colorString(error.c_str(), RED) << std::endl;
}

