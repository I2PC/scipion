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
#include "error.h"

/* Exception handling ------------------------------------------------------ */
void _Xmipp_error(const int nerr, const std::string &what,
                  const std::string &file, const long line)
{
    std::cout << nerr << ": " << what << std::endl
    << "File: " << file << " line: " << line << std::endl;
    exit(nerr);
}

// Object Constructor
XmippError::XmippError(const ErrorType nerr, const std::string &what,
                       const std::string &fileArg, const long lineArg)
{
    __errno = nerr;
    msg = what;
    file= fileArg;
    line=lineArg;

    //Store information about the stack calls
    void *array[10];
    size = backtrace(array, 10);
    strings = backtrace_symbols(array, size);
}

//Object Destructor
XmippError::~XmippError()
{
    free(strings);
}

void XmippError::printStackTrace(std::ostream& o)
{
    o << "Obtained " << size << " stack frames:" << std::endl;
    for (size_t i = 0; i < size; ++i)
        o << strings[i] << std::endl;
}

// Show message
std::ostream& operator << (std::ostream& o, XmippError& XE)
{
    o << XE.__errno << ":" << XE.getDefaultMessage() << std::endl
    << XE.msg << std::endl
    << "File: " << XE.file << " line: " << XE.line << std::endl;
    //XE.printStackTrace(o);

    return o;
}

char * XmippError::getDefaultMessage()
{
    return getDefaultMessage(__errno);
}
char * XmippError::getDefaultMessage(ErrorType e)
{
    switch (e)
    {
    case ERR_MEM_NOTENOUGH:
        return " There is not enough memory for allocation.";
    case ERR_IO_NOTEXIST:
        return " File or directory does not exists";
    case ERR_IO_NOPERM:
        return " Insufficient permissions to perform operation";
    case ERR_IO_NOREAD:
        return " Couldn't read from file";
    case ERR_IO_NOWRITE:
        return " Couldn't write to file";
    case ERR_IO_NOTFILE:
        return " It is not a file";
    case ERR_IO_NOTDIR:
        return " It is not a directory";
    case ERR_ARG_INCORRECT:
        return " Incorrect argument received";
    case ERR_ARG_MISSING:
        return " Argument missing";
    case ERR_TYPE_INCORRECT:
        return " Incorrect type received";
    case ERR_MATRIX:
        return " Matrix error";
    case ERR_MATRIX_EMPTY:
        return " The matrix is empty";
    case ERR_MATRIX_SIZE:
        return " Problem with matrix size";
    case ERR_MD:
        return " MetaData error";
    case ERR_MD_NOACTIVE:
        return " No active object in MetaData";
    case ERR_MD_NOOBJ:
        return " No exist requested object";
    case ERR_MD_BADLABEL:
        return " Unexpected label";
    case ERR_MD_SQL:
        return " Error in SQL of MetaData operations";
    case ERR_INDEX_OUTOFBOUNDS:
        return " Index out of bounds";
    case ERR_DEBUG_TEST:
        return " Just an error for debugging purpose";
    case ERR_DEBUG_IMPOSIBLE:
        return " Just for debugging: situation that can't happens";
    }
}



