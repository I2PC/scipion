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

#ifndef XMIPP_ERROR_H
#define XMIPP_ERROR_H

#include <iostream>
#include "xmipp_strings.h"

#ifdef LINUX
#include <execinfo.h>
#endif

/** @defgroup ErrorHandling Error handling
 * @ingroup DataLibrary
 *
 * The error handling is performed in two different ways depending on the
 * configuration selected for Xmipp in the file xmippConfiguration: a simple
 * error management and a second method based on C++ exceptions.
 *
 * The first method aborts the program with an error code (different for each
 * error) while the second throws an exception which might be caught by an
 * external routine or program.
 *
 * The prototype definitions in both cases are the same as they are based on
 * some macros which change with the configuration. Here goes a programming
 * example considering both implementations.
 *
 * @code
 * // Class definition
 * class ReconstructingVolume :
 * {
 * ...
 * void write(const FileName& fn) const;
 * ...
 * };
 *
 * // Class implementation
 * void ReconstructingVolume::write(const FileName& fn) const
 * {
 * ...
 * if (...)
 *     REPORT_ERROR(ERR_MULTIDIM_SIZE, "Volume too small to be stored");
 * ...
 * }
 *
 * // Use of this class in an external program
 * using std::cout;
 * using std::endl;
 * ...
 * #ifndef _NO_EXCEPTION
 * try
 * {
 *     vol_blobs.write(fn_blobs);
 * }
 * catch (XmippError XE)
 * {
 *     std::cout << XE;
 *     std::cout << "The reconstructed volume is too small to be saved in blobs";
 *     std::cout << "So, there is no blob version of it at this iteration";
 *     std::cout << "I go on processing" << std::endl;
 * }
 * #else
 * vol_blobs.write(fn_blobs);
 * #endif
 *  ...
 * @endcode
 *
 * You see that the routine implementation is the same in both cases but the
 * external program varies from one to the other as in the exception case we can
 * catch the exception and go on processing, while in the exit mode, the program
 * always is aborted. If you don't put the routine in a try-catch structure and
 * an exception is thrown then a core is generated and the program is
 * automatically aborted.
 *
 */
/** @{ */
/* Enum with errors types.
 * This enum represent the code of all possible
 * Xmipp erros that will be used to reporting errors
 * with REPORT_ERROR and EXIT_ERROR.
 * The convention for the codes is the following:
 * - All starts with prefix ERR_
 * - Follows some kind of section. Like IO for input/output errors.
 *      ERR_IO_
 *      ERR_MEM_
 *      ERR_IMG_
 * - Finally an abreviation for the error msg.
 * All error codes have a default string message
 * that can be obtained with XmippError::getDefaultMessage(ErrorType)
 */
enum ErrorType
{
    ERR_FIRST_LABEL,
    ERR_ARG_BADCMDLINE,     ///< Errors on command line parameters.
    ERR_ARG_INCORRECT,      ///< Incorrect argument received.
    ERR_ARG_MISSING,        ///< Argument missing.
    ERR_ARG_DEPENDENCE,     ///< Error with some arguments dependecies

    ERR_PROG_NOTDEF,        ///< Requiered function not implemented

    ERR_DEBUG_TEST,         ///< Just an error for debugging purpose.
    ERR_DEBUG_IMPOSIBLE,    ///< Just for debugging, situation that can't happens

    ERR_DOCFILE,            ///< Error in docfile format

    ERR_GRID,               ///< Grid general error.
    ERR_GRID_SIZE,          ///< Incorrect number of GRID volumes or shapes

    ERR_IMG_NOREAD,         ///< Cannot read image from file.
    ERR_IMG_NOWRITE,        ///< Cannot write image to file.
    ERR_IMG_UNKNOWN,        ///< Unknown image type

    ERR_INDEX_OUTOFBOUNDS,  ///< Index out of bounds.

    ERR_IO,                 ///< Input/Output general error.
    ERR_IO_NOCLOSED,        ///< File cannot be closed.
    ERR_IO_NOTEXIST,        ///< File or directory does not exists.
    ERR_IO_NOTOPEN,         ///< File cannot be open.
    ERR_IO_NOPERM,          ///< Insufficient permissions to perform operation.
    ERR_IO_NOREAD,          ///< Couldn't read from file.
    ERR_IO_NOWRITE,         ///< Couldn't write to file.
    ERR_IO_NOTFILE,         ///< It is not a file.
    ERR_IO_NOTDIR,          ///< It is not a directory.
    ERR_IO_NOPATH,          ///< Environment PATH cannot be read.
    ERR_IO_LOCKED,	    	///< Error when locking/unloking a file.
    ERR_IO_SIZE,            ///< Incorrect file size.

    ERR_MATRIX,             ///< Matrix error.
    ERR_MATRIX_DIM,         ///< Problem with matrix dimensions.
    ERR_MATRIX_EMPTY,       ///< The matrix is empty.
    ERR_MATRIX_SIZE,        ///< Problem with matrix size.

    ERR_MD,                 ///< MetaData error.
    ERR_MD_NOACTIVE,        ///< No active object in MetaData.
    ERR_MD_NOOBJ,           ///< No exist requested object.
    ERR_MD_BADLABEL,        ///< Unexpected label.
    ERR_MD_MISSINGLABEL,    ///< Missing expected label
    ERR_MD_SQL,             ///< Error in SQL of MetaData operations.
    ERR_MD_OBJECTNUMBER,    ///< Incorrect number of objects in Metadata
    ERR_MD_BADTYPE,         ///< Bad label type.
    ERR_MD_UNDEFINED,       ///< Undefined label.
    ERR_MD_BADBLOCK ,       ///< This block does not exist.

    ERR_MEM_BADREQUEST,     ///< Bad amount of memory requested.
    ERR_MEM_NOTENOUGH,      ///< There is not enough memory for allocation.
    ERR_MEM_NOTDEALLOC,     ///< Memory has not been deallocated.
    ERR_MEM_NULLPOINTER,    ///< Null pointer passed as parameter

    ERR_MMAP,               ///< Global mmap error.
    ERR_MMAP_NOTADDR,       ///< Map addressing of file has failed.

    ERR_MULTIDIM_DIM,       ///< Incorrect MultidimArray dimensions
    ERR_MULTIDIM_SIZE,      ///< Incorrect MultidimArray size
    ERR_MULTIDIM_EMPTY,     ///< MultidimArray is empty.

    ERR_NOT_IMPLEMENTED,    ///< Case or algorithm not implemented yet.

    ERR_NUMERICAL,          ///< Error related to numerical calculation.

    ERR_PARAM_INCORRECT,    ///< Parameter incorrect.
    ERR_PARAM_MISSING,      ///< Parameter missing.

    ERR_PLANS_NOCREATE,     ///< FFT Plan cannot be created.

    ERR_SELFILE,            ///< Error in docfile format

    ERR_THREADS_NOTINIT,    ///< Threads cannot be initiated.

    ERR_TYPE_INCORRECT,     ///< Incorrect type received.

    ERR_UNCLASSIFIED,       ///< Just to locate unclassified errors.

    ERR_VALUE_EMPTY,        ///< Empty value.
    ERR_VALUE_INCORRECT,    ///< Incorrect value received.
    ERR_VALUE_NOTSET,        ///< Value has not been set.

    ERR_LAST_LABEL
};


/** MACROs that proccess XMIPP exceptions
 */
#define XMIPP_TRY try{

#define XMIPP_CATCH \
	}\
    catch (XmippError &xe)\
    {\
        std::cerr << xe;\
        exit(-1);\
    }

/** Show message and exit
 *
 * This is an internal function (not to be used by programmers)
 * that shows the given message and exits with the error code.
 * This function is called when no exceptions are allowed.
 *
 */
void _Xmipp_error(const ErrorType nerr, const String& what,
                  const String &file, const long line);

/** Show message and throw exception
 *
 * This macro shows the given message and exits with the error code.
 *
 * @code
 * if (...)
 *     REPORT_ERROR(ERR_DEBUG_TEST, "Error 1");
 * @endcode
 */
#define REPORT_ERROR(nerr, ErrormMsg) throw XmippError(nerr, ErrormMsg, __FILE__, __LINE__)
/** Report error without any extra message */
//#define REPORT_ERROR(nerr) throw XmippError((ErrorType)nerr, "", __FILE__, __LINE__)

/** Exception class
 *
 * This is the class type for the errors thrown by the routines when the
 * exception handling mode is active (see Xmipp Configuration for details about
 * enabling the exception handling).
 */
class XmippError
{

private:
    //Variables to hold stack info
    char ** strings;
    size_t size;

public:
    /** Error code */
    ErrorType __errno;

    /** Message shown */
    String msg;

    /** File producing the error */
    String file;

    /** Line number */
    long line;

    /** Constructor */
    XmippError(const ErrorType nerr, const String& what,
               const String &fileArg, const long lineArg);
    /** Destructor */
    ~XmippError();

    /** Show an error */
    friend std::ostream& operator<<(std::ostream& o, XmippError& XE);

    /** Get message */
    String getMessage() const;

    /** Get Default message */
    String getDefaultMessage() const;

    /** Get default message */
    static String getDefaultMessage(ErrorType e);


};

/** Print a report warning and continue the execution.
 */
void reportWarning(const String& what);

/* @} */
#endif
