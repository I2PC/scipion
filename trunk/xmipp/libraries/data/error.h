/***************************************************************************
*
* Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
*  e-mail address 'xmipp@cnb.uam.es'
***************************************************************************/

#ifndef ERROR_H
#define ERROR_H

#include <string>
#include <iostream>

/** @defgroup ErrorHandling Error handling
 *  @ingroup DataLibrary
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
 *     REPORT_ERROR(6001, "Volume too small to be stored");
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
 * catch (Xmipp_error XE)
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
 */

/** Show message and exit
 * @ingroup ErrorHandling
 *
 * This is an internal function (not to be used by programmers)
 * that shows the given message and exits with the error code.
 * This function is called when no exceptions are allowed.
 *
 */
void _Xmipp_error(const int nerr, const std::string& what,
    const std::string &file, const long line);

/** Show message and exit
 * @ingroup ErrorHandling
 *
 * This macro shows the given message and exits with the error code.
 *
 * @code
 * if (...)
 *     EXIT_ERROR(1, "Error 1");
 * @endcode
 */
#define EXIT_ERROR(nerr, ErrormMsg) _Xmipp_error(nerr, ErrormMsg, __FILE__, __LINE__)

/** Show message and throw exception
 * @ingroup ErrorHandling
 *
 * This macro shows the given message and exits with the error code.
 *
 * @code
 * if (...)
 *     REPORT_ERROR(1, "Error 1");
 * @endcode
 */
#define REPORT_ERROR(nerr, ErrormMsg) throw Xmipp_error(nerr, ErrormMsg, __FILE__, __LINE__)

/** Exception class
 * @ingroup ErrorHandling
 *
 * This is the class type for the errors thrown by the routines when the
 * exception handling mode is active (see Xmipp Configuration for details about
 * enabling the exception handling).
 */

class Xmipp_error
{
public:
    /** Error code */
    int __errno;

    /** Message shown */
    std::string msg;

    /** File produstd::cing the error */
    std::string file;

    /** Line number */
    long line;

    Xmipp_error(const int nerr, const std::string& what,
        const std::string &fileArg, const long lineArg);
    friend std::ostream& operator<<(std::ostream& o, Xmipp_error& XE);
};

#endif
