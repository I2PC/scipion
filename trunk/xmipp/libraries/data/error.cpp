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
Xmipp_error::Xmipp_error(const int nerr, const std::string &what,
    	    	    	 const std::string &fileArg, const long lineArg)
{
    __errno = nerr;
    msg = what;
    file= fileArg;
    line=lineArg;
}

// Show message
std::ostream& operator << (std::ostream& o, Xmipp_error& XE)
{
    o << XE.__errno << ":" << XE.msg << std::endl
      << "File: " << XE.file << " line: " << XE.line << std::endl;
    return o;
}
