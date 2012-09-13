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

#include "xmipp_memory.h"
#include "xmipp_strings.h"

char*  askMemory(size_t memsize)
{
    char*  ptr = NULL;

    if ( memsize == 0 )
    {
        REPORT_ERROR(ERR_MEM_BADREQUEST, "Error in askMemory: Memory allocation size requested is zero!");
        return(NULL);
    }

    if ( ( ptr = (char *) calloc(1,memsize*sizeof(char)) ) == NULL )
    {
        REPORT_ERROR(ERR_MEM_NOTENOUGH, formatString("askMemory:: Memory allocation of %ld bytes failed",memsize));
        return(NULL);
    }

    //memset(ptr, 0, memsize);

    return(ptr);
}

int  freeMemory(void* ptr, size_t memsize)
{
    if ( ptr == NULL )
        return(0);

    if ( memsize < 1 )
    {
        return(-1);
    }

    free(ptr);
    ptr = NULL;
    return(0);
}
