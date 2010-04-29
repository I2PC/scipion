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

#ifndef _XMIPP_MEMORY
#define _XMIPP_MEMORY

#include "error.h"

/* Memory managing --------------------------------------------------------- */
///@defgroup MemoryManaging Memory management for numerical recipes
/// @ingroup DataLibrary
//@{
/** Ask memory for any type vector.
    The valid values range from v[nl] to v[nh]. If no memory is available
    an exception is thrown. NULL is returned if nh is not greater than nl*/
template <class T> void ask_Tvector(T* &v, int nl, int nh)
{
    if (nh - nl + 1 > 1)
    {
        v = (T *)malloc((unsigned)(nh - nl + 1) * sizeof(T));
        if (!v) REPORT_ERROR(1, "allocation failure in vector()");
        v -= nl;
    }
    else v = NULL;
}

/** Free memory associated to any type vector.
    After freeing v=NULL*/
template <class T> void free_Tvector(T* &v, int nl, int nh)
{
    if (v != NULL)
    {
        free((char*)(v + nl));
        v = NULL;
    }
}

/** Ask memory for any type matrix.
    The valid values range from v[nrl][ncl] to v[nrh][nch].
    If no memory is available an exception is thrown. NULL is returned if any
    nh is not greater than its nl*/
template <class T> void ask_Tmatrix(T ** &m, int nrl, int nrh,
                                    int ncl, int nch)
{
    if (nrh - nrl + 1 > 1 && nch - ncl + 1 > 1)
    {
        m = (T **) malloc((unsigned)(nrh - nrl + 1) * sizeof(T*));
        if (!m) REPORT_ERROR(1, "allocation failure 1 in matrix()");
        m -= nrl;

        for (int i = nrl; i <= nrh; i++)
        {
            m[i] = (T *) malloc((unsigned)(nch - ncl + 1) * sizeof(T));
            if (!m[i]) REPORT_ERROR(1, "allocation failure 2 in matrix()");
            m[i] -= ncl;
        }
    }
    else m = NULL;
}

/** Free memory associated to any type matrix.
    After freeing v=NULL*/
template <class T> void free_Tmatrix(T ** &m, int nrl, int nrh,
                                     int ncl, int nch)
{
    if (m != NULL)
    {
        for (int i = nrh; i >= nrl; i--) free((char*)(m[i] + ncl));
        free((char*)(m + nrl));
        m = NULL;
    }
}

/** Ask memory for any type voliume.
    The valid values range from v[nsl][nrl][ncl] to v[nsh][nrh][nch].
    If no memory is available an exception is thrown. NULL is returned if any
    nh is not greater than its nl. */
template <class T> void ask_Tvolume(T *** &m, int nsl, int nsh, int nrl,
                                    int nrh, int ncl, int nch)
{
    if (nsh - nsl + 1 > 1 && nrh - nrl + 1 > 1 && nch - ncl + 1 > 1)
    {
        m = (T ***) malloc((unsigned)(nsh - nsl + 1) * sizeof(T**));
        if (!m) REPORT_ERROR(1, "allocation failure 1 in matrix()");
        m -= nsl;

        for (int k = nsl; k <= nsh; k++)
        {
            m[k] = (T **) malloc((unsigned)(nrh - nrl + 1) * sizeof(T*));
            if (!m[k]) REPORT_ERROR(1, "allocation failure 2 in matrix()");
            m[k] -= nrl;

            for (int i = nrl; i <= nrh; i++)
            {
                m[k][i] = (T *) malloc((unsigned)(nch - ncl + 1) * sizeof(T));
                if (!m[k][i]) REPORT_ERROR(1, "allocation failure 2 in matrix()");
                m[k][i] -= ncl;
            }
        }
    }
    else m = NULL;
}

/** Free memory associated to any type volume.
    After freeing v=NULL*/
template <class T> void free_Tvolume(T *** &m, int nsl, int nsh,
                                     int nrl, int nrh, int ncl, int nch)
{
    if (m != NULL)
    {
        for (int k = nsh; k >= nsl; k--)
        {
            for (int i = nrh; i >= nrl; i--) free((char*)(m[k][i] + ncl));
            free((char*)(m[k] + nrl));
        }
        free((char*)(m + nsl));
        m = NULL;
    }
}
//@}
#endif

