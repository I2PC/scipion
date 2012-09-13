/***************************************************************************
 *
 * Authors:  Alberto Pascual Montano (pascual@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#include "xmipp_image_over.h"

// Initialise an oversampled image (ready for work) ------------------------
void ImageOver::init(int _umin, int _umax, int _uistep,
                     int _vmin, int _vmax, int _vistep,
                     int _wmin, int _wmax, int _wistep)
{
    overvmin = _vmin;
    overumin = _umin;
    overvmax = _vmax;
    overumax = _umax;
    vistep = _vistep;
    uistep = _uistep;

    overwmin = _wmin;
    overwmax = _wmax;
    wistep = _wistep;

    //   data.initZeros((_vmax-_vmin+1)*_vistep,(_umax-_umin+1)*_uistep);
    data.initZeros((_wmax - _wmin)*_wistep + 1,(_vmax - _vmin)*_vistep + 1, (_umax - _umin)*_uistep + 1);
    STARTINGZ(data) = 0;
    STARTINGY(data) = 0;
    STARTINGX(data) = 0;
    //   STARTINGY(img)=_vmin*_vistep - (_vistep-1)/2;
    //   STARTINGX(img)=_umin*_uistep - (_uistep-1)/2;
}

void ImageOver::init(MultidimArray<double> & im, int _uistep, int _vistep, int _wistep)
{
    overumin = STARTINGX(im);
    overvmin = STARTINGY(im);
    overwmin = STARTINGZ(im);
    overumax = STARTINGX(im) + XSIZE(im) - 1;
    overvmax = STARTINGY(im) + YSIZE(im) - 1;
    overwmax = STARTINGZ(im) + ZSIZE(im) - 1;
    uistep = _uistep;
    vistep = (YSIZE(im) > 1) ? ((_vistep == 0) ? uistep : _vistep) : 0;
    wistep = (ZSIZE(im) > 1) ? ((_wistep == 0) ? uistep : _wistep) : 0;

    data.initZeros(im);
    STARTINGZ(data) = 0;
    STARTINGY(data) = 0;
    STARTINGX(data) = 0;
}

// Window ------------------------------------------------------------------
void ImageOver::window(int _v0, int _u0, int _vF, int _uF)
{
    overvmin = _v0;
    overumin = _u0;
    overvmax = _vF;
    overumax = _uF;

    int newYdim = (_vF - _v0) * vistep + 1;
    int newXdim = (_uF - _u0) * uistep + 1;
    data.setXmippOrigin();
    data.selfWindow(FIRST_XMIPP_INDEX(newYdim), FIRST_XMIPP_INDEX(newXdim),
                    LAST_XMIPP_INDEX(newYdim), LAST_XMIPP_INDEX(newXdim));
    STARTINGY(data) = 0;
    STARTINGX(data) = 0;
}

// Clear -------------------------------------------------------------------
void ImageOver::clear()
{
    overvmin = overvmax = 0;
    overumin = overumax = 0;
    overwmin = overwmax = 0;
    wistep = vistep = uistep = 1;
    Image<double>::clear();
}

// Generate the normal image by averaging ----------------------------------
void ImageOver::downsample(Image< double > *I) const
{
    IMGMATRIX(*I).resize(overwmax - overwmin + 1,
                         overvmax - overvmin + 1, overumax - overumin + 1);

    double iNorm = 1./(wistep * vistep * uistep);

    for (int k = overwmin; k <= overwmax; ++k)
        for (int i = overvmin; i <= overvmax; ++i)
            for (int j = overumin; j <= overumax; ++j)
            {
                VOLVOXEL(*I, k, i, j) = 0;
                for (int w = (k - overwmin) * wistep; w < (k + 1 - overwmin)*wistep; ++w)
                    for (int v = (i - overvmin) * vistep; v < (i + 1 - overvmin)*vistep; ++v)
                        for (int u = (j - overumin) * uistep; u < (j + 1 - overumin)*uistep; ++u)
                        {
                            VOLVOXEL(*I, k, i, j) += VOLVOXEL(*this, w, v, u);
                        }
                VOLVOXEL(*I, k, i, j) *= iNorm;
            }
}
