/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
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

#include "filter.h"

#include <qinputdialog.h>

#include <data/filters.h>
#include <data/histogram.h>

/* Generic filter ---------------------------------------------------------- */

QString QtFilter::name = "Generic filter";

/* Constructor ------------------------------------------------------------- */
QtFilter::QtFilter(const Micrograph *_M)
{
    __active = false;
    __M      = _M;
}

#define I (*_img)()
/* Invert Contrast  -------------------------------------------------------- */
QString QtInvertContrastFilter::name = "Invert contrast";

void QtInvertContrastFilter::apply(Image<double> *_img)
{
    if (_img != NULL)
    {
        I *= -1;
        I += 255;
    }
}

/* Enhance Contrast  ------------------------------------------------------- */
QString QtEnhanceContrastFilter::name = "Enhance contrast";

void QtEnhanceContrastFilter::apply(Image<double> *_img)
{
    if (_img != NULL)
    {
        contrast_enhancement(_img);
    }
}

/* Substract background ---------------------------------------------------- */
QString QtSubstractBackgroundFilter::name = "Substract background";

void QtSubstractBackgroundFilter::apply(Image<double> *_img)
{
    if (_img != NULL)
    {
        substract_background_plane(IMGMATRIX(*_img));
        I.rangeAdjust(0, 255);
    }
}

/* Remove Outliers --------------------------------------------------------- */
QString QtRemoveOutlierFilter::name = "Remove Outlier values";

void QtRemoveOutlierFilter::apply(Image<double> *_img)
{
    if (_img != NULL)
    {
        histogram1D h;
        compute_hist(I, h, 100);
        double th = h.percentil(5);
        I.threshold("below", th, th);
        th = h.percentil(95);
        I.threshold("above", th, th);
        I.rangeAdjust(0, 255);
    }
}

/* LowPass ----------------------------------------------------------------- */
QString QtLowPassFilter::name = "Lowpass filter";

QtLowPassFilter::QtLowPassFilter(const Micrograph *_M): QtFilter(_M)
{
    bool ok = FALSE;
    filter.FilterShape = RAISED_COSINE;
    filter.FilterBand = LOWPASS;
    filter.raised_w = 0.02;
    filter.w1 = QInputDialog::getDouble("LowPass filter", "Digital freq.(<1/2)", 0.2);
}

void QtLowPassFilter::apply(Image<double> *_img)
{
    if (_img != NULL)
    {
        (*_img)().setXmippOrigin();
        filter.generate_mask((*_img)());
        filter.apply_mask_Space((*_img)());
        FOR_ALL_ELEMENTS_IN_ARRAY2D((*_img)())
        (*_img)()(i, j) = CLIP((*_img)()(i, j), 0, 255);
        STARTINGX((*_img)()) = 0;
        STARTINGY((*_img)()) = 0;
    }
}

/* HighPass ----------------------------------------------------------------- */
QString QtHighPassFilter::name = "Highpass filter";

QtHighPassFilter::QtHighPassFilter(const Micrograph *_M): QtFilter(_M)
{
    bool ok = FALSE;
    filter.FilterShape = RAISED_COSINE;
    filter.FilterBand = HIGHPASS;
    filter.raised_w = 0.02;
    filter.w1 = QInputDialog::getDouble("HighPass filter" , "Digital freq.(<1/2)", 0.2);
}

void QtHighPassFilter::apply(Image<double> *_img)
{
    if (_img != NULL)
    {
        (*_img)().setXmippOrigin();
        filter.generate_mask((*_img)());
        filter.apply_mask_Space((*_img)());
        (*_img)().rangeAdjust(0, 255);
        STARTINGX((*_img)()) = 0;
        STARTINGY((*_img)()) = 0;
    }
}
