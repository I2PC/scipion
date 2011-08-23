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

#include "multidim_array.h"


// Show a complex array ---------------------------------------------------
template<>
std::ostream& operator<<(std::ostream& ostrm,
                         const MultidimArray< std::complex<double> >& v)
{
    if (v.xdim == 0)
        ostrm << "NULL MultidimArray\n";
    else
        ostrm << std::endl;

    for (size_t l = 0; l < NSIZE(v); l++)
    {
        if (NSIZE(v)>1)
            ostrm << "Image No. " << l << std::endl;
        for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
        {
            if (ZSIZE(v)>1)
                ostrm << "Slice No. " << k << std::endl;
            for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
            {
                for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
                    ostrm << A3D_ELEM(v, k, i, j) << ' ';
                ostrm << std::endl;
            }
        }
    }

    return ostrm;
}

template<>
void MultidimArray< std::complex< double > >::computeDoubleMinMax(double& minval, double& maxval) const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::computeDoubleMinMax not implemented for complex.");
}
template<>
void MultidimArray< std::complex< double > >::computeDoubleMinMaxRange(double& minval, double& maxval, size_t pos, size_t size) const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::computeDoubleMinMax not implemented for complex.");
}
template<>
void MultidimArray< std::complex< double > >::rangeAdjust(std::complex< double > minF, std::complex< double > maxF)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::rangeAdjust not implemented for complex.");
}

template<>
double MultidimArray< std::complex< double > >::computeAvg() const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::computeAvg not implemented for complex.");
}

template<>
void MultidimArray<double>::computeAvgStdev(double& avg, double& stddev) const
{
    if (NZYXSIZE(*this) <= 0)
        return;

    avg = 0;
    stddev = 0;

    double* ptr=&DIRECT_MULTIDIM_ELEM(*this,0);
    size_t nmax=(nzyxdim/4)*4;

    double val;
    for (size_t n=0; n<nmax; n+=4, ptr+=4)
    {
        val=*ptr;
        avg += val;
        stddev += val * val;
        val=*(ptr+1);
        avg += val;
        stddev += val * val;
        val=*(ptr+2);
        avg += val;
        stddev += val * val;
        val=*(ptr+3);
        avg += val;
        stddev += val * val;
    }
    for (size_t n=nmax; n<nzyxdim; ++n, ptr+=1)
    {
        val=*ptr;
        avg += val;
        stddev += val * val;
    }

    avg /= NZYXSIZE(*this);

    if (NZYXSIZE(*this) > 1)
    {
        stddev = stddev / NZYXSIZE(*this) - avg * avg;
        stddev *= NZYXSIZE(*this) / (NZYXSIZE(*this) - 1);

        // Foreseeing numerical instabilities
        stddev = sqrt(fabs(stddev));
    }
    else
        stddev = 0;
}

template<>
bool operator==(const MultidimArray< std::complex< double > >& op1, const MultidimArray< std::complex< double > >& op2)
{
	double accuracy = XMIPP_EQUAL_ACCURACY;
    if (! op1.sameShape(op2) || op1.data==NULL || op2.data == NULL)
        return false;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1)
    if (   fabs(DIRECT_MULTIDIM_ELEM(op1,n).real() -
            DIRECT_MULTIDIM_ELEM(op2,n).real() > accuracy)
            ||
           fabs(DIRECT_MULTIDIM_ELEM(op1,n).imag() -
            DIRECT_MULTIDIM_ELEM(op2,n).imag() > accuracy)
       )
        return false;
    return true;
}

