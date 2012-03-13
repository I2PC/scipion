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

void MultidimArrayBase::setNdim(int Ndim)
{
    ndim = Ndim;
    nzyxdim=zyxdim*ndim;
}

/** Sets new Z dimension.
 *
 *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
 *
 */
void MultidimArrayBase::setZdim(int Zdim)
{
    zdim = Zdim;
    zyxdim=yxdim*zdim;
    nzyxdim=zyxdim*ndim;
}

/** Sets new Y dimension.
 *
 *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
 *
 */
void MultidimArrayBase::setYdim(int Ydim)
{
    ydim = Ydim;
    yxdim=(size_t)ydim*xdim;
    zyxdim=yxdim*zdim;
    nzyxdim=zyxdim*ndim;
}

/** Sets new X dimension.
  *
  *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
  *
  */
void MultidimArrayBase::setXdim(int Xdim)
{
    xdim = Xdim;
    yxdim=(size_t)ydim*xdim;
    zyxdim=yxdim*zdim;
    nzyxdim=zyxdim*ndim;
}

/** Sets new 4D dimensions.
  *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
  */
void MultidimArrayBase::setDimensions(int Xdim, int Ydim, int Zdim, int Ndim)
{
    if (((size_t)Xdim)*Ydim*Zdim*Ndim < 1)
        REPORT_ERROR(ERR_MULTIDIM_SIZE, "Dimensions' size cannot be zero nor negative.");
    ndim=Ndim;
    zdim=Zdim;
    ydim=Ydim;
    xdim=Xdim;
    yxdim=ydim*xdim;
    zyxdim=zdim*yxdim;
    nzyxdim=ndim*zyxdim;
}

/** Sets new 4D dimensions.
 *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
 */
void MultidimArrayBase::setDimensions(ArrayDim &newDim)
{
    if (newDim.ndim*newDim.zdim*newDim.ydim*newDim.xdim < 1)
        REPORT_ERROR(ERR_MULTIDIM_SIZE, "Dimensions' size cannot be zero nor negative.");
    ndim = newDim.ndim;
    zdim = newDim.zdim;
    ydim = newDim.ydim;
    xdim = newDim.xdim;

    newDim.yxdim   = yxdim   = ydim*xdim;
    newDim.zyxdim  = zyxdim  = zdim*yxdim;
    newDim.nzyxdim = nzyxdim = ndim*zyxdim;
}

/** Returns the multidimArray N,Z, Y and X dimensions.
 *
 * @code
 * V.getDimensions(Xdim, Ydim, Zdim, Ndim);
 * @endcode
 */
void MultidimArrayBase::getDimensions(int& Xdim, int& Ydim, int& Zdim, size_t &Ndim) const
{
    Xdim = xdim;
    Ydim = ydim;
    Zdim = zdim;
    Ndim = ndim;
}

void MultidimArrayBase::getDimensions(ArrayDim &adim) const
{
    adim.xdim = xdim;
    adim.ydim = ydim;
    adim.zdim = zdim;
    adim.ndim = ndim;
    adim.yxdim = yxdim;
    adim.zyxdim = zyxdim;
    adim.nzyxdim = nzyxdim;
}

ArrayDim MultidimArrayBase::getDimensions() const
{
    ArrayDim adim;
    adim.xdim = xdim;
    adim.ydim = ydim;
    adim.zdim = zdim;
    adim.ndim = ndim;
    adim.yxdim = yxdim;
    adim.zyxdim = zyxdim;
    adim.nzyxdim = nzyxdim;

    return adim;
}

/** Get dimensions.
 *
 * Returns the size of the object in a 4D vector. If the object is a matrix
 * or a vector, then the higher order dimensions will be set to 1, ie,
 * (Xdim, 1, 1) or (Xdim, Ydim, 1).
 *
 * This function is not ported to Python.
 */
void MultidimArrayBase::getDimensions(int* size) const
{
    size[0] = xdim;
    size[1] = ydim;
    size[2] = zdim;
    size[3] = ndim;
}

/** Returns the total size of the multidimArray
 *
 * @code
 * if (V.getSize() > 1) ...
 * @endcode
 */
size_t MultidimArrayBase::getSize() const
{
    return nzyxdim;
}

/** Resize the multidimarray from an ArrayDim struct
 *
 */
void MultidimArrayBase::resize(ArrayDim &adim,bool copy)
{
    resize(adim.ndim, adim.zdim, adim.ydim, adim.xdim, copy);
}

/** Copy the shape parameters
  *
  */
void MultidimArrayBase::copyShape(const MultidimArrayBase &m)
{
    ndim=m.ndim;
    zdim=m.zdim;
    ydim=m.ydim;
    xdim=m.xdim;
    yxdim=m.yxdim;
    zyxdim=m.zyxdim;
    nzyxdim=m.nzyxdim;
    zinit=m.zinit;
    yinit=m.yinit;
    xinit=m.xinit;
}

void MultidimArrayBase::setXmippOrigin()
{
    zinit = FIRST_XMIPP_INDEX(zdim);
    yinit = FIRST_XMIPP_INDEX(ydim);
    xinit = FIRST_XMIPP_INDEX(xdim);
}

void MultidimArrayBase::resetOrigin()
{
    zinit = yinit = xinit = 0;
}

void MultidimArrayBase::moveOriginTo(int k, int i, int j)
{
    zinit = k + FIRST_XMIPP_INDEX(zdim);
    yinit = i + FIRST_XMIPP_INDEX(ydim);
    xinit = j + FIRST_XMIPP_INDEX(xdim);
}

void MultidimArrayBase::moveOriginTo(int i, int j)
{
    yinit = i + FIRST_XMIPP_INDEX(ydim);
    xinit = j + FIRST_XMIPP_INDEX(xdim);
}

/** IsCorner (in 2D or 3D matrix)
         *
         * TRUE if the logical index given is a corner of the definition region of this
         * array.
         */
bool MultidimArrayBase::isCorner(const Matrix1D< double >& v) const
{

    if (v.size() < 2)
        REPORT_ERROR(ERR_MATRIX_SIZE, "isCorner: index vector has got not enough components");

    else if (ZSIZE(*this)==1)
        return ((XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this))  ||
                (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this)) ||
                (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this))  ||
                (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this)));
    else if (ZSIZE(*this)>1)
        return ((XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this)  && ZZ(v) == STARTINGZ(*this)) ||
                (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this) && ZZ(v) == STARTINGZ(*this)) ||
                (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this)  && ZZ(v) == STARTINGZ(*this))  ||
                (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this) && ZZ(v) == STARTINGZ(*this)) ||
                (XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this)  && ZZ(v) == FINISHINGZ(*this)) ||
                (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this) && ZZ(v) == FINISHINGZ(*this)) ||
                (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this)  && ZZ(v) == FINISHINGZ(*this))  ||
                (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this) && ZZ(v) == FINISHINGZ(*this)));
    else
        REPORT_ERROR(ERR_MATRIX_SIZE, formatString("isCorner: index vector has too many components. dimV= %lu matrix dim = %i", v.size(), XSIZE(*this)));
}

/** Outside
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
bool MultidimArrayBase::outside(const Matrix1D<double> &r) const
{
    if (r.size() < 1)
    {
        REPORT_ERROR(ERR_MATRIX_SIZE, "Outside: index vector has not got enough components");
    }
    else if (r.size()==1)
    {
        return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this));
    }
    else if (r.size()==2)
    {
        return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this) ||
                YY(r) < STARTINGY(*this) || YY(r) > FINISHINGY(*this));
    }
    else if (r.size()==3)
    {
        return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this) ||
                YY(r) < STARTINGY(*this) || YY(r) > FINISHINGY(*this) ||
                ZZ(r) < STARTINGZ(*this) || ZZ(r) > FINISHINGZ(*this));
    }
    else
        REPORT_ERROR(ERR_MATRIX_SIZE,"Outside: index vector has too many components");
}

void MultidimArrayBase::printShape(std::ostream& out) const
{
    if (NSIZE(*this) > 1)
        out << " Number of images = "<<NSIZE(*this);

    int dim = getDim();
    if (dim == 3)
        out<< " Size(Z,Y,X): " << ZSIZE(*this) << "x" << YSIZE(*this) << "x" << XSIZE(*this)
        << " k=[" << STARTINGZ(*this) << ".." << FINISHINGZ(*this) << "]"
        << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
        << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
    else if (dim == 2)
        out<< " Size(Y,X): " << YSIZE(*this) << "x" << XSIZE(*this)
        << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
        << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
    else if (dim == 1)
        out<< " Size(X): " << XSIZE(*this)
        << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
    else
        out << " Empty MultidimArray!";
    out<<"\n";
}

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

template<>
void MultidimArray< std::complex< double > >::getReal(MultidimArray<double> & realImg) const
{
    if (NZYXSIZE(*this) == 0)
    {
        realImg.clear();
        return;
    }

    realImg.resizeNoCopy(*this);
    double* ptr2 = MULTIDIM_ARRAY(realImg);
    double * ptr1 = (double*) MULTIDIM_ARRAY(*this);

    // Unroll the loop
    const size_t unroll=4;
    size_t nmax=(NZYXSIZE(*this)/unroll)*unroll;
    for (size_t n=0; n<nmax; n+=unroll)
    {
        DIRECT_MULTIDIM_ELEM(realImg, n)   = static_cast<double>(*ptr1++);
        ptr1++;
        DIRECT_MULTIDIM_ELEM(realImg, n+1) = static_cast<double>(*(ptr1++));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(realImg, n+2) = static_cast<double>(*(ptr1++));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(realImg, n+3) = static_cast<double>(*(ptr1++));
        ptr1++;
    }
    // Do the remaining elements
    for (size_t n=nmax; n<NZYXSIZE(*this); ++n)
    {
        DIRECT_MULTIDIM_ELEM(realImg, n) = static_cast<double>(*ptr1++);
        ptr1++;
    }

}

template<>
void MultidimArray< std::complex< double > >::getImag(MultidimArray<double> & imagImg) const
{
    if (NZYXSIZE(*this) == 0)
    {
        imagImg.clear();
        return;
    }

    imagImg.resizeNoCopy(*this);
    double* ptr2 = MULTIDIM_ARRAY(imagImg);
    double * ptr1 = (double*) MULTIDIM_ARRAY(*this);

    // Unroll the loop
    const size_t unroll=4;
    size_t nmax=(NZYXSIZE(*this)/unroll)*unroll;
    for (size_t n=0; n<nmax; n+=unroll)
    {
        DIRECT_MULTIDIM_ELEM(imagImg, n)   = static_cast<double>(*(++ptr1));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(imagImg, n+1) = static_cast<double>(*(++ptr1));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(imagImg, n+2) = static_cast<double>(*(++ptr1));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(imagImg, n+3) = static_cast<double>(*(++ptr1));
        ptr1++;
    }
    // Do the remaining elements
    for (size_t n=nmax; n<NZYXSIZE(*this); ++n)
    {
        DIRECT_MULTIDIM_ELEM(imagImg, n) = static_cast<double>(*(++ptr1));
        ptr1++;
    }

}


template<>
double MultidimArray<double>::interpolatedElement2D(double x, double y, double outside_value) const
{
    double dx0 = floor(x);
    int x0=(int)dx0;
    double fx = x - dx0;
    int x1 = x0 + 1;
    double dy0 = floor(y);
    int y0=(int)dy0;
    double fy = y - dy0;
    int y1 = y0 + 1;

    int i0=STARTINGY(*this);
    int j0=STARTINGX(*this);
    int iF=FINISHINGY(*this);
    int jF=FINISHINGX(*this);

    double d00, d10, d11, d01;
    ASSIGNVAL(d00,y0,x0);
    ASSIGNVAL(d01,y0,x1);
    ASSIGNVAL(d10,y1,x0);
    ASSIGNVAL(d11,y1,x1);

    double d0 = LIN_INTERP(fx, d00, d01);
    double d1 = LIN_INTERP(fx, d10, d11);
    return LIN_INTERP(fy, d0, d1);
}

