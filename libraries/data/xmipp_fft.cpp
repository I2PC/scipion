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

#include "xmipp_fft.h"
#include "args.h"
#include "histogram.h"
#include "external/bilib/headers/dft.h"

/* Format conversions ------------------------------------------------------ */
/** Convert whole -> half of (centro-symmetric) Fourier transforms 1D. -- */
void Whole2Half(const MultidimArray<std::complex<double> > &in,
                MultidimArray<std::complex<double> > &out)
{
    if (in.getDim() == 1)
    {
        // 1D
        int ldim = (int)(XSIZE(in) / 2) + 1;
        out.resize(ldim);
        for (int j = 0; j < ldim; j++)
            out(j) = in(j);
    }
    else if (in.getDim() == 2)
    {
        // 2D
        // This assumes squared images...
        int ldim = (int)(YSIZE(in) / 2) + 1;

        out.initZeros(ldim, XSIZE(in));
        // Fill first column only half
        for (int j = 0; j < ldim; j++)
            dAij(out, 0, j) = dAij(in, 0, j);
        // Fill rest
        for (int i = 1; i < ldim; i++)
            for (size_t j = 0; j < XSIZE(in); j++)
                dAij(out, i, j) = dAij(in, i, j);
    }
    else
        REPORT_ERROR(ERR_MULTIDIM_DIM,"ERROR: Whole2Half only implemented for 1D and 2D multidimArrays");

}

/** Convert half -> whole of (centro-symmetric) Fourier transforms 2D. -- */
void Half2Whole(const MultidimArray<std::complex<double> > &in,
                MultidimArray<std::complex<double> > &out, size_t oridim)
{
    if (in.getDim() == 1)
    {
        // 1D
        out.resizeNoCopy(oridim);
        for (size_t j = 0; j < XSIZE(in); j++)
            DIRECT_A1D_ELEM(out,j) = DIRECT_A1D_ELEM(in,j);
        for (size_t j = XSIZE(in); j < oridim; j++)
            DIRECT_A1D_ELEM(out,j) = conj(DIRECT_A1D_ELEM(in,oridim - j));
    }
    else if (in.getDim() == 2)
    {
        // 2D
        out.resizeNoCopy(oridim, XSIZE(in));

        // Old part
        for (size_t i = 0; i < YSIZE(in); i++)
            for (size_t j = 0; j < XSIZE(in); j++)
                dAij(out, i, j) = dAij(in, i, j);

        // Complete first column of old part
        for (size_t j = YSIZE(in); j < XSIZE(in); j++)
            dAij(out, 0, j) = conj(dAij(in, 0, XSIZE(in) - j));

        // New part
        for (size_t i = YSIZE(in); i < oridim; i++)
        {
            dAij(out, i, 0) = conj(dAij(in, oridim - i, 0));
            for (size_t j = 1; j < XSIZE(in); j++)
                dAij(out, i, j) = conj(dAij(in, oridim - i, XSIZE(in) - j));
        }
    }
}

/** Convert complex -> real,imag Fourier transforms 2D. -- */
void Complex2RealImag(const MultidimArray< std::complex< double > > & in,
                      MultidimArray< double > & real,
                      MultidimArray< double > & imag)
{
    real.resizeNoCopy(in);
    imag.resizeNoCopy(in);
    Complex2RealImag(MULTIDIM_ARRAY(in), MULTIDIM_ARRAY(real),
                     MULTIDIM_ARRAY(imag), MULTIDIM_SIZE(in));
}

/** Convert real,imag -> complex Fourier transforms 3D. -- */
void RealImag2Complex(const MultidimArray< double > & real,
                      const MultidimArray< double > & imag,
                      MultidimArray< std::complex< double > > & out)
{
    out.resizeNoCopy(real);
    RealImag2Complex(MULTIDIM_ARRAY(real), MULTIDIM_ARRAY(imag),
                     MULTIDIM_ARRAY(out), MULTIDIM_SIZE(real));
}

/** Direct Fourier Transform nD ------------------------------------------- */
void FourierTransform(const MultidimArray<double> &in,
                      MultidimArray< std::complex<double> > &out)
{
    if ( in.getDim() == 1 )
    {
        // 1D
        int N = XSIZE(in);
        MultidimArray<double> re(in), tmp(N), im(N), cas(N);
        out.resizeNoCopy(N);

        GetCaS(MULTIDIM_ARRAY(cas), N);
        DftRealToRealImaginary(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                               MULTIDIM_ARRAY(tmp), MULTIDIM_ARRAY(cas), N);
        RealImag2Complex(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_ARRAY(out), N);
    }
    else
    {
        // 2D and 3D
        int Status;
        MultidimArray<double> re(in), im;
        im.resizeNoCopy(in);
        out.resizeNoCopy(in);
        VolumeDftRealToRealImaginary(MULTIDIM_ARRAY(re),
                                     MULTIDIM_ARRAY(im), XSIZE(in), YSIZE(in), ZSIZE(in), &Status);
        RealImag2Complex(re,im,out);
    }

}

/** Inverse Fourier Transform nD. ----------------------------------------- */
void InverseFourierTransform(const MultidimArray< std::complex<double> > &in,
                             MultidimArray<double> &out)
{
    if ( in.getDim() == 1 )
    {
        // 1D
        int N = XSIZE(in);
        MultidimArray<double> tmp(N), im(N), cas(N);
        out.resizeNoCopy(N);

        GetCaS(MULTIDIM_ARRAY(cas), N);
        Complex2RealImag(MULTIDIM_ARRAY(in), MULTIDIM_ARRAY(out),
                         MULTIDIM_ARRAY(im), N);
        InvDftRealImaginaryToReal(MULTIDIM_ARRAY(out), MULTIDIM_ARRAY(im),
                                  MULTIDIM_ARRAY(tmp), MULTIDIM_ARRAY(cas), N);
    }
    else
    {
        // 2D and 3D
        int Status;
        MultidimArray<double> im;
        out.resizeNoCopy(in);
        im.resizeNoCopy(in);
        Complex2RealImag(in, out, im);
        VolumeInvDftRealImaginaryToReal(MULTIDIM_ARRAY(out),
                                        MULTIDIM_ARRAY(im),
                                        XSIZE(in), YSIZE(in), ZSIZE(in), &Status);
    }
}


/** Direct Fourier Transform 1D / 2D, output half of (centro-symmetric) transform ---- */
void FourierTransformHalf(const MultidimArray<double> &in,
                          MultidimArray< std::complex<double> > &out)
{

    MultidimArray<std::complex<double> > aux;
    FourierTransform(in, aux);
    Whole2Half(aux, out);
}

/** Inverse Fourier Transform 1D / 2D, input half of (centro-symmetric) transform ---- */
void InverseFourierTransformHalf(const MultidimArray< std::complex<double> > &in,
                                 MultidimArray<double> &out, int oridim)
{

    MultidimArray< std::complex<double> > aux;
    Half2Whole(in, aux, oridim);
    InverseFourierTransform(aux, out);
    out.setXmippOrigin();
}

/* Complex Fourier Transform ------------------------------------------------------ */

/** Complex Direct Fourier Transform 1D ------------------------------------------- */
void FourierTransform(const MultidimArray<std::complex<double> > &in,
                      MultidimArray< std::complex<double> > &out)
{
    // Only implemented for 1D transforms
    in.checkDimension(1);

    int N = XSIZE(in);
    MultidimArray<double> re(N), tmpre(N), tmpim(N), im(N), cas(N);
    out.resizeNoCopy(N);

    GetCaS(MULTIDIM_ARRAY(cas), N);
    Complex2RealImag(MULTIDIM_ARRAY(in), MULTIDIM_ARRAY(re),
                     MULTIDIM_ARRAY(im), N);
    DftRealImaginaryToRealImaginary(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                                    MULTIDIM_ARRAY(tmpre), MULTIDIM_ARRAY(tmpim),
                                    MULTIDIM_ARRAY(cas), N);
    RealImag2Complex(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                     MULTIDIM_ARRAY(out), N);
}

/** Complex Inverse Fourier Transform 1D. ----------------------------------------- */
void InverseFourierTransform(const MultidimArray< std::complex<double> > &in,
                             MultidimArray<std::complex<double> > &out)
{
    // Only implemented for 1D transforms
    in.checkDimension(1);

    int N = XSIZE(in);
    MultidimArray<double> tmpre(N), tmpim(N), re(N), im(N), cas(N);
    out.resizeNoCopy(N);

    GetCaS(MULTIDIM_ARRAY(cas), N);
    Complex2RealImag(MULTIDIM_ARRAY(in), MULTIDIM_ARRAY(re),
                     MULTIDIM_ARRAY(im), N);
    InvDftRealImaginaryToRealImaginary(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                                       MULTIDIM_ARRAY(tmpre), MULTIDIM_ARRAY(tmpim),
                                       MULTIDIM_ARRAY(cas), N);
    RealImag2Complex(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                     MULTIDIM_ARRAY(out), N);
}

/** Complex Fourier Transform 1D, output half of (centro-symmetric) transform ---- */
void FourierTransformHalf(const MultidimArray<std::complex<double> > &in,
                          MultidimArray< std::complex<double> > &out)
{
    // Only implemented for 1D transforms
    in.checkDimension(1);

    MultidimArray<std::complex <double> > aux;
    FourierTransform(in, aux);
    Whole2Half(aux, out);
}

/** Complex Inverse Fourier Transform 1D, input half of (centro-symmetric) transform ---- */
void InverseFourierTransformHalf(const MultidimArray< std::complex<double> > &in,
                                 MultidimArray<std::complex<double> > &out, int orixdim)
{
    // Only implemented for 1D transforms
    in.checkDimension(1);

    MultidimArray< std::complex<double> > aux;
    Half2Whole(in, aux, orixdim);
    InverseFourierTransform(aux, out);
    out.setXmippOrigin();
}

void centerFFT2(MultidimArray<double> &v)
{
    if (v.getDim() == 2)
    {
        //Just separe the even and odd dimensions case
        if (XSIZE(v) % 2 == 0 && YSIZE(v) % 2 == 0)
        {
            int xsize = XSIZE(v);
            int xhalf = xsize / 2;
            int yhalf = YSIZE(v) / 2;
            double * posA = MULTIDIM_ARRAY(v);
            double * posB = posA + xhalf;
            double * posC = posA + xsize * yhalf;
            double * posD = posC + xhalf;
            double * buffer=new double[xhalf];
            size_t bytes = xhalf * sizeof(double);

            for (int i = 0; i < yhalf; ++i,
                 posA += xsize, posB += xsize, posC += xsize, posD += xsize)
            {
                SWAP_ARRAY(posA, posD, bytes);
                SWAP_ARRAY(posB, posC, bytes);
            }
            delete []buffer;
        }
        else
        {
            //todo: implementation for the odd case needed
            CenterFFT(v, true);
        }
    }
    else
        std::cerr <<"bad dim: " << v.getDim() << std::endl;
}
/* FFT shifts ------------------------------------------------------------ */
void ShiftFFT(MultidimArray< std::complex< double > > & v,
              double xshift)
{
    v.checkDimension(1);
    double dotp, a, b, c, d, ac, bd, ab_cd;
    double xxshift = xshift / (double)XSIZE(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v)
    {
        dotp = -2 * PI * ((double)(i) * xxshift);
        a = cos(dotp);
        b = sin(dotp);
        c = DIRECT_A1D_ELEM(v,i).real();
        d = DIRECT_A1D_ELEM(v,i).imag();
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d); // (ab_cd-ac-bd = ad+bc : but needs 4 multiplications)
        DIRECT_A1D_ELEM(v,i) = std::complex<double>(ac - bd, ab_cd - ac - bd);
    }
}

void ShiftFFT(MultidimArray< std::complex< double > > & v,
              double xshift, double yshift)
{
    v.checkDimension(2);
    double dotp, a, b, c, d, ac, bd, ab_cd;
    double xxshift = xshift / (double)XSIZE(v);
    double yyshift = yshift / (double)YSIZE(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(v)
    {
        dotp = -2 * PI * ((double)(j) * xxshift + (double)(i) * yyshift);
        a = cos(dotp);
        b = sin(dotp);
        c = DIRECT_A2D_ELEM(v,i,j).real();
        d = DIRECT_A2D_ELEM(v,i,j).imag();
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d);
        DIRECT_A2D_ELEM(v,i,j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
    }
}

void ShiftFFT(MultidimArray< std::complex< double > > & v,
              double xshift, double yshift, double zshift)
{
    v.checkDimension(3);
    double dotp, a, b, c, d, ac, bd, ab_cd;
    double xxshift = -2 * PI * xshift / (double)XSIZE(v);
    double yyshift = -2 * PI * yshift / (double)YSIZE(v);
    double zzshift = -2 * PI * zshift / (double)ZSIZE(v);
    for (size_t k=0; k<ZSIZE(v); ++k)
    {
    	double zdot=(double)(k) * zzshift;
        for (size_t i=0; i<YSIZE(v); ++i)
        {
        	double zydot=zdot+(double)(i) * yyshift;
            for (size_t j=0; j<XSIZE(v); ++j)
            {
            	double *ptrv_kij=(double *)&DIRECT_A3D_ELEM(v,k,i,j);
                dotp = (double)(j) * xxshift + zydot;
                sincos(dotp,&b,&a);
                c = *ptrv_kij;
                d = *(ptrv_kij+1);
                ac = a * c;
                bd = b * d;
                ab_cd = (a + b) * (c + d);
                *ptrv_kij = ac - bd;
                *(ptrv_kij+1) = ab_cd - ac - bd;
            }
        }
    }
}

/* Position origin at center ----------------------------------------------- */
void CenterOriginFFT(MultidimArray< std::complex< double > > & v, bool forward)
{
    if ( v.getDim() == 1 )
    {
        // 1D
        double xshift = -(double)(int)(XSIZE(v) / 2);
        if (forward)
        {
            ShiftFFT(v, xshift);
            CenterFFT(v, forward);
        }
        else
        {
            CenterFFT(v, forward);
            ShiftFFT(v, -xshift);
        }
    }
    else if ( v.getDim() == 2)
    {
        // 2D

        double xshift = -(double)(int)(XSIZE(v) / 2);
        double yshift = -(double)(int)(YSIZE(v) / 2);
        if (forward)
        {
            ShiftFFT(v, xshift, yshift);
            CenterFFT(v, forward);
        }
        else
        {
            CenterFFT(v, forward);
            ShiftFFT(v, -xshift, -yshift);
        }
    }
    else if ( v.getDim() == 3)
    {
        // 3D
        double xshift = -(double)(int)(XSIZE(v) / 2);
        double yshift = -(double)(int)(YSIZE(v) / 2);
        double zshift = -(double)(int)(ZSIZE(v) / 2);
        if (forward)
        {
            ShiftFFT(v, xshift, yshift, zshift);
            CenterFFT(v, forward);
        }
        else
        {
            CenterFFT(v, forward);
            ShiftFFT(v, -xshift, -yshift, -zshift);
        }
    }
    else
        REPORT_ERROR(ERR_MULTIDIM_DIM,"CenterOriginFFT ERROR: only valis for 1D or 2D or 3D");
}

/* Xmipp image -> Xmipp PSD ------------------------------------------------ */
void xmipp2PSD(const MultidimArray<double> &input, MultidimArray<double> &output,
               bool takeLog)
{
    output = input;
    CenterFFT(output, true);
    double min_val = output.computeMax();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(output)
    {
        double pixval=A2D_ELEM(output,i,j);
        if (pixval > 0 && pixval < min_val)
            min_val = pixval;
    }
    if (takeLog)
    {
        min_val = 10 * log10(min_val);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(output)
        {
            double pixval=A2D_ELEM(output,i,j);
            if (pixval > 0)
                A2D_ELEM(output,i,j) = 10 * log10(pixval);
            else
                A2D_ELEM(output,i,j) = min_val;
        }
    }
    reject_outliers(output);
}

/* Xmipp image -> Xmipp CTF ------------------------------------------------ */
void xmipp2CTF(const MultidimArray<double> &input, MultidimArray<double> &output)
{
    output = input;
    CenterFFT(output, true);

    // Prepare PSD part
    double min_val = output(0, XSIZE(output) - 1);
    double max_val = min_val;
    bool first = true;
    int Xdim = XSIZE(output);
    int Ydim = YSIZE(output);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(output)
    {
        if ((i < Ydim / 2 && j >= Xdim / 2) || (i >= Ydim / 2 && j < Xdim / 2))
        {
            if (output(i, j) > XMIPP_EQUAL_ACCURACY &&
                (output(i, j) < min_val || first))
                min_val = output(i, j);
            if (output(i, j) > XMIPP_EQUAL_ACCURACY &&
                (output(i, j) > max_val || first))
            {
                max_val = output(i, j);
                first = false;
            }
        }
    }
    MultidimArray<double> left(YSIZE(output), XSIZE(output));
    min_val = 10 * log10(min_val);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(output)
    {
        if ((i < Ydim / 2 && j >= Xdim / 2) || (i >= Ydim / 2 && j < Xdim / 2))
        {
            if (output(i, j) > XMIPP_EQUAL_ACCURACY)
                left(i, j) = 10 * log10(output(i, j));
            else
                left(i, j) = min_val;
        }
    }
    reject_outliers(left);

    // Join both parts
    FOR_ALL_ELEMENTS_IN_ARRAY2D(output)
    if ((i < Ydim / 2 && j >= Xdim / 2) || (i >= Ydim / 2 && j < Xdim / 2))
        output(i, j) = left(i, j);
    else
        output(i, j) = ABS(output(i, j));
}
