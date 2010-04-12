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

#include "fft.h"
#include "args.h"
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
            dMij(out, 0, j) = dMij(in, 0, j);
        // Fill rest
        for (int i = 1; i < ldim; i++)
            for (int j = 0; j < XSIZE(in); j++)
                dMij(out, i, j) = dMij(in, i, j);
    }
    else
        REPORT_ERROR(1,"ERROR: Whole2Half only implemented for 1D and 2D multidimArrays");

}

/** Convert half -> whole of (centro-symmetric) Fourier transforms 2D. -- */
void Half2Whole(const MultidimArray<std::complex<double> > &in, 
		MultidimArray<std::complex<double> > &out, int oridim)
{
    if (in.getDim() == 1)
    {
        // 1D
        out.resize(oridim);
        for (int j = 0; j < XSIZE(in); j++)
            out(j) = in(j);
        for (int j = XSIZE(in); j < oridim; j++)
	out(j) = conj(in(oridim - j));
    }
    else if (in.getDim() == 2)
    {
        // 2D
        out.resize(oridim, XSIZE(in));

        // Old part
        for (int i = 0; i < YSIZE(in); i++)
            for (int j = 0; j < XSIZE(in); j++)
                dMij(out, i, j) = dMij(in, i, j);

        // Complete first column of old part
        for (int j = YSIZE(in); j < XSIZE(in); j++)
            dMij(out, 0, j) = conj(dMij(in, 0, XSIZE(in) - j));

        // New part
        for (int i = YSIZE(in); i < oridim; i++)
        {
            dMij(out, i, 0) = conj(dMij(in, oridim - i, 0));
            for (int j = 1; j < XSIZE(in); j++)
                dMij(out, i, j) = conj(dMij(in, oridim - i, XSIZE(in) - j));
        }
    }
}

/** Convert complex -> real,imag Fourier transforms 2D. -- */
void Complex2RealImag(const MultidimArray< std::complex< double > > & in,
                      MultidimArray< double > & real,
                      MultidimArray< double > & imag)
{
    real.resize(in);
    imag.resize(in);
    Complex2RealImag(MULTIDIM_ARRAY(in), MULTIDIM_ARRAY(real),
        MULTIDIM_ARRAY(imag), MULTIDIM_SIZE(in));
}

/** Convert real,imag -> complex Fourier transforms 3D. -- */
void RealImag2Complex(const MultidimArray< double > & real,
                      const MultidimArray< double > & imag,
                      MultidimArray< std::complex< double > > & out)
{
    out.resize(real);
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
        Matrix1D<double> re(in), tmp(N), im(N), cas(N);
        out.resize(N);

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
        im.resize(in);
        out.resize(in);
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
        out.resize(N);

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
        out.resize(in);
        im.resize(in);
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
void FourierTransform(const Matrix1D<std::complex<double> > &in,
		      Matrix1D< std::complex<double> > &out)
{
    int N = XSIZE(in);
    Matrix1D<double> re(N), tmpre(N), tmpim(N), im(N), cas(N);
    out.resize(N);

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
void InverseFourierTransform(const Matrix1D< std::complex<double> > &in,
			     Matrix1D<std::complex<double> > &out)
{
    int N = XSIZE(in);
    Matrix1D<double> tmpre(N), tmpim(N), re(N), im(N), cas(N);
    out.resize(N);

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
void FourierTransformHalf(const Matrix1D<std::complex<double> > &in,
                          Matrix1D< std::complex<double> > &out)
{

    Matrix1D<std::complex <double> > aux;
    FourierTransform(in, aux);
    Whole2Half(aux, out);
}

/** Complex Inverse Fourier Transform 1D, input half of (centro-symmetric) transform ---- */
void InverseFourierTransformHalf(const Matrix1D< std::complex<double> > &in,
                                 Matrix1D<std::complex<double> > &out, int orixdim)
{

    Matrix1D< std::complex<double> > aux;
    Half2Whole(in, aux, orixdim);
    InverseFourierTransform(aux, out);
    out.setXmippOrigin();
}

/* FFT shifts ------------------------------------------------------------ */
void ShiftFFT(Matrix1D< std::complex< double > > & v,
              double xshift)
{
    double dotp, a, b, c, d, ac, bd, ab_cd;
    double xxshift = xshift / (double)XSIZE(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(v)
    {
        dotp = -2 * PI * ((double)(i) * xxshift);
        a = cos(dotp);
        b = sin(dotp);
        c = DIRECT_VEC_ELEM(v,i).real();
        d = DIRECT_VEC_ELEM(v,i).imag();
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d); // (ab_cd-ac-bd = ad+bc : but needs 4 multiplications)
        DIRECT_VEC_ELEM(v,i) = std::complex<double>(ac - bd, ab_cd - ac - bd);
    }
}

void ShiftFFT(Matrix2D< std::complex< double > > & v,
              double xshift, double yshift)
{
    double dotp, a, b, c, d, ac, bd, ab_cd;
    double xxshift = xshift / (double)XSIZE(v);
    double yyshift = yshift / (double)YSIZE(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(v)
    {
        dotp = -2 * PI * ((double)(j) * xxshift + (double)(i) * yyshift);
        a = cos(dotp);
        b = sin(dotp);
        c = DIRECT_MAT_ELEM(v,i,j).real();
        d = DIRECT_MAT_ELEM(v,i,j).imag();
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d);
        DIRECT_MAT_ELEM(v,i,j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
    }
}

void ShiftFFT(Matrix3D< std::complex< double > > & v,
              double xshift, double yshift, double zshift)
{
    double dotp, a, b, c, d, ac, bd, ab_cd;
    double xxshift = xshift / (double)XSIZE(v);
    double yyshift = yshift / (double)YSIZE(v);
    double zzshift = zshift / (double)ZSIZE(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(v)
    {
        dotp = -2 * PI * ((double)(j) * xxshift + (double)(i) * yyshift + (double)(k) * zzshift);
        a = cos(dotp);
        b = sin(dotp);
        c = DIRECT_VOL_ELEM(v,k,i,j).real();
        d = DIRECT_VOL_ELEM(v,k,i,j).imag();
        ac = a * c;
        bd = b * d;
        ab_cd = (a + b) * (c + d);
        DIRECT_VOL_ELEM(v,k,i,j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
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
        REPORT_ERROR(1,"CenterOriginFFT ERROR: only valis for 1D or 2D or 3D");
}

/* Numerical derivative of a matrix ----------------------------- */
void numerical_derivative(Matrix2D<double> &M, Matrix2D<double> &D,
                          char direction, int order, int window_size, int polynomial_order)
{

    // Set D to be a copy in shape of M
    D.initZeros(M);

    Matrix1D<double> v, rotated;
    Matrix1D<double> ans; // To obtain results

    // Wrap around version of the Savitzky-Golay coefficients
    int dim = 2 * window_size + 1;
    rotated.resize(dim);

    double *pans = ans.adaptForNumericalRecipes1D();
    double *pv = v.adaptForNumericalRecipes1D();
    double *protated = rotated.adaptForNumericalRecipes1D();

    // Calculate the Savitzky-Golay filter coeficients
    savgol(protated, 2*window_size + 1, window_size,
           window_size, order, polynomial_order);

    // Savitzky-Golay filter is returned in wrap-around style, so
    // correct it to use with the convolution routine
    Matrix1D<double> coefficients(dim);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(coefficients)
    {
        int j = i + window_size;
        if (j < dim)
            coefficients(j) = rotated(i);
        else
            coefficients(j - dim) = rotated(i);
    }

    // Apply the Savitzky-Golay filter to every row or column
    if (direction == 'x')
    {
        // For every row (values in a row are values of the X direction)
        for (int i = STARTINGY(M);i <= FINISHINGY(M);i++)
        {
            M.getRow(i, v);
            series_convolution(v, coefficients, ans, false);
            ans.setRow();
            D.setRow(i, ans);
        }
    }
    else if (direction == 'y')
    {
        // For every column (values in a column are values of the Y direction)
        for (int i = STARTINGX(M);i <= FINISHINGX(M);i++)
        {
            M.getCol(i, v);
            series_convolution(v, coefficients, ans, false);
            ans.setCol();
            D.setCol(i, ans);
        }
    }
}
