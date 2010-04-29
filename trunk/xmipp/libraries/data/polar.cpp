/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
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
#include "polar.h"



void fourierTransformRings(Polar<double > & in,
                           Polar<std::complex<double> > &out,
                           Polar_fftw_plans &plans,
                           bool conjugated)
{
    Matrix1D<std::complex<double> > Fring;
    out.clear();
    for (int iring = 0; iring < in.getRingNo(); iring++)
    {

        plans.arrays[iring] = in.rings[iring];
        (plans.transformers[iring]).FourierTransform();
        (plans.transformers[iring]).getFourierAlias(Fring);
        if (conjugated)
            for (int i = 0; i < XSIZE(Fring); i++)
                Fring(i) = conj(Fring(i));
        out.rings.push_back(Fring);
    }
    out.mode = in.mode;
    out.ring_radius = in.ring_radius;
}

void inverseFourierTransformRings(Polar<std::complex<double> > & in,
                                  Polar<double > &out,
                                  Polar_fftw_plans &plans,
                                  bool conjugated)
{
    out.clear();
    for (int iring = 0; iring < in.getRingNo(); iring++)
    {
        (plans.transformers[iring]).setFourier(in.rings[iring]);
        (plans.transformers[iring]).inverseFourierTransform(); // fReal points to plans.arrays[iring]
        out.rings.push_back(plans.arrays[iring]);
    }
    out.mode = in.mode;
    out.ring_radius = in.ring_radius;
}

void rotationalCorrelation(const Polar<std::complex<double> > &M1,
                           const Polar<std::complex<double> > &M2,
                           Matrix1D<double> &angles,
                           XmippFftw &local_transformer)
{

    Matrix1D<std::complex<double> > Fsum;
    int nrings = M1.getRingNo();
    if (nrings != M2.getRingNo())
        REPORT_ERROR(1,"rotationalCorrelation: polar structures have unequal number of rings!");

    // Fsum should already be set with the right size in the local_transformer
    // (i.e. through a FourierTransform of corr)
    local_transformer.getFourierAlias(Fsum);
    Fsum.initZeros();

    // Multiply M1 and M2 over all rings and sum
    // Assume M2 is already complex conjugated!
    for (int iring = 0; iring < nrings; iring++)
    {
        double w = (2.* PI * M1.ring_radius[iring]);
        for (int i = 0; i < M1.getSampleNo(iring); i++)
            Fsum(i) += w * M1(iring,i) * M2(iring,i);
    }

    // Inverse FFT to get real-space correlations
    // The local_transformer should already have corr as setReal!!
    local_transformer.inverseFourierTransform();

    /// FIXME: the getReal should be getComplex for complex fourier transforms....
    angles.resize(XSIZE(local_transformer.getReal()));
    for (int i = 0; i < XSIZE(angles); i++)
        angles(i)=(double)i*360./XSIZE(angles);
}

// Compute the normalized Polar Fourier transform --------------------------
void normalizedPolarFourierTransform(const Matrix2D<double> &in,
                                     Polar< std::complex<double> > &out, bool flag,
                                     int first_ring, int last_ring, Polar_fftw_plans *&plans,
                                     int BsplineOrder)
{
    Polar<double> polarIn;
    if (BsplineOrder==1)
        polarIn.getPolarFromCartesianBSpline(in,first_ring,last_ring,1);
    else
    {
        Matrix2D<double> Maux;
        in.produceSplineCoefficients(Maux,3);
        polarIn.getPolarFromCartesianBSpline(Maux,first_ring,last_ring,BsplineOrder);
    }
    double mean = polarIn.computeSum(true);
    double stddev = polarIn.computeSum2(true);
    stddev = sqrt(stddev - mean * mean);
    polarIn -= mean;
    polarIn /= stddev;
    if (plans==NULL)
    {
        plans=new Polar_fftw_plans();
        polarIn.calculateFftwPlans(*plans);
    }
    fourierTransformRings(polarIn,out,*plans,flag);
}

// Best rotation -----------------------------------------------------------
double best_rotation(const Polar< std::complex<double> > &I1,
                     const Polar< std::complex<double> > &I2, XmippFftw &local_transformer)
{
    Matrix1D<double> angles;
    rotationalCorrelation(I1,I2,angles,local_transformer);

    // Compute the maximum of correlation (inside local_transformer)
    const MultidimArray<double> &corr=local_transformer.getReal();
    int imax=0;
    double maxval=DIRECT_MULTIDIM_ELEM(corr,0);
    double* ptr=NULL;
    unsigned long int n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(corr,n,ptr)
    if (*ptr > maxval)
    {
        maxval = *ptr;
        imax=n;
    }

    // Return the corresponding angle
    return angles(imax);
}

// Align rotationally ------------------------------------------------------
void alignRotationally(Matrix2D<double> &I1, Matrix2D<double> &I2,
                       int splineOrder, int wrap)
{
    I1.setXmippOrigin();
    I2.setXmippOrigin();

    Polar_fftw_plans *plans=NULL;
    Polar< std::complex<double> > polarFourierI2, polarFourierI1;
    normalizedPolarFourierTransform(I1,polarFourierI1,false,XSIZE(I1)/5,
                                    XSIZE(I1)/2,plans);
    normalizedPolarFourierTransform(I2, polarFourierI2, true, XSIZE(I2)/5,
                                    XSIZE(I2)/2,plans);

    XmippFftw local_transformer;
    Matrix1D<double> rotationalCorr;
    rotationalCorr.resize(2*polarFourierI2.getSampleNoOuterRing()-1);
    local_transformer.setReal(rotationalCorr);
    double bestRot = best_rotation(polarFourierI1,polarFourierI2,
                                   local_transformer);

    if (splineOrder==1)
        I2.selfRotate(-bestRot,wrap);
    else
        I2.selfRotateBSpline(3,-bestRot,wrap);
}
