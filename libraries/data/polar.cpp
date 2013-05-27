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

void fourierTransformRings(Polar<double> & in,
		Polar<std::complex<double> > &out, Polar_fftw_plans &plans,
		bool conjugated) {
	MultidimArray<std::complex<double> > Fring;
	out.clear();
	for (int iring = 0; iring < in.getRingNo(); iring++) {

		plans.arrays[iring] = in.rings[iring];
		(plans.transformers[iring]).FourierTransform();
		(plans.transformers[iring]).getFourierAlias(Fring);
		if (conjugated) {
			double *ptrFring_i = (double*) &DIRECT_A1D_ELEM(Fring,0);
			++ptrFring_i;
			for (size_t i = 0; i < XSIZE(Fring); ++i, ptrFring_i += 2)
				(*ptrFring_i) *= -1;
		}
		out.rings.push_back(Fring);
	}
	out.mode = in.mode;
	out.ring_radius = in.ring_radius;
}

void inverseFourierTransformRings(Polar<std::complex<double> > & in,
		Polar<double> &out, Polar_fftw_plans &plans, bool conjugated) {
	out.clear();
	for (int iring = 0; iring < in.getRingNo(); iring++) {
		(plans.transformers[iring]).setFourier(in.rings[iring]);
		(plans.transformers[iring]).inverseFourierTransform(); // fReal points to plans.arrays[iring]
		out.rings.push_back(plans.arrays[iring]);
	}
	out.mode = in.mode;
	out.ring_radius = in.ring_radius;
}

void rotationalCorrelation(const Polar<std::complex<double> > &M1,
		const Polar<std::complex<double> > &M2, MultidimArray<double> &angles,
		RotationalCorrelationAux &aux) {
	int nrings = M1.getRingNo();
	if (nrings != M2.getRingNo()) {
		char errorMsg[256];
		sprintf(
				errorMsg,
				"rotationalCorrelation: polar structures have unequal number of rings:\
    	nrings %d and M2.getRingNo %d",
				nrings, M2.getRingNo());
		std::cerr << errorMsg << std::endl;
		REPORT_ERROR(ERR_VALUE_INCORRECT, errorMsg);
	}
	// Fsum should already be set with the right size in the local_transformer
	// (i.e. through a FourierTransform of corr)
	aux.local_transformer.getFourierAlias(aux.Fsum);
	memset((double*) MULTIDIM_ARRAY(aux.Fsum), 0,
			2 * MULTIDIM_SIZE(aux.Fsum) * sizeof(double));

	// Multiply M1 and M2 over all rings and sum
	// Assume M2 is already complex conjugated!
	std::complex<double> auxC;
	for (int iring = 0; iring < nrings; iring++) {
		double w = (2. * PI * M1.ring_radius[iring]);
		int imax = M1.getSampleNo(iring);
		const MultidimArray<std::complex<double> > &M1_iring = M1.rings[iring];
		const MultidimArray<std::complex<double> > &M2_iring = M2.rings[iring];
		double *ptr1 = (double*) MULTIDIM_ARRAY(M1_iring);
		double *ptr2 = (double*) MULTIDIM_ARRAY(M2_iring);
		double *ptrFsum = (double *) MULTIDIM_ARRAY(aux.Fsum);
		for (int i = 0; i < imax; i++) {
			double a = *ptr1++;
			double b = *ptr1++;
			double c = *ptr2++;
			double d = *ptr2++;
			*(ptrFsum++) += w * (a * c - b * d);
			*(ptrFsum++) += w * (b * c + a * d);
		}
	}

	// Inverse FFT to get real-space correlations
	// The local_transformer should already have corr as setReal!!
	aux.local_transformer.inverseFourierTransform();

	angles.resize(XSIZE(aux.local_transformer.getReal()));
	double Kaux = 360. / XSIZE(angles);
	for (size_t i = 0; i < XSIZE(angles); i++)
		DIRECT_A1D_ELEM(angles,i) = (double) i * Kaux;
}

// Compute the normalized Polar Fourier transform --------------------------
void normalizedPolarFourierTransform(const MultidimArray<double> &in,
		Polar<std::complex<double> > &out, bool flag, int first_ring,
		int last_ring, Polar_fftw_plans *&plans, int BsplineOrder) {
	Polar<double> polarIn;
	if (BsplineOrder == 1)
		polarIn.getPolarFromCartesianBSpline(in, first_ring, last_ring, 1);
	else {
		MultidimArray<double> Maux;
		produceSplineCoefficients(3, Maux, in);
		polarIn.getPolarFromCartesianBSpline(Maux, first_ring, last_ring,
				BsplineOrder);
	}
	double mean, stddev;
	polarIn.computeAverageAndStddev(mean, stddev);
	polarIn.normalize(mean, stddev);
	if (plans == NULL) {
		plans = new Polar_fftw_plans();
		polarIn.calculateFftwPlans(*plans);
	}
	fourierTransformRings(polarIn, out, *plans, flag);
}

// Best rotation -----------------------------------------------------------
double best_rotation(const Polar<std::complex<double> > &I1,
		const Polar<std::complex<double> > &I2, RotationalCorrelationAux &aux) {
	MultidimArray<double> angles;
	rotationalCorrelation(I1, I2, angles, aux);

	// Compute the maximum of correlation (inside local_transformer)
	const MultidimArray<double> &corr = aux.local_transformer.getReal();
	int imax = 0;
	double maxval = DIRECT_MULTIDIM_ELEM(corr,0);
	double* ptr = NULL;
	unsigned long int n;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(corr,n,ptr)
		if (*ptr > maxval) {
			maxval = *ptr;
			imax = n;
		}

	// Return the corresponding angle
	return DIRECT_A1D_ELEM(angles,imax);
}

// Align rotationally ------------------------------------------------------
void alignRotationally(MultidimArray<double> &I1, MultidimArray<double> &I2,
		RotationalCorrelationAux &aux, int splineOrder, int wrap) {
	I1.setXmippOrigin();
	I2.setXmippOrigin();

	Polar_fftw_plans *plans = NULL;
	Polar<std::complex<double> > polarFourierI2, polarFourierI1;
	normalizedPolarFourierTransform(I1, polarFourierI1, false, XSIZE(I1) / 5,
			XSIZE(I1) / 2, plans);
	normalizedPolarFourierTransform(I2, polarFourierI2, true, XSIZE(I2) / 5,
			XSIZE(I2) / 2, plans);

	MultidimArray<double> rotationalCorr;
	rotationalCorr.resize(2 * polarFourierI2.getSampleNoOuterRing() - 1);
	aux.local_transformer.setReal(rotationalCorr);
	double bestRot = best_rotation(polarFourierI1, polarFourierI2, aux);

	MultidimArray<double> tmp = I2;
	rotate(splineOrder, I2, tmp, -bestRot, 'Z', wrap);

}

// Cartesian to polar -----------------------------------------------------
void image_convertCartesianToPolar(MultidimArray<double> &in,
		MultidimArray<double> &out, double Rmin, double Rmax, double deltaR,
		double angMin, double angMax, double deltaAng) {
	size_t NAngSteps = (size_t)floor((angMax - angMin) / deltaAng);
	size_t NRSteps = (size_t)floor((Rmax - Rmin) / deltaR);
	out.initZeros(NAngSteps, NRSteps);
	for (size_t i = 0; i < NAngSteps; ++i) {
		double s, c;
		double angle = angMin + i * deltaAng;
		sincos(angle, &s, &c);
		for (size_t j = 0; j < NRSteps; ++j) {
			double R = Rmin + j * deltaR;
			A2D_ELEM(out,i,j) = in.interpolatedElement2D(R * c, R * s);
		}
	}
}

// Cartesian to polar -----------------------------------------------------
void image_convertCartesianToPolar_ZoomAtCenter(const MultidimArray<double> &in,
		MultidimArray<double> &out, Matrix1D<double> &R, double zoomFactor,
		double Rmin, double Rmax, int NRSteps, double angMin, double angMax,
		int NAngSteps) {
	/* Octave
	 * r=0:64;plot(r,d.^(r/d),r,r,d*((r/d).^2.8));
	 */
	double deltaAng = (angMax - angMin + 1) / NAngSteps;
	double deltaR = (Rmax - Rmin + 1) / NRSteps;
	out.initZeros(NAngSteps, NRSteps);
	if (VEC_XSIZE(R) == 0) {
		R.initZeros(NRSteps);
		double Rrange = Rmax - Rmin;
		for (int j = 0; j < NRSteps; ++j)
			VEC_ELEM(R,j) = Rmin
					+ Rrange * pow(j * deltaR / Rrange, zoomFactor);
	}
	for (int i = 0; i < NAngSteps; ++i) {
		double s, c;
		double angle = angMin + i * deltaAng;
		sincos(angle, &s, &c);
		for (int j = 0; j < NRSteps; ++j) {
			double Rj = VEC_ELEM(R,j);
			A2D_ELEM(out,i,j) = in.interpolatedElement2D(Rj * c, Rj * s);
		}
	}
}

// Cartesian to cylindrical -----------------------------------------------
void volume_convertCartesianToCylindrical(const MultidimArray<double> &in,
		MultidimArray<double> &out, double Rmin, double Rmax, double deltaR,
		double angMin, double angMax, double deltaAng, Matrix1D<double> &axis) {
	size_t NAngSteps = (size_t)floor((angMax - angMin) / deltaAng);
	size_t NRSteps = (size_t)floor((Rmax - Rmin) / deltaR);
	out.initZeros(ZSIZE(in), NAngSteps, NRSteps);
	STARTINGZ(out) = STARTINGZ(in);
	STARTINGY(out) = STARTINGX(out) = 0;

	// Calculate rotation matrix
	Matrix2D<double> M;
	if (XX(axis)!=0.0 || YY(axis)!=0.0)
	{
		// 1) Calculate rotation axis v=(0,0,1) x axis
		Matrix1D<double> v(3);
		XX(v)=-YY(axis);
		YY(v)=XX(axis);

		// 2) Calculate rotation angle (0,0,1).axis=module(axis).cos(ang)
		double rotAng=acos(ZZ(axis)/axis.module());

		// 3) Calculate matrix
		rotation3DMatrix(RAD2DEG(rotAng), v, M, false);
	}

	Matrix1D<double> p(3), pc(3);
	SPEED_UP_temps012;
	for (size_t i = 0; i < NAngSteps; ++i) {
		double s, c;
		double angle = angMin + i * deltaAng;
		sincos(angle, &s, &c);
		for (int k = STARTINGZ(in); k <= FINISHINGZ(in); ++k) {
			ZZ(p)=k;
			for (size_t j = 0; j < NRSteps; ++j) {
				double R = Rmin + j * deltaR;
				YY(p)=R*s;
				XX(p)=R*c;
				if (MAT_XSIZE(M)>0)
				{
					M3x3_BY_V3x1(pc,M,p);
					A3D_ELEM(out,k,i,j) = in.interpolatedElement3D(XX(pc),YY(pc),ZZ(pc));
				}
				else
					A3D_ELEM(out,k,i,j) = in.interpolatedElement3D(XX(p),YY(p),ZZ(p));
			}
		}
	}
}

// Cartesian to cylindrical -----------------------------------------------
void volume_convertCartesianToSpherical(const MultidimArray<double> &in,
		MultidimArray<double> &out, double Rmin, double Rmax, double deltaR,
		double deltaRot, double deltaTilt) {
	size_t NRotSteps = (size_t)floor(2*PI / deltaRot);
	size_t NTiltSteps = (size_t)floor(PI / deltaRot);
	size_t NRSteps = (size_t)floor((Rmax - Rmin) / deltaR);
	out.initZeros(NRSteps, NRotSteps, NTiltSteps);
	STARTINGZ(out) = STARTINGY(out) = STARTINGX(out) = 0;

	Matrix1D<double> p(3);
	for (size_t i = 0; i < NRotSteps; ++i) {
		double srot, crot;
		double rot = i * deltaRot;
		sincos(rot, &srot, &crot);
		for (size_t j = 0; j < NRotSteps; ++j) {
			double stilt, ctilt;
			double tilt = j * deltaTilt;
			sincos(tilt, &stilt, &ctilt);
			double sc=stilt*crot;
			double ss=stilt*srot;
			for (size_t k = 0; k<NRSteps; ++k) {
				double R = Rmin + k * deltaR;
				ZZ(p)=R*ctilt;
				YY(p)=R*ss;
				XX(p)=R*sc;
				A3D_ELEM(out,k,i,j) = in.interpolatedElement3D(XX(p),YY(p),ZZ(p));
			}
		}
	}
}
