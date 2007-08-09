/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.uam.es)
 *
 *  This code is strongly based on ideas by Pawel Penczek & Zhengfan
 *  Yang as implemented in SPARX at the University of Texas - Houston 
 *  Medical School
 *
 *  see P. A. Penczek, R. Renka, and H. Schomberg,
 *      J. Opt. Soc. Am. _21_, 449 (2004)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#ifndef POLAR_H
#define POLAR_H
#include "funcs.h"
#include "matrix2d.h"
#include "gridding.h"

#define FULL_CIRCLES 0
#define HALF_CIRCLES 1
#define DONT_AVERAGE false
#define AVERAGE true
#define DONT_CONJUGATE false
#define CONJUGATE true

template<typename T>
class Polar
{
protected:
    FileName              fn_pol;       // name of the polar
public:
    int                   mode;         // Use full or half circles
    vector<double>        ring_radius;  // radius of each ring
    vector<Matrix1D<T> >  rings;        // vector with all rings
public:
    /// @defgroup PolarConstructors Polar constructors
    /// @ingroup Polars

    /** Empty constructor
     * @ingroup PolarConstructors
     *
     * An empty image with size 0x0 is created.
     *
     * @code
     * Polar P;
     * @endcode
     */
    Polar()
    {
        fn_pol = "";
	ring_radius.clear();
	rings.clear();
	mode = FULL_CIRCLES;
    }
    
    /** Copy constructor
     * @ingroup PolarConstructors
     *
     * @code
     * Polar P2(P1);
     * @endcode
     */
    Polar(const Polar& P)
    {
        rings = P.rings;
        fn_pol = P.fn_pol;
	ring_radius = P.ring_radius;
	mode = P.mode;
    }

    /// @defgroup PolarOperations Some operations
    /// @ingroup Polars
    /** Assignment
     * @ingroup PolarOperations
     */
    Polar& operator=(const Polar& P)
    {
        if (this != &P)
        {
            fn_pol = P.fn_pol;
            rings = P.rings;
	    ring_radius = P.ring_radius;
	    mode = P.mode;
        }
        return *this;
    }
  
    /** Subtract a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator-(const T val) const
    {
	Polar<T> result;
	result = *(this);
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		result(i,j) -= val;

        return result;
    }

    /** Add a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator+(const T val)
    {
	Polar<T> result;
	result = *(this);
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		result(i,j) += val;

        return result;
    }

    /** Multiply by a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator*(const T val)
    {
	Polar<T> result;
	result = *(this);
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		result(i,j) *= val;

        return result;
    }

    /** Divide by a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator/(const T val)
    {
	Polar<T> result;
	result = *(this);
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		result(i,j) /= val;

        return result;
    }

    /** Subtract a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    void operator-=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) -= val;
    }

    /** Add a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    void operator+=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) += val;
    }

    /** Multiply by a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    void operator*=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) *= val;
    }

    /** Divide by a constant pixel-by-pixel
     * @ingroup PolarOperations
     */
    void operator/=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) /= val;
    }

    /** Subtract two polars pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator-=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) -= in(i,j);
    }

     /** Add two polars pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator+=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) += in(i,j);
    }

     /** Multiply two polars pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator*=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) *= in(i,j);
    }

     /** Divide two polars pixel-by-pixel
     * @ingroup PolarOperations
     */
    Polar& operator/=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) /= in(i,j);

    }
    /** Rename polar
     * @ingroup PolarOperations
     *
     * Give a new name to the polar.
     *
     * @code
     * P.rename("new_name");
     * @endcode
     */

    void rename(FileName newName)
    {
        fn_pol = newName;
    }

    /** Empty polar
     * @ingroup PolarOperations
     *
     * This function clears the polar to an empty vector without name.
     *
     * @code
     * P.clear();
     * @endcode
     */
    void clear()
    {
        fn_pol = "";
        rings.clear();
	ring_radius.clear();
	mode = FULL_CIRCLES;
    }

    /// @defgroup PolarAccess Polar access
    /// @ingroup Polars

    /** Name access
     * @ingroup PolarAccess
     *
     * This function is used to know the name of the polar. It cannot be used to
     * assign a new one. You may use rename() for that.
     *
     * @code
     * cout << "Polar name: " << P.name() << endl;
     * @endcode
     */
    const FileName name() const
    {
        return fn_pol;
    }

    /** Number of rings access
     * @ingroup PolarAccess
     *
     * This function is used to know the number of rings in the polar. 
     *
     * @code
     * cout << "Number of rings: " << P.getRingNo() << endl;
     * @endcode
     */
    const int getRingNo() const
    {
        return rings.size();
    }

    /** Mode access
     * @ingroup PolarAccess
     *
     * This function is used to know the "mode" of the polar. 
     *
     * There are two modes: 
     * FULL_CIRCLES = 0 (used for asymmetric functions)
     * HALF_CIRCLES = 0 (used for symmetric functions, e.g. Fourier transforms)
     *
     * @code
     * cout << "Mode: " << P.getMode() << endl;
     * @endcode
     */
    const int getMode() const
    {
        return mode;
    }

    /** Number of samples in each ring access
     * @ingroup PolarAccess
     *
     * This function is used to know the number of samples in a given ring.
     *
     * @code
     * cout << "Number of samples in second ring: " << P.getSampleNo(1) << endl;
     * @endcode
     */
    const int getSampleNo(int iring) const
    {
        return XSIZE(rings[iring]);
    }

    /** The radius of each ring access
     * @ingroup PolarAccess
     *
     * This function is used to know the radius of a given ring.
     *
     * @code
     * cout << "Radius of second ring: " << P.getRadius(1) << endl;
     * @endcode
     */
    const double getRadius(int iring) const
    {
        return ring_radius[iring];
    }

    /** 1D Matrix access
     * @ingroup PolarAccess
     *
     * This operator can be used to access any ring of the polar as Matrix1D.
     *
     * @code
     * Matrix1D<double> secondring = P.getRing(1);
     * @endcode
     */
    Matrix1D< T >& getRing(int i)
    {
        return rings[i];
    }
    const  Matrix1D< T >& getRing(int i) const
    {
        return rings[i];
    }

    /** 1D Matrix access
     * @ingroup PolarAccess
     *
     * This operator can be used to set any ring of the polar with a Matrix1D.
     *
     * @code
     * Matrix1D<double> ring;
     * P.setRing(1,ring);
     * @endcode
     */
    void setRing(int i, Matrix1D< T > val) const
    {
        rings[i] = val;
    }

    /** Pixel access
     * @ingroup PolarAccess
     *
     * This operator is used to access a pixel within the polar. That
     * means, given the ring and the number of the pixel in the
     * rotational direction.
     *
     * @code
     * int ring = 1, phi = 0;
     * cout <<" first pixel of second ring= "<<P.getPixel(ring,phi)<<endl;
     * @endcode
     */
    T& getPixel(int r, int f) const
    {
        return rings[r](f);
    }

    /** Pixel access
     * @ingroup PolarAccess
     *
     * This operator is used to access a pixel within the polar. That
     * means, given the ring and the number of the pixel in the
     * rotational direction.
     *
     * @code
     * int ring = 1, phi = 0;
     * cout <<" first pixel of second ring= "<<P(ring,phi)<<endl;
     * @endcode
     */
    T& operator()(int r, int f) const
    {
	return rings[r](f);
    }

    /** Pixel access
     * @ingroup PolarAccess
     *
     * This operator is used to set a pixel within the polar. That
     * means, given the ring and the number of the pixel in the
     * rotational direction.
     *
     * @code
     * int ring = 1, phi = 0;
     * double val = 1.2;
     * P.setPixel(ring,phi,val);
     * @endcode
     */
    void setPixel(int r, int f, T val) const
    {
        rings[r](f) = val;
    }


    /** Compute sum or average of all pixels in polar rings.
     * @ingroup PolarFunctions
     *
     */
    T computeSum(bool average = DONT_AVERAGE, int mode = FULL_CIRCLES) const
    {
	T aux, sum = 0.;
	double twopi, w, N = 0;

	if (mode == FULL_CIRCLES)
	    twopi = 2.*PI;
	else if (mode == HALF_CIRCLES)
	    twopi = PI;
	else
	    REPORT_ERROR(1,"Incorrect mode for computeSum");

	for (int i = 0; i < rings.size(); i++)
	{
	    // take varying sampling into account
	    w = (twopi * ring_radius[i]) / (double) XSIZE(rings[i]);
	    for (int j = 0; j < XSIZE(rings[i]); j++)
	    {
		aux = rings[i](j);
		sum += w * aux;
		N += w;
	    }
	}
	if (N != 0. && average)
	    sum = sum / N;

	return sum;
    }
 
    /** Compute squared-sum or average squared-sum of all pixels in polar rings.
     * @ingroup PolarFunctions
     *
     */
    T computeSum2(bool average = DONT_AVERAGE, int mode = FULL_CIRCLES) const
    {
	T aux, sum2 = 0.;
	double twopi, w, N = 0;

	if (mode == FULL_CIRCLES)
	    twopi = 2.*PI;
	else if (mode == HALF_CIRCLES)
	    twopi = PI;
	else
	    REPORT_ERROR(1,"Incorrect mode for computeSum2");

	for (int i = 0; i < rings.size(); i++)
	{
	    // take varying sampling into account
	    w = (twopi * ring_radius[i]) / (double) XSIZE(rings[i]);
	    for (int j = 0; j < XSIZE(rings[i]); j++)
	    {
		aux = rings[i](j) * rings[i](j);
		sum2 += w * aux;
		N += w;
	    }
	}
	if (N != 0. && average)
	    sum2 = sum2 / N;

	// Return half the sum2 because we oversample twice...
	return sum2;
    }
 
    /** Convert cartesian Matrix2D to Polar
     * @ingroup PolarFunctions
     *
     * Use gridding for interpolation. The input Matrix2D is assumed
     * to be pre-processed for gridding purposes.
     *
     * @code
     * Polar P;
     * KaiserBessel kb;
     * Matrix2D<double> Maux;
     * produceGriddingMatrix2D(img(),Maux,kb);
     * P.getPolarFromCartesian(Maux,kb,1,15);
     * @endcode
     */
    void getPolarFromCartesian(const Matrix2D<T> &M1, KaiserBessel &kb, 
			       int first_ring, int last_ring, 
			       int ring_step = 1, double psi_step = -1., int mode1 = FULL_CIRCLES, 
			       double xoff = 0., double yoff = 0.,
			       bool const_sam = false)
    {
	int nsam, half_nsam_last;
	int nring = last_ring - first_ring + 1;
	double radius, twopi, dphi, phi; 
	double xp, yp;

	Matrix1D<T> Mring;
	rings.clear();
	ring_radius.clear();
	mode = mode1;

	// Check that no pixel falls outside the image
	if ((int)ABS(xoff) + last_ring > XSIZE(M1)/4 || 
	    (int)ABS(yoff) + last_ring > YSIZE(M1)/4)
	    REPORT_ERROR(1,"polarCoordinatesGridding: last ring falls outside image");

	if (mode == FULL_CIRCLES)
	    twopi = 2.*PI;
	else if (mode == HALF_CIRCLES)
	    twopi = PI;
	else
	    REPORT_ERROR(1,"Incorrect mode for getPolarFromCartesian");
	
	// Determine angular sampling of the last ring in degrees
	// Use half_sam_last to always have even number of sample
	// points (which is convenient for Half2Whole of FTs)
	
	if (psi_step < 0.)
	    // Oversample twice 
	    half_nsam_last = (int)( twopi * (double)last_ring );
	else 
	    // User-defined sampling
	    half_nsam_last = (int)( 0.5 * twopi / DEG2RAD(psi_step) );

	for (int iring = first_ring; iring <= last_ring; iring+= ring_step)
	{
	    radius = (double) iring;

	    if (const_sam)
		nsam = 2 * half_nsam_last;
	    else
		nsam = 2 * ( MAX(1,(int)(half_nsam_last * (double)iring / (double)last_ring)) );

	    dphi = twopi / (double)nsam;
	    Mring.resize(nsam);
	    for (int iphi = 0; iphi < nsam; iphi++)
	    {
		// polar coordinates
		phi = iphi * dphi;
		xp = sin(phi) * radius;
		yp = cos(phi) * radius;
		// Origin offsets
		xp += xoff;
		yp += yoff; 
		Mring(iphi) = (T) interpolatedElementGridding(M1,xp,yp,kb);
	    }
	    rings.push_back(Mring);
	    ring_radius.push_back(radius);
	}

    }

    /** Fourier transform all rings
     * @ingroup PolarFunctions
     *
     * 1D Fourier transform of all rings
     * Only the assymetric half is stored
     * 
     * Note that both complex and real signals can be Fourier transformed!
     *
     * @code
     * Polar<double> P;
     * Polar<complex<double> > Pf;
     * KaiserBessel kb;
     * Matrix2D<double> Maux;
     * produceGriddingMatrix2D(img(),Maux,kb);
     * P.getPolarFromCartesian(Maux,kb,1,15);
     * Pf = P.fourierTransformRings();
     * @endcode
     */
    Polar<complex<double> > fourierTransformRings(bool conjugated = DONT_CONJUGATE) const
    {
	Polar<complex<double> > out;
	Matrix1D<complex<double> > Fring;
	out.clear();
	for (int iring = 0; iring < rings.size(); iring++)
	{ 
	    FourierTransformHalf(rings[iring],Fring);
	    if (conjugated)
	    {
		for (int i = 0; i < XSIZE(Fring); i++)
		    Fring(i) = conj(Fring(i));
	    }
	    out.rings.push_back(Fring);
	}

	out.mode = mode;
	out.ring_radius = ring_radius;
	return out;
    }

};

/** @defgroup PolarRelated Polar Related functions
 * @ingroup Polar
 *
 * These functions are not methods of Polar
 */

/** Fourier-space rotational Cross-Correlation Funtion
 * @ingroup PolarRelated
 *
 *  This function returns the rotational cross-correlation
 *  function of two Polars M1 and M2 using the
 *  cross-correlation convolution theorem. 
 *
 *  M2 is assumed to be the complex conjugated.
 *
 * Note that this function can be used for real-space and
 * fourier-space correlations!!
 *
 */
template<typename T>
void rotationalCorrelation(const Polar<complex<double> > &M1,
			   const Polar<complex<double> > &M2,
			   Matrix1D<double> &angles, 
			   Matrix1D<T > &corr)
{

    Matrix1D<complex<double> > Fsum, Faux;
    complex<double> aux;

    int nrings = M1.getRingNo();
    if (nrings != M2.getRingNo())
	REPORT_ERROR(1,"getBestAnglePolarGridding: polar structures have unequal number of rings!");

    // Resize Fsum to the size of the last ring
    int nsam_last = M1.getSampleNo(nrings - 1);
    Fsum.resize(nsam_last);

    // Multiply M1 and M2 over all rings and sum
    // Assume M2 is already complex conjugated!
    for (int iring = 0; iring < nrings; iring++)
    { 
	int nsam = M1.getSampleNo(iring);
	// Penczeks weights is more-or-less 1/2piR, mine is reverse
	// because it depends whether the FFT results X or X/nsam
	double w = (2.* PI * M1.ring_radius[iring]);
	for (int i = 0; i < nsam; i++)
	{
	    aux = M1(iring,i) * M2(iring,i);
	    Fsum(i) += w * aux; 
	}
    }

    // Inverse FFT to get real-space correlations
    // As only half the FTs are stored, the original number of
    // sampling points was 2 * (nsam - 1). Note that the latter is
    // only true for even-valued whole-sizes!
    nsam_last = 2 * (nsam_last - 1);
    InverseFourierTransformHalf(Fsum,corr,nsam_last);
    corr.setStartingX(0);
    angles.resize(nsam_last);
    for (int i = 0; i < nsam_last; i++)
	angles(i)=(double)i*360./(nsam_last);

}

void inverseFourierTransformRings(const Polar<complex<double> > & in, 
				  Polar<double> & out, bool conjugated = false);
#endif
