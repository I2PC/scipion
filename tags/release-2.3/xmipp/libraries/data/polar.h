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

/// @defgroup Polar Polar coordinates 
/// @ingroup DataLibrary
//@{

/** Class for polar coodinates */
template<typename T>
class Polar
{
protected:
    FileName              fn_pol;       // name of the polar
public:
    int                   mode;         // Use full or half circles
    double                oversample;
    std::vector<double>        ring_radius;  // radius of each ring
    std::vector<Matrix1D<T> >  rings;        // vector with all rings
public:
    /** Empty constructor
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
	oversample = 1.;
    }
    
    /** Copy constructor
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
	oversample = P.oversample;
    }

    /** Assignment
     */
    Polar& operator=(const Polar& P)
    {
        if (this != &P)
        {
            fn_pol = P.fn_pol;
            rings = P.rings;
	    ring_radius = P.ring_radius;
	    mode = P.mode;
	    oversample = P.oversample;
        }
        return *this;
    }
  
    /** Subtract a constant pixel-by-pixel
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
     */
    void operator-=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) -= val;
    }

    /** Add a constant pixel-by-pixel
     */
    void operator+=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) += val;
    }

    /** Multiply by a constant pixel-by-pixel
     */
    void operator*=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) *= val;
    }

    /** Divide by a constant pixel-by-pixel
     */
    void operator/=(const T val)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) /= val;
    }

    /** Subtract two polars pixel-by-pixel
     */
    Polar& operator-=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) -= in(i,j);
    }

     /** Add two polars pixel-by-pixel
     */
    Polar& operator+=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) += in(i,j);
    }

     /** Multiply two polars pixel-by-pixel
     */
    Polar& operator*=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) *= in(i,j);
    }

     /** Divide two polars pixel-by-pixel
     */
    Polar& operator/=(const Polar<T> in)
    {
	for (int i = 0; i < rings.size(); i++)
	    for (int j = 0; j < XSIZE(rings[i]); j++)
		(*(this))(i,j) /= in(i,j);

    }

    /** Rename polar
     *
     * Give a new name to the polar.
     *
     * @code
     * P.rename("new_name");
     * @endcode
     */

    void rename(const FileName &newName)
    {
        fn_pol = newName;
    }

    /** Empty polar
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
	oversample = 1.;
    }

    /** Name access
     *
     * This function is used to know the name of the polar. It cannot be used to
     * assign a new one. You may use rename() for that.
     *
     * @code
     * std::cout << "Polar name: " << P.name() << std::endl;
     * @endcode
     */
    const FileName name() const
    {
        return fn_pol;
    }

    /** Number of rings access
     *
     * This function is used to know the number of rings in the polar. 
     *
     * @code
     * std::cout << "Number of rings: " << P.getRingNo() << std::endl;
     * @endcode
     */
    const int getRingNo() const
    {
        return rings.size();
    }

    /** Mode access
     *
     * This function is used to know the "mode" of the polar. 
     *
     * There are two modes: 
     * FULL_CIRCLES = 0 (used for asymmetric functions)
     * HALF_CIRCLES = 0 (used for symmetric functions, e.g. Fourier transforms)
     *
     * @code
     * std::cout << "Mode: " << P.getMode() << std::endl;
     * @endcode
     */
    const int getMode() const
    {
        return mode;
    }

    /** Oversample access
     *
     * This function is used to know the oversampling factor of the polar. 
     *
     * Oversample = 1 means there is no oversampling
     * Oversample > 1 means oversampling
     * Oversampling < 1 means undersampling
     *
     * @code
     * std::cout << "Oversample: " << P.getOversample() << std::endl;
     * @endcode
     */
    const double getOversample() const
    {
        return oversample;
    }

    /** Number of samples in each ring access
     *
     * This function is used to know the number of samples in a given ring.
     *
     * @code
     * std::cout << "Number of samples in second ring: " << P.getSampleNo(1) << std::endl;
     * @endcode
     */
    const int getSampleNo(int iring) const
    {
        return XSIZE(rings[iring]);
    }

    /** The radius of each ring access
     *
     * This function is used to know the radius of a given ring.
     *
     * @code
     * std::cout << "Radius of second ring: " << P.getRadius(1) << std::endl;
     * @endcode
     */
    const double getRadius(int iring) const
    {
        return ring_radius[iring];
    }

    /** 1D Matrix access
     *
     * This operator can be used to access any ring of the polar as Matrix1D.
     *
     * @code
     * Matrix1D<double> secondring = P.getRing(1);
     * @endcode
     */
    //@{
    Matrix1D< T >& getRing(int i)
    {
        return rings[i];
    }
    const  Matrix1D< T >& getRing(int i) const
    {
        return rings[i];
    }
    //@}

    /** 1D Matrix access
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
     *
     * This operator is used to access a pixel within the polar. That
     * means, given the ring and the number of the pixel in the
     * rotational direction.
     *
     * @code
     * int ring = 1, phi = 0;
     * std::cout <<" first pixel of second ring= "<<P.getPixel(ring,phi)<<std::endl;
     * @endcode
     */
    T& getPixel(int r, int f) const
    {
        return rings[r](f);
    }

    /** Pixel access
     *
     * This operator is used to access a pixel within the polar. That
     * means, given the ring and the number of the pixel in the
     * rotational direction.
     *
     * @code
     * int ring = 1, phi = 0;
     * std::cout <<" first pixel of second ring= "<<P(ring,phi)<<std::endl;
     * @endcode
     */
    T& operator()(int r, int f) const
    {
	return rings[r](f);
    }

    /** Pixel access
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

	return sum2;
    }
 
   /** Get Cartesian Coordinates of the Polar sampling
     *
     * The output of this function can be used to calculate Voronoi
     * areas, lists of neighbours etc.
     *
     * To deal with the borders of the polar structure (a maximum of)
     * "extra_shell" extra rings are calculated on the inside and outside
     * of the polar structure. 
     *
     */
    void getCartesianCoordinates(std::vector<double> &x,
				 std::vector<double> &y,
				 std::vector<T> &data,
				 const double extra_shell = GRIDDING_K/2)
    {
	double                     twopi, dphi,radius;
	int                        nsam;

	// Only for full circles for now!
	if (mode != FULL_CIRCLES)
	    REPORT_ERROR(1,"VoronoiArea only implemented for FULL_CIRCLES mode of Polar");
	else
	    twopi = 2.*PI;

	// First fill the vector with the originally sampled coordinates
	x.clear();
	y.clear();
	data.clear();
	for (int i = 0; i < rings.size(); i++)
	{
	    nsam = XSIZE(rings[i]);
	    dphi = twopi/(double)nsam;
	    radius = ring_radius[i];
	    for (int j = 0; j < nsam; j++)
	    {
		x.push_back(radius*sin(j*dphi));
		y.push_back(radius*cos(j*dphi));
		data.push_back(rings[i](j));
	    }
	}

	// Add additional points on the inside and outside of the rings
	// Add a maximum of "extra_shell" rings
	// Set data to zero here
	double first_ring  = ring_radius[0];
	double last_ring   = ring_radius[rings.size()-1];
	double outer       = last_ring + extra_shell;
	double inner       = XMIPP_MAX(0.,first_ring - extra_shell);
	for (radius = 0.; radius < outer; radius +=1.)
	{
	    if ( (radius >= inner && radius < first_ring) ||
		 ( radius <= outer && radius > last_ring) )
	    {
		nsam = 2 * (int)( 0.5 * oversample * twopi * radius );
		nsam = XMIPP_MAX(1, nsam);
		dphi = twopi / (double)nsam;
		for (int j = 0; j < nsam; j++)
		{
		    x.push_back(radius*sin(j*dphi));
		    y.push_back(radius*cos(j*dphi));
		    data.push_back(0.);
		}
	    }
	}

    }

    /** Convert cartesian Matrix2D to Polar using gridding interpolation
     *
     * The input Matrix2D is assumed to be pre-processed for gridding
     *
     * @code
     * Polar P;
     * KaiserBessel kb;
     * Matrix2D<double> Maux;
     * produceReverseGriddingMatrix2D(img(),Maux,kb);
     * P.getPolarFromCartesian(Maux,kb,1,15);
     * @endcode
     *
     */
    void getPolarFromCartesianGridding(const Matrix2D<T> &M1, KaiserBessel &kb, 
				       int first_ring, int last_ring, 
				       double xoff = 0., double yoff = 0.,
				       double oversample1 = 1., int mode1 = FULL_CIRCLES)
    {
	int nsam;
	int nring = last_ring - first_ring + 1;
	double radius, twopi, dphi, phi; 
	double xp, yp, minxp, maxxp, minyp, maxyp;

	Matrix1D<T> Mring;
	rings.clear();
	ring_radius.clear();
	mode = mode1;
	oversample = oversample1;

	if (mode == FULL_CIRCLES)
	    twopi = 2.*PI;
	else if (mode == HALF_CIRCLES)
	    twopi = PI;
	else
	    REPORT_ERROR(1,"Incorrect mode for getPolarFromCartesian");
	

	// Take 2x oversize M1 dims into account for calculating the limits 
	minxp = FIRST_XMIPP_INDEX(XSIZE(M1)/2);
	minyp = FIRST_XMIPP_INDEX(YSIZE(M1)/2);
	maxxp = LAST_XMIPP_INDEX(XSIZE(M1)/2);
	maxyp = LAST_XMIPP_INDEX(YSIZE(M1)/2);

	// Loop over all polar coordinates
	for (int iring = first_ring; iring <= last_ring; iring++)
	{
	    radius = (double) iring;
	    // Non-constant sampling!! (always even for convenient Half2Whole of FTs)
	    nsam = 2 * (int)( 0.5 * oversample * twopi * radius );
	    nsam = XMIPP_MAX(1, nsam);
	    dphi = twopi / (double)nsam;
	    Mring.resize(nsam);
	    for (int iphi = 0; iphi < nsam; iphi++)
	    {
		// from polar to original cartesian coordinates
		phi = iphi * dphi;
		xp = sin(phi) * radius;
		yp = cos(phi) * radius;

		// Origin offsets
		xp += xoff;
		yp += yoff; 

		// Wrap coordinates
                if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                    xp > maxxp + XMIPP_EQUAL_ACCURACY)
                    xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);
                if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                    yp > maxyp + XMIPP_EQUAL_ACCURACY)
                    yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);

		// Perform the convolution interpolation
		Mring(iphi) = (T) interpolatedElementReverseGridding(M1,xp,yp,kb);
	    }
	    rings.push_back(Mring);
	    ring_radius.push_back(radius);
	}

    }

    /** Convert cartesian Matrix2D to Polar using B-spline interpolation
     *
     * The input Matrix2D is assumed to be pre-processed for B-splines
     *
     * @code
     * Polar P;
     * Matrix2D<double> Maux;
     * img().produceSplineCoefficients(Maux,3);
     * P.getPolarFromCartesianBSpline(Maux,1,15);
     * @endcode
     *
     */
    void getPolarFromCartesianBSpline(const Matrix2D<T> &M1, 
				     int first_ring, int last_ring, 
				     double xoff = 0., double yoff = 0.,
				     double oversample1 = 1., int mode1 = FULL_CIRCLES)
    {
	int nsam;
	int nring = last_ring - first_ring + 1;
	double radius, twopi, dphi, phi; 
	double xp, yp, minxp, maxxp, minyp, maxyp;

	Matrix1D<T> Mring;
	rings.clear();
	ring_radius.clear();
	mode = mode1;
	oversample = oversample1;

	if (mode == FULL_CIRCLES)
	    twopi = 2.*PI;
	else if (mode == HALF_CIRCLES)
	    twopi = PI;
	else
	    REPORT_ERROR(1,"Incorrect mode for getPolarFromCartesian");
	

	// Limits of the matrix (not oversized!)
	minxp = FIRST_XMIPP_INDEX(XSIZE(M1));
	minyp = FIRST_XMIPP_INDEX(YSIZE(M1));
	maxxp = LAST_XMIPP_INDEX(XSIZE(M1));
	maxyp = LAST_XMIPP_INDEX(YSIZE(M1));

	// Loop over all polar coordinates
	for (int iring = first_ring; iring <= last_ring; iring++)
	{
	    radius = (double) iring;
	    // Non-constant sampling!! (always even for convenient Half2Whole of FTs)
	    nsam = 2 * (int)( 0.5 * oversample * twopi * radius );
	    nsam = XMIPP_MAX(1, nsam);
	    dphi = twopi / (double)nsam;
	    Mring.resize(nsam);
	    for (int iphi = 0; iphi < nsam; iphi++)
	    {
		// from polar to original cartesian coordinates
		phi = iphi * dphi;
		xp = sin(phi) * radius;
		yp = cos(phi) * radius;

		// Origin offsets
		xp += xoff;
		yp += yoff; 

		// Wrap coordinates
                if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                    xp > maxxp + XMIPP_EQUAL_ACCURACY)
                    xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);
                if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                    yp > maxyp + XMIPP_EQUAL_ACCURACY)
                    yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);

		// Perform the convolution interpolation
		Mring(iphi) = (T) M1.interpolatedElementBSpline(xp,yp,3);
	    }
	    rings.push_back(Mring);
	    ring_radius.push_back(radius);
	}

    }



    /** Fourier transform all rings
     *
     * 1D Fourier transform of all rings
     * Only the assymetric half is stored
     * 
     * Note that both complex and real signals can be Fourier transformed!
     *
     * @code
     * Polar<double> P;
     * Polar<std::complex<double> > Pf;
     * KaiserBessel kb;
     * Matrix2D<double> Maux;
     * produceReverseGriddingMatrix2D(img(),Maux,kb);
     * P.getPolarFromCartesian(Maux,kb,1,15);
     * Pf = P.fourierTransformRings();
     * @endcode
     */
    Polar<std::complex<double> > fourierTransformRings(bool conjugated = DONT_CONJUGATE) const
    {
	Polar<std::complex<double> > out;
	Matrix1D<std::complex<double> > Fring;
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

/** Fourier-space rotational Cross-Correlation Funtion
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
void rotationalCorrelation(const Polar<std::complex<double> > &M1,
			   const Polar<std::complex<double> > &M2,
			   Matrix1D<double> &angles, 
			   Matrix1D<T > &corr)
{

    Matrix1D<std::complex<double> > Fsum, Faux;
    std::complex<double> aux;

    int nrings = M1.getRingNo();
    if (nrings != M2.getRingNo())
	REPORT_ERROR(1,"rotationalCorrelation: polar structures have unequal number of rings!");

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
    STARTINGX(corr)=0;
    angles.resize(nsam_last);
    for (int i = 0; i < nsam_last; i++)
	angles(i)=(double)i*360./(nsam_last);
}

/** Inverse Fourier Transform of all rings
 *
 */
void inverseFourierTransformRings(const Polar<std::complex<double> > & in, 
				  Polar<double> & out, bool conjugated = false);

/** Convert to a single vector
 *
 * Convert complex Polar structure to a single vector
 * This may be useful for parallelization purposes.
 * 
 */
void convertPolarToSingleArray(const Polar<std::complex<double> > & in, 
			       Matrix1D<double> & out);

/** Convert to a single vector
 *
 * Convert real Polar structure to a single vector
 * This may be useful for parallelization purposes.
 * 
 */
void convertPolarToSingleArray(const Polar<double> & in, 
			       Matrix1D<double> & out);

/** Convert back from a single vector to a complex polar structure
 *
 * Convert back from a single vector to a complex polar structure.
 * This may be useful for parallelization purposes.
 * The structure of the polar should be the correct one already
 * 
 */
void convertSingleArrayToPolar(const Matrix1D<double> & in,
			       Polar<std::complex<double> > & out);

/** Convert back from a single vector to a real polar structure
 *
 * Convert back from a single vector to a real polar structure.
 * This may be useful for parallelization purposes.
 * The structure of the polar should be the correct one already
 * 
 */
void convertSingleArrayToPolar(const Matrix1D<double> & in,
			       Polar<double> & out);

/** Compute a normalized polar Fourier transform of the input image.*/
void normalizedPolarFourierTransform(const Matrix2D<double> &in,
    Polar< std::complex<double> > &out, bool flag,
    int first_ring, int last_ring);

/** Best rotation between two normalized polar Fourier transforms. */
double best_rotation(const Polar< std::complex<double> > &I1,
    const Polar< std::complex<double> > &I2);
//@}
#endif
