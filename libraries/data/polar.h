/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef POLAR_H
#define POLAR_H
#include "xmipp_funcs.h"
#include "multidim_array.h"
#include "transformations.h"
#include "xmipp_fftw.h"

#define FULL_CIRCLES 0
#define HALF_CIRCLES 1
#define DONT_CONJUGATE false
#define CONJUGATE true
#define DONT_KEEP_TRANSFORM false
#define KEEP_TRANSFORM true

/// @defgroup Polar Polar coordinates
/// @ingroup DataLibrary
//@{

/** Structure for fftw plans */
typedef struct Polar_Fftw_Plans
{
    std::vector<FourierTransformer>          transformers;
    std::vector<MultidimArray<double> >  arrays;
}
Polar_fftw_plans;

/** Class for polar coodinates */
template<typename T>
class Polar
{
protected:
    FileName                   fn_pol;       // name of the polar
public:
    int                        mode;         // Use full or half circles
    double                     oversample;
    std::vector<double>        ring_radius;  // radius of each ring
    std::vector<MultidimArray<T> >  rings;        // vector with all rings
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

    /** Destructor.
     */
    ~Polar()
    {
        rings.clear();
        ring_radius.clear();
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
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                result(i,j) -= val;

        return result;
    }

    /** Add a constant pixel-by-pixel
     */
    Polar& operator+(const T val)
    {
        Polar<T> result;
        result = *(this);
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                result(i,j) += val;

        return result;
    }

    /** Multiply by a constant pixel-by-pixel
     */
    Polar& operator*(const T val)
    {
        Polar<T> result;
        result = *(this);
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                result(i,j) *= val;

        return result;
    }

    /** Divide by a constant pixel-by-pixel
     */
    Polar& operator/(const T val)
    {
        Polar<T> result;
        result = *(this);
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                result(i,j) /= val;

        return result;
    }

    /** Subtract a constant pixel-by-pixel
     */
    void operator-=(const T val)
    {
    	size_t imax=rings.size();
        for (size_t i = 0; i < imax; i++)
        {
            MultidimArray<T> &rings_i=rings[i];
            for (size_t j = 0; j < XSIZE(rings_i); j++)
                DIRECT_A1D_ELEM(rings_i,j) -= val;
        }
    }

    /** Add a constant pixel-by-pixel
     */
    void operator+=(const T val)
    {
    	size_t imax=rings.size();
        for (size_t i = 0; i < imax; i++)
        {
            MultidimArray<T> &rings_i=rings[i];
            for (int j = 0; j < XSIZE(rings_i); j++)
                DIRECT_A1D_ELEM(rings_i,j) += val;
        }
    }

    /** Multiply by a constant pixel-by-pixel
     */
    void operator*=(const T val)
    {
    	size_t imax=rings.size();
        for (size_t i = 0; i < imax; i++)
        {
            MultidimArray<T> &rings_i=rings[i];
            for (int j = 0; j < XSIZE(rings_i); j++)
                DIRECT_A1D_ELEM(rings_i,j) *= val;
        }
    }

    /** Divide by a constant pixel-by-pixel
     */
    void operator/=(const T val)
    {
    	size_t imax=rings.size();
        double ival=1.0/val;
        for (size_t i = 0; i < imax; i++)
        {
            MultidimArray<T> &rings_i=rings[i];
            for (size_t j = 0; j < XSIZE(rings_i); j++)
                DIRECT_A1D_ELEM(rings_i,j) *= ival;
        }
    }

    /** Subtract two polars pixel-by-pixel
     */
    Polar& operator-=(const Polar<T> in)
    {
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                (*(this))(i,j) -= in(i,j);
    }

    /** Add two polars pixel-by-pixel
    */
    Polar& operator+=(const Polar<T> in)
    {
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                (*(this))(i,j) += in(i,j);
    }

    /** Multiply two polars pixel-by-pixel
    */
    Polar& operator*=(const Polar<T> in)
    {
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
                (*(this))(i,j) *= in(i,j);
    }

    /** Divide two polars pixel-by-pixel
    */
    Polar& operator/=(const Polar<T> in)
    {
        for (size_t i = 0; i < rings.size(); i++)
            for (size_t j = 0; j < XSIZE(rings[i]); j++)
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
    const FileName & name() const
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

    /** Number of samples in outer ring access
     *
     * This function is used to know the number of samples in the outer ring.
     *
     * @code
     * std::cout << "Number of samples in outer ring: " << P.getSampleNoOuterRing() << std::endl;
     * @endcode
     */
    const int getSampleNoOuterRing() const
    {
        return XSIZE(rings[rings.size()-1]);
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
     * This operator can be used to access any ring of the polar as MultidimArray.
     *
     * @code
     * MultidimArray<double> secondring = P.getRing(1);
     * @endcode
     */
    MultidimArray< T >& getRing(int i)
    {
        return rings[i];
    }
    const  MultidimArray< T >& getRing(int i) const
    {
        return rings[i];
    }

    /** 1D Matrix access
     *
     * This operator can be used to set any ring of the polar with a MultidimArray.
     *
     * @code
     * MultidimArray<double> ring;
     * P.setRing(1,ring);
     * @endcode
     */
    void setRing(int i, MultidimArray< T > val)
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
        return DIRECT_A1D_ELEM(rings[r],f);
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
    void computeAverageAndStddev(double &avg, double &stddev,
                                 int mode = FULL_CIRCLES) const
    {
        double aux, sum = 0., sum2=0.;
        double twopi, w, N = 0;

        if (mode == FULL_CIRCLES)
            twopi = 2.*PI;
        else if (mode == HALF_CIRCLES)
            twopi = PI;
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT,"Incorrect mode for computeAverageAndStddev");

        size_t imax=rings.size();
        const double *ptrRing_radius=&ring_radius[0];
        for (size_t i = 0; i < imax; ++i, ++ptrRing_radius)
        {
            // take varying sampling into account
            w = (twopi * (*ptrRing_radius)) / (double) XSIZE(rings[i]);
            const MultidimArray<T> &rings_i=rings[i];
            for (size_t j = 0; j < XSIZE(rings_i); j++)
            {
                aux = DIRECT_A1D_ELEM(rings_i,j);
                double waux=w*aux;
                sum += waux;
                sum2 += waux * aux;
                N += w;
            }
        }
        if (N>0)
        {
            sum2 = sum2 / N;
            avg = sum / N;
            stddev=sqrt(fabs(sum2-avg*avg));
        }
        else if (N != 0.)
        {
            avg = sum;
            stddev=sqrt(fabs(sum2-avg*avg));
        }
        else
            stddev=avg=0;
    }

    void normalize(double average, double stddev)
    {
    	size_t imax=rings.size();
        double istddev=1.0/stddev;
        for (size_t i = 0; i < imax; i++)
        {
            MultidimArray<T> &rings_i=rings[i];
            for (size_t j = 0; j < XSIZE(rings_i); j++)
                DIRECT_A1D_ELEM(rings_i,j) = (DIRECT_A1D_ELEM(rings_i,j)-average)*istddev;
        }
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
#define GRIDDING_K 6
    void getCartesianCoordinates(std::vector<double> &x,
                                 std::vector<double> &y,
                                 std::vector<T> &data,
                                 const double extra_shell = GRIDDING_K/2)
    {
        double                     twopi, dphi,radius;
        int                        nsam;

        // Only for full circles for now!
        if (mode != FULL_CIRCLES)
            REPORT_ERROR(ERR_VALUE_INCORRECT,"VoronoiArea only implemented for FULL_CIRCLES mode of Polar");
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
                double sine, cosine;
                sincos(j*dphi,&sine,&cosine);
                x.push_back(radius*sine);
                y.push_back(radius*cosine);
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
                    double sine, cosine;
                    sincos(j*dphi,&sine,&cosine);
                    x.push_back(radius*sine);
                    y.push_back(radius*cosine);
                    data.push_back(0.);
                }
            }
        }

    }

    /** Convert cartesian MultidimArray to Polar using B-spline interpolation
     *
     * The input MultidimArray is assumed to be pre-processed for B-splines
     *
     * @code
     * Polar P;
     * MultidimArray<double> Maux;
     * img().produceSplineCoefficients(Maux,3);
     * P.getPolarFromCartesianBSpline(Maux,1,15);
     * @endcode
     *
     */
    void getPolarFromCartesianBSpline(const MultidimArray<T> &M1,
                                      int first_ring, int last_ring, int BsplineOrder=3,
                                      double xoff = 0., double yoff = 0.,
                                      double oversample1 = 1., int mode1 = FULL_CIRCLES)
    {
        int nsam;
        double radius, twopi, dphi, phi;
        double xp, yp, minxp, maxxp, minyp, maxyp;

        MultidimArray<T> Mring;
        rings.clear();
        ring_radius.clear();
        mode = mode1;
        oversample = oversample1;

        if (mode == FULL_CIRCLES)
            twopi = 2.*PI;
        else if (mode == HALF_CIRCLES)
            twopi = PI;
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT,"Incorrect mode for getPolarFromCartesian");


        // Limits of the matrix (not oversized!)
        minxp = FIRST_XMIPP_INDEX(XSIZE(M1));
        minyp = FIRST_XMIPP_INDEX(YSIZE(M1));
        maxxp = LAST_XMIPP_INDEX(XSIZE(M1));
        maxyp = LAST_XMIPP_INDEX(YSIZE(M1));
        double minxp_e=minxp-XMIPP_EQUAL_ACCURACY;
        double minyp_e=minyp-XMIPP_EQUAL_ACCURACY;
        double maxxp_e=maxxp+XMIPP_EQUAL_ACCURACY;
        double maxyp_e=maxyp+XMIPP_EQUAL_ACCURACY;

        // Loop over all polar coordinates
        for (int iring = first_ring; iring <= last_ring; iring++)
        {
            radius = (double) iring;
            // Non-constant sampling!! (always even for convenient Half2Whole of FTs)
            nsam = 2 * (int)( 0.5 * oversample * twopi * radius );
            nsam = XMIPP_MAX(1, nsam);
            dphi = twopi / (double)nsam;
            Mring.resizeNoCopy(nsam);
            for (int iphi = 0; iphi < nsam; iphi++)
            {
                // from polar to original cartesian coordinates
                phi = iphi * dphi;
                double sine, cosine;
                sincos(phi,&sine,&cosine);
                xp = sine * radius;
                yp = cosine * radius;

                // Origin offsets
                xp += xoff;
                yp += yoff;

                // Wrap coordinates
                if (xp < minxp_e || xp > maxxp_e)
                    xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);
                if (yp < minyp_e || yp > maxyp_e)
                    yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);

                // Perform the convolution interpolation
                if (BsplineOrder==1)
                    DIRECT_A1D_ELEM(Mring,iphi) = M1.interpolatedElement2D(xp,yp);
                else
                    DIRECT_A1D_ELEM(Mring,iphi) = M1.interpolatedElementBSpline2D(xp,yp,BsplineOrder);
            }
            rings.push_back(Mring);
            ring_radius.push_back(radius);
        }
    }

    /** Precalculate a vector with FFTW plans for all rings
     *
     */
    void calculateFftwPlans(Polar_fftw_plans &out)
    {
        (out.transformers).resize(rings.size());
        (out.arrays).resize(rings.size());
        for (size_t iring = 0; iring < rings.size(); iring++)
        {
            (out.arrays)[iring] = rings[iring];
            ((out.transformers)[iring]).setReal((out.arrays)[iring]);
        }
    }

};

/** Calculate FourierTransform of all rings
 *
 *  This function returns a polar of complex<double> by calculating
 *  the FT of all rings
 *
 *  Note that the Polar_fftw_plans may be re-used for polars of the same size
 * They should initially be calculated as in the example below
 *
 * @code
 * MultidimArray<double> angles, corr;
 * Polar_fftw_plans plans;
 * Polar<std::complex<double> > F1, F2;
 *
 * M1.calculateFftwPlans(plans); // M1 is a Polar<double>
 * fourierTransformRings(M1, F1, plans, false);
 * fourierTransformRings(M1, F2, plans, true);// complex conjugated
 * @endcode
 *
 *
 */
void fourierTransformRings(Polar<double > & in,
                           Polar<std::complex<double> > &out,
                           Polar_fftw_plans &plans,
                           bool conjugated = DONT_CONJUGATE);

/** Calculate inverse FourierTransform of all rings
 *
 *  This function returns a Polar<double> by calculating
 *  the inverse FT of all rings in a Polar<std::complex<double> >
 *
 * Note that the Polar_fftw_plans may be re-used for polars of the same size
 * They should initially be calculated as in the example below
 *
 * @code
 * MultidimArray<double> angles, corr;
 * Polar_fftw_plans plans;
 * Polar<std::complex<double> > F1, F2;
 *
 * M1.calculateFftwPlans(plans); // M1 is a Polar<double>
 * fourierTransformRings(M1, F1, plans);
 * inverseFourierTransformRings(F1, M1, plans);
 * @endcode
 *
 *
 */
void inverseFourierTransformRings(Polar<std::complex<double> > & in,
                                  Polar<double > &out,
                                  Polar_fftw_plans &plans,
                                  bool conjugated = DONT_CONJUGATE);

/** Rotational Correlation Aux */
class RotationalCorrelationAux
{
public:
	FourierTransformer local_transformer;
    MultidimArray<std::complex<double> > Fsum;
};

/** Fourier-space rotational Cross-Correlation Funtion
 *
 *  This function returns the rotational cross-correlation
 *  function of two complex Polars M1 and M2 using the
 *  cross-correlation convolution theorem.
 *
 *  Note that the local_transformer should have corr (with the right size)
 *  already in its fReal, and a Fourier Transform already calculated
 *
 * @code
 * MultidimArray<double> angles, corr;
 * FourierTransformer local_transformer;
 * Polar_fftw_plans plans;
 * Polar<std::complex<double> > F1, F2;
 *
 * // M1 and M2 are already Polar<double>
 *
 * M1.calculateFftwPlans(plans);
 * fourierTransformRings(M1, F1, plans, false);
 * fourierTransformRings(M1, F2, plans, true);// complex conjugated
 * corr.resize(M1.getSampleNoOuterRing());
 * local_transformer.setReal(corr, true);
 * local_transformer.FourierTransform();
 * rotationalCorrelation(F1,F2,angles,corr,local_transformer);
 * @endcode
 *
 * Note that this function can be used for real-space and
 * fourier-space correlations!!
 *
 */
void rotationalCorrelation(const Polar<std::complex<double> > &M1,
                           const Polar<std::complex<double> > &M2,
                           MultidimArray<double> &angles,
                           RotationalCorrelationAux &aux);

/** Compute a normalized polar Fourier transform of the input image.
    If plans is NULL, they are computed and returned. */
void normalizedPolarFourierTransform(const MultidimArray<double> &in,
                                     Polar< std::complex<double> > &out, bool flag,
                                     int first_ring, int last_ring, Polar_fftw_plans *&plans,
                                     int BsplineOrder=3);

/** Best rotation between two normalized polar Fourier transforms. */
double best_rotation(const Polar< std::complex<double> > &I1,
                     const Polar< std::complex<double> > &I2, RotationalCorrelationAux &aux);

/** Align I2 rotationally to I1 */
void alignRotationally(MultidimArray<double> &I1, MultidimArray<double> &I2,
					   RotationalCorrelationAux &aux,
                       int splineOrder=1, int wrap=WRAP);

/** Produce a polar image from a cartesian image.
 * You can give the minimum and maximum radius for the interpolation, the
 * delta radius and the delta angle.
 * It is assumed that the input image has the Xmipp origin.
 * Delta Ang must be in radians. */
void image_convertCartesianToPolar(MultidimArray<double> &in, MultidimArray<double> &out,
                                   double Rmin, double Rmax, double deltaR,
                                   double angMin, double angMax, double deltaAng);

/** Produce a polar image from a cartesian image.
 * Identical to the previous function.
 * The Zoom factor must be >=1 (==1 means no zoom).
 * If R is with a 0 size, then it is assumed that this routine
 * must compute the sampling radii. If it is not 0,
 * then it is assumed that the sampling radii have already been
 * calculated.
 */
void image_convertCartesianToPolar_ZoomAtCenter(const MultidimArray<double> &in,
												MultidimArray<double> &out,
												Matrix1D<double> &R,
												double zoomFactor,
												double Rmin, double Rmax, int NRSteps,
												double angMin, double angMax, int NAngSteps);

/** Produce a cylindrical volume from a cartesian volume.
 * You can give the minimum and maximum radius for the interpolation, the
 * delta radius and the delta angle.
 * It is assumed that the input image has the Xmipp origin.
 * Delta Ang must be in radians. */
void volume_convertCartesianToCylindrical(const MultidimArray<double> &in, MultidimArray<double> &out,
                                   double Rmin, double Rmax, double deltaR,
                                   double angMin, double angMax, double deltaAng,
                                   Matrix1D<double> &axis);

/** Produce a spherical volume from a cartesian voume.
 * You can give the minimum and maximum radii for the interpolation.
 * It is assumed that the input image has the Xmipp origin.
 * Delta Ang must be in radians.
 */
void volume_convertCartesianToSpherical(const MultidimArray<double> &in,
		MultidimArray<double> &out, double Rmin, double Rmax, double deltaR=1.,
		double deltaRot=2.*PI/360., double deltaTilt=PI/180.);
//@}
#endif
