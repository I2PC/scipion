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

#ifndef _CTF_HH
#define _CTF_HH

#include "xmipp_program.h"
#include "xmipp_filename.h"
#include "metadata.h"
#include "xmipp_fft.h"
//#include "ctf_estimate_from_psd_fast.h"


const int CTF_BASIC_LABELS_SIZE = 5;
const MDLabel CTF_BASIC_LABELS[] =
	{
		MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, MDL_CTF_DEFOCUS_ANGLE, MDL_CTF_CS, MDL_CTF_Q0
	};

const int CTF_ALL_LABELS_SIZE = 13;
const MDLabel CTF_ALL_LABELS[] =
    {
        MDL_CTF_VOLTAGE,
        MDL_CTF_DEFOCUSU,
        MDL_CTF_DEFOCUSV,
        MDL_CTF_DEFOCUS_ANGLE,
        MDL_CTF_CS,
        MDL_CTF_CA,
        MDL_CTF_ENERGY_LOSS,
        MDL_CTF_LENS_STABILITY,
        MDL_CTF_CONVERGENCE_CONE,
        MDL_CTF_LONGITUDINAL_DISPLACEMENT,
        MDL_CTF_TRANSVERSAL_DISPLACEMENT,
        MDL_CTF_Q0,
        MDL_CTF_K
    };

/** Check that a metadata contains ALL CTF_BASIC_LABELS labels */
bool containsCTFBasicLabels(const MetaData &md);

/** From a give images metadata, group by CTF_BASIC_LABELS and fill a ctfs metadata.
 * Params:
 *  imgMd: input images metadata.
 *  ctfMd: the ctfs metadata that will be filled.
 *
 * Raise error if neither CTF_MODEL or ALL CTF_BASIC_LABELS are found
 * in input images metadata.
 * */
void groupCTFMetaData(const MetaData &imgMd, MetaData &ctfMd, std::vector<MDLabel> &groupbyLabels);

/**@defgroup CTFSupport CTF support classes
   @ingroup DataLibrary */
//@{
/** Precomputed values for CTF evaluation. */
class PrecomputedForCTF
{
public:
    double u_sqrt;
    double u;
    double u2;
    double u3;
    double u4;
    double ang;
    double deltaf;
};

/** CTF class.
    Here goes how to compute the radial average of a parametric CTF:

    @code
      #include <Reconstruction/CTF.hh>

      int main() {
        // Read the parametric CTF
        CTFDescription CTF;
        CTF.enable_CTF=true;
        CTF.enable_CTFnoise=false;
        CTF.read("ctf_crio.param");
        CTF.produceSideInfo();

        // Compute the CTF radial average
        double sampling_rate=3.5; // in Angstroms/pixel
        double fmax=1.0/(2.0*sampling_rate);
        for (double f=0; f<fmax; f+=fmax/100.0) {
            double CTF_at_f=0;
            double N=0;
            for(double ang=0; ang<2*PI; ang+=0.01) {
               CTF.precomputeValues(f*cos(ang),f*sin(ang));
               CTF_at_f+=CTF.getValueAt();
               N++;
            }
            std::cout << f << " " << CTF_at_f/N << std::endl;
        }
        return 0;
      }
   @endcode

   Here is another sample program for generating the CTF

   @code
#include <data/args.h>
#include <data/filters.h>
#include <reconstruction/fourier_filter.h>
#include <data/ctf.h>

#include <data/xmipp_program.h>


int main(int argc, char **argv)
{
    FileName fn_ctf, fn_root;
    int Xdim;
    try
    {
        fn_ctf=getParameter(argc,argv,"-i");
        fn_root=getParameter(argc,argv,"-o");
        Xdim=textToInteger(getParameter(argc,argv,"-xdim"));
    }
    catch (XmippError XE)
    {
        std::cerr << XE << std::endl
        << "Usage: produce_imgs \n"
        << "         -i <CTF descr file>\n"
        << "         -o <filename root>\n"
        << "         -xdim <xdim>\n"
        ;
        return 1;
    }
    std::cerr << "fn_ctf = "  << fn_ctf <<std::endl
    << "fn_root = "  << fn_root <<std::endl
    << "Xdim = "  << Xdim <<std::endl;

    try
    {
        // Read CTF model
        Image<double> img(Xdim,Xdim);
        CTFDescription ctf;
        ctf.enable_CTF=true;
        ctf.enable_CTFnoise=true;
        ctf.read(fn_ctf);
        ctf.produceSideInfo();

        double avgdef = (ctf.DeltafU + ctf.DeltafV)/2.;
        ctf.DeltafU = avgdef;
        ctf.DeltafV = avgdef;
        MultidimArray<std::complex<double> >  ctfVal;

        ctf.generateCTF(Xdim, Xdim, ctfVal);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(ctfVal)
        {
            img(i,j)=ctfVal(i,j).real();
        }
        img.write("ctf.noisy");

        ctf.clear();
        ctf.enable_CTF=true;
        ctf.enable_CTFnoise=false;
        ctf.read(fn_ctf);
        ctf.produceSideInfo();
        ctf.generateCTF(Xdim, Xdim, ctfVal);

        FOR_ALL_ELEMENTS_IN_ARRAY2D(ctfVal)
        {
            // Power spectrum density
            img(i,j)*=ctfVal(i,j).real();
        }
        img.write("ctf.noisyless");

        //      CenterFFT(Ipsd(), true); Ipsd.write(fn_root+"_psd.xmp");

    }
    catch (XmippError XE)
    {
        std::cout << XE << std::endl;
    }
    return 0;
}
   @endcode

*/

/////////////////////////////////CTF1D///////////////////////////////////////////

class CTFDescription1D
{
public:
	// Electron wavelength (Amstrongs)
	double lambda;
	// Squared frequency associated to the aperture
	// double ua2;
	// Different constants
	double K1;
	double K2;
	double K3;
	double K4;
	double K5;
	double K6;
	double K7;
	double Ksin;
	double Kcos;
	/** Standard error of defocus Gaussian function due to chromatic aberration.
		in Amstrong */
	double D;
	// Precomputed values
	PrecomputedForCTF precomputed;
	// Image of precomputed values
	std::vector<PrecomputedForCTF> precomputedImage;
	// Xdim size of the image
	int precomputedImageXdim;
public:
	/// Global gain. By default, 1
	double K;
	/// Sampling rate (A/pixel)
	double Tm;
	/// Accelerating Voltage (in KiloVolts)
	double kV;
	/// Defocus (in Angstroms). Negative values are underfocused
	double Defocus;
	/// Spherical aberration (in milimeters). Typical value 5.6
	double Cs;
	/// Chromatic aberration (in milimeters). Typical value 2
	double Ca;
	/** Mean energy loss (eV) due to interaction with sample.
		Typical value 1*/
	double espr;
	/// Objective lens stability (deltaI/I) (ppm). Typical value 1
	double ispr;
	/// Convergence cone semiangle (in mrad). Typical value 0.5
	double alpha;
	/// Longitudinal mechanical displacement (ansgtrom). Typical value 100
	double DeltaF;
	/// Transversal mechanical displacement (ansgtrom). Typical value 3
	double DeltaR;
	/// Factor for the importance of the Amplitude contrast.
	double Q0;
	/// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
	double x0;
	/// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
	double xF;
	/// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
	double y0;
	/// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
	double yF;
	/// Local CTF determination
	bool isLocalCTF;
	/// Enable CTFnoise part
	bool enable_CTFnoise;
	/// Enable CTF part
	bool enable_CTF;
	/// Global base_line
	double base_line;
	/// Gain for the gaussian term
	double gaussian_K;
	/// Gaussian width
	double sigma1;
	/// Gaussian center
	double Gc1;
	/// Gain for the square root term
	double sqrt_K;
	/// Sqrt width
	double sq;
	/// Gain for the second Gaussian term
	double gaussian_K2;
	/// Second Gaussian width
	double sigma2;
	/// Second Gaussian center
	double Gc2;
	// Background polynomial
	double bgR1, bgR2, bgR3;
	// Envelope polynomial
	double envR0, envR1, envR2;
	//Maximum frequency to estimate values
	double freq_max;
	//Extra parameters for VPP
	double phase_shift;
	double VPP_radius;

	/** Empty constructor. */
	CTFDescription1D()
	{
	  clear();
	}

	/** Read from file.
			An exception is thrown if the file cannot be open.

			If no K or sqrt_K are given then it is assumed that the user
			does not want to activate that part and the noise or the CTF
			are removed from the model unless the disable_if_not_K is set
			to false*/
	 void read(const FileName &fn, bool disable_if_not_K = true);

		/** Read from 1 row of a metadata file.

			If no K or sqrt_K are given then it is assumed that the user
			does not want to activate that part and the noise or the CTF
			are removed from the model unless the disable_if_not_K is set
			to false
			This function should be used usually if you have all ctf parameters
			in columns format metadata or after calling fillExpand having CTF_MODEL
			*/
	  void readFromMetadataRow(const MetaData &MD, size_t id, bool disable_if_not_K=true);

		/** Same as previuous reading function but providing the MDRow */
	  void readFromMdRow(const MDRow &row, bool disable_if_not_K=true);

	  /** Write current CTF model to row */
	  void setRow(MDRow &row) const;

	  /** Write to file.
			An exception is thrown if the file cannot be open.*/
	  void write(const FileName &fn);

		/// Define parameters in the command line
	  static void defineParams(XmippProgram * program);

		/// Read parameters from the command line
	  void readParams(XmippProgram * program);

		/// Show
	  friend std::ostream & operator << (std::ostream &out, const CTFDescription1D &ctf);

	  /// Clear.
	  void clear();

		/// Clear noise
	  void clearNoise();

		/// Clear pure CTF
	  void clearPureCtf();

		/** Change sampling rate.
		 * It is advisable to change the sampling rate just after reading the CTF parameters.
		 * However, unless you use precomputed values, there should not be any problem of changing
		 * it at a later stage.
		 */
	  inline void changeSamplingRate(double newTm)
	  {
		  Tm=newTm;
	  }

	  /// Produce Side information
	  void produceSideInfo();

	  /// Precompute values for a given frequency
	 void precomputeValues(double X)
	 {
		 precomputed.ang = 0;
		 precomputed.u2 = X * X;
		 precomputed.u = sqrt(precomputed.u2);
		 precomputed.u3 = precomputed.u2 * precomputed.u;
		 precomputed.u4 = precomputed.u2 * precomputed.u2;
		 precomputed.u_sqrt = sqrt(precomputed.u);

		 if (fabs(X) < XMIPP_EQUAL_ACCURACY)
			 precomputed.deltaf=0;
		 else
		 {

			 /*
			  * For a derivation of this formulae confer
			  * Principles of Electron Optics page 1380
			  * in particular term defocus and twofold axial astigmatism
			  * take into account that a1 and a2 are the coefficient
			  * of the zernike polynomials difference of defocus at 0
			  * and at 45 degrees. In this case a2=0
			  */
			 precomputed.deltaf= Defocus;//(defocus_average + defocus_deviation*cos_ellipsoid_ang_2);

		 }
	 }

	 /// Precompute values for an image
	 void precomputeValues(const MultidimArray<double> &cont_x_freq);

	 /// Precompute values for a given frequency
	 void precomputeValues(int i)
	 {
		 precomputed=precomputedImage[i+precomputedImageXdim];
		 if (precomputed.deltaf==-1)
		 {
			 precomputed.deltaf= Defocus;
		 }
	 }

	 /// Compute CTF at (U,V). Continuous frequencies
	 double getValueAt(bool show = false) const
	 {
		 double pure_CTF = enable_CTF ? getValuePureAt(show) : 0;
		 return enable_CTFnoise ? sqrt(pure_CTF*pure_CTF + getValueNoiseAt(show)) : pure_CTF;
	 }

	 /// Compute CTF damping at (U,V). Continuous frequencies
	 inline double getValueDampingAt(bool show = false) const
	 {
		 double Eespr = exp(-K3 * precomputed.u4); // OK
		 double EdeltaF = bessj0(K5 * precomputed.u2); // OK
		 double EdeltaR = SINC(precomputed.u * DeltaR); // OK
		 double aux=K7 * precomputed.u2 * precomputed.u + precomputed.deltaf * precomputed.u;
		 double Ealpha = exp(-K6 * aux * aux);
		 double E = Eespr * EdeltaF * EdeltaR * Ealpha+envR0+envR1*precomputed.u+envR2*precomputed.u2;
		 if (E < 0)
			 E	= 0;
		 if (show)
		 {
			 std::cout << "   Deltaf=" << precomputed.deltaf << std::endl;
			 std::cout << "   u,u2,u4=" << precomputed.u << " " << precomputed.u2
			 << " " << precomputed.u4 << std::endl;
			 std::cout << "   K3,Eespr=" << K3 << " " << Eespr << std::endl;
			 std::cout << "   K4,Eispr=" << K4 << " " << /*Eispr <<*/ std::endl;
			 std::cout << "   K5,EdeltaF=" << K5 << " " << EdeltaF << std::endl;
			 std::cout << "   EdeltaR=" << EdeltaR << std::endl;
			 std::cout << "   K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
			 << std::endl;
			 std::cout << "   Total atenuation(E)= " << E << std::endl;
			 std::cout << "   CTFdamp=" << E << std::endl;
		 }
		 return -K*E;
	 }

	 /// Compute CTF pure at (U,V). Continuous frequencies
	inline double getValuePureAt(bool show = false) const
	{
		double VPP;
		double check_VPP = round(VPP_radius*1000);
		if(check_VPP != 0)
			VPP = phase_shift*(1-exp(-precomputed.u2/(2*pow(VPP_radius,2.0))));
		else
			VPP = 0;
		double argument = VPP + K1 * precomputed.deltaf * precomputed.u2 + K2 *precomputed.u4;
		double sine_part, cosine_part;
		sincos(argument,&sine_part, &cosine_part); // OK
		double Eespr = exp(-K3 * precomputed.u4); // OK
		//CO: double Eispr=exp(-K4*u4); // OK
		double EdeltaF = bessj0(K5 * precomputed.u2); // OK
		double EdeltaR = SINC(precomputed.u * DeltaR); // OK
		double aux=(K7 * precomputed.u2 * precomputed.u + precomputed.deltaf * precomputed.u);
		double Ealpha = exp(-K6 * aux * aux); // OK
		// CO: double E=Eespr*Eispr*EdeltaF*EdeltaR*Ealpha;
		double E = Eespr * EdeltaF * EdeltaR * Ealpha + envR0+envR1*precomputed.u+envR2*precomputed.u2;
		if (E < 0)
			E	= 0;
		if (show)
		{
			std::cout << "   Deltaf=" << precomputed.deltaf << std::endl;
			std::cout << "   u,u2,u4=" << precomputed.u << " " << precomputed.u2
			<< " " << precomputed.u4 << std::endl;
			std::cout << "   K1,K2,argument=" << K1 << " " << K2 << " " << argument << std::endl;
			std::cout << "   K3,Eespr=" << K3 << " " << Eespr << std::endl;
			std::cout << "   K4,Eispr=" << K4 << " " << /*Eispr <<*/ std::endl;
			std::cout << "   K5,EdeltaF=" << K5 << " " << EdeltaF << std::endl;
			std::cout << "   EdeltaR=" << EdeltaR << std::endl;
			std::cout << "   K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
			<< std::endl;
			std::cout << "   Total atenuation(E)= " << E << std::endl;
			std::cout << "   K,Q0,base_line=" << K << "," << Q0 << "," << base_line << std::endl;
			std::cout << "   CTF="
			<< -K*(Ksin*sine_part - Kcos*cosine_part)*E << std::endl;
		}
		return -K*(Ksin*sine_part - Kcos*cosine_part)*E;
	}

	/// Compute CTF pure at (U,V). Continuous frequencies
	inline double getValuePureNoKAt() const
	{
		double VPP;
		double check_VPP = round(VPP_radius*1000);
		if(check_VPP != 0)
			VPP = phase_shift*(1-exp(-precomputed.u2/(2*pow(VPP_radius,2.0))));
		else
			VPP = 0;
		double argument = VPP + K1 * precomputed.deltaf * precomputed.u2 + K2 *precomputed.u4;
		double sine_part, cosine_part;
		sincos(argument,&sine_part, &cosine_part); // OK
		double Eespr = exp(-K3 * precomputed.u4); // OK
		//CO: double Eispr=exp(-K4*u4); // OK
		double EdeltaF = bessj0(K5 * precomputed.u2); // OK
		double EdeltaR = SINC(precomputed.u * DeltaR); // OK
		double aux=(K7 * precomputed.u2 * precomputed.u + precomputed.deltaf * precomputed.u);
		double Ealpha = exp(-K6 * aux * aux); // OK
		// CO: double E=Eespr*Eispr*EdeltaF*EdeltaR*Ealpha;
		double E = Eespr * EdeltaF * EdeltaR * Ealpha+envR0+envR1*precomputed.u+envR2*precomputed.u2;
		if (E < 0)
			E	= 0;
		return -(Ksin*sine_part - Kcos*cosine_part)*E;
	}

	/// Compute noise at (X,Y). Continuous frequencies, notice it is squared
	//#define DEBUG
	inline double getValueNoiseAt(bool show = false) const
	{
		double c = Gc1;
		double sigmaG1 = sigma1;

		double c2 = Gc2;
		double sigmaG2 = sigma2;

		if(show)
		{
			std::cout
			<< "   sq=" << sq << "\n"
			<< "   c=" << c << "\n"
			<< "   sigma=" << sigmaG1 << "\n"
			<< "   c2=" << c2 << "\n"
			<< "   sigma2=" << sigmaG2 << "\n"
			<< "   gaussian_K=" << gaussian_K << "\n"
			<< "   sqrt_K=" << sqrt_K << "\n"
			<< "   base_line=" << base_line << "\n";
			std::cout << "   u=" << precomputed.u << "u_sqrt=" << precomputed.u_sqrt <<") CTFnoise="
			<< base_line +
			gaussian_K*exp(-sigmaG1*(precomputed.u - c)*(precomputed.u - c)) +
			sqrt_K*exp(-sq*sqrt(precomputed.u)) -
			gaussian_K2*exp(-sigmaG2*(precomputed.u - c2)*(precomputed.u - c2)) << std::endl;
		}

		double aux=precomputed.u - c;
		double aux2=precomputed.u - c2;
		return base_line +
			   gaussian_K*exp(-sigmaG1*aux*aux) +
			   sqrt_K*exp(-sq*precomputed.u_sqrt) -
			   gaussian_K2*exp(-sigmaG2*aux2*aux2)+
			   bgR1*precomputed.u+bgR2*precomputed.u2+bgR3*precomputed.u3;
	}

	/// Compute pure CTF without damping at (U,V). Continuous frequencies
	double getValuePureWithoutDampingAt(bool show = false) const
	{
		double VPP;
		double check_VPP = round(VPP_radius*1000);
		if(check_VPP != 0)
			VPP = phase_shift*(1-exp(-precomputed.u2/(2*pow(VPP_radius,2.0))));
		else
			VPP = 0;
		double argument = VPP + K1 * precomputed.deltaf * precomputed.u2 + K2 * precomputed.u4;
		double sine_part, cosine_part;
		sincos(argument,&sine_part,&cosine_part);

		if (show)
		{
			std::cout << "   Deltaf=" << precomputed.deltaf << std::endl;
			std::cout << "   u,u2,u4=" << precomputed.u << " " << precomputed.u2
			<< " " << precomputed.u4 << std::endl;
			std::cout << "   K1,K2,sin=" << K1 << " " << K2 << " "
			<< std::endl;
			std::cout << "   Q0=" << Q0 << std::endl;
			std::cout << "   CTF without damping="
			<< -(Ksin*sine_part - Kcos*cosine_part) << std::endl;
		}
		return -(Ksin*sine_part - Kcos*cosine_part);
	}

	 /// Compute CTF pure at (U,V). Continuous frequencies
	inline double getValuePureNoPrecomputedAt(double X, bool show = false) const
	{
		double u2 = X * X;
		double u = sqrt(u2);
		double u4 = u2 * u2;
		double deltaf = Defocus;
		double VPP;
		double check_VPP = round(VPP_radius*1000);
		if(check_VPP != 0)
			VPP = phase_shift*(1-exp(-u2/(2*pow(VPP_radius,2.0))));
		else
			VPP = 0;
		double argument = VPP + K1 * deltaf * u2 + K2 * u4;
		double sine_part, cosine_part;
		sincos(argument,&sine_part, &cosine_part); // OK
		double Eespr = exp(-K3 * u4); // OK
		//CO: double Eispr=exp(-K4*u4); // OK
		double EdeltaF = bessj0(K5 * u2); // OK
		double EdeltaR = sinc(u * DeltaR); // OK
		double Ealpha = exp(-K6 * (K7 * u2 * u + deltaf * u) * (K7 * u2 * u + deltaf * u)); // OK
		// CO: double E=Eespr*Eispr*EdeltaF*EdeltaR*Ealpha;
		double E = Eespr * EdeltaF * EdeltaR * Ealpha+envR0+envR1*precomputed.u+envR2*precomputed.u2;
		if (E < 0)
			E	= 0;
		if (show)
		{
			std::cout << " Deltaf=" << deltaf << std::endl;
			std::cout << " u,u2,u4=" << u << " " << u2 << " " << u4 << std::endl;
			std::cout << " K1,K2,sin=" << K1 << " " << K2 << " "
			<< sine_part << std::endl;
			std::cout << " K3,Eespr=" << K3 << " " << Eespr << std::endl;
			//std::cout << " K4,Eispr=" << K4 << " " << /*Eispr*/ << std::endl;
			std::cout << " K5,EdeltaF=" << K5 << " " << EdeltaF << std::endl;
			std::cout << " EdeltaR=" << EdeltaR << std::endl;
			std::cout << " K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
			<< std::endl;
			std::cout << " Total atenuation(E)= " << E << std::endl;
			std::cout << " K,Q0,base_line=" << K << "," << Q0 << "," << base_line << std::endl;
			std::cout << " (X)=(" << X << ") CTF="
			<< -K*(Ksin*sine_part - Kcos*cosine_part)*E + base_line << std::endl;
		}
		return -K*(Ksin*sine_part - Kcos*cosine_part)*E;
	}


	/** Returns the continuous frequency of the zero, maximum or minimum number n in the direction u.
		u must be a unit vector, n=1,2,... Returns (-1,-1) if it is not found
		'iwhat' can be 0 (zero), 1(max), or -1 (min) */
	void lookFor(int n, const Matrix1D<double> &u, Matrix1D<double> &freq, int iwhat=0);

	/// Apply CTF to an image
	void applyCTF(MultidimArray < std::complex<double> > &FFTI, const MultidimArray<double> &I, double Ts, bool absPhase=false);

	/// Apply CTF to an image
	void applyCTF(MultidimArray <double> &I, double Ts, bool absPhase=false);

	/** Generate CTF image.
		The sample image is used only to take its dimensions. */
	template <class T1, class T2>
	void generateCTF(const MultidimArray<T1> &sample_image, MultidimArray <T2> &CTF, double Ts=-1)
	{
		if ( ZSIZE(sample_image) > 1 )
			REPORT_ERROR(ERR_MULTIDIM_DIM,"ERROR: Generate_CTF only works with 2D sample images, not 3D.");
		generateCTF(YSIZE(sample_image), XSIZE(sample_image), CTF, Ts);
		STARTINGX(CTF) = STARTINGX(sample_image);
		STARTINGY(CTF) = STARTINGY(sample_image);
	}

	/** Get profiles along a given direction.
	 * It returns the CTF profiles (BGNOISE, ENVELOPE, PSD, CTF) along the
	 * direction defined by angle. Profiles
	 * have nsamples from 0 to fmax (1/A).
	 */
	void getProfile(double fmax, int nsamples, MultidimArray<double> &profiles);

	/** Get radial average profiles.
	 * It returns the radially average CTF profiles (BGNOISE, ENVELOPE, PSD, CTF). Profiles
	 * have nsamples from 0 to fmax (1/A).
	 */
	void getAverageProfile(double fmax, int nsamples, MultidimArray<double> &profiles);

	//#define DEBUG
   /// Generate CTF image.
	template <class T>
	void generateCTF(int Ydim, int Xdim, MultidimArray < T > &CTF, double Ts=-1)
	{
		CTF.resizeNoCopy(Ydim, Xdim);
		if (Ts<0)
			Ts=Tm;
		#ifdef DEBUG
			std::cout << "CTF:\n" << *this << std::endl;
		#endif

		double iTs=1.0/Ts;
		for (int i=0; i<Ydim; ++i)
		{
			double wy;
			FFT_IDX2DIGFREQ(i, YSIZE(CTF), wy);
			double fy=wy*iTs;
			for (int j=0; j<Xdim; ++j)
			{
				double wx;
				FFT_IDX2DIGFREQ(j, XSIZE(CTF), wx);
				double fx=wx*iTs;
				precomputeValues(fx);
				A2D_ELEM(CTF, i, j) = (T) getValueAt();
				#ifdef DEBUG
						if (i == 0)
							std::cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
							<< " " << CTF(i, j) << std::endl;
				#endif
			}
		}
	}
	#undef DEBUG

	template <class T>
	void generateCTFWithoutDamping(int Ydim, int Xdim, MultidimArray < T > &CTF, double Ts=-1)
	{
		CTF.resizeNoCopy(Ydim, Xdim);
		if (Ts<0)
			Ts=Tm;
		#ifdef DEBUG
			std::cout << "CTF:\n" << *this << std::endl;
		#endif

		double iTs=1.0/Ts;
		for (int i=0; i<Ydim; ++i)
		{
			double wy;
			FFT_IDX2DIGFREQ(i, YSIZE(CTF), wy);
			double fy=wy*iTs;
			for (int j=0; j<Xdim; ++j)
			{
				double wx;
				FFT_IDX2DIGFREQ(j, XSIZE(CTF), wx);
				double fx=wx*iTs;
				precomputeValues(fx);
				A2D_ELEM(CTF, i, j) = (T) getValuePureWithoutDampingAt();
				#ifdef DEBUG
						if (i == 0)
							std::cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
							<< " " << CTF(i, j) << std::endl;
				#endif
			}
		}
	}
	#undef DEBUG

   /** Check physical meaning.
		true if the CTF parameters have physical meaning.
		Call this function after produstd::cing side information */
	bool hasPhysicalMeaning();

	/** Force physical meaning.*/
	void forcePhysicalMeaning();

};

/** Generate CTF 2D image with two CTFs.
 * The two CTFs are in fn1 and fn2. The output image is written to the file fnOut and has size Xdim x Xdim. */
void generateCTFImageWith2CTFs(const MetaData &MD1, const MetaData &MD2, int Xdim, MultidimArray<double> &imgOut);

/** compute error between two CTFS, return a single value */
double errorBetween2CTFs( MetaData &MD1,
                          MetaData &MD2,
                         size_t dim,
                         double minFreq=0.05,
                         double maxFreq=0.25);
/** Report at which resolution an astigmatism  CTF is shifted in perpendicular directions
 * by phaseRad degrees. Returns amsgtroms
 */
double errorMaxFreqCTFs( MetaData &MD1,
                         double phaseRad=HALFPI);

/** Report at which resolution these two CTF are shifted by phaseRad degrees
 *returns amsgtroms
 */
double errorMaxFreqCTFs2D( MetaData &MD1,
                         MetaData &MD2,
                         size_t xDim=256,
                         double phaseRad=HALFPI);

/** Generate an image with the PSD and the CTF
 *  Before calling the function img must have the enhanced PSD. The enhanced PSD image is modified to add the CTF. */
void generatePSDCTFImage(MultidimArray<double> &img, const MetaData &MD);


///////////////////////////// CTF2D ////////////////////////////////////////////////////

class CTFDescription: public CTFDescription1D
{
public:
    // Electron wavelength (Amstrongs)
    //double lambda;
    // Squared frequency associated to the aperture
    // double ua2;
    // Different constants
    //double K1;
    //double K2;
    //double K3;
    //double K4;
    //double K5;
    //double K6;
    //double K7;
    //double Ksin;
    //double Kcos;
    // Azimuthal angle in radians
    double rad_azimuth;
    // defocus_average = -(defocus_u + defocus_v)/2
    double defocus_average;
    // defocus_deviation = -(defocus_u - defocus_v)/2
    double defocus_deviation;
    // Gaussian angle in radians
    double rad_gaussian;
    // Second Gaussian angle in radians
    double rad_gaussian2;
    // Sqrt angle in radians
    double rad_sqrt;
    /** Standard error of defocus Gaussian function due to chromatic aberration.
        in Amstrong */
    //double D;
    // Precomputed values
    //PrecomputedForCTF precomputed;
    // Image of precomputed values
    //std::vector<PrecomputedForCTF> precomputedImage;
    // Xdim size of the image
    //int precomputedImageXdim;
public:
    /// Global gain. By default, 1
    //double K;
    /// Sampling rate (A/pixel)
    //double Tm;
    /// Accelerating Voltage (in KiloVolts)
    //double kV;
    /// Defocus in U (in Angstroms). Negative values are underfocused
    double DeltafU;
    /// Defocus in V (in Angstroms). Negative values are underfocused
    double DeltafV;
    /// Azimuthal angle (between X and U) in degrees
    double azimuthal_angle;
    // Radius of the aperture (in micras)
    // double aperture;
    /// Spherical aberration (in milimeters). Typical value 5.6
    //double Cs;
    /// Chromatic aberration (in milimeters). Typical value 2
    //double Ca;
    /** Mean energy loss (eV) due to interaction with sample.
        Typical value 1*/
    //double espr;
    /// Objective lens stability (deltaI/I) (ppm). Typical value 1
    //double ispr;
    /// Convergence cone semiangle (in mrad). Typical value 0.5
    //double alpha;
    /// Longitudinal mechanical displacement (ansgtrom). Typical value 100
    //double DeltaF;
    /// Transversal mechanical displacement (ansgtrom). Typical value 3
    //double DeltaR;
    /// Factor for the importance of the Amplitude contrast.
    //double Q0;
    /// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
    //double x0;
    /// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
    //double xF;
    /// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
    //double y0;
    /// In the case of local CTF determination x0,xF,y0,yF determines the region where the CTF is determined
    //double yF;
    /// Local CTF determination
    //bool isLocalCTF;
    /// Enable CTFnoise part
    //bool enable_CTFnoise;
    /// Enable CTF part
    //bool enable_CTF;
    /// Global base_line
    //double base_line;
    /// Gain for the gaussian term
    //double gaussian_K;
    /// Gaussian width U
    double sigmaU;
    /// Gaussian width V
    double sigmaV;
    /// Gaussian center for U
    double cU;
    /// Gaussian center for V
    double cV;
    /// Gaussian angle
    double gaussian_angle;
    /// Gain for the square root term
    //double sqrt_K;
    /// Sqrt width U
    double sqU;
    /// Sqrt width V
    double sqV;
    /// Sqrt angle
    double sqrt_angle;
    /// Gain for the second Gaussian term
    double gaussian_K2;
    /// Second Gaussian width U
    double sigmaU2;
    /// Second Gaussian width V
    double sigmaV2;
    /// Second Gaussian center for U
    double cU2;
    /// Second Gaussian center for V
    double cV2;
    /// Second Gaussian angle
    double gaussian_angle2;

    CTFDescription()
    {
    	clear();
    }

    CTFDescription(CTFDescription1D copy)
    {
    	DeltafU = DeltafV = copy.Defocus;
    	Ca = copy.Ca;
    	Cs = copy.Cs;
    	DeltaF = copy.DeltaF;
    	DeltaR = copy.DeltaR;
    	Q0 = copy.Q0;
    	cU = cV = copy.Gc1;
    	cU2 = cV2 = copy.Gc2;
    	base_line = copy.base_line;
    	gaussian_K = copy.gaussian_K;
    	gaussian_K2 = copy.gaussian_K2;
    	sigmaU = sigmaV = copy.sigma1;
    	sigmaU2 = sigmaV2 = copy.sigma2;
    	sqU = sqV = copy.sq;
    	kV = copy.kV;
    	bgR1 = copy.bgR1;
    	bgR2 = copy.bgR2;
    	bgR3 = copy.bgR3;
    	envR0 = copy.envR0;
    	envR1 = copy.envR1;
    	envR2 = copy.envR2;
    	espr = copy.espr;
    	ispr = copy.ispr;
    	alpha = copy.alpha;
    	sqrt_K = copy.sqrt_K;
    	K = copy.K;
    	Tm = copy.Tm;
    	azimuthal_angle = 0;
    	gaussian_angle = 0;
    	gaussian_angle2 = 0;
    	sqrt_angle = 0;
    	enable_CTFnoise = copy.enable_CTFnoise;
    	phase_shift = copy.phase_shift;
    	VPP_radius = copy.VPP_radius;
    }

    /** Read from file.
		An exception is thrown if the file cannot be open.

		If no K or sqrt_K are given then it is assumed that the user
		does not want to activate that part and the noise or the CTF
		are removed from the model unless the disable_if_not_K is set
		to false*/
    void read(const FileName &fn, bool disable_if_not_K = true);
    /** Read from 1 row of a metadata file.

		If no K or sqrt_K are given then it is assumed that the user
		does not want to activate that part and the noise or the CTF
		are removed from the model unless the disable_if_not_K is set
		to false
		This function should be used usually if you have all ctf parameters
		in columns format metadata or after calling fillExpand having CTF_MODEL
		*/
    void readFromMetadataRow(const MetaData &MD, size_t id, bool disable_if_not_K=true);

    /** Same as previuous reading function but providing the MDRow */
    void readFromMdRow(const MDRow &row, bool disable_if_not_K=true);

    /** Write current model to row */
    void setRow(MDRow &row) const;

    /** Write to file.
        An exception is thrown if the file cannot be open.*/
    void write(const FileName &fn);

    /// Define parameters in the command line
    static void defineParams(XmippProgram * program);

    /// Read parameters from the command line
    void readParams(XmippProgram * program);

    /// Show
    friend std::ostream & operator << (std::ostream &out, const CTFDescription &ctf);

    /// Clear.
    void clear();

    /// Clear noise
    void clearNoise();

    /// Clear pure CTF
    void clearPureCtf();

    /** Change sampling rate.
     * It is advisable to change the sampling rate just after reading the CTF parameters.
     * However, unless you use precomputed values, there should not be any problem of changing
     * it at a later stage.
     */
    inline void changeSamplingRate(double newTm)
    {
    	Tm=newTm;
    }
    /// Produce Side information
    void produceSideInfo();

    /// Precompute values for a given frequency
    void precomputeValues(double X, double Y)
    {
        precomputed.ang=atan2(Y, X);
        precomputed.u2 = X * X + Y * Y;
        precomputed.u = sqrt(precomputed.u2);
        precomputed.u3 = precomputed.u2 * precomputed.u;
        precomputed.u4 = precomputed.u2 * precomputed.u2;
        precomputed.u_sqrt = sqrt(precomputed.u);

        if (fabs(X) < XMIPP_EQUAL_ACCURACY &&
            fabs(Y) < XMIPP_EQUAL_ACCURACY)
            precomputed.deltaf=0;
        else
        {
            double ellipsoid_ang = precomputed.ang - rad_azimuth;
            /*
             * For a derivation of this formulae confer
             * Principles of Electron Optics page 1380
             * in particular term defocus and twofold axial astigmatism
             * take into account that a1 and a2 are the coefficient
             * of the zernike polynomials difference of defocus at 0
             * and at 45 degrees. In this case a2=0
             */
            double cos_ellipsoid_ang_2 = cos(2*ellipsoid_ang);
            precomputed.deltaf= (defocus_average + defocus_deviation*cos_ellipsoid_ang_2);
        }
    }

    /// Precompute values for an image
    void precomputeValues(const MultidimArray<double> &cont_x_freq,
                          const MultidimArray<double> &cont_y_freq);

    /// Precompute values for a given frequency
    void precomputeValues(int i, int j)
    {

        precomputed=precomputedImage[i*precomputedImageXdim+j];

        if (precomputed.deltaf==-1)
        {
            double ellipsoid_ang = precomputed.ang - rad_azimuth;
            double cos_ellipsoid_ang_2 = cos(2*ellipsoid_ang);
            precomputed.deltaf= (defocus_average + defocus_deviation*cos_ellipsoid_ang_2);

        }
    }

    /// Compute CTF at (U,V). Continuous frequencies
    double getValueAt(bool show = false) const
    {
        double pure_CTF = enable_CTF ? getValuePureAt(show) : 0;
        return enable_CTFnoise ? sqrt(pure_CTF*pure_CTF + getValueNoiseAt(show)) : pure_CTF;
    }

    /// Compute pure CTF without damping at (U,V). Continuous frequencies
    double getValueArgument(bool show = false) const
    {
        return K1 * precomputed.deltaf * precomputed.u2 + K2 * precomputed.u4;
    }

    /// Get Phase of the CTF
    inline double getPhaseAt() const
    {
        return K1 * precomputed.deltaf * precomputed.u2 + K2 *precomputed.u4;
    }

    /// Compute CTF pure at (U,V). Continuous frequencies
    inline double getValuePureNoPrecomputedAt(double X, double Y, bool show = false) const
    {
        double u2 = X * X + Y * Y;
        double u = sqrt(u2);
        double u4 = u2 * u2;
        //if (u2>=ua2) return 0;
        double deltaf = getDeltafNoPrecomputed(X, Y);
        double VPP;
        double check_VPP = round(VPP_radius*1000);
		if(check_VPP != 0)
			VPP = phase_shift*(1-exp(-u2/(2*pow(VPP_radius,2.0))));
		else
			VPP = 0;
        double argument = VPP + K1 * deltaf * u2 + K2 * u4;
        double sine_part, cosine_part;
        sincos(argument,&sine_part, &cosine_part); // OK
        double Eespr = exp(-K3 * u4); // OK
        //CO: double Eispr=exp(-K4*u4); // OK
        double EdeltaF = bessj0(K5 * u2); // OK
        double EdeltaR = sinc(u * DeltaR); // OK
        double Ealpha = exp(-K6 * (K7 * u2 * u + deltaf * u) * (K7 * u2 * u + deltaf * u)); // OK
        // CO: double E=Eespr*Eispr*EdeltaF*EdeltaR*Ealpha;
        double E = Eespr * EdeltaF * EdeltaR * Ealpha+envR0+envR1*precomputed.u+envR2*precomputed.u2;
        if (E < 0)
        	E	= 0;
        if (show)
        {
            std::cout << " Deltaf=" << deltaf << std::endl;
            std::cout << " u,u2,u4=" << u << " " << u2 << " " << u4 << std::endl;
            std::cout << " K1,K2,sin=" << K1 << " " << K2 << " "
            << sine_part << std::endl;
            std::cout << " K3,Eespr=" << K3 << " " << Eespr << std::endl;
            std::cout << " K4,Eispr=" << K4 << " " << /*Eispr <<*/ std::endl;
            std::cout << " K5,EdeltaF=" << K5 << " " << EdeltaF << std::endl;
            std::cout << " EdeltaR=" << EdeltaR << std::endl;
            std::cout << " K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
            << std::endl;
            std::cout << " Total atenuation(E)= " << E << std::endl;
            std::cout << " K,Q0,base_line=" << K << "," << Q0 << "," << base_line << std::endl;
            std::cout << " (X,Y)=(" << X << "," << Y << ") CTF="
            << -K*(Ksin*sine_part - Kcos*cosine_part)*E + base_line << std::endl;
        }
        return -K*(Ksin*sine_part - Kcos*cosine_part)*E;
    }

    /// Compute CTF pure at (U,V). Continuous frequencies
    inline double getValuePureNoDampingNoPrecomputedAt(double X, double Y) const
    {
        double u2 = X * X + Y * Y;
        double u = sqrt(u2);
        //if(u2 > freq_max) return 0;
        double u4 = u2 * u2;
        double deltaf = getDeltafNoPrecomputed(X, Y);
        double VPP;
        double check_VPP = round(VPP_radius*1000);
		if(check_VPP != 0)
			VPP = phase_shift*(1-exp(-u2/(2*pow(VPP_radius,2.0))));
		else
			VPP = 0;
        double argument = VPP + K1 * deltaf * u2 + K2 * u4;
        double sine_part, cosine_part;
        sincos(argument,&sine_part, &cosine_part); // OK
        return -(Ksin*sine_part - Kcos*cosine_part);
    }

    /// Deltaf at a given direction
    double getDeltafNoPrecomputed(double X, double Y) const
    {
        if (fabs(X) < XMIPP_EQUAL_ACCURACY &&
            fabs(Y) < XMIPP_EQUAL_ACCURACY)
            return 0;
        double ellipsoid_ang = atan2(Y, X) - rad_azimuth;
        double cos_ellipsoid_ang_2 = cos(2*ellipsoid_ang);
        return(defocus_average + defocus_deviation*cos_ellipsoid_ang_2);
    }

    /// Compute noise at (X,Y). Continuous frequencies, notice it is squared
    //#define DEBUG
    inline double getValueNoiseAt(bool show = false) const
    {
        double ellipsoid_ang = precomputed.ang - rad_gaussian;
        double ellipsoid_ang2 = precomputed.ang - rad_gaussian2;
        double ellipsoid_sqrt_ang = precomputed.ang - rad_sqrt;
        double cos_sqrt_ang = cos(ellipsoid_sqrt_ang);
        double cos_sqrt_ang_2 = cos_sqrt_ang*cos_sqrt_ang;
        double sin_sqrt_ang_2 = 1.0-cos_sqrt_ang_2;
        double sq = sqrt(sqU * sqU * cos_sqrt_ang_2 + sqV * sqV * sin_sqrt_ang_2);

        double cos_ang = cos(ellipsoid_ang);
        double cos_ang_2 = cos_ang*cos_ang;
        double sin_ang_2 = 1.0-cos_ang_2;
        double c = sqrt(cU * cU * cos_ang_2 + cV * cV * sin_ang_2);
        double sigma = sqrt(sigmaU * sigmaU * cos_ang_2 + sigmaV * sigmaV * sin_ang_2);

        double cos_ang2 = cos(ellipsoid_ang2);
        double cos_ang2_2 = cos_ang2*cos_ang2;
        double sin_ang2_2 = 1.0-cos_ang2_2;
        double c2 = sqrt(cU2 * cU2 * cos_ang2_2 + cV2 * cV2 * sin_ang2_2);
        double sigma2 = sqrt(sigmaU2 * sigmaU2 * cos_ang2_2 + sigmaV2 * sigmaV2 * sin_ang2_2);

        if(show)
        {
            std::cout << "   ellipsoid_ang=" << RAD2DEG(ellipsoid_ang) << std::endl
            << "   ellipsoid_sqrt_ang=" << RAD2DEG(ellipsoid_sqrt_ang) << std::endl
            << "   sq=" << sq << "\n"
            << "   c=" << c << "\n"
            << "   sigma=" << sigma << "\n"
            << "   ellipsoid_ang2=" << RAD2DEG(ellipsoid_ang2) << std::endl
            << "   cU2=" << cU2 << " cV2=" << cV2 << "\n"
            << "   cos_ang2=" << cos_ang2 << " sin_ang2_2=" << sin_ang2_2 << "\n"
            << "   c2=" << c2 << "\n"
            << "   sigmaU2=" << sigma2 << "\n"
			<< "   gaussian_K=" << gaussian_K << "\n"
			<< "   sqrt_K=" << sqrt_K << "\n"
			<< "   base_line=" << base_line << "\n";
            std::cout << "   u=" << precomputed.u << "u_sqrt=" << precomputed.u_sqrt <<") CTFnoise="
            << base_line +
            gaussian_K*exp(-sigma*(precomputed.u - c)*(precomputed.u - c)) +
            sqrt_K*exp(-sq*sqrt(precomputed.u)) -
            gaussian_K2*exp(-sigma2*(precomputed.u - c2)*(precomputed.u - c2)) << std::endl;
        }

        double aux=precomputed.u - c;
        double aux2=precomputed.u - c2;
        return base_line +
               gaussian_K*exp(-sigma*aux*aux) +
               sqrt_K*exp(-sq*precomputed.u_sqrt) -
               gaussian_K2*exp(-sigma2*aux2*aux2)+
			   bgR1*precomputed.u+bgR2*precomputed.u2+bgR3*precomputed.u3;
    }

    /** Returns the continuous frequency of the zero, maximum or minimum number n in the direction u.
        u must be a unit vector, n=1,2,... Returns (-1,-1) if it is not found
        'iwhat' can be 0 (zero), 1(max), or -1 (min) */
    void lookFor(int n, const Matrix1D<double> &u, Matrix1D<double> &freq, int iwhat=0);

    /// Apply CTF to an image
    void applyCTF(MultidimArray < std::complex<double> > &FFTI, const MultidimArray<double> &I, double Ts, bool absPhase=false);

    /// Apply CTF to an image
    void applyCTF(MultidimArray <double> &I, double Ts, bool absPhase=false);

    /** Generate CTF image.
        The sample image is used only to take its dimensions. */
    template <class T1, class T2>
    void generateCTF(const MultidimArray<T1> &sample_image, MultidimArray <T2> &CTF, double Ts=-1)
    {
        if ( ZSIZE(sample_image) > 1 )
            REPORT_ERROR(ERR_MULTIDIM_DIM,"ERROR: Generate_CTF only works with 2D sample images, not 3D.");
        generateCTF(YSIZE(sample_image), XSIZE(sample_image), CTF, Ts);
        STARTINGX(CTF) = STARTINGX(sample_image);
        STARTINGY(CTF) = STARTINGY(sample_image);
    }

    /** Get profiles along a given direction.
     * It returns the CTF profiles (BGNOISE, ENVELOPE, PSD, CTF) along the
     * direction defined by angle. Profiles
     * have nsamples from 0 to fmax (1/A).
     */
    void getProfile(double angle, double fmax, int nsamples, MultidimArray<double> &profiles);

    /** Get radial average profiles.
     * It returns the radially average CTF profiles (BGNOISE, ENVELOPE, PSD, CTF). Profiles
     * have nsamples from 0 to fmax (1/A).
     */
    void getAverageProfile(double fmax, int nsamples, MultidimArray<double> &profiles);

//#define DEBUG
   /// Generate CTF image.
    template <class T>
    void generateCTF(int Ydim, int Xdim, MultidimArray < T > &CTF, double Ts=-1)
    {
        CTF.resizeNoCopy(Ydim, Xdim);
        if (Ts<0)
        	Ts=Tm;
		#ifdef DEBUG
			std::cout << "CTF:\n" << *this << std::endl;
		#endif

        double iTs=1.0/Ts;
        for (int i=0; i<Ydim; ++i)
        {
        	double wy;
        	FFT_IDX2DIGFREQ(i, YSIZE(CTF), wy);
            double fy=wy*iTs;
        	for (int j=0; j<Xdim; ++j)
        	{
            	double wx;
            	FFT_IDX2DIGFREQ(j, XSIZE(CTF), wx);
                double fx=wx*iTs;
				precomputeValues(fx, fy);
				A2D_ELEM(CTF, i, j) = (T) getValueAt();
				#ifdef DEBUG
						if (i == 0)
							std::cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
							<< " " << CTF(i, j) << std::endl;
				#endif
        	}
        }
    }
    #undef DEBUG

    template <class T>
    void generateCTFWithoutDamping(int Ydim, int Xdim, MultidimArray < T > &CTF, double Ts=-1)
    {
        CTF.resizeNoCopy(Ydim, Xdim);
        if (Ts<0)
        	Ts=Tm;
		#ifdef DEBUG
			std::cout << "CTF:\n" << *this << std::endl;
		#endif

        double iTs=1.0/Ts;
        for (int i=0; i<Ydim; ++i)
        {
        	double wy;
        	FFT_IDX2DIGFREQ(i, YSIZE(CTF), wy);
            double fy=wy*iTs;
        	for (int j=0; j<Xdim; ++j)
        	{
            	double wx;
            	FFT_IDX2DIGFREQ(j, XSIZE(CTF), wx);
                double fx=wx*iTs;
				precomputeValues(fx, fy);
				A2D_ELEM(CTF, i, j) = (T) getValuePureWithoutDampingAt();
				#ifdef DEBUG
						if (i == 0)
							std::cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
							<< " " << CTF(i, j) << std::endl;
				#endif
        	}
        }
    }
    #undef DEBUG

    /** Check physical meaning.
        true if the CTF parameters have physical meaning.
        Call this function after produstd::cing side information */
    bool hasPhysicalMeaning();

    /** Force physical meaning.*/
    void forcePhysicalMeaning();
};

#endif
