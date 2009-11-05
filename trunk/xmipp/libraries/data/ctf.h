/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#ifndef _CTF_HH
#define _CTF_HH

#include "image.h"
#include "selfile.h"
#include <map>

/**@defgroup CTFSupport CTF support classes
   @ingroup DataLibrary */
//@{
/** CTF class.
    Here goes how to compute the radial average of a parametric CTF:

    @code
      #include <Reconstruction/CTF.hh>

      int main() {
        // Read the parametric CTF
  XmippCTF CTF;
  CTF.enable_CTF=true;
  CTF.enable_CTFnoise=false;
  CTF.read("ctf_crio.param");
  CTF.Produce_Side_Info();

        // Compute the CTF radial average
  double sampling_rate=3.5; // in Angstroms/pixel
  double fmax=1.0/(2.0*sampling_rate);
  for (double f=0; f<fmax; f+=fmax/100.0) {
     double CTF_at_f=0;
     double N=0;
     for(double ang=0; ang<2*PI; ang+=0.01) {
               CTF_at_f+=CTF.CTF_at(f*cos(ang),f*sin(ang));
        N++;
     }
     std::cout << f << " " << CTF_at_f/N << std::endl;
  }
  return 0;
      }
   @endcode

   Here is another sample program for generating the CTF, the noise
   background, the envelope and the power spectrum density

   @code
      #include <Reconstruction/CTF.hh>
      #include <XmippData/xmippFFT.hh>

      int main(int argc, char **argv) {
         FileName fn_ctf, fn_root;
         int Xdim;
         try {
            fn_ctf=getParameter(argc,argv,"-i");
            fn_root=getParameter(argc,argv,"-o");
            Xdim=textToInteger(getParameter(argc,argv,"-xdim"));
         } catch (Xmipp_error XE) {
            std::cerr << XE << std::endl
                 << "Usage: produce_imgs \n"
                 << "         -i <CTF descr file>\n"
                 << "         -o <filename root>\n"
                 << "         -xdim <xdim>\n"
            ;
            return 1;
         }

         try {
            // Read CTF model
            XmippCTF CTF;
            CTF.enable_CTF=true;
            CTF.enable_CTFnoise=true;
            CTF.read(fn_ctf);
            CTF.Produce_Side_Info();

            // Produce CTF, background and envelope
            ImageXmipp Ictf, Ibg, Ienv, Ipsd;
            Ictf().resize(Xdim,Xdim);
            Ipsd()=Ienv()=Ibg()=Ictf();

            Matrix1D<int>    idx(2);  // Indexes for Fourier plane
            Matrix1D<double> freq(2); // Frequencies for Fourier plane
            FOR_ALL_ELEMENTS_IN_MATRIX2D(Ictf()) {
               XX(idx)=j; YY(idx)=i;
               FFT_idx2digfreq(Ictf(), idx, freq);
               digfreq2contfreq(freq, freq, CTF.Tm);

               // Background
               Ibg(i,j)=CTF.CTFnoise_at(XX(freq),YY(freq));

               // Envelope
               double E=CTF.CTFdamping_at(XX(freq),YY(freq));
               Ienv(i,j)=Ibg(i,j)+E*E;

               // CTF
               Ictf(i,j)=CTF.CTFpure_at(XX(freq),YY(freq));

               // Power spectrum density
               Ipsd(i,j)=Ibg(i,j)+Ictf(i,j)*Ictf(i,j);
            }

            CenterFFT(Ibg() , true); Ibg .write(fn_root+"_bg.xmp");
            CenterFFT(Ienv(), true); Ienv.write(fn_root+"_env.xmp");
            CenterFFT(Ictf(), true); Ictf.write(fn_root+"_ctf.xmp");
            CenterFFT(Ipsd(), true); Ipsd.write(fn_root+"_psd.xmp");

         } catch (Xmipp_error XE) {
            std::cout << XE << std::endl;
         }
         return 0;
      }
   @endcode

*/
class XmippCTF
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
    // Azimuthal angle in radians
    double rad_azimuth;
    // Gaussian angle in radians
    double rad_gaussian;
    // Second Gaussian angle in radians
    double rad_gaussian2;
    // Sqrt angle in radians
    double rad_sqrt;
    /** Standard error of defocus Gaussian function due to chromatic aberration.
        in Amstrong */
    double D;
public:
    /// Global gain. By default, 1
    double K;
    /// Sampling rate (A/pixel)
    double Tm;
    /// Accelerating Voltage (in KiloVolts)
    double kV;
    /// Defocus in U (in Angstroms). Negative values are underfocused
    double DeltafU;
    /// Defocus in V (in Angstroms). Negative values are underfocused
    double DeltafV;
    /// Azimuthal angle (between X and U) in degrees
    double azimuthal_angle;
    // Radius of the aperture (in micras)
    // double aperture;
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

    /// Enable CTFnoise part
    bool enable_CTFnoise;
    /// Enable CTF part
    bool enable_CTF;
    /// Global base_line
    double base_line;
    /// Gain for the gaussian term
    double gaussian_K;
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
    double sqrt_K;
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

    /** Empty constructor. */
    XmippCTF()
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

    /** Write to file.
        An exception is thrown if the file cannot be open.*/
    void write(const FileName &fn);

    /// Usage
    void Usage();

    /// Show
    friend std::ostream & operator << (std::ostream &out, const XmippCTF &ctf);

    /// Clear.
    void clear();

    /// Clear noise
    void clear_noise();

    /// Clear pure CTF
    void clear_pure_ctf();

    /// Produce Side information
    void Produce_Side_Info();

    /// Deltaf at a given direction
    double Deltaf(double X, double Y) const
    {
        if (ABS(X) < XMIPP_EQUAL_ACCURACY &&
            ABS(Y) < XMIPP_EQUAL_ACCURACY) return 0;
        double ellipsoid_ang = atan2(Y, X) - rad_azimuth;
        double DeltafUp = DeltafU * cos(ellipsoid_ang);
        double DeltafVp = DeltafV * sin(ellipsoid_ang);
        return SGN(DeltafU)*sqrt(DeltafUp*DeltafUp + DeltafVp*DeltafVp);
    }

    double CTF_at(double X, double Y, bool show = false) const
    {
        double pure_CTF;
        if (enable_CTF)      pure_CTF = CTFpure_at(X, Y, show);
        else                 pure_CTF = 0;
        if (enable_CTFnoise) return sqrt(pure_CTF*pure_CTF + CTFnoise_at(X, Y, show));
        else                 return pure_CTF;
    }

    /// Compute CTF at (U,V). Continuous frequencies
    double CTFpure_without_damping_at(double X, double Y, bool show = false) const
    {
        double u2 = X * X + Y * Y;
        double u = sqrt(u2);
        double u4 = u2 * u2;
        // if (u2>=ua2) return 0;
        double deltaf = Deltaf(X, Y);
        double argument = K1 * deltaf * u2 + K2 * u4;
        double sine_part = sin(argument); // OK
        double cosine_part = cos(argument);
        if (show)
        {
            std::cout << "   Deltaf=" << deltaf << std::endl;
            std::cout << "   u,u2,u4=" << u << " " << u2 << " " << u4 << std::endl;
            std::cout << "   K1,K2,sin=" << K1 << " " << K2 << " "
            << sine_part << std::endl;
            std::cout << "   Q0=" << Q0 << std::endl;
            std::cout << "   (X,Y)=(" << X << "," << Y << ") CTF without damping="
            << -(sine_part + Q0*cosine_part) << std::endl;
        }
        return -(sine_part + Q0*cosine_part);
    }

    /// Compute CTF at (U,V). Continuous frequencies
    inline double CTFdamping_at(double X, double Y, bool show = false) const
    {
        double u2 = X * X + Y * Y;
        double u = sqrt(u2);
        double u4 = u2 * u2;
        double deltaf = Deltaf(X, Y);
        double Eespr = exp(-K3 * u4); // OK
        double EdeltaF = bessj0(K5 * u2); // OK
        double EdeltaR = SINC(u * DeltaR); // OK
        double Ealpha = exp(-K6 * (K7 * u2 * u + deltaf * u) * (K7 * u2 * u + deltaf * u)); // OK
        double E = Eespr * EdeltaF * EdeltaR * Ealpha;
        if (show)
        {
            std::cout << "   Deltaf=" << deltaf << std::endl;
            std::cout << "   u,u2,u4=" << u << " " << u2 << " " << u4 << std::endl;
            std::cout << "   K3,Eespr=" << K3 << " " << Eespr << std::endl;
            std::cout << "   K4,Eispr=" << K4 << " " << /*Eispr <<*/ std::endl;
            std::cout << "   K5,EdeltaF=" << K5 << " " << EdeltaF << std::endl;
            std::cout << "   EdeltaR=" << EdeltaR << std::endl;
            std::cout << "   K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
            << std::endl;
            std::cout << "   Total atenuation(E)= " << E << std::endl;
            std::cout << "   (X,Y)=(" << X << "," << Y << ") CTFdamp="
            << E << std::endl;
        }
        return -K*E;
    }

    /// Compute CTF at (U,V). Continuous frequencies
    inline double CTFpure_at(double X, double Y, bool show = false) const
    {
        double u2 = X * X + Y * Y;
        double u = sqrt(u2);
        double u4 = u2 * u2;
        // if (u2>=ua2) return 0;
        double deltaf = Deltaf(X, Y);
        double argument = K1 * deltaf * u2 + K2 * u4;
        double sine_part = sin(argument); // OK
        double cosine_part = cos(argument);
        double Eespr = exp(-K3 * u4); // OK
        //CO: double Eispr=exp(-K4*u4); // OK
        double EdeltaF = bessj0(K5 * u2); // OK
        double EdeltaR = SINC(u * DeltaR); // OK
        double Ealpha = exp(-K6 * (K7 * u2 * u + deltaf * u) * (K7 * u2 * u + deltaf * u)); // OK
        // CO: double E=Eespr*Eispr*EdeltaF*EdeltaR*Ealpha;
        double E = Eespr * EdeltaF * EdeltaR * Ealpha;
        if (show)
        {
            std::cout << "   Deltaf=" << deltaf << std::endl;
            std::cout << "   u,u2,u4=" << u << " " << u2 << " " << u4 << std::endl;
            std::cout << "   K1,K2,sin=" << K1 << " " << K2 << " "
            << sine_part << std::endl;
            std::cout << "   K3,Eespr=" << K3 << " " << Eespr << std::endl;
            std::cout << "   K4,Eispr=" << K4 << " " << /*Eispr <<*/ std::endl;
            std::cout << "   K5,EdeltaF=" << K5 << " " << EdeltaF << std::endl;
            std::cout << "   EdeltaR=" << EdeltaR << std::endl;
            std::cout << "   K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
            << std::endl;
            std::cout << "   Total atenuation(E)= " << E << std::endl;
            std::cout << "   K,Q0,base_line=" << K << "," << Q0 << "," << base_line << std::endl;
            std::cout << "   (X,Y)=(" << X << "," << Y << ") CTF="
            << -K*(sine_part + Q0*cosine_part)*E + base_line << std::endl;
        }
        return -K*(sine_part + Q0*cosine_part)*E;
    }

    /// Compute noise at (X,Y). Continuous frequencies, notice it is squared
    //#define DEBUG
    inline double CTFnoise_at(double X, double Y, bool show = false) const
    {
        double ellipsoid_ang, ellipsoid_ang2, ellipsoid_sqrt_ang;
        if (ABS(X) < XMIPP_EQUAL_ACCURACY &&
            ABS(Y) < XMIPP_EQUAL_ACCURACY)
            ellipsoid_sqrt_ang = ellipsoid_ang = ellipsoid_ang2 = 0;
        else
        {
            double yx_ang = atan2(Y, X);
            ellipsoid_ang = yx_ang - rad_gaussian;
            ellipsoid_ang2 = yx_ang - rad_gaussian2;
            ellipsoid_sqrt_ang = yx_ang - rad_sqrt;
        }
        double cos_sqrt_ang = cos(ellipsoid_sqrt_ang);
        double sin_sqrt_ang = sin(ellipsoid_sqrt_ang);
        double sqUp = sqU * cos_sqrt_ang;
        double sqVp = sqV * sin_sqrt_ang;
        double sq = sqrt(sqUp * sqUp + sqVp * sqVp);

        double cos_ang = cos(ellipsoid_ang);
        double sin_ang = sin(ellipsoid_ang);
        double cUp = cU * cos_ang;
        double cVp = cV * sin_ang;
        double c = sqrt(cUp * cUp + cVp * cVp);
        double sigmaUp = sigmaU * cos_ang;
        double sigmaVp = sigmaV * sin_ang;
        double sigma = sqrt(sigmaUp * sigmaUp + sigmaVp * sigmaVp);

        double cos_ang2 = cos(ellipsoid_ang2);
        double sin_ang2 = sin(ellipsoid_ang2);
        double cUp2 = cU2 * cos_ang2;
        double cVp2 = cV2 * sin_ang2;
        double c2 = sqrt(cUp2 * cUp2 + cVp2 * cVp2);
        double sigmaUp2 = sigmaU2 * cos_ang2;
        double sigmaVp2 = sigmaV2 * sin_ang2;
        double sigma2 = sqrt(sigmaUp2 * sigmaUp2 + sigmaVp2 * sigmaVp2);

        double w = sqrt(X * X + Y * Y);
        if (show)
        {
            std::cout << "   ellipsoid_ang=" << RAD2DEG(ellipsoid_ang) << std::endl
            << "   ellipsoid_sqrt_ang=" << RAD2DEG(ellipsoid_sqrt_ang) << std::endl
            << "   sqUp, sqVp=" << sqUp << "," << sqVp << " (" << sq << ")\n"
            << "   cUp, cVp=" << cUp << "," << cVp << " (" << c << ")\n"
            << "   sigmaUp, sigmaVp=" << sigmaUp << "," << sigmaVp << " (" << sigma << ")\n"
            << "   ellipsoid_ang2=" << RAD2DEG(ellipsoid_ang2) << std::endl
            << "   cUp2, cVp2=" << cUp2 << "," << cVp2 << " (" << c2 << ")\n"
            << "   sigmaUp2, sigmaVp2=" << sigmaUp2 << "," << sigmaVp2 << " (" << sigma2 << ")\n";
            std::cout << "   (X,Y)=(" << X << "," << Y << ") (" << w << ") CTFnoise="
            << base_line +
            gaussian_K*exp(-sigma*(w - c)*(w - c)) + sqrt_K*exp(-sq*sqrt(w)) -
            gaussian_K2*exp(-sigma2*(w - c2)*(w - c2)) << std::endl;
        }
        return base_line +
               gaussian_K*exp(-sigma*(w - c)*(w - c)) + sqrt_K*exp(-sq*sqrt(w)) -
               gaussian_K2*exp(-sigma2*(w - c2)*(w - c2));
    }

    /** Returns the continuous frequency of the zero number n in the direction u.
        u must be a unit vector, n=1,2,... Returns (-1,-1) if it is not found */
    void zero(int n, const Matrix1D<double> &u, Matrix1D<double> &freq) const;

    /// Apply CTF to an image
    void Apply_CTF(Matrix2D < std::complex<double> > &FFTI) const;

    /** Generate CTF image.
        The sample image is used only to take its dimensions. */
    template <class T>
    void Generate_CTF(const Matrix2D<T> &sample_image,
                      Matrix2D < std::complex<double> > &CTF) const
    {
        Generate_CTF(YSIZE(sample_image), XSIZE(sample_image), CTF);
        STARTINGX(CTF) = STARTINGX(sample_image);
        STARTINGY(CTF) = STARTINGY(sample_image);
    }

    /// Generate CTF image.
    void Generate_CTF(int Ydim, int Xdim,
                      Matrix2D < std::complex<double> > &CTF) const;

    /** Check physical meaning.
        true if the CTF parameters have physical meaning.
        Call this function after produstd::cing side information */
    bool physical_meaning();

    /** Force physical meaning.*/
    void force_physical_meaning();
};

/** CTF Data file.
    CTF Data files are text files with two columns: the first column
    is the name of a projection image, the second column is the
    name of the corresponding ctfparam file.
    
    Example of use:
    @code
        CTFDat ctfdat;
	ctfdat.read("ctfdat.txt");
	ctfdat.goFirstLine();
	std::cerr << "Correcting CTF phase ...\n";
	int istep = CEIL((double)ctfdat.lineNo() / 60.0);
	init_progress_bar(ctfdat.lineNo());
	int i = 0;
	while (!ctfdat.eof())
	{
            FileName fnProjection, fnCTF;
	    ctfdat.getCurrentLine(fnProjection,fnCTF);
            ...
            if (i++ % istep == 0) progress_bar(i);
	    ctfdat.nextLine();
	}
	progress_bar(ctfdat.lineNo());
    @endcode
*/
class CTFDat {
public:
    /// List with the projection files
    std::vector< FileName > fnProjectionList;

    /// List with the ctfparam files
    std::vector< FileName > fnCTFList;

    /// Iterator for going over the file
    int current;

    /** Get the CTF file for a given key.
       searchOK is true if the the key is found. If the key is not found,
       searchOK is false and the returned value is empty. */
    const FileName & getCTF(const FileName &fnProjection,
      bool& searchOK) const;

    /// Set the CTF file for a given key
    void setCTF(const FileName &fnProjection, const FileName &fnCtf);

    /// Append a CTF
    void append(const FileName &fnProjection, const FileName &fnCtf);

    /// Read the CTFDat from a file
    void read(const FileName &fnCtfdat); 

    /// Write the CTFDat to a file
    void write(const FileName &fnCtfdat) const;

    /// Move iterator to the file beginning
    void goFirstLine();

    /// Move iterator to the next image
    void nextLine();

    /// Check if at the endo of file
    bool eof() const;

    /// Returns the number of lines in the file
    int lineNo() const;

    /// Get current line
    void getCurrentLine(FileName& fnProjection, FileName& fnCtf);
    
    /// Creates a CTFDat with a selfile for which there is a single CTF
    void createFromSelfileAndSingleCTF(SelFile &SF,
       const FileName &fnCtf);
};
//@}
#endif
