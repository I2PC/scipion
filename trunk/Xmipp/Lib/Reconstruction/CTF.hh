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

#include <XmippData/xmippImages.hh>
#include <XmippInterface/xmippVTK.hh>

/**@name CTF Correction */
//@{
/** CTF class. */
class XmippCTF {
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
   // square roots of defoci
   double sqrt_DeltafU;
   double sqrt_DeltafV;
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
   /// Chromatic aberration (in milimeters). Typical value 4
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
   /// Gain for the gaussian term
   double gaussian_K;
   /// Global base_line
   double base_line;
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

   /** Empty constructor. */
   XmippCTF() {clear();}

   /** Read from file.
       An exception is thrown if the file cannot be open.

       If no K or sqrt_K are given then it is assumed that the user
       does not want to activate that part and the noise or the CTF
       are removed from the model unless the disable_if_not_K is set
       to FALSE*/
   void read(const FileName &fn, bool disable_if_not_K=TRUE) _THROW;
   
   /// Usage
   void Usage();

   /// Show
   friend ostream & operator << (ostream &out, const XmippCTF &ctf);

   /// Clear.
   void clear();

   /// Clear noise
   void clear_noise();
   
   /// Clear pure CTF
   void clear_pure_ctf();
   
   /// Produce Side information
   void Produce_Side_Info();

   /** Produce Main information.
       When producing the side information, some of the variables given
       by the user are changed in scale. This function returns them to
       the original scale so that they can be printed in the same way
       as the user is used to them */
   void Produce_Main_Info();

   /// Deltaf at a given direction
   double Deltaf(double X, double Y) const {
      if (ABS(X)<XMIPP_EQUAL_ACCURACY &&
          ABS(Y)<XMIPP_EQUAL_ACCURACY) return 0;
      double ellipsoid_ang=atan2(Y,X)-rad_azimuth;
      double DeltafUp=DeltafU*cos(ellipsoid_ang);
      double DeltafVp=DeltafV*sin(ellipsoid_ang);
      return SGN(DeltafU)*sqrt(DeltafUp*DeltafUp+DeltafVp*DeltafVp);
   }

   double CTF_at(double X, double Y, bool show=FALSE) const {
      double pure_CTF;
      if (enable_CTF)      pure_CTF=CTFpure_at(X,Y,show);
      else                 pure_CTF=0;
      if (enable_CTFnoise) return sqrt(pure_CTF*pure_CTF+CTFnoise_at(X,Y,show));
      else                 return pure_CTF;
   }

   /// Compute CTF at (U,V). Continuous frequencies
   double CTFpure_without_damping_at(double X, double Y, bool show=FALSE) const {
      double u2=X*X+Y*Y;
      double u=sqrt(u2);
      double u4=u2*u2;
      // if (u2>=ua2) return 0;
      double deltaf=Deltaf(X,Y);
      double argument=K1*deltaf*u2+K2*u4;
      double sine_part=sin(argument); // OK
      double cosine_part=cos(argument);
      if (show) {
         cout << "   Deltaf=" << deltaf << endl;
         cout << "   u,u2,u4=" << u << " " << u2 << " " << u4 << endl;
         cout << "   K1,K2,sin=" << K1 << " " << K2 << " "
              << sine_part << endl;
	 cout << "   Q0=" << Q0 << endl;
         cout << "   (X,Y)=(" << X << "," << Y << ") CTF without damping="
              << -(sine_part+Q0*cosine_part) << endl;
      }
      return -(sine_part+Q0*cosine_part);
   }

   /// Compute CTF at (U,V). Continuous frequencies
   double CTFpure_at(double X, double Y, bool show=FALSE) const {
      double u2=X*X+Y*Y;
      double u=sqrt(u2);
      double u4=u2*u2;
      // if (u2>=ua2) return 0;
      double deltaf=Deltaf(X,Y);
      double argument=K1*deltaf*u2+K2*u4;
      double sine_part=sin(argument); // OK
      double cosine_part=cos(argument);
      double Eespr=exp(-K3*u4); // OK
      //CO: double Eispr=exp(-K4*u4); // OK
      double EdeltaF=bessj0(K5*u2); // OK
      double EdeltaR=SINC(u*DeltaR); // OK
      double Ealpha=exp(-K6*(K7*u2*u+deltaf*u)*(K7*u2*u+deltaf*u)); // OK
      // CO: double E=Eespr*Eispr*EdeltaF*EdeltaR*Ealpha;
      double E=Eespr*EdeltaF*EdeltaR*Ealpha;
      if (show) {
         cout << "   Deltaf=" << deltaf << endl;
         cout << "   u,u2,u4=" << u << " " << u2 << " " << u4 << endl;
         cout << "   K1,K2,sin=" << K1 << " " << K2 << " "
              << sine_part << endl;
         cout << "   K3,Eespr=" << K3 << " " << Eespr << endl;
         cout << "   K4,Eispr=" << K4 << " " << /*Eispr <<*/ endl;
         cout << "   K5,EdeltaF=" << K5 << " " << EdeltaF << endl;
         cout << "   EdeltaR=" << EdeltaR << endl;
         cout << "   K6,K7,Ealpha=" << K6 << " " << K7 << " " << Ealpha
              << endl;
         cout << "   Total atenuation(E)= " << E << endl;
	 cout << "   K,Q0,base_line=" << K << "," << Q0 << "," << base_line << endl;
         cout << "   (X,Y)=(" << X << "," << Y << ") CTF="
              << -K*(sine_part+Q0*cosine_part)*E+base_line << endl;
      }
      return -K*(sine_part+Q0*cosine_part)*E;
   }

   /** Compute CTF gaussian at (U,V). Continuous frequencies.
       This is the pure gaussian, without any scaling factor. Its range
       is between 0 and 1.*/
   double CTFgaussian_at(double X, double Y, bool show=false) const {
      double ellipsoid_ang;
      if (gaussian_K==0.0) return 0.0;
      if (ABS(X)<XMIPP_EQUAL_ACCURACY &&
          ABS(Y)<XMIPP_EQUAL_ACCURACY) ellipsoid_ang=0;
      else ellipsoid_ang=atan2(Y,X)-rad_gaussian;
      double cUp=cU*cos(ellipsoid_ang);
      double cVp=cV*sin(ellipsoid_ang);
      double c=sqrt(cUp*cUp+cVp*cVp);
      double sigmaUp=sigmaU*cos(ellipsoid_ang);
      double sigmaVp=sigmaV*sin(ellipsoid_ang);
      double sigma=sqrt(sigmaUp*sigmaUp+sigmaVp*sigmaVp);
      double w=sqrt(X*X+Y*Y);
      if (show) {
         cout << "   ellipsoid_ang=" << RAD2DEG(ellipsoid_ang) << endl
              << "   cUp, cVp=" << cUp << "," << cVp << " (" << c << ")\n"
              << "   sigmaUp, sigmaVp=" << sigmaUp << "," << sigmaVp << " (" << sigma << ")\n";
         cout << "   (X,Y)=(" << X << "," << Y << ") (" << w << ") Gaussian="
              << exp(-sigma*(w-c)*(w-c)) << endl;
      }
      return exp(-sigma*(w-c)*(w-c));	
   }

   /// Compute noise at (X,Y). Continuous frequencies, notice it is squared
   //#define DEBUG
   double CTFnoise_at(double X, double Y, bool show=false) const {
      double ellipsoid_ang;
      if (ABS(X)<XMIPP_EQUAL_ACCURACY &&
          ABS(Y)<XMIPP_EQUAL_ACCURACY) ellipsoid_ang=0;
      else ellipsoid_ang=atan2(Y,X)-rad_gaussian;
      double sqUp=sqU*cos(ellipsoid_ang);
      double sqVp=sqV*sin(ellipsoid_ang);
      double sq=sqrt(sqUp*sqUp+sqVp*sqVp);
      double cUp=cU*cos(ellipsoid_ang);
      double cVp=cV*sin(ellipsoid_ang);
      double c=sqrt(cUp*cUp+cVp*cVp);
      double sigmaUp=sigmaU*cos(ellipsoid_ang);
      double sigmaVp=sigmaV*sin(ellipsoid_ang);
      double sigma=sqrt(sigmaUp*sigmaUp+sigmaVp*sigmaVp);
      double w=sqrt(X*X+Y*Y);
      if (show) {
         cout << "   ellipsoid_ang=" << RAD2DEG(ellipsoid_ang) << endl
              << "   sqUp, sqVp=" << sqUp << "," << sqVp << " (" << sq << ")\n"
              << "   cUp, cVp=" << cUp << "," << cVp << " (" << c << ")\n"
              << "   sigmaUp, sigmaVp=" << sigmaUp << "," << sigmaVp << " (" << sigma << ")\n";
         cout << "   (X,Y)=(" << X << "," << Y << ") (" << w << ") CTFnoise="
              << base_line+
             gaussian_K*exp(-sigma*(w-c)*(w-c))+sqrt_K*exp(-sq*sqrt(w)) << endl;
      }
      return base_line+
             gaussian_K*exp(-sigma*(w-c)*(w-c))+sqrt_K*exp(-sq*sqrt(w));	
   }

   /** Returns the frequency of the zero number n in the direction u.
       u must be a unit vector, n=1,2,... Returns (-1,-1) if it is not found */
   void zero(int n, const matrix1D<double> &u, matrix1D<double> &freq);

   /// Apply CTF to an image
   void Apply_CTF(vtkImageData * &FFTI) const;

   /// Generate CTF image
   void Generate_CTF(vtkImageData * FFTI, vtkImageData *&CTF) const;

   /** Check physical meaning.
       TRUE if the CTF parameters have physical meaning.
       Call this function after producing side information */
   bool physical_meaning();
};
//@}
#endif
