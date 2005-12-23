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

#include "../CTF.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFFT.hh>

/* Read -------------------------------------------------------------------- */
void XmippCTF::read(const FileName &fn, bool disable_if_not_K) {
   FILE *fh_param;
   if ((fh_param = fopen(fn.c_str(), "r")) == NULL)
      REPORT_ERROR(1,
         (string)"XmippCTF::read: There is a problem "
         "opening the file "+fn);

   try {
      Tm=AtoF(get_param(fh_param,"sampling_rate",0,"1"));
      if (enable_CTF) {
         DeltafU=AtoF(get_param(fh_param,"defocusU",0,"0"));
         if (check_param(fh_param,"defocusV"))
              DeltafV=AtoF(get_param(fh_param,"defocusV",0));
         else DeltafV=DeltafU;
         azimuthal_angle=AtoF(get_param(fh_param,"azimuthal_angle",0,"0"));
         kV=AtoF(get_param(fh_param,"voltage",0,"100"));
         Cs=AtoF(get_param(fh_param,"spherical_aberration",0,"0"));
         Ca=AtoF(get_param(fh_param,"chromatic_aberration",0,"0"));
         espr=AtoF(get_param(fh_param,"energy_loss",0,"0"));
         ispr=AtoF(get_param(fh_param,"lens_stability",0,"0"));
         alpha=AtoF(get_param(fh_param,"convergence_cone",0,"0"));
         DeltaF=AtoF(get_param(fh_param,"longitudinal_displace",0,"0"));
         DeltaR=AtoF(get_param(fh_param,"transversal_displace",0,"0"));
         Q0=AtoF(get_param(fh_param,"Q0",0,"0"));
         K=AtoF(get_param(fh_param,"K",0,"0"));
         if (K==0 && disable_if_not_K) enable_CTF=false;
      }
      
      if (enable_CTFnoise) {
	 gaussian_K    =AtoF(get_param(fh_param,"gaussian_K",0,"0"));
	 base_line     =AtoF(get_param(fh_param,"base_line",0,"0"));
	 sigmaU        =AtoF(get_param(fh_param,"sigmaU",0,"0"));
	 if (check_param(fh_param,"sigmaV"))
            sigmaV     =AtoF(get_param(fh_param,"sigmaV",0));
	 else sigmaV   =sigmaU;
	 cU            =AtoF(get_param(fh_param,"cU",0,"0"));
	 if (check_param(fh_param,"cV"))
            cV         =AtoF(get_param(fh_param,"cV",0));
	 else cV       =cU;
	 gaussian_angle=AtoF(get_param(fh_param,"gaussian_angle",0,"0"));
	 sqU           =AtoF(get_param(fh_param,"sqU",0,"0"));
	 if (check_param(fh_param,"sqV"))
            sqV        =AtoF(get_param(fh_param,"sqV",0));
	 else sqV      =sqU;
	 sqrt_angle=AtoF(get_param(fh_param,"sqrt_angle",0,"0"));
	 sqrt_K        =AtoF(get_param(fh_param,"sqrt_K",0,"0"));
         if (gaussian_K==0 && sqrt_K==0 && base_line==0  && disable_if_not_K)
            enable_CTFnoise=false;
      }
      
   } catch (Xmipp_error XE) {
      cout << XE << endl;
      REPORT_ERROR(1,(string)"There is an error reading "+fn);
   }
   fclose(fh_param);
}
   
/* Write ------------------------------------------------------------------- */
void XmippCTF::write(const FileName &fn) {
   ofstream fh_param;
   fh_param.open(fn.c_str());
   if (!fh_param)
      REPORT_ERROR(1,(string)"Xmipp_CTF::write: Cannot open "+fn+
         " for output");
   fh_param << *this << endl;
   fh_param.close();
}

/* Usage ------------------------------------------------------------------- */
void XmippCTF::Usage() {
   cerr << "  [defocusU=<DeltafU>]              : Defocus in Angstroms (Ex: -800)\n"
        << "  [defocusV=<DeltafV=DeltafU>]      : If astigmatism\n"
        << "  [azimuthal_angle=<ang=0>]         : Angle between X and U (degrees)\n"
        << "  [sampling_rate=<Tm=1>]            : Angstroms/pixel\n"
        << "  [voltage=<kV=100>]                : Accelerating voltage (kV)\n"
        << "  [spherical_aberration=<Cs=0>]     : Milimiters. Ex: 5.6\n"
        << "  [chromatic_aberration=<Ca=0>]     : Milimiters. Ex: 4\n"
        << "  [energy_loss=<espr=0>]            : eV. Ex: 1\n"
        << "  [lens_stability=<ispr=0>]         : ppm. Ex: 1\n"
        << "  [convergence_cone=<alpha=0>]      : mrad. Ex: 0.5\n"
        << "  [longitudinal_displace=<DeltaF=0>]: Angstrom. Ex: 100\n"
        << "  [transversal_displace=<DeltaR=0>] : Angstrom. Ex: 3\n"
        << "  [Q0=<Q0=0>]                       : Percentage of cosine\n"
        << "  [K=<K=1>]                         : Global gain\n"
	<< endl
        << "  [gaussian_K=<K=0>]                : Gaussian gain\n"
        << "  [sigmaU=<s=0>]                    : gaussian width\n"
        << "  [sigmaV=<s=0>]                    : if astigmatism\n"
        << "  [cU=<s=0>]                        : gaussian center (in cont. freq)\n"
        << "  [cV=<s=0>]                        : if astigmatism\n"
        << "  [gaussian_angle=<ang=0>]          : Angle between X and U (degrees)\n"
        << "  [sqrt_K=<K=0>]                    : Square root gain\n"
        << "  [sqU=<sqU=0>]                     : Square root width\n"
        << "  [sqV=<sqV=0>]                     : if astigmatism\n"
        << "  [sqrt_angle=<ang=0>]              : Angle between X and U (degrees)\n"
        << "  [base_line=<b=0>]                 : Global base line\n"
   ;
}

/* Show -------------------------------------------------------------------- */
ostream & operator << (ostream &out, const XmippCTF &ctf) {
   if (ctf.enable_CTF) {
      out
       /* << "   Aperture:      " << ctf.aperture        << endl */
	  << "sampling_rate=        " << ctf.Tm		  << endl
	  << "voltage=              " << ctf.kV		  << endl
	  << "defocusU=             " << ctf.DeltafU	  << endl
	  << "defocusV=             " << ctf.DeltafV	  << endl
	  << "azimuthal_angle=      " << ctf.azimuthal_angle << endl
	  << "spherical_aberration= " << ctf.Cs              << endl
	  << "chromatic_aberration= " << ctf.Ca              << endl
	  << "energy_loss=          " << ctf.espr            << endl
	  << "lens_stability=       " << ctf.ispr            << endl
	  << "convergence_cone=     " << ctf.alpha           << endl
	  << "longitudinal_displace=" << ctf.DeltaF          << endl
	  << "transversal_displace= " << ctf.DeltaR          << endl
	  << "Q0=                   " << ctf.Q0              << endl
	  << "K=                    " << ctf.K               << endl
      ;
   }
   if (ctf.enable_CTFnoise) {
      out << "gaussian_K=           " << ctf.gaussian_K << endl
	  << "sigmaU=               " << ctf.sigmaU	<< endl
	  << "sigmaV=               " << ctf.sigmaV	<< endl
	  << "cU=                   " << ctf.cU 	<< endl
	  << "cV=                   " << ctf.cV 	<< endl
	  << "gaussian_angle=       " << ctf.gaussian_angle << endl
	  << "sqrt_K=               " << ctf.sqrt_K	<< endl
	  << "sqU=                  " << ctf.sqU	<< endl
	  << "sqV=                  " << ctf.sqV	<< endl
          << "sqrt_angle=           " << ctf.sqrt_angle << endl
	  << "base_line=            " << ctf.base_line  << endl
      ;
   }
   return out;
}

/* Default values ---------------------------------------------------------- */
void XmippCTF::clear() {
   enable_CTF=true;
   enable_CTFnoise=false;
   clear_noise();
   clear_pure_ctf();
}

void XmippCTF::clear_noise() {
   cU=cV=sigmaU=sigmaV=gaussian_angle=gaussian_K=0;
   sqU=sqV=sqrt_K=sqrt_angle=0;
   base_line=0;
}

void XmippCTF::clear_pure_ctf() {
   // aperture=10;
   enable_CTF=true;
   enable_CTFnoise=false;
   Tm=2;
   kV=100;
   DeltafU=DeltafV=azimuthal_angle=0;
   Cs=Ca=espr=ispr=alpha=DeltaF=DeltaR=0;
   K=1;
   Q0=0;
}

/* Produce Side Information ------------------------------------------------ */
void XmippCTF::Produce_Side_Info() {
   // Change units
   // aperture *= 1e4; 
   double local_alpha=alpha/1000;
   double local_Cs = Cs*1e7;
   double local_Ca = Ca*1e7;
   double local_kV = kV*1e3;
   double local_ispr= ispr*1e6;
   rad_azimuth=DEG2RAD(azimuthal_angle);
   rad_gaussian=DEG2RAD(gaussian_angle);
   rad_sqrt=DEG2RAD(sqrt_angle);

   // lambda=h/sqrt(2*m*e*kV)
   //    h: Planck constant
   //    m: electron mass
   //    e: electron charge
   // lambda=0.387832/sqrt(kV*(1.+0.000978466*kV)); // Hewz: Angstroms
   lambda=12.3/sqrt(local_kV*(1.+local_kV*1e-6)); // ICE
   
   // Phase shift for spherical aberration
   // X(u)=-PI*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
   // ICE: X(u)=-PI/2*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
   //          = K1*deltaf(u)*u^2         +K2*u^4
   K1=PI/2*2*lambda;
   K2=PI/2*local_Cs*lambda*lambda*lambda;

   // Envelope
   // D(u)=Ed(u)*Ealpha(u)
   // Ed(u)=exp(-1/2*PI^2*lambda^2*D^2*u^4)
   // Ealpha(u)=exp(-PI^2*alpha^2*u^2*(Cs*lambda^2*u^2+Deltaf(u))^2)
   // ICE: Eespr(u)=exp(-(1/4*PI*Ca*lambda*espr/kV)^2*u^4/log2)
   // ICE: Eispr(u)=exp(-(1/2*PI*Ca*lambda*ispr)^2*u^4/log2)
   // ICE: EdeltaF(u)=bessj0(PI*DeltaF*lambda*u^2)
   // ICE: EdeltaR(u)=sinc(u*DeltaR)
   // ICE: Ealpha(u)=exp(-PI^2*alpha^2*(Cs*lambda^2*u^3+Deltaf(u)*u)^2)
   // CO: K3=pow(0.25*PI*Ca*lambda*(espr/kV,2)/log(2); Both combines in new K3
   // CO: K4=pow(0.5*PI*Ca*lambda*ispr,2)/log(2);
   K3=pow(0.25*PI*local_Ca*lambda*(espr/kV+2*local_ispr),2)/log(2.0);
   K5=PI*DeltaF*lambda;
   K6=PI*PI*alpha*alpha;
   K7=local_Cs*lambda*lambda;
}

/* Zero -------------------------------------------------------------------- */
void XmippCTF::zero(int n, const matrix1D<double> &u, matrix1D<double> &freq) {
   double wmax=1/(2*Tm);
   double wstep=wmax/300;
   int sign_changes=0;
   double last_ctf=CTF_at(0,0), ctf;
   double w;
   for (w=0; w<=wmax; w+=wstep) {
      V2_BY_CT(freq,u,w);
      ctf=CTFpure_at(XX(freq),YY(freq));
      if (SGN(ctf)!=SGN(last_ctf)) {
         sign_changes++;
	 if (sign_changes==n) break;
      }
      last_ctf=ctf;
   }
   if (sign_changes!=n) {
      VECTOR_R2(freq,-1,-1);
   } else {
      // Compute more accurate zero
      w+=ctf*wstep/(last_ctf-ctf);
      V2_BY_CT(freq,u,w);
   }
}

/* Apply the CTF to an image ----------------------------------------------- */
void XmippCTF::Apply_CTF(matrix2D < complex<double> > &FFTI) const {
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(FFTI) {
      XX(idx)=j; YY(idx)=i;
      FFT_idx2digfreq(FFTI, idx, freq);
      double ctf=CTF_at(XX(freq),YY(freq));
      FFTI(i,j)*=ctf;
   }
}

/* Generate CTF Image ------------------------------------------------------ */
//#define DEBUG
void XmippCTF::Generate_CTF(int Ydim, int Xdim,
   matrix2D < complex<double> > &CTF) const {
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   CTF.resize(Ydim, Xdim);
   #ifdef DEBUG
      cout << "CTF:\n" << *this << endl;
   #endif
   FOR_ALL_ELEMENTS_IN_MATRIX2D(CTF) {
      XX(idx)=j; YY(idx)=i;
      FFT_idx2digfreq(CTF, idx, freq);
      digfreq2contfreq(freq, freq, Tm);
      CTF(i,j)=CTF_at(XX(freq),YY(freq));
      #ifdef DEBUG
         if (i==0)
            cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
                 << " " << CTF(i,j) << endl;
      #endif
   }
}
#undef DEBUG

/* Physical meaning -------------------------------------------------------- */
//#define DEBUG
bool XmippCTF::physical_meaning() {
   bool retval;
   if (enable_CTF) {
      retval=
          K>=0       && base_line>=0  &&
          kV>=50     && kV<=1000      &&
	  espr>=0    && espr<=20      &&
	  ispr>=0    && ispr<=20      &&
	  Cs>=0      && Cs<=20        &&
	  Ca>=0      && Ca<=20        &&
	  alpha>=0   && alpha<=5      &&
	  DeltaF>=0  && DeltaF<=1000  &&
	  DeltaR>=0  && DeltaR<=100   &&
	  Q0>=-0.40  && Q0<=0         &&
          DeltafU<=0 && DeltafV<=0    &&
	  CTF_at(0,0)>=0;
      #ifdef DEBUG
         cout << *this << endl;
         cout << "K>=0       && base_line>=0  " << (K>=0       && base_line>=0) << endl
      	      << "kV>=50     && kV<=1000      " << (kV>=50     && kV<=1000)     << endl
	      << "espr>=0    && espr<=20      " << (espr>=0    && espr<=20)     << endl
	      << "ispr>=0    && ispr<=20      " << (ispr>=0    && ispr<=20)     << endl
	      << "Cs>=0      && Cs<=20        " << (Cs>=0      && Cs<=20)       << endl
	      << "Ca>=0      && Ca<=20        " << (Ca>=0      && Ca<=20)       << endl
	      << "alpha>=0   && alpha<=5      " << (alpha>=0   && alpha<=5)     << endl
	      << "DeltaF>=0  && DeltaF<=1000  " << (DeltaF>=0  && DeltaF<=1000) << endl
	      << "DeltaR>=0  && DeltaR<=100   " << (DeltaR>=0  && DeltaR<=100)  << endl
	      << "Q0>=-0.40  && Q0<=0	      " << (Q0>=-0.40  && Q0<=0)        << endl
      	      << "DeltafU<=0 && DeltafV<=0    " << (DeltafU<=0 && DeltafV<=0)   << endl
	      << "CTF_at(0,0)>=0	      " << (CTF_at(0,0)>=0)	        << endl
         ;
         cout << "CTF_at(0,0)=" << CTF_at(0,0,true) << endl;
      #endif
   } else retval=true;
   bool retval2;
   if (enable_CTFnoise) {
      double min_sigma=MIN(sigmaU,sigmaV);
      double min_c=MIN(cU,cV);
      retval2=
	  gaussian_K>=0     &&
	  base_line>=0      &&
	  sigmaU>=0         && sigmaV>=0     	   &&
	  sigmaU<=100e3     && sigmaV<=100e3       &&
	  cU>=0             && cV>=0         	   &&
	  sqU>=0            && sqV>=0        	   &&
	  sqrt_K>=0         &&
	  gaussian_angle>=0 && gaussian_angle<=90 &&
	  sqrt_angle>=0     && sqrt_angle<=90
      ;
      if (min_sigma>0)   retval2=retval2 && ABS(sigmaU-sigmaV)/min_sigma<3;
      if (min_c>0)       retval2=retval2 && ABS(cU-cV)/min_c<3;
      if (gaussian_K!=0) retval2=retval2 && (cU*Tm>=0.01) && (cV*Tm>=0.01);
      #ifdef DEBUG
         cout << *this << endl;
         cout << "gaussian_K>=0     &&  		  " << (gaussian_K>=0	  )			 << endl
      	      << "base_line>=0      &&  		  " << (base_line>=0	  )			 << endl
	      << "sigmaU>=0	    && sigmaV>=0	  " << (sigmaU>=0	  && sigmaV>=0) 	 << endl
	      << "sigmaU<=100e3     && sigmaV<=100e3      " << (sigmaU<=100e3     && sigmaV<=100e3)      << endl
	      << "cU>=0 	    && cV>=0		  " << (cU>=0		  && cV>=0 )		 << endl
	      << "sqU>=0	    && sqV>=0		  " << (sqU>=0  	  && sqV>=0)		 << endl
	      << "sqrt_K>=0	    &&  		  " << (sqrt_K>=0)				 << endl
	      << "gaussian_angle>=0 && gaussian_angle<=90 " << (gaussian_angle>=0 && gaussian_angle<=90) << endl
	      << "sqrt_angle>=0     && sqrt_angle<=90     " << (sqrt_angle>=0     && sqrt_angle<=90)     << endl
	      << "ABS(sigmaU-sigmaV)/min_sigma<3          " << (ABS(sigmaU-sigmaV)/MAX(min_sigma,1)<3)   << endl
	      << "ABS(cU-cV)/min_c<3                      " << (ABS(cU-cV)/MAX(min_c,1)<3)               << endl
	      << "(cU*Tm>=0.01) && (cV*Tm>=0.01)          " << ((cU*Tm>=0.01) && (cV*Tm>=0.01))          << endl
         ;
      #endif
   } else retval2=true;
   #ifdef DEBUG
       cout << "Retval= " << retval << " retval2= " << retval2 << endl;
   #endif
   return retval && retval2;
}
#undef DEBUG
