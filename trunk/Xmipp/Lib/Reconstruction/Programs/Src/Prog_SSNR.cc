/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2003)
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

#include "../Prog_SSNR.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippProjection.hh>
#include <XmippData/xmippFFT.hh>
#include <Reconstruction/projection.hh>

// Read parameters from command line ---------------------------------------
void Prog_SSNR_prm::read(int argc, char **argv) {
   if (!check_param(argc,argv,"-radial_avg")) {
      fn_S=get_param(argc,argv,"-S");
      fn_N=get_param(argc,argv,"-N");
      fn_Ssel=get_param(argc,argv,"-selS");
      fn_Nsel=get_param(argc,argv,"-selN");
      generate_images=check_param(argc,argv,"-generate_images");
   } else 
      fn_VSSNR=get_param(argc,argv,"-VSSNR");
   ring_width=AtoF(get_param(argc,argv,"-ring","4"));
   Tm=AtoF(get_param(argc,argv,"-sampling_rate","1"));
   fn_out=get_param(argc,argv,"-o","");
}

// Show parameters ---------------------------------------------------------
ostream & operator << (ostream &out, const Prog_SSNR_prm &prm) {
   out << "Signal:         " << prm.fn_S       << endl
       << "Noise:          " << prm.fn_N       << endl
       << "Signal selfile: " << prm.fn_Ssel    << endl
       << "Noise  selfile: " << prm.fn_Nsel    << endl
       << "Volumetric SSNR:" << prm.fn_VSSNR   << endl
       << "Output:         " << prm.fn_out     << endl
       << "Ring width:     " << prm.ring_width << endl
       << "Sampling rate:  " << prm.Tm         << endl
       << "Generate images:" << prm.generate_images << endl;
   ;
   return out;
}

// Usage -------------------------------------------------------------------
void Prog_SSNR_prm::usage() const {
   cerr << " Estimation from images -----------------------------------------\n"
        << "SSNR\n"
        << "   -S <Volume>           : Signal volume\n"
        << "   -N <Volume>           : Noise volume\n"
        << "   -selS <Selfile>       : Selfile with experimental images\n"
        << "   -selN <Selfile>       : Selfile with noise images\n"
        << "  [-o <Text file=\"\">]    : Output file\n"
        << "  [-ring <w=4>]          : Ring width for the SSNR computation\n"
        << "  [-sampling_rate <Tm=1>]: Sampling rate A/pixel\n"
        << "  [-generate_images]     : generate SSNR images\n"
        << " Estimation by radial averaging ---------------------------------\n"
        << " SSNR -radial_avg\n"
        << "   -VSSNR <fn_vol>       : Volume with the Volumetric SSNR\n"
        << "  [-ring <w=4>]          : Ring width for the SSNR computation\n"
        << "  [-sampling_rate <Tm=1>]: Sampling rate A/pixel\n"
        << "  [-o <Text file=\"\">]    : Output file\n"
   ; 
}

// Produce side Info -------------------------------------------------------
void Prog_SSNR_prm::produce_side_info() _THROW {
   if (fn_VSSNR=="") {
      S.read(fn_S); S().set_Xmipp_origin();
      N.read(fn_N); N().set_Xmipp_origin();
      if (!S().same_shape(N()))
         REPORT_ERROR(1,
            "SSNR: Signal and Noise volumes are not of the same size");

      SF_S.read(fn_Ssel);
      SF_N.read(fn_Nsel);
      int sYdim, sXdim; SF_S.ImgSize(sYdim, sXdim);
      int nYdim, nXdim; SF_N.ImgSize(nYdim, nXdim);
      if (sYdim!=nYdim || sYdim!=YSIZE(S()) ||
          sXdim!=nXdim || sXdim!=XSIZE(S()))
         REPORT_ERROR(1,
            "SSNR: conflict among the projection/projection sizes "
            "or projection/volume");

      if (SF_S.ImgNo()!=SF_N.ImgNo())
         REPORT_ERROR(1,
            "SSNR: the number of projections in both selfiles is different");
   } else {
      VSSNR.read(fn_VSSNR);
      VSSNR().set_Xmipp_origin();
   }
}

// Estimate SSNR 1D --------------------------------------------------------
void Prog_SSNR_prm::Estimate_SSNR_1D(
   matrix2D<double> &output) {
   matrix2D<double> output_S, output_N;

   // Compute the uncorrected SSNR for the signal
   Compute_Uncorrected_SSNR_1D(S(), SF_S, ring_width, Tm, output_S, false);
   
   // Compute the uncorrected SSNR for the noise
   Compute_Uncorrected_SSNR_1D(N(), SF_N, ring_width, Tm, output_N, true);

   // Correct the SSNR and produce output
   output.resize(YSIZE(output_S), 9);
   for (int i=0; i<YSIZE(output_S); i++) {
      double w;
      FFT_IDX2DIGFREQ(i,XSIZE(S()),w);
      output(i,0)=i;
      output(i,1)=w*1/Tm;
      double SSNR=output_S(i,2)/output_N(i,2);
      if (SSNR>1) output(i,2)=10*log10(SSNR-1); // Corrected SSNR
      else        output(i,2)=-1000;            // In fact it should be -inf
      output(i,3)=output_S(i,2);
      output(i,4)=output_S(i,3);
      output(i,5)=output_S(i,4);
      output(i,6)=output_N(i,2);
      output(i,7)=output_N(i,3);
      output(i,8)=output_N(i,4);
   }
}

// Estimate SSNR 2D --------------------------------------------------------
void Prog_SSNR_prm::Estimate_SSNR_2D() {
   cerr << "Computing the SSNR 2D ...\n";
   init_progress_bar(SF_S.ImgNo());
   int imgno=0;
   while (!SF_S.eof()) {
      ImageXmipp Is, In;
      Is.read(SF_S.NextImg()); Is().set_Xmipp_origin();
      In.read(SF_N.NextImg()); In().set_Xmipp_origin();
      
      Projection Iths, Ithn;
      project_Volume(S(), Iths, YSIZE(Is()), XSIZE(Is()),
         Is.rot(), Is.tilt(), Is.psi());
      project_Volume(N(), Ithn, YSIZE(Is()), XSIZE(Is()),
         Is.rot(), Is.tilt(), Is.psi());
         
      Is()-=Iths();
      In()-=Ithn();

      matrix2D< complex<double> > FFT_Is;   FourierTransform(Is  (),FFT_Is  );
      matrix2D< complex<double> > FFT_Iths; FourierTransform(Iths(),FFT_Iths);
      matrix2D< complex<double> > FFT_In;   FourierTransform(In  (),FFT_In  );
      matrix2D< complex<double> > FFT_Ithn; FourierTransform(Ithn(),FFT_Ithn);

      // Compute the amplitudes
      ImageXmipp S2s; FFT_magnitude(FFT_Iths,S2s()); S2s()*=S2s();
      ImageXmipp N2s; FFT_magnitude(FFT_Is  ,N2s()); N2s()*=N2s();
      ImageXmipp S2n; FFT_magnitude(FFT_Ithn,S2n()); S2n()*=S2n();
      ImageXmipp N2n; FFT_magnitude(FFT_In  ,N2n()); N2n()*=N2n();

      // Compute the SSNR image
      ImageXmipp SSNR2D; SSNR2D().init_zeros(S2s());
      FOR_ALL_ELEMENTS_IN_MATRIX2D(S2s()) {
         double ISSNR=0, alpha=0, SSNR=0;
         if (N2s(i,j)>XMIPP_EQUAL_ACCURACY) ISSNR=S2s(i,j)/N2s(i,j);
         if (N2n(i,j)>1e-3) alpha=S2n(i,j)/N2n(i,j);
         if (alpha   >XMIPP_EQUAL_ACCURACY) SSNR=MAX(ISSNR/alpha-1,0);
         if (SSNR    >XMIPP_EQUAL_ACCURACY) SSNR2D(i,j)=10*log10(SSNR+1);
      }
      CenterFFT(SSNR2D(),true);
      
      // Set angles
      SSNR2D.rot ()=Is.rot ();
      SSNR2D.tilt()=Is.tilt();
      SSNR2D.psi ()=Is.psi ();

      // Save image
      FileName fn_img_out=fn_out+ItoA(Is.name().get_number(),5)+".xmp";
      SSNR2D.write(fn_img_out);

      // Finished with this image
      if (++imgno%50==0) progress_bar(imgno);
   }
}

// Estimate SSNR 1D --------------------------------------------------------
void Prog_SSNR_prm::Radial_average(
   matrix2D<double> &output) {

   // Compute the radial average ...........................................
   matrix1D<double> VSSNR_avg((int)(XSIZE(VSSNR())/2-ring_width));
   matrix1D<double> K1D(VSSNR_avg);
   matrix1D<int>    idx(3);
   matrix1D<double> freq(3);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(VSSNR()) {
      VECTOR_R3(idx,j,i,k);
      FFT_idx2digfreq(VSSNR(), idx, freq);
      if (XX(freq)<0) continue;

      // Look for the index corresponding to this frequency
      double w=freq.module();
      double widx=w*XSIZE(VSSNR());
      if (widx>=XSIZE(VSSNR_avg)) continue;
      int l0=MAX(0,CEIL(widx-ring_width));
      int lF=FLOOR(widx);

      double VSSNRkij=pow(10,VSSNR(k,i,j)/10)-1;

      for (int l=l0; l<=lF; l++) {
         VSSNR_avg(l)+=VSSNRkij;
         K1D(l)++;
      }
   }

   FOR_ALL_ELEMENTS_IN_MATRIX1D(VSSNR_avg)
      if (K1D(i)!=0) VSSNR_avg(i)/=K1D(i);

   // Produce output .......................................................
   output.resize(XSIZE(VSSNR_avg), 3);
   for (int i=0; i<XSIZE(VSSNR_avg); i++) {
      double w;
      FFT_IDX2DIGFREQ(i,XSIZE(VSSNR()),w);
      output(i,0)=i;
      output(i,1)=w*1/Tm;
      double SSNR=VSSNR_avg(i);
      if (SSNR>1) output(i,2)=10*log10(SSNR-1); // Corrected SSNR
      else        output(i,2)=-1000;            // In fact it should be -inf
   }
}

// Compute uncorrected SSNR ------------------------------------------------
//#define DEBUG
void Compute_Uncorrected_SSNR_1D(matrix3D<double> &V,
   SelFile &SF, double ring_width, double Tm,
   matrix2D<double> &output, bool treat_as_noise) {

   matrix1D<double> S21D((int)(XSIZE(V)/2-ring_width)),
                    N21D((int)(XSIZE(V)/2-ring_width)),
                    K1D ((int)(XSIZE(V)/2-ring_width)),
                    SSNR1D;

   cerr << "Computing the SSNR ...\n";
   init_progress_bar(SF.ImgNo());
   int imgno=0;

   bool first=true;
   while (!SF.eof()) {
      // Read experimental image
      ImageXmipp Iexp;
      Iexp.read(SF.NextImg());
      Iexp().set_Xmipp_origin();
      
      // Generate the projection of the volume in the same direction
      Projection Ith; 
      project_Volume(V, Ith, YSIZE(Iexp()), XSIZE(Iexp()),
         Iexp.rot(), Iexp.tilt(), Iexp.psi());
      
      // Compute the difference
      if (!treat_as_noise) Iexp()-=Ith();

      // Compute the FFT of both images
      matrix2D< complex<double> > FFT_Iexp; FourierTransform(Iexp(),FFT_Iexp);
      matrix2D< complex<double> > FFT_Ith ; FourierTransform(Ith (),FFT_Ith );
      
      // Compute the amplitudes
      ImageXmipp S2; FFT_magnitude(FFT_Ith,S2()); S2()*=S2();
      ImageXmipp N2; FFT_magnitude(FFT_Iexp,N2()); N2()*=N2();
      STARTINGX(S2())=STARTINGY(S2())=0;
      STARTINGX(N2())=STARTINGY(N2())=0;

      // Average over rings
      matrix1D<int>    idx(2);
      matrix1D<double> freq(2);

      FOR_ALL_ELEMENTS_IN_MATRIX2D(S2()) {
         XX(idx)=j; YY(idx)=i;
         FFT_idx2digfreq(FFT_Iexp, idx, freq);
         if (XX(freq)<0) continue;

         // Look for the index corresponding to this frequency
         double w=freq.module();
         double widx=w*XSIZE(V);
         if (widx>=XSIZE(S21D)) continue;
         int l0=MAX(0,CEIL(widx-ring_width));
         int lF=FLOOR(widx);

         double signal=S2(i,j);
         double noise =N2(i,j);
         //cout << "Processing " << j << "," << i << " --> "
         //     << freq.transpose() << " " << w << " " << widx << " --> "
         //     << l0 << "," << lF
         //     << " S2(i,j)=" << S2(i,j) << " N2(i,j)=" << N2(i,j) << endl;

         for (int l=l0; l<=lF; l++) {
            S21D(l)+=signal;
            N21D(l)+=noise;
            K1D(l)++;
         }
      }

      #ifdef DEBUG
         cout << "Image: " << Iexp.name() << endl;
         S2()=10*log10(S2()); S2.write("PPPS2.xmp");
         N2()=10*log10(N2()); N2.write("PPPN2.xmp");
         Iexp.write("PPPdiff.xmp");
         matrix2D<double> output(XSIZE(S21D),5);
         int imax=0;
         FOR_ALL_ELEMENTS_IN_MATRIX1D(S21D) {
            double w;
            FFT_IDX2DIGFREQ(i,XSIZE(V),w);
            if (w<0) {imax=i; break;}
            output(i,0)=i;
            output(i,1)=w*1/Tm;
            output(i,2)=S21D(i)/N21D(i);
            output(i,3)=10*log10(S21D(i)/imgno);
            output(i,4)=10*log10(N21D(i)/imgno);
            imax++;
         }
         output.resize(imax,5);
         output.write("PPPs2n2.txt");
         char c; cout << "Press any key\n"; cin >> c;
      #endif
      
      // Finished with this image
      if (++imgno%50==0) progress_bar(imgno);
   }
   imgno=SF.ImgNo();
   progress_bar(imgno);
   
   // Compute the SSNR
   SSNR1D=S21D/N21D;
   S21D/=K1D;
   N21D/=K1D;
   
   output.resize(XSIZE(SSNR1D),5);
   int imax=0;
   FOR_ALL_ELEMENTS_IN_MATRIX1D(SSNR1D) {
      double w;
      FFT_IDX2DIGFREQ(i,XSIZE(V),w);
      if (w<0) {imax=i; break;}
      output(i,0)=i;
      output(i,1)=w*1/Tm;
      output(i,2)=SSNR1D(i);
      output(i,3)=10*log10(S21D(i)/imgno);
      output(i,4)=10*log10(N21D(i)/imgno);
      imax++;
   }
   output.resize(imax,5);
}
#undef DEBUG

// Compute uncorrected SSNR ------------------------------------------------
void Compute_SSNR_images(matrix3D<double> &S, matrix3D<double> &N,
   SelFile &SF_S, SelFile &SF_N, FileName &fn_root,
   SelFile &SF_out) {

   cerr << "Computing the SSNR images ...\n";
   init_progress_bar(SF_S.ImgNo());
   int imgno=0;

   bool first=true;
   while (!SF_S.eof()) {
      // Read experimental image
      ImageXmipp Iexp;
      Iexp.read(SF_S.NextImg());
      Iexp().set_Xmipp_origin();
      
      // Generate the projection of the volume in the same direction
      Projection Ith; 
      project_Volume(S, Ith, YSIZE(Iexp()), XSIZE(Iexp()),
         Iexp.rot(), Iexp.tilt(), Iexp.psi());
      
      // Compute the difference
      Iexp()-=Ith();

      // Compute the FFT of both images
      matrix2D< complex<double> > FFT_Iexp; FourierTransform(Iexp(),FFT_Iexp);
      matrix2D< complex<double> > FFT_Ith ; FourierTransform(Ith (),FFT_Ith );
      
      // Compute the amplitudes
      ImageXmipp S2; FFT_magnitude(FFT_Ith,S2()); S2()*=S2();
      ImageXmipp N2; FFT_magnitude(FFT_Iexp,N2()); N2()*=N2();
      STARTINGX(S2())=STARTINGY(S2())=0;
      STARTINGX(N2())=STARTINGY(N2())=0;

      // Finished with this image
      if (++imgno%50==0) progress_bar(imgno);
   }
   imgno=SF_S.ImgNo();
   progress_bar(imgno);
}

// Main --------------------------------------------------------------------
void ROUT_SSNR(Prog_SSNR_prm &prm, matrix2D<double> &output) {
   cout << prm;
   prm.produce_side_info();
   if (prm.fn_VSSNR=="") {
      if (!prm.generate_images) {
         prm.Estimate_SSNR_1D(output);
         if (prm.fn_out!="") output.write(prm.fn_out);
         else                output.write(prm.fn_S.insert_before_extension("_SSNR"));
      } else prm.Estimate_SSNR_2D();
   } else {
      prm.Radial_average(output);
      if (prm.fn_out!="") output.write(prm.fn_out);
      else                output.write(prm.fn_VSSNR.insert_before_extension("_radial_avg"));
   }
}
