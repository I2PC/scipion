/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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
   fn_S=get_param(argc,argv,"-S");
   fn_N=get_param(argc,argv,"-N");
   fn_Ssel=get_param(argc,argv,"-selS");
   fn_Nsel=get_param(argc,argv,"-selN");
   fn_out=get_param(argc,argv,"-o","");
   ring_width=AtoF(get_param(argc,argv,"-ring","4"));
   Tm=AtoF(get_param(argc,argv,"-sampling_rate","1"));
}

// Show parameters ---------------------------------------------------------
ostream & operator << (ostream &out, const Prog_SSNR_prm &prm) {
   out << "Signal:         " << prm.fn_S       << endl
       << "Noise:          " << prm.fn_N       << endl
       << "Signal selfile: " << prm.fn_Ssel    << endl
       << "Noise  selfile: " << prm.fn_Nsel    << endl
       << "Output:         " << prm.fn_out     << endl
       << "Ring width:     " << prm.ring_width << endl
       << "Sampling rate:  " << prm.Tm         << endl
   
   ;
   return out;
}

// Usage -------------------------------------------------------------------
void Prog_SSNR_prm::usage() const {
   cerr << "SSNR\n"
        << "   -S <Volume>           : Signal volume\n"
        << "   -N <Volume>           : Noise volume\n"
        << "   -selS <Selfile>       : Selfile with experimental images\n"
        << "   -selN <Selfile>       : Selfile with noise images\n"
        << "  [-o <Text file=\"\">]    : Output file\n"
        << "  [-ring_width <w=4>]    : Ring width for the SSNR computation\n"
        << "  [-sampling_rate <Tm=1>]: Sampling rate A/pixel\n"
   ; 
}

// Produce side Info -------------------------------------------------------
void Prog_SSNR_prm::produce_side_info() _THROW {
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
}

// Estimate SSNR 1D --------------------------------------------------------
void Prog_SSNR_prm::Estimate_SSNR_1D(
   matrix2D<double> &output) {
   matrix2D<double> output_S, output_N;

   // Compute the uncorrected SSNR for the signal
   Compute_Uncorrected_SSNR_1D(S(), SF_S, ring_width, Tm, output_S);

   // Compute the uncorrected SSNR for the noise
   Compute_Uncorrected_SSNR_1D(N(), SF_N, ring_width, Tm, output_N);

   // Correct the SSNR and produce output
   output.resize(YSIZE(output_S), 9);
   for (int i=0; i<YSIZE(output_S); i++) {
      double w;
      FFT_IDX2DIGFREQ(i,XSIZE(S()),w);
      output(i,0)=i;
      output(i,1)=w*1/Tm;
      output(i,2)=10*log10(output_S(i,2)/output_N(i,2)-1); // Corrected SSNR
      output(i,3)=output_S(i,2);
      output(i,4)=output_S(i,3);
      output(i,5)=output_S(i,4);
      output(i,6)=output_N(i,2);
      output(i,7)=output_N(i,3);
      output(i,8)=output_N(i,4);
   }
}

// Compute uncorrected SSNR ------------------------------------------------
//#define DEBUG
void Compute_Uncorrected_SSNR_1D(matrix3D<double> &V,
   SelFile &SF, double ring_width, double Tm,
   matrix2D<double> &output) {

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
      Iexp()-=Ith();

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
   prm.produce_side_info();
   prm.Estimate_SSNR_1D(output);
   if (prm.fn_out!="") output.write(prm.fn_out);
   else                output.write(prm.fn_S.insert_before_extension("_SSNR"));
}
