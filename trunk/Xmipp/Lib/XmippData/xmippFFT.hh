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
/*****************************************************************************/
/* Fourier Transforms                                                        */
/*****************************************************************************/
#ifndef _XMIPP_FFT_HH
   #define _XMIPP_FFT_HH

#include <complex>
#include "xmippMatrices1D.hh"
#include "xmippMatrices2D.hh"
#include "xmippMatrices3D.hh"
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippFuncs.hh>

/**@name Fourier Transforms */
//@{
/**@name Index <--> Frequency, Continuous <--> Discrete*/
//@{
   /** Index to frequency.
       Given an index and a size of the FFT, this function returns the
       corresponding digital frequency (-1/2 to 1/2). */
   #define FFT_IDX2DIGFREQ(idx,size,freq) \
       freq=((double)((idx)<((size)>>1))?(idx):-(size)+(idx))/(double)(size);

   /** Frequency to index (int).
       Given a frequency and a size of the FFT, this macro returns the
       corresponding integer index. */
   #define DIGFREQ2FFT_IDX(freq,size,idx) {\
       (idx)=(int)(ROUND((size)*(freq))); if ((idx)<0) (idx)+=(int)(size);}

   /** Frequency to index (double).
       Given a frequency and a size of the FFT, this macro returns the
       corresponding double index. */
   #define DIGFREQ2FFT_IDX_DOUBLE(freq,size,idx) {\
       (idx)=((size)*(freq)); if ((idx)<0) (idx)+=(size);}

   /** Index to frequency.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   template <class T>
   void FFT_idx2digfreq(T &v, const matrix1D<int> &idx,
      matrix1D<double> &freq) {
         if (XSIZE(idx)<1 || XSIZE(idx)>3)
            REPORT_ERROR(1,"FFT_idx2digfreq: Index is not of the correct size");
         freq.resize(XSIZE(idx));

         int size[3]; v.get_size(size);
         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
            FFT_IDX2DIGFREQ(VEC_ELEM(idx,i),size[i],VEC_ELEM(freq,i));
   }

   /** Frequency to index.
       This function can be used with vectors of any size (1,2,3).
       The Digital spectrum is limited between -1/2 and 1/2.
       If the vector has got more than 3 coordinates, then an exception
       is thrown*/
   template <class T>
   void digfreq2FFT_idx(T &v, const matrix1D<double> &freq,
      matrix1D<int> &idx) {
         if (XSIZE(freq)<1 || XSIZE(freq)>3)
            REPORT_ERROR(1,"digfreq2FFT_idx: freq is not of the correct size");
         idx.resize(XSIZE(freq));

         int size[3]; v.get_size(size);
         FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
            DIGFREQ2FFT_IDX(VEC_ELEM(freq,i),size[i],VEC_ELEM(idx,i));
   }


   /** Digital to Continuous frequency.
       The pixel size must be given in Amstrongs. The digital frequency is
       between [-1/2,1/2].*/
   inline void digfreq2contfreq(const matrix1D<double> &digfreq,
      matrix1D<double> &contfreq, double pixel_size)
      {contfreq=digfreq/pixel_size;}
   
   /** Continuous to Digital frequency.
       The pixel size must be given in Amstrongs. The digital frequency is
       between [-1/2,1/2].*/
   inline void contfreq2digfreq(const matrix1D<double> &contfreq,
      matrix1D<double> &digfreq, double pixel_size)
      {digfreq=contfreq*pixel_size;}
//@}

/**@name Format conversions.*/
//@{
    /** Conversion from whole -> half 2D.*/
    void Whole2Half(const matrix2D<complex<double> > &in, 
		    matrix2D<complex<double> > &out );

    /** Conversion from half -> whole 2D.*/
    void Half2Whole(const matrix2D<complex<double> > &in, 
		    matrix2D<complex<double> > &out, int oriydim);

//@}

/**@name Fourier Transforms */
//@{
    /** Direct Fourier Transform 1D.*/
    void FourierTransform(const matrix1D<double> &in,
       matrix1D< complex<double> > &out);

    /** Direct Fourier Transform 2D.*/
    void FourierTransform(const matrix2D<double> &in,
       matrix2D< complex<double> > &out);

    /** Direct Fourier Transform 3D.*/
    void FourierTransform(const matrix3D<double> &in,
       matrix3D< complex<double> > &out);

    /** Inverse Fourier Transform 1D.*/
    void InverseFourierTransform(const matrix1D< complex<double> > &in,
       matrix1D<double> &out);

    /** Inverse Fourier Transform 2D.*/
    void InverseFourierTransform(const matrix2D< complex<double> > &in,
       matrix2D<double> &out);

    /** Inverse Fourier Transform 3D.*/
    void InverseFourierTransform(const matrix3D< complex<double> > &in,
       matrix3D<double> &out);

    /** Direct Fourier Transform 2D, output half of (centro-symmetric) transform */
    void FourierTransformHalf(const matrix2D<double> &in,
       matrix2D< complex<double> > &out);

    /** Inverse Fourier Transform 2D, input half of (centro-symmetric) transform */
    void InverseFourierTransformHalf(const matrix2D< complex<double> > &in,
       matrix2D<double> &out, int oriydim);

//@}

/**@name Operations with the Fourier Transforms */
//@{
    /** CenterFFT 1D. */
    template <class T>
    void CenterFFT(matrix1D<T> &v, bool forward) {
         matrix1D<T> aux;
         int l, shift;

         l=XSIZE(v); aux.resize(l); shift=(int)(l/2);
         if (!forward) shift=-shift;
         // Shift the input in an auxiliar vector
         for (int i=0; i<l; i++) {
            int ip=i+shift; if (ip<0) ip+=l; else if (ip>=l) ip-=l;
            aux(ip)=DIRECT_VEC_ELEM(v,i);
         }
         // Copy the vector
         for (int i=0; i<l; i++)
            DIRECT_VEC_ELEM(v,i)=DIRECT_VEC_ELEM(aux,i);
    }

    /** CenterFFT 2D. */
    template <class T>
    void CenterFFT(matrix2D<T> &v, bool forward) {
         matrix1D<T> aux;
         int l, shift;

         // Shift in the X direction
         l=XSIZE(v); aux.resize(l); shift=(int)(l/2);
         if (!forward) shift=-shift;
         for (int i=0; i<YSIZE(v); i++) {
            // Shift the input in an auxiliar vector
            for (int j=0; j<l; j++) {
               int jp=j+shift; if (jp<0) jp+=l; else if (jp>=l) jp-=l;
               aux(jp)=DIRECT_MAT_ELEM(v,i,j);
            }
            // Copy the vector
            for (int j=0; j<l; j++)
               DIRECT_MAT_ELEM(v,i,j)=DIRECT_VEC_ELEM(aux,j);
         }

         // Shift in the Y direction
         l=YSIZE(v); aux.resize(l); shift=(int)(l/2);
         if (!forward) shift=-shift;
         for (int j=0; j<XSIZE(v); j++) {
            // Shift the input in an auxiliar vector
            for (int i=0; i<l; i++) {
               int ip=i+shift; if (ip<0) ip+=l; else if (ip>=l) ip-=l;
               aux(ip)=DIRECT_MAT_ELEM(v,i,j);
            }
            // Copy the vector
            for (int i=0; i<l; i++)
               DIRECT_MAT_ELEM(v,i,j)=DIRECT_VEC_ELEM(aux,i);
         }
    }

    /** CenterFFT 3D. */
    template <class T>
    void CenterFFT(matrix3D<T> &v, bool forward) {
         matrix1D<T> aux;
         int l, shift;

         // Shift in the X direction
         l=XSIZE(v); aux.resize(l); shift=(int)(l/2);
         if (!forward) shift=-shift;
         for (int k=0; k<ZSIZE(v); k++)
            for (int i=0; i<YSIZE(v); i++) {
               // Shift the input in an auxiliar vector
               for (int j=0; j<l; j++) {
                  int jp=j+shift; if (jp<0) jp+=l; else if (jp>=l) jp-=l;
                  aux(jp)=DIRECT_VOL_ELEM(v,k,i,j);
               }
               // Copy the vector
               for (int j=0; j<l; j++)
                  DIRECT_VOL_ELEM(v,k,i,j)=DIRECT_VEC_ELEM(aux,j);
            }

         // Shift in the Y direction
         l=YSIZE(v); aux.resize(l); shift=(int)(l/2);
         if (!forward) shift=-shift;
         for (int k=0; k<ZSIZE(v); k++)
            for (int j=0; j<XSIZE(v); j++) {
               // Shift the input in an auxiliar vector
               for (int i=0; i<l; i++) {
                  int ip=i+shift; if (ip<0) ip+=l; else if (ip>=l) ip-=l;
                  aux(ip)=DIRECT_VOL_ELEM(v,k,i,j);
               }
               // Copy the vector
               for (int i=0; i<l; i++)
                  DIRECT_VOL_ELEM(v,k,i,j)=DIRECT_VEC_ELEM(aux,i);
            }

         // Shift in the Z direction
         l=ZSIZE(v); aux.resize(l); shift=(int)(l/2);
         if (!forward) shift=-shift;
         for (int i=0; i<YSIZE(v); i++)
            for (int j=0; j<XSIZE(v); j++) {
               // Shift the input in an auxiliar vector
               for (int k=0; k<l; k++) {
                  int kp=k+shift; if (kp<0) kp+=l; else if (kp>=l) kp-=l;
                  aux(kp)=DIRECT_VOL_ELEM(v,k,i,j);
               }
               // Copy the vector
               for (int k=0; k<l; k++)
                  DIRECT_VOL_ELEM(v,k,i,j)=DIRECT_VEC_ELEM(aux,k);
            }
    }

    /** FFT Magnitude 1D. */
    void FFT_magnitude(const matrix1D< complex<double> > &v,
       matrix1D<double> &mag);

    /** FFT Magnitude 2D. */
    void FFT_magnitude(const matrix2D< complex<double> > &v,
       matrix2D<double> &mag);

    /** FFT Magnitude 3D. */
    void FFT_magnitude(const matrix3D< complex<double> > &v,
       matrix3D<double> &mag);

    /** FFT Phase 1D. */
    void FFT_phase(const matrix1D< complex<double> > &v,
       matrix1D<double> &phase);

    /** FFT Phase 2D. */
    void FFT_phase(const matrix2D< complex<double> > &v,
       matrix2D<double> &phase);

    /** FFT Phase 3D. */
    void FFT_phase(const matrix3D< complex<double> > &v,
       matrix3D<double> &phase);

    /** Autocorrelation function of an Xmipp matrix.
        Fast calcuation of the autocorrelation matrix of a given one using
        Fast Fourier Transform. (Using the correlation theorem) */
  template <class T>
     void auto_correlation_matrix(matrix2D<T> const &Img,
        matrix2D<double> &R) {
         // Compute the Fourier Transform
         matrix2D< complex<double> > FFT1; FourierTransform(Img, FFT1);

         // Multiply FFT1 * FFT1'
         FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
            FFT1(i,j)*=conj(FFT1(i,j));

         // Invert the product, in order to obtain the correlation image
         InverseFourierTransform(FFT1,R);

         // Center the resulting image to obtain a centered autocorrelation
         CenterFFT(R,true);
     }

   /** Autocorrelation function of an Xmipp matrix.
       Fast calcuation of the correlation matrix on two matrices using
       Fast Fourier Transform. (Using the correlation theorem).
       The output matrix must be already resized*/
   template <class T>
     void correlation_matrix(matrix2D<T> const &m1,matrix2D<T> const &m2,
        matrix2D<double> &R) {
         // Compute the Fourier Transforms
         matrix2D< complex<double> > FFT1; FourierTransform(m1, FFT1);
         matrix2D< complex<double> > FFT2; FourierTransform(m2, FFT2);

         // Multiply FFT1 * FFT2'
         FOR_ALL_ELEMENTS_IN_MATRIX2D(FFT1)
            FFT1(i,j)*=conj(FFT2(i,j));

         // Invert the product, in order to obtain the correlation image
         InverseFourierTransform(FFT1,R);

         // Center the resulting image to obtain a centered autocorrelation
         CenterFFT(R,true);
     }

   /** Fourier-Ring-Correlation between two 2D-matrices using Fast Fourier Transform.*/
   template <class T>
     void fourier_ring_correlation(matrix2D<T> const &m1,matrix2D<T> const &m2, double sampling_rate,
        matrix1D<double> &freq, matrix1D<double> &frc, matrix1D<double> &frc_noise) {
        if (!m1.same_shape(m2))
           REPORT_ERROR(1,"Fourier_ring_correlation: matrices have different shapes!");

        matrix2D<T> aux(m1);
        matrix1D<int> origin(3),radial_count;
        matrix1D<double> tmp1,tmp2;
        matrix2D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
        matrix2D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
        matrix2D <double> realFT1; realFT1.resize(FT1);
        int dim=(int)FT1.RowNo()/2;
        origin.init_zeros();

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
          dMij(aux,i,j)=abs(dMij(FT1,i,j))*abs(dMij(FT1,i,j));
        }
        tmp1.init_zeros();
        radial_average(aux,origin,tmp1,radial_count,true);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
          dMij(aux,i,j)=abs(dMij(FT2,i,j))*abs(dMij(FT2,i,j));
        }
        tmp2.init_zeros();
        radial_average(aux,origin,tmp2,radial_count,true);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(FT1) {
          dMij(realFT1,i,j)=real(conj(dMij(FT1,i,j))*dMij(FT2,i,j));
        }
        frc.init_zeros();
        radial_average(realFT1,origin,frc,radial_count,true);
        frc.resize(dim);
        frc_noise.resize(dim);
        freq.resize(dim);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
          int j=i;
          VEC_ELEM(freq,i)=(double)j/(dim*2*sampling_rate);
          VEC_ELEM(frc,i)=VEC_ELEM(frc,i)/sqrt(VEC_ELEM(tmp1,i)*VEC_ELEM(tmp2,i));
          VEC_ELEM(frc_noise,i)=2/sqrt((double)VEC_ELEM(radial_count,i));
        }
     }

   /** Fourier-Ring-Correlation between two 3D-matrices using Fast Fourier Transform.*/
   template <class T>
     void fourier_ring_correlation(matrix3D<T> const &m1,matrix3D<T> const &m2, double sampling_rate,
        matrix1D<double> &freq, matrix1D<double> &frc, matrix1D<double> &frc_noise) {
        if (!m1.same_shape(m2))
          REPORT_ERROR(1,"Fourier_ring_correlation: matrices have different shapes!");

        matrix3D<T> aux(m1);
        matrix1D<int> origin(3),radial_count;
        matrix1D<double> tmp1,tmp2;
        matrix3D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
        matrix3D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
        matrix3D <double> realFT1; realFT1.resize(FT1);
        int dim=(int)FT1.RowNo()/2;
        origin.init_zeros();

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
          dVkij(aux,k,i,j)=abs(dVkij(FT1,k,i,j))*abs(dVkij(FT1,k,i,j));
        }
        tmp1.init_zeros();
        radial_average(aux,origin,tmp1,radial_count,true);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
          dVkij(aux,k,i,j)=abs(dVkij(FT2,k,i,j))*abs(dVkij(FT2,k,i,j));
        }
        tmp2.init_zeros();
        radial_average(aux,origin,tmp2,radial_count,true);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1) {
          dVkij(realFT1,k,i,j)=real(conj(dVkij(FT1,k,i,j))*dVkij(FT2,k,i,j));
        }
        frc.init_zeros();
        radial_average(realFT1,origin,frc,radial_count,true);
        frc.resize(dim);
        frc_noise.resize(dim);
        freq.resize(dim);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
          int j=i;
          VEC_ELEM(freq,i)=(double)j/(dim*2*sampling_rate);
          VEC_ELEM(frc,i)=VEC_ELEM(frc,i)/sqrt(VEC_ELEM(tmp1,i)*VEC_ELEM(tmp2,i));
          VEC_ELEM(frc_noise,i)=2/sqrt((double)VEC_ELEM(radial_count,i));
        }
     }

   /** Differential Phase Residualbetween two 2D-matrices using Fast Fourier Transform.*/
   template <class T>
     void differential_phase_residual(matrix2D<T> const &m1,matrix2D<T> const &m2, double sampling_rate,
        matrix1D<double> &freq, matrix1D<double> &dpr) {
        if (!m1.same_shape(m2))
          REPORT_ERROR(1,"Differential phase residual: matrices have different shapes!");

        matrix2D<T> aux(m1);
        matrix1D<int> origin(3),radial_count;
        matrix1D<double> tmp1,tmp2;
        matrix1D<complex <double> > tmp3;
        matrix2D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
        matrix2D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
        int dim=(int)FT1.RowNo()/2;

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
          dMij(aux,i,j)=abs(dMij(FT1,i,j))+abs(dMij(FT2,i,j));
        }
        tmp1.init_zeros();
        radial_average(aux,origin,tmp1,radial_count,true);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(aux) {
          dMij(aux,i,j)*=realWRAP(RAD2DEG(
                           (atan2(dMij(FT1,i,j).imag(),dMij(FT1,i,j).real())) -
                           (atan2(dMij(FT2,i,j).imag(),dMij(FT2,i,j).real()))) 
				        ,-180,180);
          dMij(aux,i,j)*=realWRAP(RAD2DEG(
                           (atan2(dMij(FT1,i,j).imag(),dMij(FT1,i,j).real())) -
                           (atan2(dMij(FT2,i,j).imag(),dMij(FT2,i,j).real()))) 
				        ,-180,180);
        }
        tmp2.init_zeros();
        radial_average(aux,origin,tmp2,radial_count,true);
        dpr=tmp2/tmp1;
        dpr=SQRTnD(dpr);
        dpr.resize(dim);
        freq.resize(dim);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
          int j=i;
          VEC_ELEM(freq,i)=(double)j/(dim*2*sampling_rate);
        }
     }

   /** Differential Phase Residualbetween two 3D-matrices using Fast Fourier Transform.*/
   template <class T>
     void differential_phase_residual(matrix3D<T> const &m1,matrix3D<T> const &m2, double sampling_rate,
        matrix1D<double> &freq, matrix1D<double> &dpr) {
        if (!m1.same_shape(m2))
          REPORT_ERROR(1,"Differential phase residual: matrices have different shapes!");
        matrix3D<T> aux(m1);
        matrix1D<int> origin(3),radial_count;
        matrix1D<double> tmp1,tmp2;
        matrix1D<complex <double> > tmp3;
        matrix3D<complex <double> > FT1; FourierTransform(m1,FT1); CenterFFT(FT1,true);
        matrix3D<complex <double> > FT2; FourierTransform(m2,FT2); CenterFFT(FT2,true);
        int dim=(int)FT1.RowNo()/2;

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
          dVkij(aux,k,i,j)=abs(dVkij(FT1,k,i,j))+abs(dVkij(FT2,k,i,j));
        }
        tmp1.init_zeros();
        radial_average(aux,origin,tmp1,radial_count,true);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(aux) {
          dVkij(aux,k,i,j)*=realWRAP(RAD2DEG(
                           (atan2(dVkij(FT1,k,i,j).imag(),dVkij(FT1,k,i,j).real())) -
                           (atan2(dVkij(FT2,k,i,j).imag(),dVkij(FT2,k,i,j).real()))) 
				        ,-180,180);
          dVkij(aux,k,i,j)*=realWRAP(RAD2DEG(
                           (atan2(dVkij(FT1,k,i,j).imag(),dVkij(FT1,k,i,j).real())) -
                           (atan2(dVkij(FT2,k,i,j).imag(),dVkij(FT2,k,i,j).real()))) 
				        ,-180,180);
        }
        tmp2.init_zeros();
        radial_average(aux,origin,tmp2,radial_count,true);
        dpr=tmp2/tmp1;
        dpr=SQRTnD(dpr);
        dpr.resize(dim);
        freq.resize(dim);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(freq) {
          int j=i;
          VEC_ELEM(freq,i)=(double)j/(dim*2*sampling_rate);
        }
     }

   /** Signal to noise ratio for 2D */
   #define NGRUPOS 20
   template <class T>
     void my_ssnr(matrix2D<T> const &AverageImage, SelFile  &SF_sel, double sampling_rate,
		 matrix1D<double> &freq, matrix1D<double> &ssnr,
		 matrix1D<double> &pixel,bool apply_geo) {
        Image       Iaverage_sub_group,Id;
        double      dummy;
        SelFile  SF_tmp, SF_vector[NGRUPOS];
        SelLine  line;
        int N=NGRUPOS;

        if(SF_sel.ImgNo(SelLine::ACTIVE) < NGRUPOS)
            REPORT_ERROR(1,(string)"FFT::my_ssnr: my_ssnr I need more than " +
                         ItoA(SF_sel.ImgNo(SelLine::ACTIVE)) +
                          " images in file " + SF_sel.name());

        int contImg=SF_sel.ImgNo();
        //  SF_tmp=SF_sel.randomize();
        SF_tmp=SF_sel;
        int Nsub=(int)((double)contImg/N);//floor()

        //assing images
        for (int i=0; i<N; i++) {
          SF_vector[i].clear();
          SF_tmp.go_beginning();
          SF_tmp.jump_lines(Nsub*i);
          for (int nn=0; nn<Nsub; nn++) {
            SF_vector[i].insert(SF_tmp.current());
            SF_tmp.NextImg();
          }
        }

        //all the images has been assigned?
        if (contImg%N!=0){
           for(int i=(N-(contImg%N));i<N;i++){
	       SF_vector[i].insert(SF_tmp.current());
	       SF_tmp.NextImg();}

        }

        //alloc memory for results      
        int dim=(int)AverageImage.RowNo()/2;
        ssnr.resize(dim);
        pixel.resize(dim);
        freq.resize(dim); 

        //compute average and fft
        //fft radial averaged image
        matrix2D<complex <double> > FTAverageImage,
                                    FTaverageSubGroup;
        FourierTransform(AverageImage,FTAverageImage);
        CenterFFT(FTAverageImage,true);
        for (int i=0; i<N; i++) {
           SF_vector[i].get_statistics(Iaverage_sub_group,Id,dummy,dummy,apply_geo);
           FourierTransform(Iaverage_sub_group(),FTaverageSubGroup);
           CenterFFT(FTaverageSubGroup,true);   
           //calculate ssNR
           my_ssnr_step(FTAverageImage,FTaverageSubGroup,
                        ssnr,pixel,SF_vector[i].ImgNo()); 

         }
         int dim2=(int)AverageImage.RowNo();
         FOR_ALL_ELEMENTS_IN_MATRIX1D(ssnr) {
             VEC_ELEM(ssnr,i)=(contImg*(NGRUPOS-1))/(VEC_ELEM(ssnr,i)-1);
             VEC_ELEM(freq,i)=((double)i+0.5)/(dim2*sampling_rate);
             if(VEC_ELEM(ssnr,i)>NGRUPOS)
                VEC_ELEM(ssnr,i) = NGRUPOS;
         }
     }

   /** Signal to noise ratio for 2D, process one image. freq and ssnr
       must have correct size */
   void my_ssnr_step(matrix2D< complex<double> > const &AverageImage,
                     matrix2D< complex<double> > const &FTaverageSubGroup,
	             matrix1D<double> &ssnr,matrix1D<double> &pixel, int n);

   /** Series convolution function. Gives the convolution of two series
       given as Xmipp Vectors. Result is stored in result vector.
       Fast calcuation of the convolution result using
       Fast Fourier Transform. If FullConvolution, set by default to
       FALSE, is TRUE the full convolution series is returned. Otherwise
       the convolution vector refers only to the valid values, whose
       number is the greater dimension of the two series.
       Note: Complex numbers are allowed */
   template <class T>
      void series_convolution(matrix1D<T> &series1,matrix1D<T> &series2,
         matrix1D<T> &result, bool FullConvolution=false) {
         // Store dimension of series
         int dim1=series1.get_dim();
         int dim2=series2.get_dim();

         // Resize series to the size of the resulting series
         // (Zeros are stored in the expanded values)
         series1.resize(dim1+dim2-1);
         series2.resize(series1);
         result .resize(series1);

          // Fourier Transform the two series
          matrix1D<complex <double> > FFT1; FourierTransform(series1, FFT1);
          matrix1D<complex <double> > FFT2; FourierTransform(series2, FFT2);

          // Multiply the vectors element by element to do the 
          // convolution in the Fourier space.
          FFT1*=FFT2;

          // Recover the convolution result by inverse FFT
          InverseFourierTransform(FFT1,result);

          // Restore the dimensions
          series1.resize(dim1);
          series2.resize(dim2);

          /* If the full convolution is required, nothing more remains to be done.
             Otherwise, if the valid values are required, return the central ones. */
          if (!FullConvolution) {
             // First get the maximum dimension of the series,
             // which is the dimension of the result
             int dim=MAX(dim1,dim2);

             // Determine the number of values to discard
             int discard=result.get_dim()-dim;

             // Divide it by two as we have to discard them in both sides of
             // the vector
             discard/=2;  // Integer division is intended

             // Copy required values (simple displacement of values)
             for (int i=STARTINGX(result); i<STARTINGX(result)+dim; i++)
                result(i)=result(i+discard);

             // and finally resize to discard not copied values
             result.resize(dim);
          }
      }

   /** Numerical_derivative: This function computes the numerical derivative
       of a matrix in Y direction (rows) or X direction (columns) of a given
       matrix, using a Savitzky-Golay filter on every row or column, and then
       convolving.Input matrix is M, result is stored in D, direction can have
       values of 'x' or 'y'. Order is the derivative order.
       Window size and polynomial order are parameters for the
       Savitzky-Golay filter that define the number of points forward
       (+window_size)
       and backward (-window_size) considered to calculate the filter, and the degree of
       the polynomial to interpolate these values, respectively. Default values are
       window_size=2 and polynomial_order=4, which are equivalent to a
       5-point algorithm to calculate the derivate and give good results.
       But they can be changed provided that polynomial_order <= 2*window_size.
       As one can expect, the values of the matrix in a border of size
       window_size are not accurate ones, as there aren't enough points to perform
       the estimation. As a rule of thumb, the greater the window, the more the
       filtering and the less the precission of the derivatives, and the
       greater the order of the polynomial, the greater the precission.  */
    void numerical_derivative(matrix2D<double> &M,matrix2D<double> &D,
       char direction, int order, int window_size=2, int polynomial_order=4);
//@}

//@}
#endif
