/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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

#include <data/args.h>
#include <data/geometry.h>
#include <stdio.h>
#include "micrograph_phase_flipping.h"

void Prog_micrograph_phase_flipping::show(void)
{
    std::cout
    << "input_micrograph:      " << fn_in    << std::endl
    << "output_micrograph:     " << fn_out << std::endl
    << "ctf_param_file:        " << fnt_ctf << std::endl
    ;
}
/********
In order to make a clean implementation 4 arrays are needed
1) micrograph
2) init data for Fourier transform
3) output data for fourier transform
4) phase corrected micrograph

The micrograph (that may be in any kind of variable) will be access througt 
mmap
before allocating space for the xmipp image we will free fIn

*/
void Prog_micrograph_phase_flipping::run(void)
{
    //#define DEBUGTIME
    #ifdef  DEBUGTIME
    #include <ctime>
    int clo = clock();
    #endif
    // Read input micrograph
    M_in.open_micrograph(fn_in,reversed);
    int Ydim, Xdim;
    M_in.size(Xdim, Ydim); //micrograph size
    //out image with phase flipped
    #ifdef  DEBUGTIME
    printf ("open micrograph %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    //image dimension
    int fNdim = 2;
    int * fN ;
    fN = new int[fNdim];
    //image size
    fN[0] = Xdim;
    fN[1] = Ydim;
    
    //get access to output image 1D array
    //init fourier transform object
    bool inplace=false;
    xmippFftw      myfft(fNdim, // transform dimension 
                         fN, // size of each dimensio
                         inplace, //inplace transformation (only for 1D)
                         NULL); // if NULL fftw will alloc input 
                                // and output memory
    #ifdef  DEBUGTIME
    printf ("create fft %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    //create plan (fftw stuff)
    myfft.Init("ES",FFTW_FORWARD,false);
    //compute fft size
    //int fTotalSize;
    //fTotalSize=1.;
    //for (int i=0; i<fNdim; i++){
    //   fTotalSize*=fN[i];
    //}
    #ifdef  DEBUGTIME
    printf ("create plan %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    // copy data to input matrix allocated in 
    // xmippFftw myfft(fNdim, fN, inplace,Mout);
    // Note that we are using the xmippImage array as input array
    for(int j=0;j<Xdim;j++)
        for(int i=0;i<Ydim;i++)
            {
            myfft.fIn[i + Ydim * j]=(double)M_in(j,i) / myfft.fTotalSize;//x,y
            }
    M_in.close_micrograph();
        
    #ifdef  DEBUGTIME
    printf ("copy data from micrograph to array %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    //do the forward transform.
    myfft.Transform();
    #ifdef  DEBUGTIME
    printf ("do fft %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    // Compute CTF
    XmippCTF ctf;
    ctf.clear();
    ctf.read(fnt_ctf);
    ctf.enable_CTF = true;
    ctf.enable_CTFnoise = false;
    ctf.Produce_Side_Info();

    Matrix1D<int>    idx(2);  // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    //ROB add std::
    std::complex<double> * cfOut, * cfIn;
    cfOut = (std::complex<double> *)myfft.fOut; //fOut is double *
    cfIn  = (std::complex<double> *)myfft.fIn;  //fIn is double *
    //complex<double> * aux_cfOut, * aux_cfIn;
    int xsize = Xdim;
    int ysize = (int) Ydim/2 +1;

    #ifdef  DEBUGTIME
    printf ("Some book keeping %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    
    ////////////////////you must shift the center
    ////////////////////
    int ii=0;
    for (int j=0;j<xsize/2;j++)
        for (int i=0;i<ysize;i++)
           {
           //XX(idx)=j; YY(idx)=i;
           FFT_IDX2DIGFREQ(j, Xdim, VEC_ELEM(freq, 0));
           FFT_IDX2DIGFREQ(i, Ydim, VEC_ELEM(freq, 1));
           //FFT_idx2digfreq(M_out(), idx, freq); 
           digfreq2contfreq(freq, freq, ctf.Tm);
           if(ctf.CTFpureNodumping_at(XX(freq),YY(freq))<0)
           {
               cfIn[(j*ysize+i)] = (-1.)*(cfOut[(j*ysize+i)]);//*(complex<double>)(+1,0);
           }
           else
           {
               cfIn[(j*ysize+i)] = (cfOut[(j*ysize+i)]);//*(complex<double>)(+1,0);
           }
//           if (idx.module()> 0)//
//                cfIn[ii]  =0.;  //
           //ii++;
           }
    for (int j=xsize/2;j<xsize;j++)
        for (int i=0;i<ysize;i++)
           {
           //XX(idx)=xsize-j; YY(idx)=i;
           FFT_IDX2DIGFREQ(xsize-j, Xdim, VEC_ELEM(freq, 0));
           FFT_IDX2DIGFREQ(i, Ydim, VEC_ELEM(freq, 1));
           //FFT_idx2digfreq(M_out(), idx, freq); 
           digfreq2contfreq(freq, freq, ctf.Tm);
           if(ctf.CTFpureNodumping_at(XX(freq),YY(freq))<0)
           {
               cfIn[(j*ysize+i)] = (-1.)*(cfOut[(j*ysize+i)]);//*(complex<double>)(+1,0);
           }
           else
           {
               cfIn[(j*ysize+i)] = (cfOut[(j*ysize+i)]);//*(complex<double>)(+1,0);
           }
        }
    //#define DEBUG5
    #ifdef  DEBUG5
    ImageXmipp     M_out_aux(Ydim,Xdim);
    for (int j=0;j<xsize;j++)
        for (int i=0;i<ysize;i++)
           {
           M_out_aux(i,j)=log(sqrt((myfft.fOut[(j*ysize+i)*2]) *
                                   (myfft.fOut[(j*ysize+i)*2]) +
                                   (myfft.fOut[(j*ysize+i)*2+1])*
                                   (myfft.fOut[(j*ysize+i)*2+1])))
                          ;
           }
    //remove dc term
    M_out_aux(0,0)=0.;
       
    M_out_aux.write("mypower.spi");
    for (int i=0;i<ysize;i++)
        for (int j=0;j<xsize/2;j++)
           {
           //XX(idx)=j; YY(idx)=i;
           FFT_IDX2DIGFREQ(j, Xdim, VEC_ELEM(freq, 0));
           FFT_IDX2DIGFREQ(i, Ydim, VEC_ELEM(freq, 1));
           //FFT_idx2digfreq(M_out(), idx, freq); 
           digfreq2contfreq(freq, freq, ctf.Tm);
           M_out_aux(i,j) = ctf.CTFpureNodumping_at(XX(freq),YY(freq));
           }
    for (int i=0;i<ysize;i++)
        for (int j=xsize/2;j<xsize;j++)
           {
           //XX(idx)=xsize-j; YY(idx)=i;
           FFT_IDX2DIGFREQ(xsize-j, Xdim, VEC_ELEM(freq, 0));
           FFT_IDX2DIGFREQ(i, Ydim, VEC_ELEM(freq, 1));
           //FFT_idx2digfreq(M_out(), idx, freq); 
           digfreq2contfreq(freq, freq, ctf.Tm);
           M_out_aux(i,j) = ctf.CTFpureNodumping_at(XX(freq),YY(freq));
           }
    M_out_aux.write("myctf.spi");
    #endif


    #ifdef  DEBUGTIME
    printf ("invert phase %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif

    //invert the transform
    myfft.Init("ES",FFTW_BACKWARD,false);
    myfft.Transform();
    #ifdef  DEBUGTIME
    printf ("inverse transform %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    myfft.delete_fIn();
    ImageXmipp     M_out(Ydim,Xdim);
    #ifdef  DEBUGTIME
    printf ("alloc memory for output image %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif


    // write image as spider

    double * imgOut;
    imgOut = (double *)MULTIDIM_ARRAY(M_out());
    //may be a for will be saffer
    //M_out().data = imgOut;
    ii=0;
    
    for(int j=0;j<Xdim;j++)
        for(int i=0;i<Ydim;i++)
            {
            imgOut[j+i*Xdim]=myfft.fOut[i + Ydim * j];
            }

    #ifdef  DEBUGTIME
    printf ("write new micrograph %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    // Close
    M_out.write(fn_out);
}
