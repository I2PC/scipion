/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

/* ------------------------------------------------------------------------- */
/* Includes                                                                  */
/* ------------------------------------------------------------------------- */
#include <reconstruction/precompute_sampling.h>
#include <data/fftw.h>
#include <data/error.h>
#include <iostream>

/* ------------------------------------------------------------------------- */
/* Program                                                                   */
/* ------------------------------------------------------------------------- */
#define MYSIZE 10000000
int main(int argc, char *argv[])
{
    #define DEBUGTIME
    #ifdef  DEBUGTIME
    #include <ctime>
    int clo = clock();
    #endif

    bool inplace = false; 
    VolumeXmipp img;
    img.read("kk.vol");
    //img().transpose();
    #ifdef  DEBUGTIME
    printf ("image read in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    
    #define DIMENSION 2
    xmippFftw fft_img(img());
    #ifdef  DEBUGTIME
    printf ("image transformed in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
#ifdef NEVERDEFINE
    std::complex<double> * IMG;
    IMG = (std::complex<double> *)fft_img.fOut;
    
    int ysize = (int)(((double)YSIZE(img()) / 2.) +1.);
    int xsize = XSIZE(img());
    
    img().initZeros();
    long int ii=0;
    //xmipp y contigous, fftw x contigous
    for (int j=0;j<xsize;j++)
        for (int i=0;i<ysize;i++)
        {
           DIRECT_MAT_ELEM(img(),i,j) = abs(IMG[ii]);
           ii++;
        }
    img.write("pp.xmp");
#endif
   fft_img.img_bandpass_filter(.1,0.);
#ifdef NEVERDEFINE
    std::complex<double> * IMG;
    IMG = (std::complex<double> *)fft_img.fOut;
    
    int ysize = (int)(((double)YSIZE(img()) / 2.) +1.);
    int xsize = XSIZE(img());
    
    img().initZeros();
    long int ii=0;
    //xmipp y contigous, fftw x contigous
    for (int j=0;j<xsize;j++)
        for (int i=0;i<ysize;i++)
        {
           DIRECT_MAT_ELEM(img(),i,j) = abs(IMG[ii]);
           ii++;
        }
    img.write("pp.xmp");
#endif
    #ifdef  DEBUGTIME
    printf ("image filtered in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    
    double * tmp;
    tmp = fft_img.fIn;
    fft_img.fIn=fft_img.fOut;
    fft_img.fOut=tmp;
    
    fft_img.Init("ES",FFTW_BACKWARD);
    #ifdef  DEBUGTIME
    printf ("backwards init in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    
    fft_img.Transform();
    #ifdef  DEBUGTIME
    printf ("image backtransformed in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
/*
    double * IMG;
    IMG = fft_img.fOut;

    int xsize = (int)(((double)XSIZE(img()) / 2.) +1.);
    int ysize = YSIZE(img());
    
    img().initZeros();
    long int ii=0;
    for (int i=0;i<ysize;i++)
        for (int j=0;j<xsize;j++)
        {
           DIRECT_MAT_ELEM(img(),i,j) = abs(IMG[ii]);
           ii++;
        }
*/
    img.write("kk.xmp");
    fft_img.fOut=NULL;//do not free, matrix2d destructor will do it
/**
For processing a second image read again and call to Transform, no need to create another fourier
object
*/
}
