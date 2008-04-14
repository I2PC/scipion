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
    ImageXmipp img;
    img.read("proj00001.xmp");

    #ifdef  DEBUGTIME
    printf ("image read in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    
    #define DIMENSION 2
    xmippFftw fft_img(img(), false, true);
    fft_img.Init("ES",FFTW_FORWARD,false);
    #ifdef  DEBUGTIME
    printf ("Plan created in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
    fft_img.CenterRealDataBeforeTransform();

    fft_img.Transform();
    #ifdef  DEBUGTIME
    printf ("image transformed in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif

    complex<double> * IMG;
    IMG = (complex<double> *)fft_img.fOut;

    double * img_ptr;
    img_ptr=MULTIDIM_ARRAY(img());
    int xsize = (int)(((double)XSIZE(img()) / 2.) +1.);
    int ysize = YSIZE(img());
    
    int ii=0;
    for(int i=0; i<img().xdim; i++)
       for(int j=0; j<img().ydim; j++)
          img_ptr[ii++]=0.;

    ii=0;
    for (int i=0;i<ysize;i++)
        for (int j=0;j<xsize;j++)
           {
           DIRECT_MAT_ELEM(img(),i,j) = abs(IMG[ii]);
           ii++;
           }
    #ifdef  DEBUGTIME
    printf ("computed abs value in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif
       
    img.write("proj00001_fftwmod.xmp");
    #ifdef  DEBUGTIME
    printf ("write in  %f mseconds\n", (double)(clock()-clo)/1000 );
    clo = clock();
    #endif

/**
For processing a second image read again and call to Transform, no need to create another fourier
object
*/
}
