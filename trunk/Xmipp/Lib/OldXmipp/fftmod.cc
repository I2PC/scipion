/***************************************************************************
 *
 * Authors:     Irene Martinez
 *              Roberto Marabini
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


  /* **********************************************************************

        This file contains the routine to calculate the fast Fourier 
        transform magnitude of an image. The image is in a "rows x cols"
	matrix of floats, and it will be overwritten by the result of:

        final_pixel = log( || F.F.T.( ini_pixel ) || +1 )

        Where :
                final_pixel    : It's the final pixel value.

                ini_pixel      : It's the initial pixel value.

                log            : It's the natural logarithm.

                F.F.T          : It's the fast Fourier transform .

                || F.F.T.() || : The fast Fourier transform magnitude.

   ************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <spider.h>
#include <groe.h>

#define Re(i,j)  fImg[i][2*(j)]
#define Im(i,j)  fImg[i][2*(j)+1]

/**************************************************************************/

int fast_fourier_trans_img(float **imgModulo,int rows,int cols)
{
int i,j, cols2;
float **fImg;
float par=-(1.);

 if ( (fImg=(float **)imalloc(rows,cols,FOURIER)) ==NULL)
  {        
     fprintf(stderr,"\nCan't allocate memory \n");
     return(0);
  }

 for( i=0; i <rows; i++)
  {                          /* multiplico por 1 y -1 alternativamente */
    for( j=0; j<cols ; j++)
     {
        par *= -1.0;
        fImg[i][j] = par * imgModulo[i][j];
     }
    par *= -1.0;
  }

                 /*  calculo de los valores de la transfomada de fourier  */

 if ( image_fft(fImg,rows,cols,DIRECT) == ERROR)
  {
     fprintf(stderr,"\n FFT coefficients error");
     return(0);
  }
            /* Calculo del modulo y la fase de  la transformada    */

 cols2=(int) (cols/2);

 for (i=0; i<rows ; i++)
    for (j=0; j<cols;  j++)
      {
        if (j<=cols2 )
            imgModulo[i][j] = fmodulo( Re(i,j),Im(i,j));
        else  if (i==0)
            imgModulo[0][j] = fmodulo(Re(0,cols-j),Im(0,cols-j));
        else
            imgModulo[i][j] = fmodulo(Re(rows-i,cols-j),Im(rows-i,cols-j));
      }

 return(1);
}

/****************************************************************************/

float fmodulo(float x,float y)
{
  double res;

  res = log ( (double)(1.0 + sqrt(x*x+y*y) ));
  return( (float) res);

}  /* end of fmodulo   */

/****************************************************************************/
/****************************************************************************/



