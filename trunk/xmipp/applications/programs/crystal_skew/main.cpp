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
// This program "unskrew" a volume produced by the MRC package
// (An converted to Spider format)

#include <data/volume.h>
#include <data/image.h>
#include <data/args.h>

#include <cstdio>

struct TransfMat {
    float Sx,           /* Scale: Scale factor for X axis. */
          Hy,    	/* Shear: transform along X axis. */
          Tx;    	/* Translation: Offset in X axis. */
    };                  /* Structure for Transformation Matrix. */

void Usage(char *argv[]);
void CompTransfMatrix(struct TransfMat *trmatrix, float gamma, int ny);
float Interpo(float x, int x1, int x2, float y1, float y2);

int main (int argc, char **argv) {
   FileName		fn_out, fn_input;
   float 		gamma;
   int 			minz, maxz;
   VolumeXmipp     	volume_in;
   VolumeXmipp     	volume_out;
   int 			nx, ny, nz, newnx;
   struct TransfMat 	transformatrix;      /* Transformation Matrix. */
   float            	sheartransY;         /* Shear and Translation factors
                                                for row <y>. */
   int 			skewnx;              /* Nx of skewed & scaled image. */
   int			z,y;
   int			i;

   // Get command line parameters ------------------------------------------
   try {
      fn_input     = get_param(argc,argv,"-i","");
      fn_out       = get_param(argc,argv,"-o","");
      gamma        = AtoF(get_param(argc,argv,"-gamma"));
      minz        = AtoI(get_param(argc,argv,"-zmin"));
      maxz        = AtoI(get_param(argc,argv,"-zmax"));
   } catch (Xmipp_error XE) {Usage(argv);}
   try {
   // get volume dimensions ---------------------------------------------------
   if (Is_ImageXmipp(fn_input)) {
       cout << "\nBy the moment I only process volumes"<<endl;
       exit(1);
       }
   else if (Is_VolumeXmipp(fn_input)) {
      cout << "\n Reading file: " << fn_input << endl;
      volume_in.read(fn_input);
      }
   else {
      cout << "\n Can not read file: " << fn_input << endl;
   }

   /* Checking out the value of gamma. */
   if(gamma == 90){
       printf("\n Skew angle = 90 degrees. No skewing applied.\n\n");
       exit(1);
       }
   else if(gamma >= 180 || gamma <= 0){
       printf("\n Skew angle must be in ]0,180[.\n\n");
       exit(1);
       }

   // Computing Transformation Matrix.
   nz=ZSIZE(VOLMATRIX(volume_in));
   ny=YSIZE(VOLMATRIX(volume_in));
   nx=XSIZE(VOLMATRIX(volume_in));

   CompTransfMatrix(&transformatrix, gamma, ny);

   /* Checking out the zmax/zmin. */
   if( minz<0 || minz >= maxz || maxz>nz ){
       printf("\n minz or maxz values are wrong. No skewing applied.\n\n");
       exit(1);
       }

   // Computing Dimensions of new image.
   newnx = (int)(nx * transformatrix.Sx + ny * fabsf(transformatrix.Hy));
   if(newnx%2) newnx++;                  /* Even number of columns is wanted. */

   // output volume
   VOLMATRIX(volume_out).resize(maxz-minz+1,ny,newnx);
   cout << " Skewing file with gama= " << gamma << endl;

   for(z=minz;z<=(maxz); z++){
      for(y=0,sheartransY = transformatrix.Tx,
          skewnx=(int) floor(nx *transformatrix.Sx);
                  y<ny; y++, sheartransY += transformatrix.Hy)
	  {/*skew line */
	     register int skewx, xa, xb;
	     int skewx0, skewx1, skewxn;
	     float x, *ptrline;
	     float Sx=transformatrix.Sx;

	     /* Limits in output line. */
	     skewx1 = (int)ceil(sheartransY);
	     skewxn = (int)floor((nx - 1) * Sx + sheartransY);
	     /* Skewing line. */
	     for(skewx=skewx1; skewx<=skewxn; skewx++){
		 x = (skewx - sheartransY) / Sx;
		 xa = (int)floor(x); xb = (int)ceil(x);
		 if(xa == xb)
        	     VOLVOXEL(volume_out,z-minz,y,skewx) = VOLVOXEL(volume_in,z,y,xa);
		 else
        	     VOLVOXEL(volume_out,z-minz,y,skewx) = Interpo(x, xa, xb,
	                                        	     VOLVOXEL(volume_in,z,y,xa),
							     VOLVOXEL(volume_in,z,y,xb));
		 }

	     /* Wrap around. */
	     skewx0 = (int)ceil(sheartransY - Sx);
	     if(skewx0 >= 0) {
		for(skewx=skewx0; skewx<skewx1; skewx++){
		   x = (skewx - sheartransY) / Sx;
		   xa = (int)floor(x); xb = (int)ceil(x);
		   if(xa == xb)
		     VOLVOXEL(volume_out,z-minz,y,skewx) = VOLVOXEL(volume_in,z,y,xa);
		   else
		     VOLVOXEL(volume_out,z-minz,y,skewx) = Interpo(x, xa, xb,
	  					      VOLVOXEL(volume_in,z,y,nx-1),
	        				      VOLVOXEL(volume_in,z,y,0));
		}
		if(skewx0)
		   for(i=0;i<skewx0;i++)
	              VOLVOXEL(volume_out,z-minz,y,i)=
			   VOLVOXEL(volume_out,z-minz,y,i+skewnx);
	     }

	     if(skewxn + 1 < newnx){
  		for(i=0;i<(newnx - skewxn - 1);i++)
                   VOLVOXEL(volume_out,z-minz,y,i+skewxn + 1)=
          		   VOLVOXEL(volume_out,z-minz,y,i+(skewxn + 1 - skewnx));

	     }

	  }/*SkewLine*/
    }/* for z end */
      cout << "\n Writting file: " << fn_out <<" ..."<< endl;
      volume_out.write(fn_out);
      cout << "\n End writting file: " << fn_out <<" ..."<< endl;

   } catch (Xmipp_error XE) {cout << XE;}
   exit(0);
}
/*----------------------------------------------------------------------------*/

void Usage (char *argv[]) {
    cout << "Purpose:\n";
    cout << "    Skew (undeform) a volume\n";

    cout << "Usage:" << argv[0] <<" -i filename -o filename -gamma skew_amgle_degrees -zmin Zmin -zmax Zmax" << endl << endl;
    cout << endl;

    exit(1);

}
/*----------------------------------------------------------------------------*/

void CompTransfMatrix(struct TransfMat *trmatrix, float gamma, int ny)
{
/* Transformation matrix. */
trmatrix->Hy = (float) tan(PI/2 - DEG2RAD(gamma));  /* Hy = tan(90 - gamma). */
trmatrix->Sx = (float) sqrt(1.0 + trmatrix->Hy * trmatrix->Hy);
trmatrix->Tx = (gamma < 90.0)? 0.0 : - ny * trmatrix->Hy;
}
/*----------------------------------------------------------------------------*/

float Interpo(float x, int x1, int x2, float y1, float y2)
{
return (y1 + (x - x1) * (y2 - y1)/(x2 - x1));
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Skew {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Skew/Help/skew.html";
      help="Applies a skew to the slected Z planes in a volume";
      OPEN MENU menu_skew;
      COMMAND LINES {
         +usual: xmipp_skew -i $FILE_IN -o $FILE_OUT
                    -gamma $GAMMA -zmin $ZMIN -zmax $ZMAX
      }
      PARAMETER DEFINITIONS {
         $FILE_IN  {label="Input volume"; type=file existing;}
         $FILE_OUT {label="Output volume"; type=file; by default=$FILE_IN;}
         $GAMMA {label="Skewing angle"; type=float;}
         $ZMIN {label="Minimum Z"; type=integer;}
         $ZMAX {label="Maximum Z"; type=integer;}
      }
   }

   MENU menu_skew {
      "I/O Parameters"
      $FILE_IN
      $FILE_OUT
      "Skewing parameters"
      $GAMMA
      $ZMIN
      $ZMAX
   }
*/
