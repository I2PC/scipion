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

#ifndef _XMIPP_ROTATIONALSPECTRUM_HH
#  define _XMIPP_ROTATIONALSPECTRUM_HH

#include "matrix2d.h"

/**@name Rotational spectrum */
//@{

/** Cylindrical_Wave_Decomposition class. */
class Cylindrical_Wave_Decomposition {
public:
   /// In: Minimum harmonics
   int numin;
   /// In: Maximum harmonics
   int numax;
   /// In: Center of symmetry (x)
   double x0;
   /// In: Center of symmetry (x)
   double y0;
   /// In: Minimum integration radius
   double r1;
   /// In: Maximum integration radius
   double r2;
   /// In: Integration increment
   double r3;
   /// Out: Ir
   int ir;
   /// Out: Ampcos
   matrix1D<double> out_ampcos;
   /// Out: Ampsin
   matrix1D<double> out_ampsin;

   /// Show this object
   friend ostream & operator << (ostream &_out,
      const Cylindrical_Wave_Decomposition &_cwd);

   /// Interpolate image value (bilinear)
   double interpolate(matrix2D<double> &img, double y,double x);

   /// Compute the Cylindrical Wave decomposition of an image
   void compute_cwd(matrix2D<double> &img);
};

/** Rotational spectrum.
    Example of use:
    \begin{verbatim}
	 int main(int argc, char **argv) {
	    ImageXmipp I(argv[1]);

	    int rl=0;
	    int rh=22;
	    int dr=1;

	    Rotational_Spectrum spt;
	    spt.rl=rl;
	    spt.rh=rh;
	    spt.dr=dr;
	    spt.numin=1;
	    spt.numax=15;
	    spt.x0=(double)XSIZE(I())/2;
	    spt.y0=(double)YSIZE(I())/2;
	    spt.compute_rotational_spectrum (I(),rl,rh,dr,rh-rl);
	    cout << spt.rot_spectrum << endl;
	    return 0;
	 }
    \end{verbatim}
*/
class Rotational_Spectrum {
public:
   /** In: Ir.
       This value is given by the Cylindrical Wave Decomposition.
       It is not necessary to fill this field if the rotational spectrum
       is computed on an image.*/
   int ir;
   /// In: Minimum harmonics
   int numin;
   /// In: Maximum harmonics
   int numax;
   /// In: Center of symmetry (x)
   double x0;
   /// In: Center of symmetry (y)
   double y0;
   /// In: Minimum integration radius
   double rl;
   /// In: Maximum integration radius
   double rh;
   /// In: Integration increment
   double dr;
   /// Out: Rotational spectrum
   matrix1D<double> rot_spectrum;
   /// Show
   friend ostream & operator << (ostream &_out,
      const Rotational_Spectrum &_spt);
   /** Compute rotational spectrum of an image.
       xr1 is the minimum integration radius. xr2 the maximum integration
       radius. xdr the increment, and xr the length of integration.
       Usually, xr=xr2-xr1. */
   void compute_rotational_spectrum (matrix2D<double> &img,
      double xr1, double xr2, double xdr, double xr);
   /// Compute rotational spectrum using CWD
   void compute_rotational_spectrum (Cylindrical_Wave_Decomposition &cwd,
      double xr1, double xr2, double xdr, double xr);
   /// Read parameters from command line
   void read(int argc, char **argv);
   /// Usage
   void usage();
};
//@}
#endif
