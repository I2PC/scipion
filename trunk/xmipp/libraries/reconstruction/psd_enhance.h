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
#ifndef _PROG_ENHANCE_PSD_HH
#  define _PROG_ENHANCE_PSD_HH

#include <data/progs.h>

/**@name Enhance PSD program */
//@{
/* Enhance PSD Program Parameters ------------------------------------------ */
/** Parameter class for the project program */
class Prog_Enhance_PSD_Parameters: public Prog_parameters {
public:
   /// Center PSD before working
   bool center;

   /// Take log10 before working
   bool take_log;

   /// Bandpass filter low frequency (in Fourier space, max 0.5)
   double filter_w1;

   /// Bandpass filter high frequency (in Fourier space, max 0.5)
   double filter_w2;

   /// Decay width (raised cosine)
   double decay_width;

   /// Lower frequency for the mask (in Fourier space, max 0.5)
   double mask_w1;

   /// Higher frequency for the mask (in Fourier space, max 0.5)
   double mask_w2;
public:
   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message. */
   void read(int argc, char **argv);

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage();

   /** Show parameters. */
   void show();

   /** Apply to a single PSD.
       The steps are basically: outlier removal, band pass filtration, masking
       and normalization. */
   void apply(matrix2D<double> &PSD);
};
//@}
#endif
