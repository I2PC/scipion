/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.csic.es (2008)
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

#ifndef _ANGULAR_CLASS_AVERAGE_H
#define _ANGULAR_CLASS_AVERAGE_H

#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/filters.h>
#include <data/mask.h>

/**@defgroup ClassAverage Create class averages from projection
   matching docfiles
   @ingroup ReconsLibraryPrograms */
//@{
/** angular_class_average parameters. */
class Prog_angular_class_average_prm {

public:

    /** Input and library docfiles */
    DocFile          DF, DFlib;
    /** Output rootnames */
    FileName         fn_out, fn_out1, fn_out2;
    /** Column numbers */
    int              col_rot, col_tilt, col_psi, col_xshift, col_yshift, col_mirror, col_select;          
    /** Upper and lower selection limits */
    double           limit0, limitF;
    /** Flags wether to use limit0 and limitF selection */
    bool             do_limit0, do_limitF;
    /** Flag whether to apply mirror operations */
    bool             do_mirrors;
    /** Flag whether also to write out class averages of random halves of the data */
    bool             do_split;
    /** One empty image with correct dimensions */
    ImageXmipp       Iempty;

public:
  /// Read arguments from command line
  void read(int argc, char **argv);

  /// Show
  void show();

  /// Usage
  void usage();

  /** Make shiftmask and calculate nr_psi */
  void produceSideInfo();

  /** Process a single class */
  void processOneClass(int &dirno, 
		       double &lib_rot, 
		       double &lib_tilt);

};				
//@}
#endif
