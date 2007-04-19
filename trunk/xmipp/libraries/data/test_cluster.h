/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2004)
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


#ifndef _PROG_TEST_CLUSTER
   #define _PROG_TEST_CLUSTER

#include "funcs.h"
#include "mask.h"
#include "selfile.h"
#include "matrix2d.h"

/**@name Test cluster */
//@{
/** Test cluster parameters */
class Test_cluster_parameters {
public:
   typedef enum {EUCLIDEAN,
         CORRELATION,
	 MAHALANOBIS} Distance_type;

   /// Selfile with the cluster images
   FileName fn_selfile;

   /// Filename of the output histogram
   FileName fn_out;

   /// Mask to apply, maybe none
   Mask_Params mask;

   /** Distance.
       Valid distances are EUCLIDEAN, CORRELATION, MAHALANOBIS. */
   Distance_type distance;

public:
   // Selfile with all the cluster
   SelFile SF_in;

   // Covariance matrix
   matrix2D<double> covariance;

   // Mean vector
   matrix1D<double> mean;

public:
   /// Empty constructor
   Test_cluster_parameters();

   /// Read parameters from command line
   void read(int argc, char **argv);

   /** Produce side info.
       The mask is created.*/
   void produce_side_info();

   /// Show parameters. This function calls show_specific
   void show();

   /// Usage. This function calls usage_specific
   void usage();

   /// Build covariance matrix
   void build_covariance_matrix();

   /// Effectively do the work
   void test_cluster();
};
//@}
#endif
