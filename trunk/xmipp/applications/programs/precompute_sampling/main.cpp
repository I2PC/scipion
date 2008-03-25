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
#include <data/error.h>
#include <iostream>

/* ------------------------------------------------------------------------- */
/* Program                                                                   */
/* ------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    int test_case=0;
    test_case = textToInteger(getParameter(argc, argv, "-test_case", "0"));
    
    XmippSampling mysampling;
    mysampling.read_sampling_file("img");
    std::ofstream filestr; 
    filestr.open ("precompute.bild");

    // big sphere in the center
    filestr    << ".color white" 
	   << std::endl
	   << ".sphere 0 0 0 .95"
	   << std::endl
	   ;
    // initial point   
    filestr    << ".color yellow" 
	   << std::endl
	   << ".sphere "
       << mysampling.no_redundant_sampling_points_vector[test_case].transpose()
       << " .021"
	   << std::endl
	   ;

    for (int i = 0; i < mysampling.my_neighbors[test_case].size(); i++)
        {
        if(mysampling.my_neighbors_psi[test_case][i]==0)
            filestr    << ".color red" 
	                   << std::endl
	                   << ".sphere "
                       << mysampling.no_redundant_sampling_points_vector[
                                    mysampling.my_neighbors[test_case][i]].transpose()
                       << " .02"
	                   << std::endl
	                   ;
        else           
            filestr    << ".color green" 
	                   << std::endl
	                   << ".sphere "
                       << mysampling.no_redundant_sampling_points_vector[
                                    mysampling.my_neighbors[test_case][i]].transpose()
                       << " .02"
	                   << std::endl
	                   ;
        }
    filestr.close();
}
