/***************************************************************************
 *
 * Authors:     Roberto Marabini
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
/* This file contains functions related to the SAMPLING Transform */

#ifndef _SAMPLING1_HH
#define _SAMPLING1_HH
#include <vector>
#include <data/docfile.h>
#include <data/macros.h>
#include <data/matrix1d.h>
#include <reconstruction/symmetries.h>
#include <data/geometry.h>

#define cte_w 1.107149
/**@name Sampling */
//@{
/** Routines with sampling the direction Sphere

       A triangular grid based on an icosahedron was first introduced in a
meteorological model by Sadourny et al. (1968) and Williamson (1969). The
approach outlined here, especially the code implementation, is based on the work
of Baumgardner (1995).
 */

class XmippSampling
{
public:
    struct  vertices
    {
        double rot;
        double tilt;
        double psi;
    };
    typedef vector<vertices> Vect_angles;
    /** Geographical co-ordinates of the home vertices of the 10 diamonds
     as angles*/
    Vect_angles vertices_angles;

    /** Geographical co-ordinates of the home vertices of the 10 diamonds
    as angles*/
    vector <Matrix1D<double> > vertices_vectors;


    /** sampling rate in radians */
    double sampling_rate_rad;

    /** number of samples */
    int number_of_samples;

    /** neighborhood  in radians */
    double neighborhood_radius_rad;

    /** cosine of neighborhood  s */
    double cos_neighborhood_radius;

    /** vector with neighbors */
    vector<vector<int> >  my_neighbors;
    /** vector with neighbors psi*/
    vector<vector<double> > my_neighbors_psi;

    /** vector with sampling points described by vectors */
    vector <Matrix1D<double> > sampling_points_vector;
    /** vector with sampling points described by angles */
    vector <Matrix1D<double> > sampling_points_angles;

    /** vector with sampling points described by vectors, only store
        the non redundant part */
    vector <Matrix1D<double> > no_redundant_sampling_points_vector;
    /** vector with sampling points described by vectors, only store
        the non redundant part */
    vector <Matrix1D<double> > no_redundant_sampling_points_angles;

    /** Default constructor. sampling in degrees*/
    XmippSampling();

    /** symmetry file */
    FileName symmetry_file;

    /** symmetry information **/
    SymList  SL;

    /* Docfile with proyection directions as coordinates */
    DocFile           DFvectors;
    
    /* Docfile with proyection directions as euler angles */
    DocFile           DFangles;

    /** Compuute edge sampling points
        if you are looking only for directtions set only_half_sphere = true
    */
    void Compute_sampling_points(bool only_half_sphere = true,
                                 double max_tilt= +91.,
                                 double min_tilt= -91.);
    /** fill edge */
    void fill_edge(Matrix1D<double> starting_point,
                   Matrix1D<double> ending_point,
                   vector <Matrix1D<double> > &edge_vector,
                   bool FLAG_END
                  );
    /** fill distance */
    void fill_distance(Matrix1D<double> starting_point,
                       Matrix1D<double> ending_point,
                       vector <Matrix1D<double> > &edge_vector,
                       int number,
                       bool only_half_spheree,
                       double min_z= -10.,
                       double max_z= +10.
                      );
    /** set sampling rate */
    void SetSampling(double sampling);

    /** set neighborhood distance */
    void SetNeighborhoodRadius(double neighborhood);

    /* eliminate redundant points,
        symmetry group, symmetry order */

    void remove_redundant_points(string symmetry, int sym_order);

    /// Usage
    //void Usage();

    /* sorting criteria for euler angles */
    int sort_func(Matrix1D<double> & a, Matrix1D<double> & b);

    /** create symmetry file from introduced symmetry
        see  SymList class */
    void create_sym_file(string symmetry, int sym_order);

    /** save assymetric unit sampling in a doc file */
    void create_asym_unit_file(FileName docfilename);

    /** for each point i in the assymetric sampling unit cell
    compute the neighbors inside the assymetric unit cell,
    save not only the neighbors but the angle psi
    */

    void compute_neighbors(void);
};
//@}
#endif
