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
#include <iterator>

#define cte_w 1.107149
/**@defgroup SphereSampling sampling (Sampling the projection sphere)
   @ingroup ReconsLibraryPrograms */
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
    typedef std::vector<vertices> Vect_angles;
    /** Geographical co-ordinates of the home vertices of the 10 diamonds
     as angles*/
    Vect_angles vertices_angles;

    /** Geographical co-ordinates of the home vertices of the 10 diamonds
    as angles*/
    std::vector <Matrix1D<double> > vertices_vectors;


    /** sampling rate in radians */
    double sampling_rate_rad;

    /** sampling rate for the unit vectors */
    double sampling_noise;

    /** number of samples */
    int number_of_samples;

    /** neighborhood  in radians */
    double neighborhood_radius_rad;

    /** cosine of neighborhood  s */
    double cos_neighborhood_radius;

    /** vector with neighbors */
    std::vector<std::vector<int> >  my_neighbors;
    /** vector with experimental images per sampling point */
    std::vector<std::vector<int> >  my_exp_img_per_sampling_point;
#ifdef MYPSI    
    /** vector with neighbors psi*/
    std::vector<std::vector<double> > my_neighbors_psi;
#endif
    /** vector with angular distance between points (dot_product)*/
    std::vector<std::vector<double> > my_cross_correlation;

    /** vector with sampling points described by vectors */
    std::vector <Matrix1D<double> > sampling_points_vector;
    /** vector with sampling points described by angles */
    std::vector <Matrix1D<double> > sampling_points_angles;
    /** vector with symmetry matrices */
    std::vector <Matrix2D<double> > R_repository;
    std::vector <Matrix2D<double> > L_repository;
    /** vector with product of experimental images and L and R */
    std::vector <Matrix1D<double> > exp_data_projection_direction_by_L_R;
#ifdef MYPSI    
    /** vector with psi resulting from product of experimental images and L and R */
    std::vector <double > exp_data_projection_direction_by_L_R_psi;
#endif
    
    /** vector with sampling points described by vectors, only store
        the non redundant part */
    std::vector <Matrix1D<double> > no_redundant_sampling_points_vector;
    /** vector with sampling points described by angles, only store
        the non redundant part */
    std::vector <Matrix1D<double> > no_redundant_sampling_points_angles;

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
                   std::vector <Matrix1D<double> > &edge_vector,
                   bool FLAG_END
                  );
    /** fill distance */
    void fill_distance(Matrix1D<double> starting_point,
                       Matrix1D<double> ending_point,
                       std::vector <Matrix1D<double> > &edge_vector,
                       int number,
                       bool only_half_spheree,
                       double min_z= -10.,
                       double max_z= +10.
                      );
    /** fill R and L Repository (vector with symmetry matrices) */
    void fill_L_R_repository(void);
                      
    /** set sampling rate */
    void SetSampling(double sampling);

    /** set sampling noise for projection vectors create in the unit
        sphere */
    void SetNoise(double deviation, int my_seed=-1);

    /** set neighborhood distance */
    void SetNeighborhoodRadius(double neighborhood);

    /* eliminate redundant points,
        symmetry group, symmetry order */
    void remove_redundant_points(const int symmetry, int sym_order);

    /* eliminate redundant points,
        symmetry group, symmetry order 
        This function first calls remove_redundant_points,
        and then checks each point versus all others to calculate an angular distance
        If this distance is less than 0.8 times the angular sampling, the point is deleted
    */
    void remove_redundant_points_exhaustive(const int symmetry, 
                                            int sym_order, 
                                            bool only_half_sphere, 
                                            double max_ang);

    /// Usage
    //void Usage();

    /* sorting criteria for euler angles */
    int sort_func(Matrix1D<double> & a, Matrix1D<double> & b);

    /** create symmetry file from introduced symmetry
        see  SymList class */
    void create_sym_file(FileName simFp,int  symmetry, int sym_order);

    /** save assymetric unit sampling in a doc file */
    void create_asym_unit_file(const FileName& docfilename);

    /** for each point i in the assymetric sampling unit cell
    compute the neighbors inside the assymetric unit cell,
    save not only the neighbors but the angle psi
    */

    void compute_neighbors(bool only_winner=false);
   /** Save neighbors in a propietary ascii file. The structure is as
       follows 
      [vnum]
      [size1]
      [vec1_neighbors]
#ifdef MYPSI    
      [vec1_psi]
#endif
      [vec1_crosscorrelation]
      [size2]
      [vec2_neighbors]
#ifdef MYPSI    
      [vec2_psi]
#endif
      [vec2_crosscorrelation]
      ...
      [sizen]
      [vecn_neighbors]
#ifdef MYPSI    
      [vecn_psi]
#endif
      [vecn_crosscorrelation]
      
      for the neighbors and
      X1_angle  Y1_angle  Z1_angle
      X1_vector Y1_vector Z1_vector
      X2_angle  Y2_angle  Z2_angle
      X2_vector Y2_vector Z2_vector
      ...
      Xn_angle  Yn_angle  Zn_angle
      Xn_vector Yn_vector Zn_vector
      for the sampling points

      where vnum is the number of vectors, sizen is the number of elements in
      that vector and vecn is the elements. 
   
   */
   void save_sampling_file(FileName outfilename,bool write_vectors=true);
   /** Read neighbors in a propietary ascii file. The structure is as
       follows 
      [vnum]
      [size1]
      [vec1]
      [size2]
      [vec2]
      ...
      [sizen]
      [vecn]
      
      for the neighbors and
      X1 Y1 Z1
      X2 Y2 Z2
      ...
      Xn Yn Zn
      for the sampling points

      where vnum is the number of vectors, sizen is the number of elements in
     that vector and vecn is the elements. 
   
   */
   void read_sampling_file(FileName infilename,bool read_vectors=true);
   
   /** remove all those points that are further away from experimental data
       than neighborhood_radius_rad */
   void remove_points_far_away_from_experimental_data();

   /** Find the closest sampling point for a docfile of experimental projections*/   
   void find_closest_sampling_point(DocFile &DFi,
                                    FileName output_file_root);
   void find_closest_sampling_point(FileName FnexperimentalImages,
                                    FileName output_file_root);

   /**for each sampling point find the experimental images
      closer to that point than to any other */
   void  find_closest_experimental_point();

   /** Precalculate exp_data by symmetry matrices (speeds up calculations)*/
   void fill_exp_data_projection_direction_by_L_R(DocFile &DFi);
   void fill_exp_data_projection_direction_by_L_R(FileName FnexperimentalImages);
};
//@}
#endif
