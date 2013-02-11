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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
/* This file contains functions related to the SAMPLING Transform */

#ifndef _SAMPLING1_HH
#define _SAMPLING1_HH
#include <vector>
#include <iterator>

#include "matrix1d.h"
#include "xmipp_filename.h"
#include "metadata.h"
#include "symmetries.h"

#define cte_w 1.107149
/**@defgroup SphereSampling sampling (Sampling the projection sphere)
   @ingroup DataLibrary */
//@{
/** Routines with sampling the direction Sphere
    A triangular grid based on an icosahedron was first introduced in a
    meteorological model by Sadourny et al. (1968) and Williamson (1969). The
    approach outlined here, especially the code implementation, is based on the work
    of Baumgardner (1995).
*/
class Sampling
{
public:
    /** Type to define a vertex */
    struct Vertex
    {
        double rot;
        double tilt;
        double psi;
    };

    /** Typename to contain a list of vertices */
    typedef std::vector<Vertex> Vect_angles;

    /** Geographical co-ordinates of the home Vertex of the 10 diamonds
     as angles*/
    Vect_angles vertices_angles;

    /** Geographical co-ordinates of the home Vertex of the 10 diamonds
    as vectors*/
    std::vector <Matrix1D<double> > vertices_vectors;

    /** sampling rate in radians */
    double sampling_rate_rad;

    /** sampling rate for the unit vectors */
    double sampling_noise;

    /** number of samples */
    size_t number_of_samples;

    /** number of samples in the assymetric unit*/
    size_t numberSamplesAsymmetricUnit;

    /** neighborhood  in radians */
    double neighborhood_radius_rad;

    /** cosine of neighborhood  s */
    double cos_neighborhood_radius;

    /** vector with neighbors */
    std::vector<std::vector<size_t> >  my_neighbors;

    /** vector with experimental images per sampling point */
    std::vector<std::vector<size_t> >  my_exp_img_per_sampling_point;

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
    /** vector with product of experimental images and L and R */
    std::vector <FileName > exp_data_fileNames;
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
    /** vector with the indexes of each sampling point */
    std::vector <size_t> no_redundant_sampling_points_index;

    /** Verbose */
    int verbose;

    /** Default constructor. sampling in degrees*/
    Sampling();

    /** 'is equal to' (equality).*/
    bool operator==(const Sampling& op) const;

    /** symmetry file */
    FileName symmetry_file;

    /** symmetry information **/
    SymList  SL;

    /** Compute edge sampling points
        if you are looking only for directtions set only_half_sphere = true
    */
    void computeSamplingPoints(bool only_half_sphere = true,
                                 double max_tilt= +91.,
                                 double min_tilt= -91.);
    /** fill edge */
    void fillEdge(Matrix1D<double> starting_point,
                   Matrix1D<double> ending_point,
                   std::vector <Matrix1D<double> > &edge_vector,
                   bool FLAG_END
                  );
    /** fill distance */
    void fillDistance(Matrix1D<double> starting_point,
                       Matrix1D<double> ending_point,
                       std::vector <Matrix1D<double> > &edge_vector,
                       int number,
                       bool only_half_spheree,
                       double min_z= -10.,
                       double max_z= +10.
                      );
    /** fill R and L Repository (vector with symmetry matrices) */
    void fillLRRepository(void);

    /** set sampling rate */
    void setSampling(double sampling);

    /** set sampling noise for projection vectors create in the unit
        sphere */
    void setNoise(double deviation, int my_seed=-1);

    /** set neighborhood distance */
    void setNeighborhoodRadius(double neighborhood);

    /* eliminate redundant points,
        symmetry group, symmetry order */
    void removeRedundantPoints(const int symmetry, int sym_order);

    /* eliminate redundant points,
        symmetry group, symmetry order 
        This function first calls removeRedundantPoints,
        and then checks each point versus all others to calculate an angular distance
        If this distance is less than 0.8 times the angular sampling, the point is deleted
    */
    void removeRedundantPointsExhaustive(const int symmetry,
                                            int sym_order,
                                            bool only_half_sphere,
                                            double max_ang);

    /* sorting criteria for euler angles */
    int sortFunc(Matrix1D<double> & a, Matrix1D<double> & b);

    /** create symmetry file from introduced symmetry
        see  SymList class */
    void createSymFile(FileName simFp,int  symmetry, int sym_order);

    /** save assymetric unit sampling in a doc file */
    void createAsymUnitFile(const FileName& docfilename);

    /** for each point i in the assymetric sampling unit cell
    compute the neighbors inside the assymetric unit cell,
    save not only the neighbors but the angle psi
    */
    void computeNeighbors(bool only_winner=false);

    /** Save neighbors as threee metadtada blocks
     * 1) header with sampling rate and angular distance
     * 2) one row for ach experimentald data with neighbours
     * 3) sampling points
    */
    void saveSamplingFile(const FileName &fn_base, bool write_vectors = true, bool write_sampling_sphere = false);

    /** Read neighbors i
    */
    void readSamplingFile(const FileName &infilename,bool read_vectors=true, bool write_sampling_sphere = false);

    /** remove all those points that are further away from experimental data
        than neighborhood_radius_rad */
    void removePointsFarAwayFromExperimentalData();

    /** Find the closest sampling point for a docfile of experimental projections*/
    void findClosestSamplingPoint(MetaData &DFi,
                                     const FileName &output_file_root);

    /** Find the closest sampling point for a docfile of experimental projections*/
    void findClosestSamplingPoint(const FileName &FnexperimentalImages,
                                     const FileName &output_file_root);

    /**for each sampling point find the experimental images
       closer to that point than to any other */
    void findClosestExperimentalPoint();

    /** Precalculate exp_data by symmetry matrices (speeds up calculations)*/
    void fillExpDataProjectionDirectionByLR(MetaData &DFi);

    /** Precalculate exp_data by symmetry matrices (speeds up calculations)*/
    void fillExpDataProjectionDirectionByLR(const FileName &FnexperimentalImages);
};
//@}
#endif
