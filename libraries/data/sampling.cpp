/***************************************************************************
 *
 * Authors:    Roberto Marabini
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
#include "sampling.h"
#include "matrix2d.h"

/* Default Constructor */
Sampling::Sampling()
{
    Vertex aux;
    aux.rot =   -PI / 5.;
    aux.tilt = PI / 2. - cte_w  ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot =    PI / 5.;
    aux.tilt = PI / 2. - cte_w  ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = 3. * PI / 5.;
    aux.tilt = PI / 2. - cte_w  ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = 5. * PI / 5.;
    aux.tilt = PI / 2. - cte_w  ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = -3. * PI / 5.;
    aux.tilt = PI / 2. - cte_w  ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot =    0. / 5.;
    aux.tilt = -cte_w + PI / 2. ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = 2. * PI / 5.;
    aux.tilt = -cte_w + PI / 2. ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = 4. * PI / 5.;
    aux.tilt = -cte_w + PI / 2. ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = -4. * PI / 5.;
    aux.tilt = -cte_w + PI / 2. ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);
    aux.rot = -2. * PI / 5.;
    aux.tilt = -cte_w + PI / 2. ;
    aux.psi =  0.;
    vertices_angles.push_back(aux);

    vertices_vectors.push_back(vectorR3(0., 0., 1.));
    vertices_vectors.push_back(vectorR3(0.723606900230461, -0.525731185781806, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.723606900230461, 0.525731185781806, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.276393239417711, 0.850650928976665, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.8944273172062, 0., 0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.276393239417711, -0.850650928976665, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.8944273172062, 0., -0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.276393242471372, 0.850650927984471, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.723606898343194, 0.525731188379405, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.723606898343194, -0.525731188379405, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.276393242471372, -0.850650927984471, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(0., 0., -1.));

    sampling_noise=0.0;
    cos_neighborhood_radius=-1.01;

    numberSamplesAsymmetricUnit=-1;
    exp_data_fileNames.clear();

    verbose=1;
    //#define DEBUG1
#ifdef  DEBUG1

    for (int i = 0;
         i < vertices_vectors.size();
         i++)
        std::cout  <<  vertices_vectors[].transpose()  << std::endl;
#endif
#undef DEBUG1
}

bool Sampling::operator==(const Sampling& op) const
{

    return (XMIPP_EQUAL_REAL(sampling_rate_rad, op.sampling_rate_rad) &&
            XMIPP_EQUAL_REAL(cos_neighborhood_radius, op.cos_neighborhood_radius) &&
            no_redundant_sampling_points_angles == op.no_redundant_sampling_points_angles &&
            no_redundant_sampling_points_vector == op.no_redundant_sampling_points_vector &&
            no_redundant_sampling_points_index == op.no_redundant_sampling_points_index &&
            my_neighbors == op.my_neighbors //&&
            //exp_data_fileNames == op.exp_data_fileNames
           );

}

void Sampling::setSampling(double sampling)
{
    sampling_rate_rad = DEG2RAD(sampling);
    number_of_samples = ROUND(cte_w / sampling_rate_rad)+1;
    if (number_of_samples < 3)
    {
        std::cerr << "maximum value of angular sampling rate is "
        << cte_w*0.5*180./PI
        << std::endl;
        exit(1);
    }
}

void Sampling::setNoise(double noise_deviation, int my_seed)
{
    sampling_noise = noise_deviation;
    init_random_generator(my_seed);
}

void Sampling::setNeighborhoodRadius(double neighborhood)
{
    if(neighborhood<0)
        cos_neighborhood_radius=-1.01;
    else if(neighborhood>180.001)
        REPORT_ERROR(ERR_ARG_INCORRECT,"Neighborhood can not be greater than 180. \n \
                     Use any negative value to cover the whole space ");
    else
    {
        neighborhood_radius_rad = DEG2RAD(neighborhood);
        cos_neighborhood_radius = cos(neighborhood_radius_rad);
    }
}

/* Compute edge sampling points using Baumgardner  1995 */
void Sampling::computeSamplingPoints(bool only_half_sphere,
                                     double max_tilt,
                                     double min_tilt)
{
    /** vector to decimate the triangles */
    std::vector <Matrix1D<double> > edge_vector_start;
    /** vector to decimate the triangles */
    std::vector <Matrix1D<double> > edge_vector_end;
    // I need 10 auxiliary vector for edges
    Matrix1D<double> starting_point, ending_point;
    double max_z;
    double min_z;
    sampling_points_angles.clear();
    sampling_points_vector.clear();

    /* this is wrong, revert to previous code ROB
    if(min_tilt <= 0.)
        min_z= -1.;
    else
        min_z=-fabs(cos(PI * min_tilt / 180.));

    if(max_tilt >= 180.)
        max_z=1.;
    else if (max_tilt<=90)
    	max_z=-fabs(cos(PI * max_tilt / 180.));
    else
        max_z=fabs(cos(PI * max_tilt / 180.));
    */
    if(max_tilt >= 90.)
        max_z=10.;
    else
        max_z=sin(PI * max_tilt / 180.);
    if(min_tilt <= -90.)
        min_z= -10.;
    else
        min_z=sin(PI * min_tilt / 180.);

    //01a
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[1];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[6];
    ending_point = vertices_vectors[1];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //01b
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[2];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[6];
    ending_point = vertices_vectors[2];
    fillEdge(starting_point, ending_point, edge_vector_end, true);
    //02a
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[2];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[7];
    ending_point = vertices_vectors[2];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //02b
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[3];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[7];
    ending_point = vertices_vectors[3];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //03a
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[3];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[8];
    ending_point = vertices_vectors[3];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //03b
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[4];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[8];
    ending_point = vertices_vectors[4];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //04a
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[4];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[9];
    ending_point = vertices_vectors[4];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //04b
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[5];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[9];
    ending_point = vertices_vectors[5];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //05a
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[5];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[10];
    ending_point = vertices_vectors[5];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //05b
    starting_point = vertices_vectors[0];
    ending_point = vertices_vectors[1];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[10];
    ending_point = vertices_vectors[1];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //06a
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[10];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[5];
    ending_point = vertices_vectors[10];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //06b
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[9];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[5];
    ending_point = vertices_vectors[9];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //07a
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[9];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[4];
    ending_point = vertices_vectors[9];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //07b
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[8];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[4];
    ending_point = vertices_vectors[8];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //08a
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[8];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[3];
    ending_point = vertices_vectors[8];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //08b
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[7];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[3];
    ending_point = vertices_vectors[7];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //09a
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[7];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[2];
    ending_point = vertices_vectors[7];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //09b
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[6];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[2];
    ending_point = vertices_vectors[6];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //10a
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[6];
    fillEdge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[1];
    ending_point = vertices_vectors[6];
    fillEdge(starting_point, ending_point, edge_vector_start, true);
    //10b
    starting_point = vertices_vectors[11];
    ending_point = vertices_vectors[10];
    fillEdge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[1];
    ending_point = vertices_vectors[10];
    fillEdge(starting_point, ending_point, edge_vector_end, true);

    //#define DEBUG2
#ifdef  DEBUG2

    for (int i = 0;
         i < edge_vector_start.size();
         i++)
    {
        std::cout  <<  ".sphere " << edge_vector_start[i]  << " 1 1 " << std::endl;
        std::cout  <<  ".sphere " << edge_vector_end[i]  << " 1 2 " << std::endl;
    }
    //std::cout  <<  ending_point.transpose()    << " 1.1 1.5 " << std::endl;
#endif
#undef DEBUG2
    // add main corners that are not part of the basic octahedrons
    {
        int i=11;
        if (   (only_half_sphere && ZZ(vertices_vectors[i]) < 0.0)
               || ZZ(vertices_vectors[i]) < min_z
               || ZZ(vertices_vectors[i]) > max_z
           )
            //if (only_half_sphere && ZZ(vertices_vectors[i]) < 0.0)
            ;
        else
            sampling_points_vector.push_back(vertices_vectors[i]);
    }
    {
        int i=0;
        if (   (only_half_sphere && ZZ(vertices_vectors[i]) < 0.0)
               || ZZ(vertices_vectors[i]) < min_z
               || ZZ(vertices_vectors[i]) > max_z
           )
            //if (only_half_sphere && ZZ(vertices_vectors[i]) < 0.0)
            ;
        else
            sampling_points_vector.push_back(vertices_vectors[i]);
    }

    // add edges
    for (size_t i = 0;
         i < edge_vector_start.size();
         i++)
    {
        if (i < number_of_samples * 10 - 15)
        {
            if (   (only_half_sphere && ZZ(edge_vector_start[i]) < 0.0)
                   || ZZ(edge_vector_start[i]) < min_z
                   || ZZ(edge_vector_start[i]) > max_z
               )
                //if (only_half_sphere && ZZ(edge_vector_start[i]) < 0.0)
                continue;
            else
                sampling_points_vector.push_back(edge_vector_start[i]);
        }
        else
        {
            if (   (only_half_sphere && ZZ(edge_vector_end[i]) < 0.0)
                   || ZZ(edge_vector_end[i]) < min_z
                   || ZZ(edge_vector_end[i]) > max_z
               )
                //if (only_half_sphere && ZZ(edge_vector_end[i]) < 0.0)
                continue;
            else
                sampling_points_vector.push_back(edge_vector_end[i]);
        }
    }
    //#define DEBUG3
#ifdef  DEBUG3
    std::ofstream filestr;
    filestr.open ("debug3_1.bild");
    for (int i = 0;
         i < sampling_points_vector.size();
         i++)
    {
        filestr  <<  ".color 1 0 0" << std::endl
        <<  ".sphere " << sampling_points_vector[i] <<
        " .025" << std::endl;
        sampling_points_vector[i];//.add(0.0, sampling_noise, "gaussian");
        sampling_points_vector[i].selfNormalize();
        filestr  <<  ".color 0 0 1" << std::endl
        <<  ".sphere " << sampling_points_vector[i]  <<
        " .027" << std::endl;
    }
    filestr.close();

#endif
#undef DEBUG3

    // add in between points
    int j = 0;
    bool j_flag = false;
    for (size_t i = 0;
         i < edge_vector_start.size();
         i++)
    {
        if ((j % (number_of_samples - 1)) == 0 && j != 0)
        {
            j = 0;
            j_flag = true;
        }
        if ((j % (number_of_samples - 2)) == 0 && j != 0  && j_flag == true)
        {
            j = 0;
            j_flag = false;
        }
        fillDistance(edge_vector_start[i],
                     edge_vector_end[i],
                     sampling_points_vector,
                     (j + 1) % number_of_samples,
                     only_half_sphere,min_z,max_z);
        j++;
    }
    //#define DEBUG3
#ifdef  DEBUG3
    for (int i = 0;
         i < sampling_points_vector.size();
         i++)
    {
        std::cout  <<  ".color 1 0 0" << std::endl
        <<  ".sphere " << sampling_points_vector[i].transpose()  <<
        " .025" << std::endl;
    }
#endif

    //noisify angles
    if(sampling_noise!=0.0)
    {
        for (size_t n = 0;
             n < sampling_points_vector.size();
             n++)
        {
            FOR_ALL_ELEMENTS_IN_MATRIX1D(sampling_points_vector[n])
            (sampling_points_vector[n])(i) += rnd_gaus(0., sampling_noise);
            sampling_points_vector[n].selfNormalize();
        }
    }

    //#define DEBUG3
#ifdef  DEBUG3
    for (int i = 0;
         i < sampling_points_vector.size();
         i++)
    {
        std::cout  <<  ".color 0 1 0" << std::endl
        <<  ".sphere " << sampling_points_vector[i].transpose()  <<
        " .03" << std::endl;
    }
#endif
#undef DEBUG3

    // store sampling points as angles
    Matrix1D<double> aux(3), aux1(3);
    ZZ(aux) = 0.;
    for (size_t i = 0;
         i < sampling_points_vector.size();
         i++)
    {
        XX(aux) = atan2(YY(sampling_points_vector[i]),
                        XX(sampling_points_vector[i]));
        YY(aux) = acos(ZZ(sampling_points_vector[i]));
        if (YY(aux) < 0.)
            YY(aux) += PI;
        sampling_points_angles.push_back(aux*180. / PI);
    }
#ifdef NEVERDEFINED
    //sort points
    int k;
    int bound = sampling_points_angles.size() - 1;
    int t;
    int last_swap;
    Matrix1D<double> aux_angle(3);
    Matrix1D<double> aux_vector(3);

    while (bound)
    {
        last_swap = 0;
        for (k = 0; k < bound; k++)
        {
            aux = sampling_points_angles[k]; /* aux is a maximum of A[0]..A[k] */
            aux1 = sampling_points_vector[k]; /* aux is a maximum of A[0]..A[k] */
            if (sortFunc(aux, sampling_points_angles[k+1]))
            {
                sampling_points_angles[k] = sampling_points_angles[k+1];
                //                sampling_points_angles[k] = sampling_points_angles[k+1];
                sampling_points_angles[k+1] = aux; /*swap*/
                sampling_points_vector[k] = sampling_points_vector[k+1];
                sampling_points_vector[k] = sampling_points_vector[k+1];
                sampling_points_vector[k+1] = aux1; /*swap*/
                last_swap = k; /* mark the last swap position */
            }//if
        }//for
        bound = last_swap; /*elements after bound already sorted */
    }//while


#endif
    //#define DEBUG3
#ifdef  DEBUG3
    std::ofstream filestr;
    filestr.open ("sampling_file.bild");

    filestr    << ".color yellow"
    << std::endl
    ;
    for (int i = 0;
         i < sampling_points_vector.size();
         i++)
    {
        filestr  <<  ".sphere " << sampling_points_vector[i]  << " 0.018 " << std::endl;
    }
    filestr.close();
#endif
#undef DEBUG3

    numberSamplesAsymmetricUnit=sampling_points_vector.size();
}

// return 1 if a should go first 0 is equal -1 if before
int Sampling::sortFunc(Matrix1D<double> &t, Matrix1D<double> &a)
{
    if (YY(t) - 0.000001 > YY(a))
    {
        return 1;
    }
    else if (YY(t) + 0.000001 < YY(a))
    {
        return 0;
    };
    if (XX(t)  > XX(a))
    {
        return 1;
    }
    else
    {
        return 0;
    };
}

void Sampling::fillEdge(Matrix1D<double> starting_point,
                        Matrix1D<double> ending_point,
                        std::vector <Matrix1D<double> > & edge_vector,
                        bool END_FLAG
                       )
{
    Matrix1D<double> v_aux(3);

    double alpha;
    double beta;
    double gamma;
    // skip first corener, already computed;
    double upsilon = acos(dotProduct(starting_point, ending_point));
    for (size_t i1 = 1; i1 < number_of_samples; i1++)
    {
        gamma  = (double)i1 / (number_of_samples - 1);
        alpha  = sin((1. - gamma) * upsilon) / (sin(upsilon));
        beta   = sin(gamma * upsilon) / sin(upsilon);
        v_aux = alpha * starting_point + beta * ending_point;
        v_aux.selfNormalize();
        if (beta > 0.9999 && END_FLAG)
            continue;
        edge_vector.push_back(v_aux);
    }
}
void Sampling::fillDistance(Matrix1D<double> starting_point,
                            Matrix1D<double> ending_point,
                            std::vector <Matrix1D<double> > &
                            sampling_points_vector,
                            int my_number_of_samples,
                            bool only_half_sphere,
                            double min_z,
                            double max_z
                           )
{
    Matrix1D<double> v_aux(3);

    double alpha;
    double beta;
    double gamma;
    // skip first corener, already computed;
    double upsilon = acos(dotProduct(starting_point, ending_point));

    for (int i1 = 1; i1 < my_number_of_samples; i1++)
    {
        gamma  = (double)i1 / (my_number_of_samples);
        alpha  = sin((1. - gamma) * upsilon) / (sin(upsilon));
        beta   = sin(gamma * upsilon) / sin(upsilon);
        v_aux = alpha * starting_point + beta * ending_point;
        v_aux.selfNormalize();
        if (   (only_half_sphere && ZZ(v_aux) < 0.0)
               || ZZ(v_aux) < min_z
               || ZZ(v_aux) > max_z
           )
            //if (only_half_sphere && ZZ(v_aux) < 0.0)
            continue;
        else
            sampling_points_vector.push_back(v_aux);
    }
    /*
        //Remove points not in tilt range, check only the edges
        for (int i = 0;
             i < auxCounter;
             i++)
        {
            if (  ZZ(sampling_points_vector[i]) < min_z
               || ZZ(sampling_points_vector[i]) > max_z
               )
            sampling_points_vector[i].remove();
        }
    */
}

#define CLEAR_VECTORS() \
no_redundant_sampling_points_vector.clear();\
no_redundant_sampling_points_angles.clear();\
no_redundant_sampling_points_index.clear()

#define CREATE_INDEXES() \
    size_t __size = no_redundant_sampling_points_angles.size();\
    no_redundant_sampling_points_index.resize(__size, 0);\
    size_t * __ptrIndex = &(no_redundant_sampling_points_index[0]);\
    for (size_t i = 1; i < __size; ++i)\
      __ptrIndex[i] = i

void Sampling::removeRedundantPoints(const int symmetry, int sym_order)
{
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  aux(3, 3);
    Matrix1D<double>  row1(3), row2(3);

    //int j_end=0;
    Matrix1D<double>  row(3);

    CLEAR_VECTORS();

    if (symmetry == pg_CN)
    {//OK
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= (-180. / sym_order) &&
                XX(sampling_points_angles[i]) <= (180. / sym_order))
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry == pg_CI  ||
             symmetry == pg_CS )
    {//OK
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (YY(sampling_points_angles[i]) <= 90)
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNV )
    {//OK
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >=    0. / sym_order &&
                XX(sampling_points_angles[i]) <=  180. / sym_order)
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNH )
    {//OK
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= -180. / sym_order &&
                XX(sampling_points_angles[i]) <=  180. / sym_order &&
                YY(sampling_points_angles[i]) <=    90.
               )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_SN )
    {//OK
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= -180.*2. / sym_order &&
                XX(sampling_points_angles[i]) <=  180.*2. / sym_order &&
                YY(sampling_points_angles[i]) <=    90.
               )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DN )
    {
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= -180. / (sym_order)  + 90. &&
                XX(sampling_points_angles[i]) <=  180. / (sym_order)  + 90. &&
                YY(sampling_points_angles[i]) <=    90.
               )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNV )
    {
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >=    0. + 90. &&
                XX(sampling_points_angles[i]) <=  180. / (sym_order) +90. &&
                YY(sampling_points_angles[i]) <=    90.
               )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNH )
    {
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= 0. &&
                XX(sampling_points_angles[i]) <=  180. / (sym_order) &&
                YY(sampling_points_angles[i]) <=   90.
               )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_T )
    {//OK
        Matrix1D<double>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(-0.942809, 0., 0.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_3_fold_axis_3(3);
        _3_fold_axis_2_by_3_fold_axis_3 = vectorR3(0.471405, 0.272165, 0.7698);
        _3_fold_axis_2_by_3_fold_axis_3.selfNormalize();
        Matrix1D<double>  _3_fold_axis_3_by_3_fold_axis_1(3);
        _3_fold_axis_3_by_3_fold_axis_1 = vectorR3(0.471404, 0.816497, 0.);
        _3_fold_axis_3_by_3_fold_axis_1.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if ((XX(sampling_points_angles[i]) >=     90. &&
                XX(sampling_points_angles[i]) <=   150. )||
                XX(sampling_points_angles[i]) ==     0
               )
                if (
                    dotProduct(sampling_points_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(sampling_points_vector[i], _3_fold_axis_2_by_3_fold_axis_3) >= 0 &&
                    dotProduct(sampling_points_vector[i], _3_fold_axis_3_by_3_fold_axis_1) >= 0
                )
                {
                    no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                    no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_TD )
    {//OK
        Matrix1D<double>  _2_fold_axis_1_by_3_fold_axis_2(3);
        _2_fold_axis_1_by_3_fold_axis_2 = vectorR3(-0.942809, 0., 0.);
        _2_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_3_fold_axis_5(3);
        _3_fold_axis_2_by_3_fold_axis_5 = vectorR3(0.471405, 0.272165, 0.7698);
        _3_fold_axis_2_by_3_fold_axis_5.selfNormalize();
        Matrix1D<double>  _3_fold_axis_5_by_2_fold_axis_1(3);
        _3_fold_axis_5_by_2_fold_axis_1 = vectorR3(0., 0.471405, -0.666667);
        _3_fold_axis_5_by_2_fold_axis_1.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            //           if ( XX(sampling_points_angles[i])>=     120. &&
            //                 XX(sampling_points_angles[i])<=   150. ||
            //                 XX(sampling_points_angles[i])==     0
            //              )
            if (
                dotProduct(sampling_points_vector[i], _2_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_2_by_3_fold_axis_5) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_5_by_2_fold_axis_1) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_TH )
    {//OK
        Matrix1D<double>  _3_fold_axis_1_by_2_fold_axis_1(3);
        _3_fold_axis_1_by_2_fold_axis_1 = vectorR3(-0.816496, 0., 0.);
        _3_fold_axis_1_by_2_fold_axis_1.selfNormalize();
        Matrix1D<double>  _2_fold_axis_1_by_2_fold_axis_2(3);
        _2_fold_axis_1_by_2_fold_axis_2 = vectorR3(0.707107, 0.408248, -0.57735);
        _2_fold_axis_1_by_2_fold_axis_2.selfNormalize();
        Matrix1D<double>  _2_fold_axis_2_by_3_fold_axis_1(3);
        _2_fold_axis_2_by_3_fold_axis_1 = vectorR3(-0.408248, -0.707107, 0.);
        _2_fold_axis_2_by_3_fold_axis_1.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            //           if ( XX(sampling_points_angles[i])>=     120. &&
            //                 XX(sampling_points_angles[i])<=   150. ||
            //                 XX(sampling_points_angles[i])==     0
            //              )
            if (
                dotProduct(sampling_points_vector[i], _3_fold_axis_1_by_2_fold_axis_1) >= 0 &&
                dotProduct(sampling_points_vector[i], _2_fold_axis_1_by_2_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _2_fold_axis_2_by_3_fold_axis_1) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_O )
    {//OK
        Matrix1D<double>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(0., -1., 1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3(1., 1., 0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<double>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3(-1., 1., 0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if ((XX(sampling_points_angles[i]) >=   45. &&
                 XX(sampling_points_angles[i]) <=  135. &&
                 YY(sampling_points_angles[i]) <=  90.) ||
                XX(sampling_points_angles[i]) ==  0.
               )
                if (
                    dotProduct(sampling_points_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(sampling_points_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(sampling_points_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                )
                {
                    no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                    no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_OH )
    {//OK
        Matrix1D<double>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(0., -1., 1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3(1., 1., 0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<double>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3(-1., 1., 0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >=   90. &&
                XX(sampling_points_angles[i]) <=  135. &&
                YY(sampling_points_angles[i]) <=  90.)
                if (
                    dotProduct(sampling_points_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(sampling_points_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(sampling_points_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                )
                {
                    no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                    no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I || symmetry  == pg_I2)
    {//OK
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(-0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3(0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I1)
    {//OK
        Matrix2D<double>  A(3, 3);
        Euler_angles2matrix(0, 90, 0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();

        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3)
    {//OK
        Matrix2D<double>  A(3, 3);
        Euler_angles2matrix(0,31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();

        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4)
    {//OK
        Matrix2D<double>  A(3, 3);
        Euler_angles2matrix(0,-31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();

        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(sampling_points_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
                dotProduct(sampling_points_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) <= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5)
    {//OK
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else if (symmetry  == pg_IH || symmetry  == pg_I2H)
    {//OK
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(-0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3(0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  vectorR3(1.,0.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I1H)
    {//OK
        Matrix2D<double>  A(3, 3);
        Euler_angles2matrix(0, 90, 0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                        -0.8090170074556163,
                                        0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  A * vectorR3(1.,0.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(sampling_points_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3H)
    {//OK
        Matrix2D<double>  A(3, 3);
        Euler_angles2matrix(0,31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3(0.,1.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   >= 0 &&
                dotProduct(sampling_points_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   >= 0 &&
                dotProduct(sampling_points_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(sampling_points_vector[i],
                           _3_fold_axis_by_2_fold_axis)     >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4H)
    {//OK
        Matrix2D<double>  A(3, 3);
        Euler_angles2matrix(0,-31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3(0.,1.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (size_t i = 0; i < sampling_points_angles.size(); i++)
        {
            if (
                dotProduct(sampling_points_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(sampling_points_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
                dotProduct(sampling_points_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) <= 0 &&
                dotProduct(sampling_points_vector[i],
                           _3_fold_axis_by_2_fold_axis)     >= 0
            )
            {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5H)
    {//OK
        std::cerr << "ERROR: pg_I5H Symmetry not implemented" << std::endl;
        exit(0);
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }

    CREATE_INDEXES();
    numberSamplesAsymmetricUnit=no_redundant_sampling_points_vector.size();

}

void Sampling::removeRedundantPointsExhaustive(const int symmetry,
        int sym_order,
        bool only_half_sphere,
        double max_ang)
{
    // Maximum distance
    double cos_max_ang = cos(DEG2RAD(max_ang));
    double my_dotProduct;
    //int j_end=0;
    Matrix1D<double>  direction(3), direction1(3);

    // First call to conventional removeRedundantPoints
    removeRedundantPoints(symmetry, sym_order);
//    for (int isym = 0; isym < no_redundant_sampling_points_vector.size(); isym++)
//    	std::cout << "sampling:" << no_redundant_sampling_points_vector[isym];
    std::vector <Matrix1D<double> > old_vector = no_redundant_sampling_points_vector;
    std::vector <Matrix1D<double> > old_angles = no_redundant_sampling_points_angles;

    CLEAR_VECTORS();

    // Precalculate symmetry matrices
    fillLRRepository();

    // Then check all points versus each other
    for (size_t i = 0; i < old_angles.size(); i++)
    {
        //direction1=(old_vector[i]).transpose();
        direction1=old_vector[i];
        bool uniq = true;
        for (size_t j = 0; j < R_repository.size(); j++)
        {
            for (size_t k = 0; k < no_redundant_sampling_points_vector.size(); k++)
            {
                direction =  L_repository[j] *
                             (no_redundant_sampling_points_vector[k].transpose() *
                              R_repository[j]).transpose();
                //Calculate distance
                my_dotProduct = dotProduct(direction,direction1);
                if (only_half_sphere)
                    my_dotProduct = ABS(my_dotProduct);

                if (my_dotProduct > cos_max_ang)
                {
                    uniq = false;
                    break;
                }
            }// for k
            if (!uniq)
                break;
        } // for j
        if (uniq)
        {
            no_redundant_sampling_points_vector.push_back(old_vector[i]);
            no_redundant_sampling_points_angles.push_back(old_angles[i]);
        }
    } // for i

    CREATE_INDEXES();

}


//THIs FUNCTION IS NOW OBSOLETE
//SINCE readSymmetryFile does not longer need a file
//use symmetry functions instead
/* Create symmetry file----------------------------------------------------- */
void Sampling::createSymFile(FileName simFp,int symmetry, int sym_order)
{
    symmetry_file = simFp + ".sym";
    std::ofstream SymFile;
    SymFile.open(symmetry_file.c_str(), std::ios::out);
    if (symmetry == pg_CN)
    {
        SymFile << "rot_axis " << sym_order << " 0 0 1";
    }
    else if (symmetry == pg_CI)
    {
        SymFile << "inversion ";
    }
    else if (symmetry == pg_CS)
    {
        SymFile << "mirror_plane 0 0 1";
    }
    else if (symmetry == pg_CNV)
    {
        SymFile << "rot_axis " << sym_order << " 0 0 1";
        SymFile << std::endl;
        SymFile << "mirror_plane 0 1 0";
    }
    else if (symmetry == pg_CNH)
    {
        SymFile << "rot_axis " << sym_order << " 0 0 1";
        SymFile << std::endl;
        SymFile << "mirror_plane 0 0 -1";
    }
    else if (symmetry == pg_SN)
    {
        int order = sym_order / 2;
        SymFile << "rot_axis " << order << " 0 0 1";
        SymFile << std::endl;
        SymFile << "inversion ";
    }
    else if (symmetry == pg_DN)
    {
        SymFile << "rot_axis " << sym_order << " 0 0 1";
        SymFile << std::endl;
        SymFile << "rot_axis " << "2" << " 0 1 0";
    }
    else if (symmetry == pg_DNV)
    {
        SymFile << "rot_axis " << sym_order << " 0 0 1";
        SymFile << std::endl;
        SymFile << "rot_axis " << "2" << " 0 1 0";
        SymFile << std::endl;
        SymFile << "mirror_plane 0 1 0";
    }
    else if (symmetry == pg_DNH)
    {
        SymFile << "rot_axis " << sym_order << " 0 0 1";
        SymFile << std::endl;
        SymFile << "rot_axis " << "2" << " 1 0 0";
        SymFile << std::endl;
        SymFile << "mirror_plane 0 0 1";
    }
    else if (symmetry == pg_T)
    {
        SymFile << "rot_axis " << "3" << "  0. 0. 1.";
        SymFile << std::endl;
        SymFile << "rot_axis " << "2" << " 0. 0.816496 0.577350";
    }
    else if (symmetry == pg_TD)
    {
        SymFile << "rot_axis " << "3" << "  0. 0. 1.";
        SymFile << std::endl;
        SymFile << "rot_axis " << "2" << " 0. 0.816496 0.577350";
        SymFile << std::endl;
        SymFile << "mirror_plane 1.4142136 2.4494897 0.0000000";
    }
    else if (symmetry == pg_TH)
    {
        SymFile << "rot_axis " << "3" << "  0. 0. 1.";
        SymFile << std::endl;
        SymFile << "rot_axis " << "2" << " 0. -0.816496 -0.577350";
        SymFile << std::endl;
        SymFile << "inversion";
    }
    else if (symmetry == pg_O)
    {
        SymFile << "rot_axis " << "3" << "  .5773502  .5773502 .5773502";
        SymFile << std::endl;
        SymFile << "rot_axis " << "4" << " 0 0 1";
    }
    else if (symmetry == pg_OH)
    {
        SymFile << "rot_axis " << "3" << "  .5773502  .5773502 .5773502";
        SymFile << std::endl;
        SymFile << "rot_axis " << "4" << " 0 0 1";
        SymFile << std::endl;
        SymFile << "mirror_plane 0 1 1";
    }
    else if (symmetry == pg_I)
    {
        SymFile << "rot_axis 2  0             0          1";
        SymFile << std::endl;
        SymFile << "rot_axis 5 -1.618033989  -1           0";
        SymFile << std::endl;
        SymFile << "rot_axis 3 -0.53934467   -1.4120227   0";
    }
    else if (symmetry == pg_IH)
    {
        SymFile << "rot_axis 2  0             0          1";
        SymFile << std::endl;
        SymFile << "rot_axis 5 -1.618033989  -1           0";
        SymFile << std::endl;
        SymFile << "rot_axis 3 -0.53934467   -1.4120227   0";
        SymFile << std::endl;
        SymFile << "mirror_plane 1 0 0";
    }
    else
    {
        std::cerr << "ERROR:: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }
    SymFile.close();
    SL.readSymmetryFile(symmetry_file);

}
void Sampling::createAsymUnitFile(const FileName &docfilename)
{
    MetaData DF;
    FileName tmp_filename;
    //#define CHIMERA
#ifdef CHIMERA

    std::ofstream filestr;
    filestr.open ("create_asym_unit_file.bild");
    filestr    << ".color white"
    << std::endl
    << ".sphere 0 0 0 .95"
    << std::endl
    ;
    filestr    << ".color green"
    << std::endl
    ;
#endif

    MDRow row;
    for (size_t i = 0; i < no_redundant_sampling_points_vector.size(); i++)
    {
#ifdef CHIMERA
        filestr  << ".sphere "
        << no_redundant_sampling_points_vector[i]
        << " 0.018"
        << std::endl
        ;
#endif

        row.setValue(MDL_NEIGHBOR, no_redundant_sampling_points_index[i]);
        row.setValue(MDL_ANGLE_ROT,XX(no_redundant_sampling_points_angles[i]));
        row.setValue(MDL_ANGLE_TILT,YY(no_redundant_sampling_points_angles[i]));
        row.setValue(MDL_ANGLE_PSI,ZZ(no_redundant_sampling_points_angles[i]));
        row.setValue(MDL_X,XX(no_redundant_sampling_points_vector[i]));
        row.setValue(MDL_Y,YY(no_redundant_sampling_points_vector[i]));
        row.setValue(MDL_Z,ZZ(no_redundant_sampling_points_vector[i]));
        DF.addRow(row);
    }
#ifdef CHIMERA
    filestr.close();
#endif
    #undef CHIMERA

    tmp_filename = docfilename + "_angles.doc";
    DF.write(tmp_filename);
}

#define FN_SAMPLING_NEI(base) formatString("neighbors@%s_sampling.xmd", base.c_str())
#define FN_SAMPLING_PROJ(base) formatString("projectionDirections@%s_sampling.xmd", base.c_str())
#define FN_SAMPLING_EXTRA(base) formatString("extra@%s_sampling.xmd", base.c_str())
#define FN_SAMPLING_SPHERE(base) formatString("projectionDirectionsSphere@%s_sampling.xmd", base.c_str())

void Sampling::saveSamplingFile(const FileName &fn_base, bool write_vectors, bool write_sampling_sphere)
{
    MetaData md;
    MDRow row;

    row.setValue(MDL_SAMPLINGRATE, sampling_rate_rad);
    row.setValue(MDL_NEIGHBORHOOD_RADIUS, cos_neighborhood_radius);
    row.setValue(MDL_POINTSASYMETRICUNIT,numberSamplesAsymmetricUnit);
    md.setComment("data_extra -> sampling description;"\
                  " data_neighbors --> List with order of each"\
                  "experimental images and its neighbors"\
                 );
    md.setColumnFormat(false);
    md.addRow(row);
    md.write(FN_SAMPLING_EXTRA(fn_base), MD_OVERWRITE);

    md.clear();
    md.setColumnFormat(true);
    row.clear();
    size_t size = my_neighbors.size();
    //Write first block with experimental images order and its neighbors
    bool writeFileName = !exp_data_fileNames.empty();
    for(size_t i = 0; i < size; ++i)
    {
        row.setValue(MDL_ORDER,i+FIRST_IMAGE);
        if (writeFileName)
            row.setValue(MDL_IMAGE,exp_data_fileNames[i]);
        row.setValue(MDL_NEIGHBORS, my_neighbors[i]);
        md.addRow(row);
    }

    md.write(FN_SAMPLING_NEI(fn_base), MD_APPEND);

    //Write projection directions
    md.clear();
    row.clear();
    size = no_redundant_sampling_points_index.size();

    for (size_t i = 0; i < size; ++i)
    {
        Matrix1D<double> &angles = no_redundant_sampling_points_angles[i];
        row.setValue(MDL_NEIGHBOR, no_redundant_sampling_points_index[i]);
        row.setValue(MDL_ANGLE_ROT, XX(angles));
        row.setValue(MDL_ANGLE_TILT, YY(angles));
        row.setValue(MDL_ANGLE_PSI, ZZ(angles));

        if (write_vectors)
        {
            Matrix1D<double> &vectors = no_redundant_sampling_points_vector[i];
            row.setValue(MDL_X, XX(vectors));
            row.setValue(MDL_Y, YY(vectors));
            row.setValue(MDL_Z, ZZ(vectors));
        }
        md.addRow(row);
    }
    md.write(FN_SAMPLING_PROJ(fn_base), MD_APPEND);

    if (write_sampling_sphere)
    {
        md.clear();
        row.clear();
        size = sampling_points_angles.size();

        for (size_t i = 0; i < size; ++i)
        {
            row.setValue(MDL_NEIGHBOR, no_redundant_sampling_points_index[i]);
            Matrix1D<double> &angles = sampling_points_angles[i];
            row.setValue(MDL_ANGLE_ROT, XX(angles));
            row.setValue(MDL_ANGLE_TILT, YY(angles));
            row.setValue(MDL_ANGLE_PSI, ZZ(angles));

            if (write_vectors)
            {
                Matrix1D<double> &vectors = sampling_points_vector[i];
                row.setValue(MDL_X, XX(vectors));
                row.setValue(MDL_Y, YY(vectors));
                row.setValue(MDL_Z, ZZ(vectors));
            }
            md.addRow(row);
        }
        md.write(FN_SAMPLING_SPHERE(fn_base), MD_APPEND);

    }

}


void Sampling::readSamplingFile(const FileName &fn_base, bool read_vectors, bool write_sampling_sphere)
{
    //Read extra info
    MetaData md(FN_SAMPLING_EXTRA(fn_base));
    size_t id = md.firstObject();
    md.getValue(MDL_SAMPLINGRATE, sampling_rate_rad, id);
    md.getValue(MDL_NEIGHBORHOOD_RADIUS, cos_neighborhood_radius, id);
    md.getValue(MDL_POINTSASYMETRICUNIT,numberSamplesAsymmetricUnit,id);

    //Read neighbors
    md.read(FN_SAMPLING_NEI(fn_base));
    my_neighbors.resize(md.size());
    bool readFileName = md.containsLabel(MDL_IMAGE);
    if (readFileName)
    {
        exp_data_fileNames.clear();
        exp_data_fileNames.resize(md.size());
    }
    int ii=0;
    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        if (readFileName)
            md.getValue(MDL_IMAGE,exp_data_fileNames[ii++], __iter.objId);
        md.getValue(MDL_NEIGHBORS, my_neighbors[__iter.objIndex], __iter.objId);
    }

    //Read projection directions
    md.read(FN_SAMPLING_PROJ(fn_base));
    size_t size = md.size();
    no_redundant_sampling_points_index.resize(size);
    no_redundant_sampling_points_angles.resize(size);
    if (read_vectors)
        no_redundant_sampling_points_vector.resize(size);

    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        size_t &i = __iter.objIndex;
        size_t &id = __iter.objId;

        Matrix1D<double> &angles = no_redundant_sampling_points_angles[i];
        angles.resizeNoCopy(3);
        md.getValue(MDL_NEIGHBOR, no_redundant_sampling_points_index[i], id);
        md.getValue(MDL_ANGLE_ROT, XX(angles), id);
        md.getValue(MDL_ANGLE_TILT, YY(angles), id);
        md.getValue(MDL_ANGLE_PSI, ZZ(angles), id);

        if (read_vectors)
        {
            Matrix1D<double> &vectors = no_redundant_sampling_points_vector[i];
            vectors.resizeNoCopy(3);
            md.getValue(MDL_X, XX(vectors), id);
            md.getValue(MDL_Y, YY(vectors), id);
            md.getValue(MDL_Z, ZZ(vectors), id);
        }
    }

    if (write_sampling_sphere)
    {
        //Read projection directions
        md.read(FN_SAMPLING_SPHERE(fn_base));
        size_t size = md.size();
        sampling_points_angles.resize(size);
        if (read_vectors)
            sampling_points_vector.resize(size);

        FOR_ALL_OBJECTS_IN_METADATA(md)
        {
            size_t &i = __iter.objIndex;
            size_t &id = __iter.objId;

            Matrix1D<double> &angles = sampling_points_angles[i];
            angles.resizeNoCopy(3);
            md.getValue(MDL_ANGLE_ROT, XX(angles), id);
            md.getValue(MDL_ANGLE_TILT, YY(angles), id);
            md.getValue(MDL_ANGLE_PSI, ZZ(angles), id);

            if (read_vectors)
            {
                Matrix1D<double> &vectors = sampling_points_vector[i];
                vectors.resizeNoCopy(3);
                md.getValue(MDL_X, XX(vectors), id);
                md.getValue(MDL_Y, YY(vectors), id);
                md.getValue(MDL_Z, ZZ(vectors), id);
            }
        }

    }
}

void Sampling::computeNeighbors(bool only_winner)
{
    double my_dotProduct;
    double winner_dotProduct;
    Matrix1D<double>  row(3);
    Matrix2D<double>  L(4, 4), R(4, 4);
    std::vector<size_t>  aux_neighbors;
    std::vector<double> aux_neighbors_psi;
    std::vector <Matrix1D<double> > exp_data_projection_direction;
    Matrix1D<double>  direction(3);
    bool new_reference=true;
    my_neighbors.clear();
#ifdef MYPSI

    my_neighbors_psi.clear();
#endif

    // calculate some sizes only once
    size_t exp_data_projection_direction_by_L_R_size = exp_data_projection_direction_by_L_R.size();
    size_t no_redundant_sampling_points_vector_size = no_redundant_sampling_points_vector.size();

    if (verbose)
    {
    	std::cout << "Find valid sampling points based on the neighborhood" <<std::endl;
    	init_progress_bar(exp_data_projection_direction_by_L_R_size);
    }
    size_t ratio = exp_data_projection_direction_by_L_R_size / 60;
    ratio = XMIPP_MAX(ratio, 1);

    for(size_t j = 0; j < exp_data_projection_direction_by_L_R_size;)
    {
        if ((j%ratio) == 0 && verbose)
            progress_bar(j);

#ifdef MYPSI

        aux_neighbors_psi.clear();
#endif

        aux_neighbors.clear();
        size_t * aux_neighborsArray = NULL;
        for (size_t k = 0; k < R_repository.size(); k++,j++)
        {
            winner_dotProduct = -1.;
            for (size_t i = 0; i < no_redundant_sampling_points_vector_size; ++i)
            {
                my_dotProduct = dotProduct(no_redundant_sampling_points_vector[i],
                                           exp_data_projection_direction_by_L_R[j]);

                if (my_dotProduct > cos_neighborhood_radius)
                {
                    if(aux_neighbors.size()==0)
                    {
                        aux_neighbors.push_back(no_redundant_sampling_points_index[i]);
                        winner_dotProduct=my_dotProduct;

#ifdef MYPSI

                        aux_neighbors_psi.push_back(exp_data_projection_direction_by_L_R_psi[j]);
#endif

                    }
                    else
                    {
                        new_reference = true;
                        if(only_winner)
                        {
                            //std::cerr << "DEBUG_JM: In ONLY_WINNER" <<std::endl;
                            if(winner_dotProduct<my_dotProduct)
                            {
                                if(winner_dotProduct!=-1)
                                    aux_neighbors.pop_back();
#ifdef MYPSI

                                if(winner_dotProduct!=-1)
                                    aux_neighbors_psi.pop_back();
#endif

                                winner_dotProduct=my_dotProduct;
                            }
                            else
                            {
                                new_reference=false;
                            }
                        }
                        else
                        {
                            //precalculate size saves time here but
                            //not in the whole loop
                            aux_neighborsArray =  &aux_neighbors[0];
                            size_t _size = aux_neighbors.size();
                            //std::cerr << "DEBUG_JM: no_redundant_sampling_points_index[i]: " << no_redundant_sampling_points_index[i] << std::endl;
                            // for (size_t kkk = 0; kkk < aux_neighbors.size(); ++kkk)
                            //   std::cerr << aux_neighbors[kkk] << " ";
                            for( size_t l=0;l<  _size;l++)
                            {
                                //if (aux_neighbors[l]==i)
                                if (aux_neighborsArray[l]==no_redundant_sampling_points_index[i])
                                {
                                    new_reference=false;
                                    break;
                                }
                            }
                            //std::cerr << "DEBUG_JM: new_reference: " << new_reference << std::endl;
                        }
                        if (new_reference)
                        {
                            //std::cerr << formatString("DEBUG_JM: j %lu k %lu i %lu ", j, k, i) << std::endl;
                            aux_neighbors.push_back(no_redundant_sampling_points_index[i]);

                            //for (size_t kkk = 0; kkk < aux_neighbors.size(); ++kkk)
                            //  std::cerr << aux_neighbors[kkk] << " ";
                            // std::cerr << std::endl;

#ifdef MYPSI

                            aux_neighbors_psi.push_back(exp_data_projection_direction_by_L_R_psi[j]);
#endif

                        }
                    }
                    //same sampling point should appear only once
                    //note that psi recorded here may be different from psi
                    //recorded in _closest_sampling_points because
                    //may refer to a different sampling point
                    //in fact every point is degenerated
                }
            }//for i;
        }//for k
        my_neighbors.push_back(aux_neighbors);
#ifdef MYPSI

        my_neighbors_psi.push_back(aux_neighbors_psi);
#endif

    }//for j
    if (verbose)
    	progress_bar(exp_data_projection_direction_by_L_R_size);

    //#define DEBUG
#ifdef DEBUG

    for (int i=0;i< my_neighbors.size();i++)
        for (int j=0;j< my_neighbors[i].size();j++)
            std::cerr << "image:" << i << " "<< my_neighbors[i][j]<<std::endl;
    exit(1);
#endif
#undef DEBUG

    //#define CHIMERA
#ifdef CHIMERA

    std::ofstream filestr;
    filestr.open ("compute_neighbors.bild");
    filestr    << ".color white"
    << std::endl
    << ".sphere 0 0 0 .95"
    << std::endl
    ;
    int exp_image=1;
    filestr    <<  ".color yellow" << std::endl
    <<  ".sphere "   << exp_data_projection_direction_by_L_R[exp_image*R_repository.size()]
    <<  " .021"      << std::endl;
    for(int i=(exp_image*R_repository.size());
        i< (exp_image+1)*R_repository.size();
        i++)
    {
        filestr    <<  ".color red" << std::endl
        <<  ".sphere "   << exp_data_projection_direction_by_L_R[i]
        <<  " .017"      << std::endl;
    }
    double blue;
    for(int i=0;
        i< my_neighbors[exp_image].size();
        i++)
    {
#ifdef MYPSI
        blue = (my_neighbors_psi[exp_image][i]+180.)/360.;
#else

        blue = 1.;
#endif

        filestr    <<  ".color 0 0 " << blue  << std::endl
        <<  ".sphere "   <<
        no_redundant_sampling_points_vector[my_neighbors[exp_image][i]]
        <<  " .019"      << std::endl;
        //std::cerr << my_neighbors_psi[exp_image][i] << std::endl;
    }
    filestr.close();

#endif
    #undef CHIMERA

}
/** Remove all those points no closer than neighborhood_radius_rad
    */
#define REMOVE_LAST(vector) vector[i] = vector[my_end]; vector.pop_back();

void Sampling::removePointsFarAwayFromExperimentalData()
{
    double my_dotProduct;
    Matrix1D<double>  row(3),direction(3);
    Matrix2D<double>  L(4, 4), R(4, 4);

    size_t my_end = no_redundant_sampling_points_vector.size() - 1;

    for (size_t i = 0; i <= my_end; i++)
    {
        bool my_delete=true;
        for (size_t j=0; my_delete && j< exp_data_projection_direction_by_L_R.size();j++)
        {
            my_dotProduct = dotProduct(no_redundant_sampling_points_vector[i],
                                       exp_data_projection_direction_by_L_R[j]);
            if (my_dotProduct > cos_neighborhood_radius)
                my_delete=false;
        }//for j
        if(my_delete)
        {
            REMOVE_LAST(no_redundant_sampling_points_vector);
            REMOVE_LAST(no_redundant_sampling_points_angles);
            REMOVE_LAST(no_redundant_sampling_points_index);

            --my_end;
            --i;//since a point has been swaped we should repeat the same index
        }// if(my_delete)
    }//for i end
    //#define CHIMERA
#ifdef CHIMERA
    std::ofstream filestr;
    filestr.open ("remove_points_far_away_from_experimental_data.bild");
    filestr    << ".color white"
    << std::endl
    << ".sphere 0 0 0 .95"
    << std::endl
    ;
    filestr    << ".color green"
    << std::endl
    ;
    //green neighbours
    for (int i = 0;
         i < no_redundant_sampling_points_vector.size();
         i++)
    {
        filestr    <<  ".color green" << std::endl
        <<  ".sphere " << no_redundant_sampling_points_vector[i].transpose()  <<
        " .018" << std::endl;
    }
    filestr.close();

#endif
    #undef CHIMERA
}
void Sampling::findClosestSamplingPoint(const FileName &FnexperimentalImages,
                                        const FileName &output_file_root)
{
    //read input files
    MetaData DFi;
    DFi.read(FnexperimentalImages);//experimental points
    findClosestSamplingPoint(DFi,output_file_root);

}
void Sampling::findClosestSamplingPoint(MetaData &DFi,
                                        const FileName &output_file_root)
{
    double my_dotProduct,my_dotProduct_aux;
    Matrix1D<double>  row(3),direction(3);
    Matrix1D<double> docline;
    docline.initZeros(7);//three original angles, one winnir, new angles
    Matrix2D<double>  L(4, 4), R(4, 4);
    int winner_sampling=-1;
#if defined(CHIMERA) || defined(MYPSI)
    int winner_exp_L_R=-1;
#endif

    MetaData DFo;
    size_t id;

    DFo.setComment("Original rot, tilt, psi, Xoff, Yoff are stored as comments");

    //#define DEBUG3
#ifdef  DEBUG3

    std::ofstream filestr;
    filestr.open ("find_closest_sampling_point.bild");
    int exp_image=1;
#endif

    MDIterator iter(DFi);
    for(size_t i=0;i< exp_data_projection_direction_by_L_R.size();)
    {
        my_dotProduct=-2;
        for (size_t k = 0; k < R_repository.size(); k++,i++)
        {
#ifdef  DEBUG3
            //experimental points plus symmetry
            if( i>(exp_image*R_repository.size()-1) && i< ((exp_image+1)*R_repository.size()))
            {
                filestr    <<  ".color red" << std::endl
                <<  ".sphere "   << exp_data_projection_direction_by_L_R[i]
                <<  " .019"      << std::endl;
            }
#endif
            for(size_t j=0;j< no_redundant_sampling_points_vector.size();j++)
            {
                my_dotProduct_aux =
                    dotProduct(exp_data_projection_direction_by_L_R[i],
                               no_redundant_sampling_points_vector[j]);

                if ( my_dotProduct_aux > my_dotProduct)
                {
                    my_dotProduct = my_dotProduct_aux;
                    winner_sampling = j;
#if defined(CHIMERA) || defined(MYPSI)
                    winner_exp_L_R  = i;
#endif
                }
            }//for j
        }//for k
#ifdef  DEBUG3
        if( i==  ((exp_image+1)*R_repository.size()) )
        {
            filestr    <<  ".color yellow" << std::endl
            <<  ".sphere "   << no_redundant_sampling_points_vector[winner_sampling]
            <<  " .020"      << std::endl;
        }
#endif
        //add winner to the DOC fILE
        std::string fnImg, comment;
        double aux;
        DFi.getValue(MDL_IMAGE, fnImg, iter.objId);
        DFi.getValue(MDL_ANGLE_ROT,aux, iter.objId);
        comment+=floatToString(aux)+" ";
        DFi.getValue(MDL_ANGLE_TILT,aux, iter.objId);
        comment+=floatToString(aux)+" ";
        DFi.getValue(MDL_ANGLE_PSI,aux, iter.objId);
        comment+=floatToString(aux)+" ";
        DFi.getValue(MDL_SHIFT_X,aux, iter.objId);
        comment+=floatToString(aux)+" ";
        DFi.getValue(MDL_SHIFT_Y,aux, iter.objId);
        comment+=floatToString(aux);
        id = DFo.addObject();
        DFo.setValue(MDL_COMMENT,comment, id);
        DFo.setValue(MDL_IMAGE,fnImg, id);
        DFo.setValue(MDL_REF, winner_sampling, id);
#ifdef MYPSI

        DFo.set(6, exp_data_projection_direction_by_L_R_psi[winner_exp_L_R]);
#endif

        DFo.setValue(MDL_NEIGHBOR, no_redundant_sampling_points_index[winner_sampling], id);
        DFo.setValue(MDL_ANGLE_ROT,XX(no_redundant_sampling_points_angles[winner_sampling]), id);
        DFo.setValue(MDL_ANGLE_TILT,YY(no_redundant_sampling_points_angles[winner_sampling]), id);
        DFo.setValue(MDL_ANGLE_PSI,ZZ(no_redundant_sampling_points_angles[winner_sampling]), id);

        iter.moveNext();
    }//for i
    if (output_file_root.size() > 0)
        DFo.write(output_file_root+ "_closest_sampling_points.doc");
#ifdef  DEBUG3

    filestr.close();
#endif
#undef DEBUG3
}

void Sampling::findClosestExperimentalPoint()
{
    double my_dotProduct,my_dotProduct_aux;
    Matrix1D<double>  row(3),direction(3);
    int winner_sampling=-1;
    int winner_exp=-1;
    //#define CHIMERA
#ifdef CHIMERA

    std::vector<std::vector<int> >  aux_vec;
    aux_vec.resize(no_redundant_sampling_points_vector.size());
#endif

    std::vector<std::vector<size_t> >  aux_my_exp_img_per_sampling_point;

    //resize vector
    aux_my_exp_img_per_sampling_point.resize(
        no_redundant_sampling_points_vector.size());

    for(size_t i=0,l=0;i< exp_data_projection_direction_by_L_R.size();l++)
    {
        my_dotProduct=-2;
        for (size_t k = 0; k < R_repository.size(); k++,i++)
        {
            for(size_t j=0;j< no_redundant_sampling_points_vector.size();j++)
            {
                my_dotProduct_aux =
                    dotProduct(exp_data_projection_direction_by_L_R[i],
                               no_redundant_sampling_points_vector[j]);

                if ( my_dotProduct_aux > my_dotProduct)
                {
                    my_dotProduct = my_dotProduct_aux;
                    winner_sampling = j;
#ifdef CHIMERA

                    winner_exp_L_R  = i;
#endif

                    winner_exp = l;
                }
            }//for j
        }//for k
        aux_my_exp_img_per_sampling_point[winner_sampling].push_back(winner_exp);
#ifdef CHIMERA

        aux_vec[winner_sampling].push_back(winner_exp_L_R);
#endif

    }//for i aux_my_exp_img_per_sampling_point
    for(size_t i=0;i< aux_my_exp_img_per_sampling_point.size();i++)
        if(aux_my_exp_img_per_sampling_point[i].size()!=0)
            my_exp_img_per_sampling_point.push_back(aux_my_exp_img_per_sampling_point[i]);
#ifdef CHIMERA

    std::ofstream filestr;
    filestr.open ("find_closest_experimental_point.bild");
    filestr    << ".color white"
    << std::endl
    << ".sphere 0 0 0 .95"
    << std::endl
    ;
    filestr    << ".color red"
    << std::endl
    ;
    for (int i = 0;
         i < no_redundant_sampling_points_vector.size();
         i++)
    {
        filestr    <<  ".sphere " << no_redundant_sampling_points_vector[i].transpose()  <<
        " .018" << std::endl;
    }
    int my_sample_point=5;
    filestr    << ".color green"
    << std::endl
    ;
    int ii;
    for (int i = 0;
         i < my_exp_img_per_sampling_point[my_sample_point].size();
         i++)
    {
        ii=aux_vec[my_sample_point][i];
        filestr    << ".sphere "
        << exp_data_projection_direction_by_L_R[ii].transpose()
        << " .017" << std::endl;
    }
    filestr.close();

#endif
    #undef CHIMERA
    //#define DEBUG4
#ifdef DEBUG4

    std::ofstream filestr;
    filestr.open ("find_closest_experimental_point.txt");

    for (int i = 0;
         i < my_exp_img_per_sampling_point.size();
         i++)
    {   //for each sampling point write its experimental images
        filestr << i << std::endl;
        for (int j = 0;
             j < my_exp_img_per_sampling_point[i].size();
             j++)
        {
            filestr    << my_exp_img_per_sampling_point[i][j]
            << " " ;
        }
        filestr << std::endl;
    }
    filestr.close();

#endif
    #undef DEBUG4
}

void Sampling::fillLRRepository(void)
{
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    //NEXT 2 ROB
    R_repository.clear();
    L_repository.clear();
    R_repository.push_back(Identity);
    L_repository.push_back(Identity);
    for (int isym = 0; isym < SL.symsNo(); isym++)
    {
        SL.getMatrices(isym, L, R);
        R.resize(3, 3);
        L.resize(3, 3);
        R_repository.push_back(R);
        L_repository.push_back(L);
    }
    //#define DEBUG3
#ifdef  DEBUG3
    for (int isym = 0; isym < R_repository.size(); isym++)
    {
        std::cout << R_repository[isym];
        std::cout << L_repository[isym];
    }
#endif
#undef DEBUG3
}

void Sampling::fillExpDataProjectionDirectionByLR(
    const FileName &FnexperimentalImages)
{
    //read input files
    MetaData DFi;
    DFi.read(FnexperimentalImages);//experimental points
    fillExpDataProjectionDirectionByLR(DFi);
}

void Sampling::fillExpDataProjectionDirectionByLR(MetaData &DFi)
{
    std::vector <Matrix1D<double> > exp_data_projection_direction;
    Matrix1D<double>  direction(3);
    DFi.firstObject();
    //#define CHIMERA
#ifdef CHIMERA

    std::ofstream filestr;
    filestr.open ("exp_data_projection_direction_by_L_R.bild");
    filestr    << ".color white"
    << std::endl
    << ".sphere 0 0 0 .95"
    << std::endl
    ;
    filestr    << ".color green"
    << std::endl
    ;
#endif

    double img_tilt,img_rot,img_psi;
    FileName imgName;
    exp_data_fileNames.clear();
    FOR_ALL_OBJECTS_IN_METADATA(DFi)
    {
        DFi.getValue(MDL_ANGLE_ROT,img_rot,__iter.objId);
        DFi.getValue(MDL_ANGLE_TILT,img_tilt,__iter.objId);
        DFi.getValue(MDL_ANGLE_PSI,img_psi,__iter.objId);
        Euler_direction(img_rot, img_tilt, img_psi, direction);
        exp_data_projection_direction.push_back(direction);
        DFi.getValue(MDL_IMAGE,imgName,__iter.objId);
        exp_data_fileNames.push_back(imgName);
    }

    exp_data_projection_direction_by_L_R.clear();
#ifdef MYPSI

    exp_data_projection_direction_by_L_R_psi.clear();
#endif

    for (size_t i = 0; i < exp_data_projection_direction.size(); i++)
        for (size_t j = 0; j < R_repository.size(); j++)
        {
            direction =  L_repository[j] *
                         (exp_data_projection_direction[i].transpose() *
                          R_repository[j]).transpose();
            exp_data_projection_direction_by_L_R.push_back(direction);
#ifdef MYPSI

            Euler_apply_transf(L_repository[j],
                               R_repository[j], img_rot,
                               img_tilt,
                               img_psi,
                               rotp,
                               tiltp,
                               psip);
            exp_data_projection_direction_by_L_R_psi.push_back(psip);
#endif
            //#define CHIMERA
#ifdef CHIMERA

            filestr << ".sphere " << direction.transpose()
            << " 0.02" << std::endl
            ;
#endif

        }
    //#define CHIMERA
#ifdef CHIMERA
    filestr.close();
#endif
    #undef CHIMERA
}
