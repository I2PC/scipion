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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#include "sampling.h"

/* Default Constructor */
XmippSampling::XmippSampling()
{

    vertices aux;
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
    
    //#define DEBUG1
#ifdef  DEBUG1
    for (int i = 0;
         i < vertices_vectors.size();
         i++)
        std::cout  <<  vertices_vectors[].transpose()  << std::endl;
#endif
#undef DEBUG1
}

void XmippSampling::SetSampling(double sampling)
{
    sampling_rate_rad = DEG2RAD(sampling);
    number_of_samples = ROUND(cte_w / sampling_rate_rad);
}

void XmippSampling::SetNoise(double noise_deviation, int my_seed)
{
    sampling_noise = noise_deviation;
    init_random_generator(my_seed);
}

void XmippSampling::SetNeighborhoodRadius(double neighborhood)
{
    neighborhood_radius_rad = DEG2RAD(neighborhood);
    cos_neighborhood_radius = cos(neighborhood_radius_rad);
}

/* Compute edge sampling points using Baumgardner  1995 */
void XmippSampling::Compute_sampling_points(bool only_half_sphere,
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

    if(max_tilt >= 90.)
       max_z=10.;
    else
       max_z=sin(PI * max_tilt / 180.);   
    if(min_tilt <= -90.)
       min_z= -10.;
    else
       min_z=sin(PI * min_tilt / 180.);  

    //01a
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[1];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[6]; ending_point = vertices_vectors[1];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //01b
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[2];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[6]; ending_point = vertices_vectors[2];
    fill_edge(starting_point, ending_point, edge_vector_end, true);
    //02a
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[2];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[7]; ending_point = vertices_vectors[2];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //02b
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[3];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[7]; ending_point = vertices_vectors[3];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //03a
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[3];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[8]; ending_point = vertices_vectors[3];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //03b
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[4];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[8]; ending_point = vertices_vectors[4];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //04a
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[4];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[9]; ending_point = vertices_vectors[4];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //04b
    starting_point = vertices_vectors[0]; ending_point = vertices_vectors[5];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[9]; ending_point = vertices_vectors[5];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //05a
    starting_point = vertices_vectors[0];  ending_point = vertices_vectors[5];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[10]; ending_point = vertices_vectors[5];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //05b
    starting_point = vertices_vectors[0];  ending_point = vertices_vectors[1];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[10]; ending_point = vertices_vectors[1];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //06a
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[10];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[5];  ending_point = vertices_vectors[10];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //06b
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[9];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[5];  ending_point = vertices_vectors[9];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //07a
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[9];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[4];  ending_point = vertices_vectors[9];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //07b
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[8];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[4];  ending_point = vertices_vectors[8];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //08a
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[8];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[3];  ending_point = vertices_vectors[8];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //08b
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[7];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[3];  ending_point = vertices_vectors[7];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //09a
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[7];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[2];  ending_point = vertices_vectors[7];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //09b
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[6];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[2];  ending_point = vertices_vectors[6];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //10a
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[6];
    fill_edge(starting_point, ending_point, edge_vector_start, false);
    starting_point = vertices_vectors[1]; ending_point = vertices_vectors[6];
    fill_edge(starting_point, ending_point, edge_vector_start, true);
    //10b
    starting_point = vertices_vectors[11]; ending_point = vertices_vectors[10];
    fill_edge(starting_point, ending_point, edge_vector_end, false);
    starting_point = vertices_vectors[1];  ending_point = vertices_vectors[10];
    fill_edge(starting_point, ending_point, edge_vector_end, true);

    //#define DEBUG2
#ifdef  DEBUG2
    for (int i = 0;
         i < edge_vector_start.size();
         i++)
    {
        std::cout  <<  edge_vector_start[i].transpose()  << " 1 1 " << std::endl;
        std::cout  <<  edge_vector_end[i].transpose()  << " 1 2 " << std::endl;
    }
    //std::cout  <<  ending_point.transpose()    << " 1.1 1.5 " << std::endl;
#endif
#undef DEBUG2
/*

    // add  main corners
    for (int i = 0;
         i < vertices_vectors.size();
         i++)
    {
        if (   (only_half_sphere && ZZ(vertices_vectors[i]) < 0.0)
           || ZZ(vertices_vectors[i]) < min_z
           || ZZ(vertices_vectors[i]) > max_z
           )
        //if (only_half_sphere && ZZ(vertices_vectors[i]) < 0.0)
            continue;
        else
            sampling_points_vector.push_back(vertices_vectors[i]);
    }
*/
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
    for (int i = 0;
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
     for (int i = 0;
	  i < sampling_points_vector.size();
	  i++)
     {
         std::cout  <<  ".color 1 0 0" << std::endl
                    <<  ".sphere " << sampling_points_vector[i].transpose()  << 
	               " .025" << std::endl;
	 sampling_points_vector[i].add_noise(0.0, sampling_noise, "gaussian");
	 sampling_points_vector[i].selfNormalize();
         std::cout  <<  ".color 0 0 1" << std::endl
                    <<  ".sphere " << sampling_points_vector[i].transpose()  << 
	               " .027" << std::endl;
     }
    
#endif
#undef DEBUG3

    // add in between points
    int j = 0;
    bool j_flag = false;
    for (int i = 0;
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
        fill_distance(edge_vector_start[i],
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
    if(sampling_noise!=0.0){ 
	 for (int i = 0;
              i < sampling_points_vector.size();
              i++)
	 {
             sampling_points_vector[i].add_noise(0.0, sampling_noise, "gaussian");
	     sampling_points_vector[i].selfNormalize();
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
    double rot, tilt, psi;
    for (int i = 0;
         i < sampling_points_vector.size();
         i++)
    {
        XX(aux) = atan2(YY(sampling_points_vector[i]),
                        XX(sampling_points_vector[i]));
        YY(aux) = acos(ZZ(sampling_points_vector[i]));
        if (YY(aux) < 0.) YY(aux) += PI;
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
            if (sort_func(aux, sampling_points_angles[k+1]))
            {
                sampling_points_angles[k] = sampling_points_angles[k+1];
                sampling_points_angles[k] = sampling_points_angles[k+1];
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
    for (int i = 0;
         i < sampling_points_angles.size();
         i++)
    {
        std::cout  <<  sampling_points_angles[i].transpose()  << " 1 1 " << std::endl;
    }
#endif
#undef DEBUG3

}

// return 1 if a should go first 0 is equal -1 if before
int XmippSampling::sort_func(Matrix1D<double> &t, Matrix1D<double> &a)
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

void XmippSampling::fill_edge(Matrix1D<double> starting_point,
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
    for (int i1 = 1; i1 < number_of_samples; i1++)
    {
        gamma  = (double)i1 / (number_of_samples - 1);
        alpha  = sin((1. - gamma) * upsilon) / (sin(upsilon));
        beta   = sin(gamma * upsilon) / sin(upsilon);
        v_aux = alpha * starting_point + beta * ending_point;
        v_aux.selfNormalize();
        if (beta > 0.9999 && END_FLAG) continue;
        edge_vector.push_back(v_aux);
    }
}
void XmippSampling::fill_distance(Matrix1D<double> starting_point,
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
/*
    //we need to remove projections with tilt angle not
    //in between min_tilt-max_tilt
    int auxCounter;
    auxCounter=sampling_points_vector.size();
*/
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

void XmippSampling::remove_redundant_points(const int symmetry,
        int sym_order)
{
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  aux(3, 3);
    Matrix1D<double>  row1(3), row2(3);
    bool valid = true;
    double rot, tilt, psi = 0;
    double rotp, tiltp, psip = 0.;
    no_redundant_sampling_points_vector.clear();
    double aux1, aux2;
    bool match = false;

    //int j_end=0;
    Matrix1D<double>  row(3);

    no_redundant_sampling_points_vector.clear();
    no_redundant_sampling_points_angles.clear();
    double my_dotProduct;
    if (symmetry == pg_CN)
    {//OK
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= -180. / (sym_order) &&
                XX(sampling_points_angles[i]) <=  180. / (sym_order) &&
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >=    0. &&
                XX(sampling_points_angles[i]) <=  180. / (sym_order) &&
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >=     90. &&
                XX(sampling_points_angles[i]) <=   150. ||
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
        for (int i = 0; i < sampling_points_angles.size(); i++)
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
    else if (symmetry  == pg_I )
    {//OK
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        for (int i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >= -31.7174745559 &&
                XX(sampling_points_angles[i]) <=  31.7174745559 &&
                YY(sampling_points_angles[i]) <=  90 &&
                YY(sampling_points_angles[i]) >=  69.094843368)
                if (
                    dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                    dotProduct(sampling_points_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
                {
                    no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                    no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_IH )
    {//OK
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        for (int i = 0; i < sampling_points_angles.size(); i++)
        {
            if (XX(sampling_points_angles[i]) >=  0. &&
                XX(sampling_points_angles[i]) <=  31.7174745559 &&
                YY(sampling_points_angles[i]) <=  90 &&
                YY(sampling_points_angles[i]) >=  69.094843368)
                if (
                    dotProduct(sampling_points_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0
                )
                {
                    no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                    no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
        }// for i
    }
    else
    {
        std::cerr << "eRROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }

}
//THIs FUNCTION IS NOW OBSOLETE SIN 
//SINCE read_sym_file does not longer need a file
//use symmetry functions instead
/* Create symmetry file----------------------------------------------------- */
void XmippSampling::create_sym_file(FileName simFp,int symmetry, int sym_order)
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
        SymFile << "mirror_plane 0 0 1";
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
    SL.read_sym_file(symmetry_file);

}
void XmippSampling::create_asym_unit_file(const FileName &docfilename)
{
    DocFile DFvectors, DFangles;
    DFvectors.clear();
    DFangles.clear();
    FileName tmp_filename;
    for (int i = 0; i < no_redundant_sampling_points_vector.size(); i++)
    {
        DFvectors.insert_data_line(no_redundant_sampling_points_vector[i]);
        DFangles.insert_data_line(no_redundant_sampling_points_angles[i]);
    }
    tmp_filename = docfilename + "_vectors.doc";
    DFvectors.write(tmp_filename);
    tmp_filename = docfilename + "_angles.doc";
    DFangles.write(tmp_filename);
}

void XmippSampling::compute_neighbors(void)
{
    double rot,  tilt,  psi;
    double rotp, tiltp, psip;
    double my_dotProduct;
    Matrix1D<double>  row(3);
    Matrix2D<double>  L(4, 4), R(4, 4);
    std::vector<int>  aux_neighbors;
    std::vector<double> aux_neighbors_psi;
    for (int i = 0; i < no_redundant_sampling_points_vector.size(); i++)
    {
        aux_neighbors_psi.clear();
        aux_neighbors.clear();
        //#define DEBUGH
        //#define DEBUGHH
        for (int j = 0; j < no_redundant_sampling_points_vector.size(); j++)
        {
#ifdef DEBUGH
            std::cerr << "i,j " << i << " " << j << std::endl;
#endif
            if (i == j) continue;
            //check identity
            my_dotProduct = dotProduct(no_redundant_sampling_points_vector[i],
                                         no_redundant_sampling_points_vector[j]);

            if (my_dotProduct > cos_neighborhood_radius)
            {
                aux_neighbors.push_back(j);
                aux_neighbors_psi.push_back(0.0);
            }
            else
            {
                for (int isym = 0; isym < SL.SymsNo(); isym++)
                {
#ifdef DEBUGH
                    std::cerr << "            isym " << isym << std::endl;
#endif
                    SL.get_matrices(isym, L, R);
                    R.resize(3, 3);
                    row =  no_redundant_sampling_points_vector[j].transpose() * R;
                    my_dotProduct = dotProduct(no_redundant_sampling_points_vector[i],
                                                 row);

#ifdef DEBUGHH
                    std::cerr << "            row out " << row << std::endl;
#endif
                    if (my_dotProduct > cos_neighborhood_radius)
                    {
                        aux_neighbors.push_back(j);
                        rot = XX(no_redundant_sampling_points_angles[j]);
                        tilt = YY(no_redundant_sampling_points_angles[j]);
                        psi = ZZ(no_redundant_sampling_points_angles[j]);
                        L.resize(3, 3); // Erase last row and column
                        Euler_apply_transf(L, R, rot, tilt, psi, rotp, tiltp, psip);
                        aux_neighbors_psi.push_back(psip);
#ifdef DEBUGHH
                        Euler_direction(rotp, tiltp, psip, row);
                        std::cerr << "row in " << row << std::endl;
#endif
                        break;
                    }
                }//for (int isym=0
            }
        }//for (int j=0;
        my_neighbors.push_back(aux_neighbors);
        my_neighbors_psi.push_back(aux_neighbors_psi);
    }//for i
#undef DEBUGH
}
