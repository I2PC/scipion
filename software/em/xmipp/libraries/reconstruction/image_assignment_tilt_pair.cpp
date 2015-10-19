/***************************************************************************
 * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
 *
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

#include "image_assignment_tilt_pair.h"
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time
#include <algorithm>
#include <iostream>
#include <data/xmipp_image.h>
#include <data/micrograph.h>
#include <delaunay/delaunay.h>
#include <delaunay/dcel.h>
#include <geometry.h>
//#include <relion-1.3/include/relion-1.3/src/matrix1d.h>
#include <data/matrix1d.h>

void ProgassignmentTiltPair::readParams()
{
	fnuntilt = getParam("--untiltcoor");
	fntilt = getParam("--tiltcoor");
	fnmic = getParam("--tiltmicsize");
	fndir = getParam("--odir");
	mshift = getDoubleParam("--maxshift");
	particle_size = getDoubleParam("--particlesize");
	thr = getDoubleParam("--threshold");
}

void ProgassignmentTiltPair::defineParams()
{
	//usage
	addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
	//params

	addParamsLine("  [--untiltcoor <md_file=\"\">]    : Untilt coordinates");
	addParamsLine("  [--tiltcoor <md_file=\"\">]    : Tilt coordinates");
	addParamsLine("  [--tiltmicsize <img_file=\"\">]    : Tilt micrography");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--maxshift <s=1000>]   : Maximum shift");
	addParamsLine("  [--particlesize <p=100>]   : Particle size");
	addParamsLine("  [--threshold <d=0.3>]      : If the distance between two points is lesser than threshold*particlesize, they will be the same point");
}

double mean(Matrix1D<double> v)
{
	int Num_elem = v.vdim;
	double sum_v = v.sum();
	return sum_v/Num_elem;
}

void ProgassignmentTiltPair::run()
{
	std::cout << "Starting..." << std::endl;
	//LOAD METADATA and TILTPAIRS
	MetaData md_untilt, md_tilt, mduntilt, mdtilt;;
	size_t Ndim, Zdim, Ydim , Xdim, objId, breakassignment = 0;


	getImageSize(fnmic, Xdim, Ydim, Zdim, Ndim);
	/*M_in.read(fnmic);
	getImageSize(M_in,Xdim,Ydim,Zdim,Ndim);*/

	md_untilt.read(fnuntilt);
	md_tilt.read(fntilt);


	int x, y, len_u=0;
	//storing points and creating Delaunay triangulations
	struct Delaunay_T delaunay_untilt;
	Matrix1D<double> ux(md_untilt.size()), uy(md_untilt.size());
	init_Delaunay( &delaunay_untilt, md_untilt.size());
	FOR_ALL_OBJECTS_IN_METADATA(md_untilt)
	{

		md_untilt.getValue(MDL_XCOOR, x, __iter.objId);
		md_untilt.getValue(MDL_YCOOR, y, __iter.objId);

		VEC_ELEM(ux,len_u) = x;
		VEC_ELEM(uy,len_u) = y;

		insert_Point( &delaunay_untilt, x, y);
		len_u = len_u +1;
	}
	std::cout << "len_u = " << len_u << std::endl;


	create_Delaunay_Triangulation( &delaunay_untilt, 0);
	//printf("CREATED1d\n");
	int tri_number_untilt = get_Number_Real_Faces(delaunay_untilt.dcel);

	struct Delaunay_T delaunay_tilt;
	init_Delaunay( &delaunay_tilt, md_tilt.size());
	FOR_ALL_OBJECTS_IN_METADATA(md_tilt)
	{
		md_tilt.getValue(MDL_XCOOR, x,__iter.objId);
		md_tilt.getValue(MDL_YCOOR, y,__iter.objId);
		insert_Point( &delaunay_tilt, x, y);
	}

	write_DCEL(delaunay_tilt.dcel, 0, "tiltdata.txt");
	create_Delaunay_Triangulation( &delaunay_tilt, 1);
	int tri_number_tilt = get_Number_Real_Faces(delaunay_tilt.dcel);

	// Triangle areas
	struct Point_T p, q, r, u1, u2, u3, t1, t2, t3;
	MultidimArray<double> trig_untilt_area(tri_number_untilt+1), trig_tilt_area(tri_number_tilt+1), sortedArea;
	Matrix1D<double> trig_untilt_area_1D(tri_number_untilt+1), trig_tilt_area_1D(tri_number_tilt+1);

	for (int i=1; i<tri_number_untilt+1 ;i++)
	{
		get_Face_Points(&delaunay_untilt, i, &p, &q, &r); //i is the triangle number
		A1D_ELEM(trig_untilt_area,i-1) = triangle_area(p.x, p.y, q.x, q.y, r.x, r.y);
	}

	for (int i=1; i<tri_number_tilt+1 ;i++)
	{
		get_Face_Points(&delaunay_tilt, i, &p, &q, &r);
		A1D_ELEM(trig_tilt_area,i-1) = triangle_area(p.x, p.y, q.x, q.y, r.x, r.y);
	}
	double t_mean_area=trig_tilt_area.computeAvg();
	double u_mean_area=trig_untilt_area.computeAvg();

	trig_untilt_area.sort(sortedArea);

	//std:cout << "AVERAGE AREA = " << trig_untilt_area << std::endl;


	double threshold_area = sortedArea( round((tri_number_untilt+1)*0.2) );

	//std:cout << "THRESHOLD AREA = " << threshold_area << std::endl;


	int bestInliers=0;
	Matrix2D<double> A_matrix, A_coarse;

	Matrix1D<double> T, T_coarse, def_T, t_test(2), u(2), dist_vec, dist_vec2;
	Matrix2D<double> invW;
	double trace_A, det_A, sqrt_aux, Eig_A1, Eig_A2, discriminant_A, dist, estimator;
	struct Point_T t_dist, t_closest;
	double cos_tilt = t_mean_area/u_mean_area;

	//std::cout << "cos tilt = " << cos_tilt << std::endl;
	//std::cout << "cos tilt2 = " << 1/cos_tilt << std::endl;
	//std::cout << "tilt = " << acos(cos_tilt)*180/PI << std::endl;

	A_matrix.initZeros(2,2);

	if (verbose==1)
	{
		std::cerr << "Exploring triangle matching" << std::endl;
		init_progress_bar(tri_number_untilt);
	}
	std::cout << "Untilt triangles = " << tri_number_untilt << std::endl;
	std::cout << "Tilt triangles = " << tri_number_tilt << std::endl;

	///////////// COARSE PHASE///////////////
	if (verbose==1)
		std::cerr << "Coarse Phase" << std::endl;

	for (int k=1; k<tri_number_untilt; k++)
	{
		std::cout << "Iteration k = " << k << std::endl;
		if (trig_untilt_area(k) < threshold_area)
			continue;

		for (int j=843; j<tri_number_tilt; j++)
		{
			std::cout << "Iteration j = " << j << std::endl;
			if (trig_untilt_area(k) > trig_tilt_area(j))
			{
				get_Face_Points(&delaunay_untilt, k, &u1, &u2, &u3);
				get_Face_Points(&delaunay_tilt, j, &t1, &t2, &t3);

				for (int i=0; i<6; i++)
				{
					if (i==0){
						def_affinity(u1.x, u1.y, u2.x, u2.y, u3.x, u3.y, t1.x, t1.y, t2.x, t2.y, t3.x, t3.y, A_matrix, T, invW);}
					if (i==1)
						def_affinity(u1.x, u1.y, u3.x, u3.y, u2.x, u2.y, t1.x, t1.y, t3.x, t3.y, t2.x, t2.y, A_matrix, T, invW);
					if (i==2)
						def_affinity(u2.x, u2.y, u1.x, u1.y, u3.x, u3.y, t2.x, t2.y, t1.x, t1.y, t3.x, t3.y, A_matrix, T, invW);
					if (i==3)
						def_affinity(u2.x, u2.y, u3.x, u3.y, u1.x, u1.y, t2.x, t2.y, t3.x, t3.y, t1.x, t1.y, A_matrix, T, invW);
					if (i==4)
						def_affinity(u3.x, u3.y, u1.x, u1.y, u2.x, u2.y, t3.x, t3.y, t1.x, t1.y, t2.x, t2.y, A_matrix, T, invW);
					if (i==5)
						def_affinity(u3.x, u3.y, u2.x, u2.y, u1.x, u1.y, t3.x, t3.y, t2.x, t2.y, t1.x, t1.y, A_matrix, T, invW);

					trace_A = MAT_ELEM(A_matrix, 0, 0) + MAT_ELEM(A_matrix, 1, 1);
					det_A = MAT_ELEM(A_matrix, 0, 0)*MAT_ELEM(A_matrix, 1, 1) - MAT_ELEM(A_matrix, 0, 1)*MAT_ELEM(A_matrix, 1, 0);

//					if ( fabs(det_A - cos_tilt)>0.1)
//						continue;

					discriminant_A = 0.25*(trace_A)*(trace_A) - det_A;

					if (discriminant_A < DBL_EPSILON)
						continue;

					sqrt_aux = sqrt(discriminant_A);
					Eig_A1 = trace_A/2 + sqrt_aux;
					Eig_A2 = trace_A/2 - sqrt_aux;

					if (Eig_A1<0.34 || Eig_A2<0.34 || fabs(Eig_A1-1)>0.05)
						continue;
					dist_vec.initConstant(len_u, -1);
					dist_vec2.initConstant(len_u, 10);
					std::cout << "len_u = " << len_u << std::endl;
					for (int tt=223; tt<len_u; tt++)
					{
						VEC_ELEM(u,0) = VEC_ELEM(ux,tt);
						VEC_ELEM(u,1) = VEC_ELEM(uy,tt);
						t_test = A_matrix*u + T;

						if (VEC_ELEM(t_test,0)<0 || VEC_ELEM(t_test,0)>Xdim || VEC_ELEM(t_test,1)<0 || VEC_ELEM(t_test,1)>Ydim)
							continue;

						t_dist.x = VEC_ELEM(t_test,0);
						t_dist.y = VEC_ELEM(t_test,1);
						std::cout << "Checkpoint 2 " << std::endl;
						std::cout << "Coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
						std::cout << "Closest " << select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) << std::endl;
						if (!select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) )
						{
							breakassignment = 1;
							//write_DCEL(delaunay_tilt.dcel, 0, "tiltdata_dcel.txt");
							std::cout << "Coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
							break;
						}
						std::cout << "Checkpoint 3 " << std::endl;
						VEC_ELEM(dist_vec,tt) = dist;
						std::cout << "dist = " << VEC_ELEM(dist_vec,tt) << std::endl;
						std::cout << "iteration tt = " << tt << std::endl;

					}

					int counter = 0;
					int inliers2 = 0;
					double dist2 = 0;

					for (int kk=0; kk<len_u; kk++)
					{
						if (VEC_ELEM(dist_vec,kk)>-1)
						{
							counter = counter + 1;
							if (VEC_ELEM(dist_vec,kk) < thr*particle_size)
							{
								dist2 = VEC_ELEM(dist_vec,kk) + dist2;
								inliers2 = inliers2 + 1;
							}
						}
					}

					if (counter == 0)
						continue;
					std::cout << "ultimo for " << std::endl;
					if (inliers2 > bestInliers)
					{
						estimator = dist2/inliers2;
						bestInliers = inliers2;
						A_coarse = A_matrix;
						T_coarse = T;
					}
					if ((inliers2 == bestInliers) && (dist2/inliers2 < estimator) )
					{
						estimator = dist2/inliers2;
						A_coarse = A_matrix;
						T_coarse = T;
					}
				}
			}
		//std::cout << "Iteration = " << k << ", " << j << std::endl;
		}
		if (verbose==1 && k%100==0)
			progress_bar(k);
	}


	std::cout << "Coarse Inliers = " << bestInliers << std::endl;

	////////////// REFINEMENT PHASE///////////////

	if (bestInliers >3)
	{
		if (verbose==1)
		std::cerr << "Refinement Phase" << std::endl;

		dist_vec.initConstant(len_u, -1);
		Matrix1D<double> t_cloy(bestInliers), u_clox(bestInliers), u_cloy(bestInliers);
		Matrix2D<double> invref, invW2;
		PseudoInverseHelper pseudoInverter;
		Matrix1D<double> t_clo(2*bestInliers);
		invref.initZeros(2*bestInliers,6);
		invW2.initZeros(2*bestInliers,6);
		int count = 0;

		Matrix1D<double> X;
		PseudoInverseHelper h;
		for (int k=0; k<len_u; k++)
		{
			VEC_ELEM(u,0) = VEC_ELEM(ux,k);
			VEC_ELEM(u,1) = VEC_ELEM(uy,k);
			t_test = A_coarse*u + T_coarse;
			t_dist.x = VEC_ELEM(t_test,0);
			t_dist.y = VEC_ELEM(t_test,1);
			//std::cout << "Checkpoint 2 " << std::endl;
			//std::cout << "Coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
			if (!select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) )
			{
				breakassignment = 1;
				std::cout << "Coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
				break;
			}
			//std::cout << "dist = " << dist << std::endl;
			if (dist<thr*particle_size)
			{
				VEC_ELEM(t_clo,count) = t_closest.x;
				VEC_ELEM(t_clo,count+bestInliers) = t_closest.y;
				MAT_ELEM(invW2,count,0) = VEC_ELEM(u,0);
				MAT_ELEM(invW2,count,1) = VEC_ELEM(u,1);
				MAT_ELEM(invW2,count,4) = 1;

				MAT_ELEM(invW2, count+bestInliers, 2) = VEC_ELEM(u,0);
				MAT_ELEM(invW2, count+bestInliers, 3) = VEC_ELEM(u,1);
				MAT_ELEM(invW2, count+bestInliers, 5) = 1;

				count = count + 1;
			}
			//std::cout << "vector = " << t_clo << std::endl;

		}
		if (breakassignment)
		{
			std::cerr << "ERROR IN TRIANGULATION OR CLOSEST NEIGHBOUR" << std::endl;
		}

		Matrix2D<double> def_A;
		def_A.initZeros(2,2);
		def_T.initZeros(2);


		if (count >= bestInliers)
		{
			//std::cout << "Count > bestInliers" << std::endl;
			h.A=invW2;

			h.b=t_clo;
			solveLinearSystem(h,X);
			MAT_ELEM(def_A,0,0) = VEC_ELEM(X,0);
			MAT_ELEM(def_A,0,1) = VEC_ELEM(X,1);
			MAT_ELEM(def_A,1,0) = VEC_ELEM(X,2);
			MAT_ELEM(def_A,1,1) = VEC_ELEM(X,3);
			VEC_ELEM(def_T,0) = VEC_ELEM(X,4);
			VEC_ELEM(def_T,1) = VEC_ELEM(X,5);
			std::cout << "Best fitting is gotten using refinement phase" << std::endl;
			std::cout << "Refinement Inliers = " << count << std::endl;
			std::cout << "Coarse Inliers = " << bestInliers << std::endl;
			std::cout << "A" << def_A << std::endl;
			std::cout << "---------------------------" << std::endl;
			std::cout << "T" << def_T << std::endl;
		}
		else
		{
			std::cout << "Best fitting is gotten using coarse phase" << std::endl;
			std::cout << "Coarse Inliers = " << bestInliers << std::endl;
			std::cout << "Refinement Inliers = " << count << std::endl;
			def_A = A_coarse;
			def_T = T_coarse;
			std::cout << "A" << def_A << std::endl;
			std::cout << "---------------------------" << std::endl;
			std::cout << "T" << def_T << std::endl;
		}

		//std::cout << "Checkpoint 1" << std::endl;
		//////WRITING TILT PAIRS/////////
		//std::cout << "checkpoint 1" << std::endl;
		if (!breakassignment)
		{
		//std::cout << "Checkpoint 2" << std::endl;
		for (int k=0; k<len_u; k++)
		{
			VEC_ELEM(u,0) = VEC_ELEM(ux,k);
			VEC_ELEM(u,1) = VEC_ELEM(uy,k);
			t_test = def_A*u + def_T;
			t_dist.x = VEC_ELEM(t_test,0);
			t_dist.y = VEC_ELEM(t_test,1);
			//std::cout << "Checkpoint 1" << std::endl;
			//std::cout << "Coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
			if (!select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) )
			{
				breakassignment = 1;
				break;
			}
			//std::cout << "Checkpoint 2" << std::endl;

			if (dist<thr*particle_size)
			{
				objId = mduntilt.addObject();
				//mduntilt.setValue(MDL_XCOOR,3,objId);
				mduntilt.setValue(MDL_XCOOR,(int) VEC_ELEM(u,0),objId);
				//std::cout << mduntilt;
				mduntilt.setValue(MDL_YCOOR,(int) VEC_ELEM(u,1),objId);
				//mduntilt.setValue(MDL_YCOOR,2,objId);

				objId = mdtilt.addObject();
				mdtilt.setValue(MDL_XCOOR,(int) t_closest.x,objId);
				mdtilt.setValue(MDL_YCOOR,(int) t_closest.y,objId);


			}
		}
		}
		if (breakassignment)
			std::cerr << "ERROR IN TRIANGULATION OR CLOSEST NEIGHBOUR" << std::endl;

	mduntilt.write((String)"particles@"+fndir+'/'+fnuntilt.getBaseName() + ".pos" );
	mdtilt.write((String)"particles@"+fndir+'/'+fntilt.getBaseName() + ".pos" );
	}
	else
	{
		std::cerr << " Matching not found" << std::endl;
	}
	delete_Delaunay(&delaunay_untilt);
	delete_Delaunay(&delaunay_tilt);

}
