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

#include "assignation_tilt_pair.h"
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

void ProgassignationTiltPair::readParams()
{
	fnuntilt = getParam("--untiltcoor");
	fntilt = getParam("--tiltcoor");
	fnmic = getParam("--tiltmic");
	fndir = getParam("--odir");
	mshift = getDoubleParam("--maxshift");
	particle_size = getDoubleParam("--particlesize");
	thr = getDoubleParam("--threshold");
}

void ProgassignationTiltPair::defineParams()
{
	//TODO change the values <---->
	//usage
	addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
	//params

	addParamsLine("  [--untiltcoor <md_file=\"\">]    : Untilt coordinates");
	addParamsLine("  [--tiltcoor <md_file=\"\">]    : Tilt coordinates");
	addParamsLine("  [--tiltmic <md_file=\"\">]    : Tilt micrography");
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

void ProgassignationTiltPair::run()
{



	//LOAD METADATA and TILTPAIRS

	MetaData md_untilt, md_tilt, M_in;
	size_t Ndim, Zdim, Ydim , Xdim;

	M_in.read(fnmic);
	getImageSize(M_in,Xdim,Ydim,Zdim,Ndim);

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



	int tri_number_untilt = get_Number_Real_Faces(delaunay_untilt.dcel);


	struct Delaunay_T delaunay_tilt;
	init_Delaunay( &delaunay_tilt, md_tilt.size());
	FOR_ALL_OBJECTS_IN_METADATA(md_tilt)
	{
		md_tilt.getValue(MDL_XCOOR, x,__iter.objId);
		md_tilt.getValue(MDL_YCOOR, y,__iter.objId);
		insert_Point( &delaunay_tilt, x, y);
	}

	//write_DCEL(delaunay_tilt.dcel, 0, "tiltdata.txt");
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


	double threshold_area = sortedArea( round((tri_number_untilt+1)*0.2) );

	int bestInliers=0;
	Matrix2D<double> A, def_A;

	Matrix1D<double> T, def_T, t_test(2), u(2), dist_vec, dist_vec2;
	Matrix2D<double> invW;
	double trace_A, det_A, sqrt_aux, Eig_A1, Eig_A2, discriminant_A, dist, estimator;
	struct Point_T t_dist, t_closest;
	double cos_tilt = t_mean_area/u_mean_area;


	std::cout << "cos_tilt" << cos_tilt << std::endl;
	std::cout << "------------------------------"  << std::endl;


	A.initZeros(2,2);

	///////////// COARSE PHASE///////////////
	for (int k=5; k<tri_number_untilt; k++)
	{
		//std::cout << k << "/" << tri_number_untilt << std::endl; //<< "  J = " << tri_number_tilt << std::endl;
		if (trig_untilt_area(k) < threshold_area)
			continue;
		//std::cin.get();
		for (int j=0; j<tri_number_tilt; j++)
		{
			//std::cout << "j= " << j << std::endl;
			if (trig_untilt_area(k) > trig_tilt_area(j))
			{
				get_Face_Points(&delaunay_untilt, k, &u1, &u2, &u3);
				get_Face_Points(&delaunay_tilt, j, &t1, &t2, &t3);

				//std::cout << u1.x << " " << u1.y << " " << u2.x << " " << u2.y << " " << u3.x << " " << u3.y << std::endl;

				for (int i=0; i<6; i++)
				{
					//std::cout << "i = " << i << std::endl;
					if (i==0){
						//std::cout << u1.x << " " << u1.y << " " <<  u2.x << " " << u2.y << " " << u3.x << " " << u3.y << std::endl;
						def_affinity(u1.x, u1.y, u2.x, u2.y, u3.x, u3.y, t1.x, t1.y, t2.x, t2.y, t3.x, t3.y, A, T, invW);}
						//std::cout << "case0" << std::endl;}
					if (i==1)
						def_affinity(u1.x, u1.y, u3.x, u3.y, u2.x, u2.y, t1.x, t1.y, t3.x, t3.y, t2.x, t2.y, A, T, invW);
						//std::cout << "case1" << std::endl;}
					if (i==2)
						def_affinity(u2.x, u2.y, u1.x, u1.y, u3.x, u3.y, t2.x, t2.y, t1.x, t1.y, t3.x, t3.y, A, T, invW);
						//std::cout << "case2" << std::endl;}
					if (i==3)
						def_affinity(u2.x, u2.y, u3.x, u3.y, u1.x, u1.y, t2.x, t2.y, t3.x, t3.y, t1.x, t1.y, A, T, invW);
						//std::cout << "case3" << std::endl;}
					if (i==4)
						def_affinity(u3.x, u3.y, u1.x, u1.y, u2.x, u2.y, t3.x, t3.y, t1.x, t1.y, t2.x, t2.y, A, T, invW);
						//std::cout << "case4" << std::endl;}
					if (i==5)
						def_affinity(u3.x, u3.y, u2.x, u2.y, u1.x, u1.y, t3.x, t3.y, t2.x, t2.y, t1.x, t1.y, A, T, invW);
						//std::cout << "case5" << std::endl;}


					//std::cout << "Matrix" << A << std::endl;
					//std::cout << k << " " << j << std::endl;

					trace_A = MAT_ELEM(A, 0, 0) + MAT_ELEM(A, 1, 1);
					//std::cout << "traza = " << trace_A << std::endl;
					det_A = MAT_ELEM(A, 0, 0)*MAT_ELEM(A, 1, 1) - MAT_ELEM(A, 0, 1)*MAT_ELEM(A, 1, 0);
					//std::cout << "determinant = " << det_A << std::endl;

					if ( fabs(det_A - cos_tilt)>0.1)
						continue;

					discriminant_A = 0.25*(trace_A)*(trace_A) - det_A;
					//std::cout << "discriminant" << discriminant_A << std::endl;
					if (discriminant_A < DBL_EPSILON)
						continue;

					sqrt_aux = sqrt(discriminant_A);
					Eig_A1 = trace_A/2 + sqrt_aux;
					Eig_A2 = trace_A/2 - sqrt_aux;
					//std::cout << "Eig1" << Eig_A1 << std::endl;
					//std::cout << "Eig2" << Eig_A2 << std::endl;

					if (Eig_A1<0.34 || Eig_A2<0.34 || fabs(Eig_A1-1)>0.05)
						continue;
					dist_vec.initConstant(len_u, -1);
					dist_vec2.initConstant(len_u, 10);
					//std::cout << "ux" << VEC_ELEM(ux,400) << std::endl;
					//std::cout << "len_u = " << len_u << std::endl;
					for (int tt=0; tt<len_u; tt++)

					{
						VEC_ELEM(u,0) = VEC_ELEM(ux,tt);
						VEC_ELEM(u,1) = VEC_ELEM(uy,tt);
						t_test = A*u + T;
						//std::cout << "u = (" << VEC_ELEM(u,0) << "," << VEC_ELEM(u,1) << ")" << std::endl;
						//std::cout << "t = (" << t_dist.x << "," << t_dist.y << ")" << std::endl;
						// Condition I: Affinity must be inside the micrography
						//std::cout << "Xdim Y dim" << Xdim << " " << Ydim << std::endl;

						if (VEC_ELEM(t_test,0)<0 || VEC_ELEM(t_test,0)>Xdim || VEC_ELEM(t_test,1)<0 || VEC_ELEM(t_test,1)>Ydim)
							continue;


						t_dist.x = VEC_ELEM(t_test,0);
						t_dist.y = VEC_ELEM(t_test,1);
						//std::cout << "t = (" << t_dist.x << "," << t_dist.y << ")" << std::endl;
						//std::cout << "Entro1" << tt << std::endl;
						select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist);
						//std::cout << "Entro " << tt << std::endl;
						//std::cout << "dist" << dist << std::endl;
						//std::cout << "coordenada x" << t_closest.x << " coordenada y " << t_closest.y << std::endl;

						VEC_ELEM(dist_vec,tt) = dist;
						//exit (0);

					}

					int counter = 0;
					int inliers2 = 0;
					double dist2 = 0;
					//std::cout << "len_u = " << len_u << std::endl;
					for (int kk=0; kk<len_u; kk++)
					{
						//std::cout << "VEC_ELEM = " << VEC_ELEM(dist_vec,kk) << std::endl;
						if (VEC_ELEM(dist_vec,kk)>-1)
						{
							counter = counter + 1;
							//std::cout << "Pasa " << "inliers = " << inliers2 << "   " << "dist = " << VEC_ELEM(dist_vec,kk) << std::endl;
							if (VEC_ELEM(dist_vec,kk) < thr*particle_size)
							{
								dist2 = VEC_ELEM(dist_vec,kk) + dist2;
								inliers2 = inliers2 + 1;
								//std::cout << "inliers = " << inliers2 << std::endl;
							}
						}
					}
					if (counter == 0)
						continue;

					//mean(dist_vec2);
					if (inliers2 > bestInliers)
					{
						//estimator = mean(dist_vec2);
						estimator = dist2/inliers2;
						bestInliers = inliers2;
						std::cout << k << "/" << tri_number_untilt << std::endl;
						std::cout << "Inliers = " << bestInliers << std::endl;
						std::cout << "Estimator = " << estimator << std::endl;
						std::cout << "--------------------------------------------" << std::endl;
						def_A = A;
						def_T = T;
					}
					if ((inliers2 == bestInliers) && (dist2/inliers2 < estimator) )
					{
						estimator = dist2/inliers2;
						std::cout << k << "/" << tri_number_untilt << std::endl;
						std::cout << "Inliers = " << bestInliers << std::endl;
						std::cout << "Estimator = " << estimator << std::endl;
						std::cout << "--------------------------------------------" << std::endl;
						def_A = A;
						def_T = T;
					}
				}
			}
		}
	}

#ifdef NEVERDEFINED
	////////////// REFINEMENT PHASE///////////////
	dist_vec.initConstant(len_u, -1);
	Matrix1D<double> t_clo(2*bestInliers), t_cloy(bestInliers), u_clox(bestInliers), u_cloy(bestInliers);
	Matrix2D<double> invref;
	invref.initZeros(2*bestInliers,6);
	int count = 0;
	for (int k=0; k<len_u; k++)
	{
		VEC_ELEM(u,0) = VEC_ELEM(ux,k);
		VEC_ELEM(u,1) = VEC_ELEM(uy,k);
		t_test = def_A*u + def_T;

		t_dist.x = VEC_ELEM(t_test,0);
		t_dist.y = VEC_ELEM(t_test,1);

		select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist);

		if (dist<thr)
		{
			VEC_ELEM(t_clo,count) = t_closest.x;
			VEC_ELEM(t_clo,count+bestInliers+1) = t_closest.y;
			MAT_ELEM(invW,count,0) = VEC_ELEM(u,0);
			MAT_ELEM(invW,count,1) = VEC_ELEM(u,1);
			MAT_ELEM(invW,count,4) = 1;
			MAT_ELEM(invW, count+bestInliers+1, 2) = VEC_ELEM(u,0);
			MAT_ELEM(invW, count+bestInliers+1, 3) = VEC_ELEM(u,1);
			MAT_ELEM(invW, count+bestInliers+1, 5) = 1;

			count = count + 1;
		}
	}



	solveLinearSystem(PseudoInverseHelper &h, Matrix1D<double> &result);








#endif
//	delete_Delaunay( &delaunay_untilt);
//	delete_Delaunay( &delaunay_tilt);

}


