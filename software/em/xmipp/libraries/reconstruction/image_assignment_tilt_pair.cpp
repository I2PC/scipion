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
//#include <data/multidim_array.h>

void ProgassignmentTiltPair::readParams()
{
	fnuntilt = getParam("--untiltcoor");
	fntilt = getParam("--tiltcoor");
	fnmic = getParam("--tiltmicsize");
	fndir = getParam("--odir");
	mshift = getDoubleParam("--maxshift");
	tiltest = getDoubleParam("--tiltangle");
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
	addParamsLine("  [--tiltangle <s=-1>]   : Tilt angle estimation, the method will look for the assignment in the interval of [tiltangle-15º, til_angle+15º]");
	addParamsLine("  [--particlesize <p=100>]   : Particle size");
	addParamsLine("  [--threshold <d=0.3>]      : If the distance between two points is lesser than threshold*particlesize, they will be the same point");
}

void ProgassignmentTiltPair::search_affine_transform(float u1x, float u1y, float u2x, float u2y, float u3x, float u3y, float t1x,
		float t1y, float t2x, float t2y, float t3x, float t3y,
		Matrix1D<double> ux, Matrix1D<double> uy, size_t Xdim, size_t Ydim, struct Delaunay_T &delaunay_tilt,
		int &bestInliers, Matrix2D<double> &A_coarse, Matrix1D<double> &T_coarse)
{
	double estimator, dist;
	struct Point_T t_dist, t_closest;
	Matrix1D<double> T;
	Matrix2D<double> invW;
	Matrix2D<double> A_matrix;

	A_matrix.initZeros(2,2);

	for (int i=0; i<6; i++)
	{
		if (i==0){
			def_affinity(u1x, u1y, u2x, u2y, u3x, u3y, t1x, t1y, t2x, t2y, t3x, t3y, A_matrix, T, invW);}
		if (i==1)
			def_affinity(u1x, u1y, u2x, u2y, u3x, u3y, t1x, t1y, t3x, t3y, t2x, t2y, A_matrix, T, invW);
		if (i==2)
			def_affinity(u1x, u1y, u2x, u2y, u3x, u3y, t2x, t2y, t1x, t1y, t3x, t3y, A_matrix, T, invW);
		if (i==3)
			def_affinity(u1x, u1y, u2x, u2y, u3x, u3y, t2x, t2y, t3x, t3y, t1x, t1y, A_matrix, T, invW);
		if (i==4)
			def_affinity(u1x, u1y, u2x, u2y, u3x, u3y, t3x, t3y, t1x, t1y, t2x, t2y, A_matrix, T, invW);
		if (i==5)
			def_affinity(u1x, u1y, u2x, u2y, u3x, u3y, t3x, t3y, t2x, t2y, t1x, t1y, A_matrix, T, invW);

		double trace_A = MAT_ELEM(A_matrix, 0, 0) + MAT_ELEM(A_matrix, 1, 1);
		double det_A = MAT_ELEM(A_matrix, 0, 0)*MAT_ELEM(A_matrix, 1, 1) - MAT_ELEM(A_matrix, 0, 1)*MAT_ELEM(A_matrix, 1, 0);


//		if ( fabs(det_A - cos_tilt_max)>0.1)
//			continue;

		double discriminant_A = 0.25*(trace_A)*(trace_A) - det_A;

		if (discriminant_A < DBL_EPSILON)
			continue;

		double sqrt_aux = sqrt(discriminant_A);
		double Eig_A1 = trace_A/2 + sqrt_aux;
		double Eig_A2 = trace_A/2 - sqrt_aux;

		if (Eig_A1<0.34 || Eig_A2<0.34 || fabs(Eig_A1-1)>0.05)
			continue;

		Matrix1D<double> u(2), dist_vec, dist_vec2, t_test(2);
		size_t len_u=VEC_XSIZE(ux);

		dist_vec.initConstant(len_u, -1);
		dist_vec2.initConstant(len_u, 10);

		for (size_t tt=0; tt<len_u; tt++)
		{
			VEC_ELEM(u,0) = VEC_ELEM(ux,tt);
			VEC_ELEM(u,1) = VEC_ELEM(uy,tt);
			t_test = A_matrix*u + T;


			if (VEC_ELEM(t_test,0)<0 || VEC_ELEM(t_test,0)>Xdim || VEC_ELEM(t_test,1)<0 || VEC_ELEM(t_test,1)>Ydim)
				continue;

			t_dist.x = VEC_ELEM(t_test,0);
			t_dist.y = VEC_ELEM(t_test,1);

			if (!select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) )
			{

				write_DCEL(delaunay_tilt.dcel, 0, "tiltdata_dcel.txt");
				std::cerr << "WARNING IN TRIANGULATION OR CLOSEST NEIGHBOUR" << std::endl;
				printf("Maybe the warning involves the tilt coordinates ( %f , %f ) \n", t_dist.x, t_dist.y);
				continue;
			}
			VEC_ELEM(dist_vec,tt) = dist;
		}

		int counter = 0;
		size_t inliers2 = 0;
		double dist2 = 0;

		for (size_t kk=0; kk<len_u; kk++)
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

		if (inliers2 > bestInliers)
		{
			estimator = dist2/inliers2;
			bestInliers = inliers2;
			A_coarse = A_matrix;
			T_coarse = T;
			std::cout << bestInliers << std::endl;
		}
		else if ((inliers2 == bestInliers) && (dist2/inliers2 < estimator) )
		{
			estimator = dist2/inliers2;
			A_coarse = A_matrix;
			T_coarse = T;
		}
	}
}

void ProgassignmentTiltPair::run()
{
	std::cout << "Starting..." << std::endl;

	//LOAD METADATA and TRIANGULATIONS
	MetaData md_untilt, md_tilt, mduntilt, mdtilt;
	size_t Ndim, Zdim, Ydim , Xdim, objId;

	getImageSize(fnmic, Xdim, Ydim, Zdim, Ndim);

	bool exist1, exist2;
	bool flag_assign = true;

	exist1 = fnuntilt.existsTrim();
	exist2 = fntilt.existsTrim();



//	std::cout << "fnuntilt = " << fnuntilt << std::endl;
//	std::cout << "exist 1 = " << exist1 << std::endl;
//	std::cout << "fntilt = " << fntilt << std::endl;
//	std::cout << "exist 2 = " << exist2 << std::endl;
//	std::cout << "exist = " << (exist1 && exist2) << std::endl;


	if ((exist1 && exist2) == true)
	{
	md_untilt.read(fnuntilt);
	md_tilt.read(fntilt);

	int x, y, len_u=0, len_t=0;

	//storing untilted points and creating Delaunay triangulation
	struct Delaunay_T delaunay_untilt;
	Matrix1D<double> ux(md_untilt.size()), uy(md_untilt.size()), tx(md_tilt.size()), ty(md_tilt.size());
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
	create_Delaunay_Triangulation( &delaunay_untilt, 0);
	int tri_number_untilt = get_Number_Real_Faces(delaunay_untilt.dcel);

	//storing tilted points and creating Delaunay triangulation
	struct Delaunay_T delaunay_tilt;
	init_Delaunay( &delaunay_tilt, md_tilt.size());
	FOR_ALL_OBJECTS_IN_METADATA(md_tilt)
	{
		md_tilt.getValue(MDL_XCOOR, x,__iter.objId);
		md_tilt.getValue(MDL_YCOOR, y,__iter.objId);

		VEC_ELEM(tx,len_t) = x;
		VEC_ELEM(ty,len_t) = y;

		insert_Point( &delaunay_tilt, x, y);
		len_t = len_t +1;
	}
	int thrs;
	if (len_u<len_t)
	{
		 thrs = len_u;
	}
	else
	{
		thrs = len_t;
	}

	if (len_u == 0 || len_t == 0)
	{
		flag_assign = false;
	}

	if (flag_assign == true)
	{
	create_Delaunay_Triangulation( &delaunay_tilt, 1);
	int tri_number_tilt = get_Number_Real_Faces(delaunay_tilt.dcel);
	//////////////DATA LOAD AND TRIANGULATIONS FINISHED

	//DETERMINING TRIANGLES AREAS
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

	trig_untilt_area.sort(sortedArea);

	double threshold_area = sortedArea( round((tri_number_untilt+1)*0.2) );
	//////////////TRIANGLES ARES FINISHED

	//DEFINING PARAMETERS FOR SEARCH_AFFINE_TRANSFORM
	int bestInliers=0, bestInliers_con=0;
	Matrix2D<double> A_coarse;
	Matrix1D<double> T_coarse, def_T, t_test(2), u(2), dist_vec, dist_vec2;

	double trace_A, det_A, sqrt_aux, Eig_A1, Eig_A2, discriminant_A, dist, estimator;
	struct Point_T t_dist, t_closest;
	//////////// DEFINING PARAMETERS FOR SEARCH_AFFINE_TRANSFORM

	//TILT ANGLE ESTIMATION
	double tilt_max = (tiltest + 15)*PI/180, tilt_min = (tiltest - 15)*PI/180, cos_tilt_max, cos_tilt_min;

	if (tiltest == -1)
	{
		cos_tilt_max = 0;
		cos_tilt_min = 1;
	}
	else
	{
		cos_tilt_max= cos(tilt_max);
		cos_tilt_min= cos(tilt_min);
	}
	////////////TILT ANGLE ESTIMATION


	if (verbose==1)
	{
		std::cerr << "Exploring triangle matching" << std::endl;
		init_progress_bar(tri_number_untilt);
	}


	/////////////////////////////////////// COARSE PHASE ///////////////////////////////////////
	if (verbose==1)
		std::cerr << "Coarse Phase" << std::endl;

	for (int k=0; k<tri_number_untilt; k++)
	{
		if (trig_untilt_area(k) < threshold_area)
			continue;

		for (int j=0; j<tri_number_tilt; j++)
		{
			//std::cout << "Iteration k = " << k << "     Iteration j = " << j << std::endl;
			if ( (trig_untilt_area(k) < trig_tilt_area(j)) )
				continue;

			if ( (trig_untilt_area(k)*cos_tilt_min < trig_tilt_area(j)) || (trig_untilt_area(k)*cos_tilt_max > trig_tilt_area(j)) )
				continue;

			get_Face_Points(&delaunay_untilt, k, &u1, &u2, &u3);
			get_Face_Points(&delaunay_tilt, j, &t1, &t2, &t3);

			search_affine_transform(u1.x, u1.y, u2.x, u2.y, u3.x, u3.y, t1.x,
					t1.y, t2.x, t2.y, t3.x, t3.y,
					ux, uy, Xdim, Ydim, delaunay_tilt, bestInliers,
					A_coarse, T_coarse);
		}
		if (verbose==1 && k%100==0)
			progress_bar(k);
	}
	std::cout << "Coarse Inliers = " << bestInliers << std::endl;


	///////////////////////////////////// REFINEMENT PHASE /////////////////////////////////////
		if (verbose==1)
		std::cerr << "Refinement Phase" << std::endl;

		dist_vec.initConstant(len_u, -1);
		Matrix1D<double> t_clo(2*bestInliers), X;
		Matrix2D<double> invref, invW2;
		//PseudoInverseHelper pseudoInverter;
		invref.initZeros(2*bestInliers,6);
		invW2.initZeros(2*bestInliers,6);
		int count = 0;

		PseudoInverseHelper h;
		for (int k=0; k<len_u; k++)
		{
			VEC_ELEM(u,0) = VEC_ELEM(ux,k);
			VEC_ELEM(u,1) = VEC_ELEM(uy,k);
			t_test = A_coarse*u + T_coarse;
			t_dist.x = VEC_ELEM(t_test,0);
			t_dist.y = VEC_ELEM(t_test,1);
			if (!select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) )
			{
				std::cerr << "WARNING IN TRIANGULATION OR CLOSEST NEIGHBOUR (in Coarse Phase)" << std::endl;
				std::cout << "Maybe the warning involves the tilt coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
				continue;
			}
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

		}


		Matrix2D<double> def_A, A_con;
		Matrix1D<double> T_con;
		T_con.initZeros(2);
		A_con.initZeros(2,2);
		def_A.initZeros(2,2);
		def_T.initZeros(2);
		bool flag;


		if ((count >= bestInliers) && (count >= 0.2*thrs))
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
			flag = true;
		}
		if (((count <= bestInliers) && (bestInliers >= 0.2*thrs)) && !flag)
		{
			std::cout << "Best fitting is gotten using coarse phase" << std::endl;
			std::cout << "Coarse Inliers = " << bestInliers << std::endl;
			std::cout << "Refinement Inliers = " << count << std::endl;
			def_A = A_coarse;
			def_T = T_coarse;
			std::cout << "A" << def_A << std::endl;
			std::cout << "---------------------------" << std::endl;
			std::cout << "T" << def_T << std::endl;
			flag = true;
		}
		if (((count <=bestInliers || count >=bestInliers) && ((count < 0.2*thrs) || (bestInliers < 0.2*thrs))) && !flag)
		{
			std::cerr << " Matching not found" << std::endl;
			std::cerr << " Starting contingency method" << std::endl;
			std::cerr << " This phase can take a long time" << std::endl;
			//Defining triangles in untilted space
			double window=4*particle_size;
			double window_t=window *(1-0.34);  //0.34 (approx 0.3) is cos(70), tilts higher than 70 degrees are not considered
			std::vector<Matrix1D<double> > allTrianglesU;
			Matrix1D<double> triangleu(7);
			for (int k = 0; k< len_u; k++)
			{
				VEC_ELEM(triangleu,0)=VEC_ELEM(ux,k);
				VEC_ELEM(triangleu,1)=VEC_ELEM(uy,k);
				for (int j=k+1; j<len_u; j++)
				{
					if ( ( (VEC_ELEM(ux,j)>VEC_ELEM(ux,k)+window) || (VEC_ELEM(ux,j)<VEC_ELEM(ux,k)-window)) ||
					     ( (VEC_ELEM(uy,j)>VEC_ELEM(uy,k)+window) || (VEC_ELEM(uy,j)<VEC_ELEM(uy,k)-window)) )
						continue;
					VEC_ELEM(triangleu,2)=VEC_ELEM(ux,j);
					VEC_ELEM(triangleu,3)=VEC_ELEM(uy,j);
					for (int l=j+1; l<len_u; l++)
					{
						if ( ( (VEC_ELEM(ux,l)>VEC_ELEM(ux,k)+window) || (VEC_ELEM(ux,l)<VEC_ELEM(ux,k)-window)) ||
						     ( (VEC_ELEM(uy,l)>VEC_ELEM(uy,k)+window) || (VEC_ELEM(uy,l)<VEC_ELEM(uy,k)-window)) )
							continue;

						VEC_ELEM(triangleu,6) = triangle_area(VEC_ELEM(ux,k), VEC_ELEM(uy,k), VEC_ELEM(ux,j), VEC_ELEM(uy,j), VEC_ELEM(ux,l), VEC_ELEM(uy,l));
						if (VEC_ELEM(triangleu,6) < threshold_area)
							continue;

						VEC_ELEM(triangleu,4) = VEC_ELEM(ux,l);
						VEC_ELEM(triangleu,5) = VEC_ELEM(uy,l);
						allTrianglesU.push_back(triangleu);
					}
				}
			}
			//Defining triangles in tilted space
			std::vector<Matrix1D<double> > allTrianglesT;
			Matrix1D<double> trianglet(7);
			for (int k = 0; k< len_t; k++)
			{
				VEC_ELEM(trianglet,0)=VEC_ELEM(tx,k);
				VEC_ELEM(trianglet,1)=VEC_ELEM(ty,k);
				for (int j=k+1; j<len_t; j++)
				{
					if ( ( (VEC_ELEM(tx,j)>VEC_ELEM(tx,k)+window_t) || (VEC_ELEM(tx,j)<VEC_ELEM(tx,k)-window_t)) ||
						 ( (VEC_ELEM(ty,j)>VEC_ELEM(ty,k)+window_t) || (VEC_ELEM(ty,j)<VEC_ELEM(ty,k)-window_t)) )
						continue;
					VEC_ELEM(trianglet,2)=VEC_ELEM(tx,j);
					VEC_ELEM(trianglet,3)=VEC_ELEM(ty,j);
					for (int l=j+1; l<len_t; l++)
					{
						if ( ( (VEC_ELEM(tx,l)>VEC_ELEM(tx,k)+window_t) || (VEC_ELEM(tx,l)<VEC_ELEM(tx,k)-window_t)) ||
							 ( (VEC_ELEM(ty,l)>VEC_ELEM(ty,k)+window_t) || (VEC_ELEM(ty,l)<VEC_ELEM(ty,k)-window_t)) )
							continue;
						VEC_ELEM(trianglet,6) = triangle_area(VEC_ELEM(tx,k), VEC_ELEM(ty,k), VEC_ELEM(tx,j), VEC_ELEM(ty,j), VEC_ELEM(tx,l), VEC_ELEM(ty,l));
						VEC_ELEM(trianglet,4) = VEC_ELEM(tx,l);
						VEC_ELEM(trianglet,5) = VEC_ELEM(ty,l);
						allTrianglesT.push_back(trianglet);
					}
				}
			}

			//Searching the affine application

			for (size_t ku=0; ku<allTrianglesU.size(); ku++)
			{

				const Matrix1D<double> &triangleU=allTrianglesU[ku];

				for (size_t kt=0; kt<allTrianglesT.size(); kt++)
				{
					const Matrix1D<double> &triangleT=allTrianglesT[kt];

					if (VEC_ELEM(triangleU,6)<VEC_ELEM(triangleT,6))
						continue;

					search_affine_transform(VEC_ELEM(triangleU,0), VEC_ELEM(triangleU,1), VEC_ELEM(triangleU,2), VEC_ELEM(triangleU,3),
											VEC_ELEM(triangleU,4), VEC_ELEM(triangleU,5), VEC_ELEM(triangleT,0), VEC_ELEM(triangleT,1),
											VEC_ELEM(triangleT,2), VEC_ELEM(triangleT,3), VEC_ELEM(triangleT,4), VEC_ELEM(triangleT,5),
											ux, uy, Xdim, Ydim, delaunay_tilt, bestInliers_con,
											A_con, T_con);
				}
			}
			std::cout << "Contingency Inliers = " << bestInliers_con << std::endl;
			if ((bestInliers_con > bestInliers) && (bestInliers_con>count))
			{
				std::cout << "Best fitting is gotten using contingency phase" << std::endl;
				std::cout << "Contingency Inliers = " << bestInliers_con << std::endl;
				std::cout << "Refinement Inliers = " << count << std::endl;
				def_A = A_con;
				def_T = T_con;
				std::cout << "A" << def_A << std::endl;
				std::cout << "---------------------------" << std::endl;
				std::cout << "T" << def_T << std::endl;
			}
//			else
//			{
//				def_A = A_coarse;
//				def_T = T_coarse;
//			}
		}

		///////////////////////////////////// WRITING TILT PAIRS /////////////////////////////////////

		for (int k=0; k<len_u; k++)
		{
			VEC_ELEM(u,0) = VEC_ELEM(ux,k);
			VEC_ELEM(u,1) = VEC_ELEM(uy,k);
			t_test = def_A*u + def_T;
			t_dist.x = VEC_ELEM(t_test,0);
			t_dist.y = VEC_ELEM(t_test,1);
			if (!select_Closest_Point(&delaunay_tilt, &t_dist, &t_closest, &dist) )
			{
				std::cerr << "WARNING IN TRIANGULATION OR CLOSEST NEIGHBOUR (in writing phase)" << std::endl;
				std::cout << "Maybe the warning involves the tilt coordinates " << "(" << t_dist.x << ", " << t_dist.y << ")" << std::endl;
				continue;
			}

			if (dist<thr*particle_size)
			{
				objId = mduntilt.addObject();
				mduntilt.setValue(MDL_XCOOR,(int) VEC_ELEM(u,0),objId);
				mduntilt.setValue(MDL_YCOOR,(int) VEC_ELEM(u,1),objId);

				objId = mdtilt.addObject();
				mdtilt.setValue(MDL_XCOOR,(int) t_closest.x,objId);
				mdtilt.setValue(MDL_YCOOR,(int) t_closest.y,objId);
			}
		}
		delete_Delaunay(&delaunay_untilt);
		delete_Delaunay(&delaunay_tilt);
	}
	}
	else
	{
		std::cout << "WARNING: input pos. file does not exist, therefore output files are empty" << std::endl;
	}
	if (flag_assign == false)
	{
		std::cout << "WARNING: input pos file is empty." << std::endl;
	}

	mduntilt.write((String)"particles@"+fndir+'/'+fnuntilt.getBaseName() + ".pos" );
	mdtilt.write((String)"particles@"+fndir+'/'+fntilt.getBaseName() + ".pos" );
	/////////////////////////////////////////////////////////////////////////////////////////////////

}
