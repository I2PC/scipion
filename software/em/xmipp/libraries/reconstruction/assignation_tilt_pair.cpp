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

void ProgassignationTiltPair::readParams()
{
    fnuntilt = getParam("--untiltcoor");
    fntilt = getParam("--tiltcoor");
    fndir = getParam("--odir");
    mshift = getIntParam("--maxshift");
    particle_size = getIntParam("--particlesize");
    radius = getIntParam("--threshold");
}

void ProgassignationTiltPair::defineParams()
{
	//TODO change the values <---->
    //usage
    addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
    //params
    addParamsLine("  [--untiltcoor <md_file=\"\">]    : Volume to validate");
    addParamsLine("  [--tiltcoor <md_file=\"\">]    : Volume to validate");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--maxshift <s=1000>]   : Maximum shift");
    addParamsLine("  [--particlesize <p=100>]   : Particle size");
    addParamsLine("  [--threshold <d=0.3>]      : if the distance between two points is lesser than threshold*particlesize, the they will be the same point");
}

void ProgassignationTiltPair::run()
{

	//LOAD METADATA and TILTPAIRS
	MetaData md_untilt, md_tilt;
    md_untilt.read(fnuntilt);
    md_tilt.read(fntilt);

    int x, y;
    //storing points and creating Delaunay triangulations
    struct Delaunay_T delaunay_untilt;
    init_Delaunay( &delaunay_untilt, md_untilt.size()-1);
    FOR_ALL_OBJECTS_IN_METADATA(md_untilt)
    {
        md_untilt.getValue(MDL_XCOOR, x, __iter.objId);
        md_untilt.getValue(MDL_YCOOR, y, __iter.objId);
		insert_Point( &delaunay_untilt, x, y);
    }
    int tri_number_untilt = get_Number_Real_Faces(delaunay_untilt.dcel);
    create_Delaunay_Triangulation( &delaunay_untilt);

    struct Delaunay_T delaunay_tilt;
    init_Delaunay( &delaunay_tilt, md_tilt.size()-1);
    FOR_ALL_OBJECTS_IN_METADATA(md_tilt)
    {
        md_tilt.getValue(MDL_XCOOR, x,__iter.objId);
        md_tilt.getValue(MDL_YCOOR, y,__iter.objId);
        insert_Point( &delaunay_tilt, x, y);
    }
    int tri_number_tilt = get_Number_Real_Faces(delaunay_tilt.dcel);
    create_Delaunay_Triangulation( &delaunay_tilt);


    // Triangle areas
    struct Point_T p, q, r, u1, u2, u3, t1, t2, t3;
	MultidimArray<double> trig_untilt_area(tri_number_untilt+1), trig_tilt_area(tri_number_tilt+1), sortedArea;

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
	Matrix2D<double> A;
	Matrix1D<double> T;
	Matrix2D<double> invW;

	for (int k=0; k<tri_number_untilt; k++)
	{
		if (trig_untilt_area(k) < threshold_area)
			continue;
		for (int j=0; j<tri_number_untilt; j++)
		{
			if (trig_untilt_area(k) > trig_tilt_area(j))
			{
				get_Face_Points(&delaunay_untilt, k, &u1, &u2, &u3);
				get_Face_Points(&delaunay_tilt, j, &t1, &t2, &t3);
				for (int i=0; i<6; i++)
					{
						switch (i)
						{
							case 0:
							def_affinity(u1.x, u1.y, u2.x, u2.y, u3.x, u3.y, t1.x, t1.y, t2.x, t2.y, t3.x, t3.y, A, T, invW);
							case 1:
							def_affinity(u1.x, u1.y, u3.x, u3.y, u2.x, u2.y, t1.x, t1.y, t3.x, t3.y, t2.x, t2.y, A, T, invW);
							case 2:
							def_affinity(u2.x, u2.y, u1.x, u1.y, u3.x, u3.y, t2.x, t2.y, t1.x, t1.y, t3.x, t3.y, A, T, invW);
							case 3:
							def_affinity(u2.x, u2.y, u3.x, u3.y, u1.x, u1.y, t2.x, t2.y, t3.x, t3.y, t1.x, t1.y, A, T, invW);
							case 4:
							def_affinity(u3.x, u3.y, u1.x, u1.y, u2.x, u2.y, t3.x, t3.y, t1.x, t1.y, t2.x, t2.y, A, T, invW);
							case 5:
							def_affinity(u3.x, u3.y, u2.x, u2.y, u1.x, u1.y, t3.x, t3.y, t2.x, t2.y, t1.x, t1.y, A, T, invW);
						}
						//***************
					}
				std::cout << "caca" << std::endl;
			}
		}
	}




#ifdef NEVERDEFINED





	print_Point(&p);
	print_Point(&q);
	print_Point(&r);
	print_DCEL(delaunay_tilt.dcel);

	double dist;
	select_Closest_Point(&delaunay_tilt, &p, &q, &dist); // p es el punto de netrada y q el mas cercano dist su distancia.
														 //Si no encuentra ningun punto o falla entonces devuelve un cero.




    //Coordenadas triangulo untilt ==>
    //Coordenadas triangulo untilt ==>
    //Numero de triangulos en untilt ==> num_tri_u
    //Numero de triangulos en tilt ==> num_tri_t
    //////////////////////////



    //TODO change size
    //double arr_u[1000], arr_t[1000];
    std::vector<double> arr_u[num_tri_u], arr_t[num_tri_t];

    for (int k = 1; k < num_tri_u; k++)
    {
    	//DELAUNAY UNTILT TRIANGLE INDICES
    	triangle_area(1, 0, 2 ,1, 5, 7, &arr_u[k]);
    }

    for (int k = 1; k < num_tri_t; k++)
    {
    	//DELAUNAY TILT TRIANGLE INDICES
    	triangle_area(1, 0, 2 ,1, 5, 7, & arr_t[k]);
    }

    /////////// MEAN//////////
    std::vector<double> v;
    std::vector<double> diff(v.size());
    double mean;
    transform(v.begin(), v.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
    //////////////////////////

    ////////PERMUTATIONS/////////
    //Se pueden poner todas a mano, en vez de hacer una funcion que las calcule
    /////////////////////////////

#endif
    delete_Delaunay( &delaunay_untilt);
    delete_Delaunay( &delaunay_tilt);

}


