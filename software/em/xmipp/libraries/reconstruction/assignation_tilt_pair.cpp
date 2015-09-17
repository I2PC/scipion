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
    addParamsLine("  [--maxshift <md_file=\"\">]   : Maximum shift");
    addParamsLine("  [--particlesize <p=100>]   : Particle size");
    addParamsLine("  [--threshold <d=0.3>]      : Sampling rate in A/px");
}


void ProgassignationTiltPair::triangle_area(double x1, double y1,
		double x2, double y2, double x3, double y3, double &trigarea)
{
    trigarea = ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))/2.0;
    (trigarea > 0.0) ? trigarea : -trigarea;
}


void ProgassignationTiltPair::run()
{

	//LOAD METADATA and TILTPAIRS
	MetaData md_untilt, md_tilt;
    md_untilt.read(fnuntilt);
    md_tilt.read(fntilt);

    //Metadadata size
    int sis_un, sis_t;
    sis_un = md_untilt.size();
    sis_t = md_tilt.size();

    int x_untilt[sis_un], y_untilt[sis_un], x_tilt[sis_t], y_tilt[sis_t];
    int x_aux, y_aux;

    int counter_u=0, counter_t=0;

    //metadata coordinates are stored in arrays x_untilt, y_untilt
    FOR_ALL_OBJECTS_IN_METADATA(md_untilt)
    {
        md_untilt.getValue(MDL_XCOOR, x_aux,__iter.objId);
        md_untilt.getValue(MDL_YCOOR, y_aux,__iter.objId);
        x_untilt[counter_u] = x_aux;
        y_untilt[counter_u] = y_aux;
        counter_u = counter_u + 1;
    }

    FOR_ALL_OBJECTS_IN_METADATA(md_tilt)
    {
        md_tilt.getValue(MDL_XCOOR, x_aux,__iter.objId);
        md_tilt.getValue(MDL_YCOOR, y_aux,__iter.objId);
        x_tilt[counter_t] = x_aux;
        y_tilt[counter_t] = y_aux;
        counter_t = counter_t + 1;
    }

    //DELAUNAY TRIANGULATION//
    struct Delaunay_T delaunay_untilt;
    struct Delaunay_T delaunay_tilt;

    double x, y;
    int tri_number_untilt, tri_number_tilt;
    //Untilt triangulation
    if (init_Delaunay( &delaunay_untilt, counter_u-1) == 1)
    {
    	for (int i=0; i<counter_u-1 ;i++)
    	{
    		x = x_untilt[i];
    		y = y_untilt[i];
    		insert_Point( &delaunay_untilt, x, y);
    	}
    	tri_number_untilt = get_Number_Real_Faces( delaunay_untilt.dcel);
    	create_Delaunay_Triangulation( &delaunay_untilt);
    }

    //Tilt triangulation
    if (init_Delaunay( &delaunay_tilt, counter_t-1) == 1)
    {
    	for (int i=0; i<counter_t-1 ;i++)
    	{
    		x = (double) x_tilt[i];
    		y = (double) y_tilt[i];

    		insert_Point( &delaunay_tilt, x, y);
    	}
    	create_Delaunay_Triangulation( &delaunay_tilt);
    	tri_number_tilt = get_Number_Real_Faces(delaunay_tilt.dcel); //Number of triangles in delaunay decomposition
    }

	struct Point_T p, q, r;

	for (int i=1; i<tri_number_untilt+1 ;i++)
	{
		get_Face_Points(&delaunay_tilt, 1, &p, &q, &r); //10 es el "numero" de triangulo
		triangle_area(p.x, p.y,
				double x2, double y2, double x3, double y3, double &trigarea)

	}


	print_Point(&p);
	print_Point(&q);
	print_Point(&r);
	print_DCEL(delaunay_tilt.dcel);

	double dist;
	select_Closest_Point(&delaunay_tilt, &p, &q, &dist); // p es el punto de netrada y q el mas cercano dist su distancia.
														 //Si no encuentra ningun punto o falla entonces devuelve un cero.


    #ifdef NEVERDEFINED

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


