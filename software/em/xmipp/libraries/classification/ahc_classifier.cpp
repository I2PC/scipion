/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#include <alglib/src/dataanalysis.h>

#include "ahc_classifier.h"


void AHCClassifier::clusterData(const Matrix2D<double> &X, int numberOfClusters, int distance, int linkageType)
{
	alglib::integer_1d_array cidx;
    try
    {
		// See example at http://www.alglib.net/translator/man/manual.cpp.html#example_clst_ahc
		alglib::real_2d_array xy;
		xy.setcontent(MAT_YSIZE(X),MAT_XSIZE(X),MATRIX2D_ARRAY(X));

		alglib::clusterizerstate s;
		alglib::ahcreport rep;

		// Run Hierarchical Clustering
		clusterizercreate(s);
		clusterizersetahcalgo(s, linkageType);
		clusterizersetpoints(s, xy, distance);
		clusterizerrunahc(s, rep);

		// Now get K clusters
		alglib::integer_1d_array cz;
		clusterizergetkclusters(rep, numberOfClusters, cidx, cz);
    }
    catch (alglib::ap_error e)
    {
        REPORT_ERROR(ERR_UNCLASSIFIED,e.msg);
    }

    // Collect results
    clusterAssigned.resizeNoCopy(MAT_YSIZE(X));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(clusterAssigned)
    	VEC_ELEM(clusterAssigned,i)=cidx[i];

    std::vector<int> dummy;
    cluster.clear();
    for (int k=0; k<numberOfClusters; ++k)
    	cluster.push_back(dummy);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(clusterAssigned)
    	cluster[VEC_ELEM(clusterAssigned,i)].push_back(i);
}

void AHCClassifier::clusterWithDistance(const Matrix2D<double> &D, int numberOfClusters, int linkageType)
{
	alglib::integer_1d_array cidx;
    try
    {
		// See example at https://www.tol-project.org/svn/tolp/OfficialTolArchiveNetwork/AlgLib/CppTools/source/alglib/manual.cpp.html#example_clst_kclusters
		alglib::real_2d_array d;
		d.setcontent(MAT_YSIZE(D),MAT_XSIZE(D),MATRIX2D_ARRAY(D));

		alglib::clusterizerstate s;
		alglib::ahcreport rep;

		// Run Hierarchical Clustering
		clusterizercreate(s);
		clusterizersetahcalgo(s, linkageType);
	    clusterizersetdistances(s, d, true);
		clusterizerrunahc(s, rep);

		// Now get K clusters
		alglib::integer_1d_array cz;
		clusterizergetkclusters(rep, numberOfClusters, cidx, cz);
    }
    catch (alglib::ap_error e)
    {
        REPORT_ERROR(ERR_UNCLASSIFIED,e.msg);
    }

    // Collect results
    clusterAssigned.resizeNoCopy(MAT_YSIZE(D));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(clusterAssigned)
    	VEC_ELEM(clusterAssigned,i)=cidx[i];

    std::vector<int> dummy;
    cluster.clear();
    for (int k=0; k<numberOfClusters; ++k)
    	cluster.push_back(dummy);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(clusterAssigned)
    	cluster[VEC_ELEM(clusterAssigned,i)].push_back(i);
}

