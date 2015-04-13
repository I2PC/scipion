/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#include "classify_compare_classes.h"
#include <data/metadata_extension.h>

// Evaluate classes program -------------------------------------------
void ProgCompareClass::readParams()
{
	fnClass1=getParam("--i1");
	fnClass2=getParam("--i2");
	fnOut=getParam("-o");
	append=checkParam("--append");
}

void ProgCompareClass::defineParams()
{
    addUsageLine("Compare two classifications");
    addUsageLine("+The program provides information about which class of classification 1 corresponds to which class of classification 2");
    addUsageLine("+");
    addParamsLine(" --i1 <infile1>             : Metadata with the classification 1");
    addParamsLine(" --i2 <infile2>             : Metadata with the classification 2");
    addParamsLine(" -o <outfile>              : Output text file");
    addParamsLine(" [--append]                : Append text to output");
    addExampleLine("xmipp_classify_compare_classes --i1 2D/CL2D/run_001/results_level_00_classes.xmd --i2 2D/CL2D/run_001/results_level_01_classes.xmd -o comparison_level_00_01.txt");
}

void ProgCompareClass::show()
{
	if (!verbose)
		return;
	std::cout
	<< "Input1: " << fnClass1 << std::endl
	<< "Input2: " << fnClass2 << std::endl
	<< "Output: " << fnOut << std::endl;
}

void ProgCompareClass::run()
{
	MetaData MD1(formatString("classes@%s",fnClass1.c_str()));
	MetaData MD2(formatString("classes@%s",fnClass2.c_str()));
        std::vector<int> ref1, ref2;
        int aux;
        FOR_ALL_OBJECTS_IN_METADATA(MD1)
        {
            MD1.getValue(MDL_REF,aux,__iter.objId);
            ref1.push_back(aux);
        }
        FOR_ALL_OBJECTS_IN_METADATA(MD2)
        {
            MD2.getValue(MDL_REF,aux,__iter.objId);
            ref2.push_back(aux);
        }

	Matrix2D<int> comparisonMatrix(MD1.size(),MD2.size());
	Matrix1D<int> MD1classSize(MAT_YSIZE(comparisonMatrix)), MD2classSize(MAT_XSIZE(comparisonMatrix));

	// Read the size of the individual classes
	MetaData MDclass1, MDclass2;
	FOR_ALL_ELEMENTS_IN_MATRIX1D(MD1classSize)
	{
		MDclass1.read(formatString("class%06d_images@%s",ref1[i],fnClass1.c_str()));
		VEC_ELEM(MD1classSize,i)=MDclass1.size();
	}

	FOR_ALL_ELEMENTS_IN_MATRIX1D(MD2classSize)
	{
		MDclass2.read(formatString("class%06d_images@%s",ref2[i],fnClass2.c_str()));
		VEC_ELEM(MD2classSize,i)=MDclass2.size();
	}

	// Now compare the two classifications
	for (size_t i=0; i<MAT_YSIZE(comparisonMatrix); i++)
	{
		MDclass1.read(formatString("class%06d_images@%s",ref1[i],fnClass1.c_str()));
		for (size_t j=0; j<MAT_XSIZE(comparisonMatrix); j++)
		{
			MDclass2.read(formatString("class%06d_images@%s",ref2[j],fnClass2.c_str()));
			MDclass2.intersection(MDclass1,MDL_IMAGE);
			MAT_ELEM(comparisonMatrix,i,j)=MDclass2.size();
		}
	}

	// Report analysis
	std::ofstream fhOut;
	if (append)
	{
		fhOut.open(fnOut.c_str(),std::ios::app);
		fhOut << "\n\n------------------------------------------------------------------------\n";
	}
	else
		fhOut.open(fnOut.c_str());
	if (!fhOut)
		REPORT_ERROR(ERR_IO_NOWRITE,fnOut);
	fhOut << "Comparison of " << fnClass1 << " and " << fnClass2 << std::endl;
	fhOut << "Analysis of " << fnClass1 << " =======================\n";
	for (size_t i=0; i<MAT_YSIZE(comparisonMatrix); i++)
	{
		fhOut << "Class " << formatString("class%06d_images@%s",ref1[i],fnClass1.c_str()) << ": " << VEC_ELEM(MD1classSize,i) << " images\n";
		for (size_t j=0; j<MAT_XSIZE(comparisonMatrix); j++)
			if (MAT_ELEM(comparisonMatrix,i,j)>0)
				fhOut << "   " << 100.0*MAT_ELEM(comparisonMatrix,i,j)/VEC_ELEM(MD1classSize,i) << "% are in class " << formatString("class%06d_images@%s",j+1,fnClass2.c_str()) << std::endl;
	}
	fhOut << "\n\nAnalysis of " << fnClass2 << " =======================\n";
	for (size_t j=0; j<MAT_XSIZE(comparisonMatrix); j++)
	{
		fhOut << "Class " << formatString("class%06d_images@%s",ref2[j],fnClass2.c_str()) << ": " << VEC_ELEM(MD2classSize,j) << " images\n";
		for (size_t i=0; i<MAT_YSIZE(comparisonMatrix); i++)
			if (MAT_ELEM(comparisonMatrix,i,j)>0)
				fhOut << "   " << 100.0*MAT_ELEM(comparisonMatrix,i,j)/VEC_ELEM(MD2classSize,j) << "% are in class " << formatString("class%06d_images@%s",i+1,fnClass1.c_str()) << std::endl;
	}
	fhOut.close();
}
