/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

// Translated from MATLAB code by Yoel Shkolnisky

#include <mpi.h>
#include "mpi_performance_test.h"
#include <data/mask.h>
#include <data/metadata_extension.h>


// Empty constructor =======================================================
ProgPerformanceTest::ProgPerformanceTest(int argc, char **argv)
{
    node=new MpiNode(argc,argv);
    if (!node->isMaster())
        verbose=0;
}

// MPI destructor
ProgPerformanceTest::~ProgPerformanceTest()
{
    delete node;
}

// Read arguments ==========================================================
void ProgPerformanceTest::readParams()
{
    fnIn = getParam("-i");
}

// Show ====================================================================
void ProgPerformanceTest::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input:               " << fnIn << std::endl
    ;
}

// usage ===================================================================
void ProgPerformanceTest::defineParams()
{
    addUsageLine("Makes a rotational invariant representation of the image collection");
    addParamsLine("    -i <selfile>               : Selfile with experimental images");
    addExampleLine("mpirun -np 4 `which xmipp_mpi_image_rotational_pca` -i images.stk --oroot images_eigen --thr 4");
}

// Produce side info =====================================================
void ProgPerformanceTest::produceSideInfo()
{
    TimeStamp t0;
    annotate_time(&t0);
    MetaData MDin(fnIn);
    print_elapsed_time(t0,false);
}

// Run ====================================================================
void ProgPerformanceTest::run()
{
    show();
    if (!system("hostname"))
    	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
    produceSideInfo();
}
