/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 * 				Javier Vargas (jvargas@cnb.csic.es)
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

#ifndef _PROG_IMAGE_SORT
#define _PROG_IMAGE_SORT

#include <data/xmipp_program.h>
#include <data/basic_pca.h>
#include <data/histogram.h>
#include <reconstruction/fringe_processing.h>
#include <data/filters.h>
#include <math.h>


/**@defgroup ImageSort Image sort by statistics
   @ingroup ReconsLibrary */
//@{
class ProgSortByStatistics: public XmippProgram
{
public:
    FileName fn, fn_out, fn_train;
    bool addToInput;

public:
    double cutoff, per;

    std::vector<PCAMahalanobisAnalyzer> pcaAnalyzer;
    MetaData SF, SFtrain;

public:
    void clear();

    void readParams();

    void defineParams();

    void processInputPrepare(MetaData &SF);

    //Method to obtain a zScore of each extracted particle. The process is based on a segmentation approach using
    //the Spiral Phase Transform (SPT), then the respective binary image is fitted to an ellipse and we obtain the major
    //and minor axis. Also a other coarser descriptor based on the image histogram is used to detect outliers. Additionaly,
    //the centroid of the normalized image (SPT) is used as a descriptor but with a very low weight in the zScore
    //JavierVargas : September 2012 BCU
    void processInprocessInputPrepareSPTH(MetaData &SF);

    void run();
};
//@}
#endif

