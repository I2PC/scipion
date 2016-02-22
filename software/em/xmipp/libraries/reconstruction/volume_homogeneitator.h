/***************************************************************************
 * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
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


#include <data/xmipp_program.h>
#include <data/projection.h>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/gpu/gpu.hpp"


/**@defgroup Volume Homogeneitator
   @ingroup ReconsLibrary */
//@{
class ProgVolumeHomogeneitator: public XmippProgram
{
public:
    FileName fnVol, fnRef, fnSetOfImgIn, fnSetOfImgOut;

public:
    void readParams();

    void defineParams();

    // Converts a XMIPP MultidimArray to OpenCV matrix
    void xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat);

    // Converts an OpenCV matrix to XMIPP MultidimArray
    void opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray);

    //Method to get two same volume but from different class and correcting all images of one of the volume with respect
    //to the another one as a reference, using optical flow algorithm. This is to later merging the corrected images
    //to the images of the reference map to reconstruct a volume with better resolution
    //JavierVargas : February 2016 BCU
    void run();
};
//@}


