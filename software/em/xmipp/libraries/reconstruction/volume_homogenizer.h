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
#include <data/multidim_array.h>
#include <reconstruction/fourier_filter.h>
#include "opencv2/core/core.hpp"
#include "opencv2/video/video.hpp"



/**@defgroup Volume Homogenizer
   @ingroup ReconsLibrary */
//@{
class ProgVolumeHomogenizer: public XmippProgram
{
public:
    FileName fnVol, fnRef, fnSetOfImgIn, fnSetOfImgOut;

    FileName fnTestOut1, fnTestOut2, fnTestOut3;

    int winSize;

    double cutFreq;

    MetaData mdPartialParticles;

    size_t rank, Nprocessors;

public:

    ProgVolumeHomogenizer();

    void readParams();

    void defineParams();

    // Converts a XMIPP MultidimArray to OpenCV matrix
    void xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat);

    // Converts an OpenCV float matrix to an OpenCV Uint8 matrix
    void convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat);

    // Converts an OpenCV matrix to XMIPP MultidimArray
    void opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray);

    //Method to get two same volume but from different class with different conformation and correcting all images of one of the volume with respect
    //to the another one as a reference, using optical flow algorithm. This is to later merging the corrected images
    //to the images of the reference map to reconstruct a volume with improved resolution
    void run();

private:

    /// Gather alignment
    virtual void gatherResults() {}

    /// Synchronize with other processors
    virtual void synchronize() {}
};
//@}
