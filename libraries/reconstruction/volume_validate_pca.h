/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#ifndef VOLUME_VALIDATE_PCA_H_
#define VOLUME_VALIDATE_PCA_H_


#endif /* VOLUME_VALIDATE_PCA_H_ */

#include <data/xmipp_program.h>
#include <data/metadata.h>
#include "reconstruct_fourier.h"
#include "angular_project_library.h"
#include "dimred/pca.h"

/**@defgroup VolumeValidatePCA Validation of volume consistency with respect to the provided classes
   @ingroup ReconsLibrary */
//@{

/** PCA validation volume parameters. */
class ProgVolumeValidationPCA: public XmippProgram
{
public:
    /** Filenames */
    FileName fnClasses, fnOut,fnSym, fnAngles;

    /** Number of intermediate volumes to generate*/
    int NVols;

    /** Number of classes to generate the intermediate volumes*/
    int NClasses;

private:
    size_t xdim, ydim, zdim, ndim;
    FileName fnVol;

public: // Internal members
    MetaData mdClasses, mdReconstruction;

  public:
    /// Read arguments from command line
    void readParams();

    /// Read arguments from command line
    void defineParams();

    /** Show. */
    void show();

    /** Run. */
    void run();

    /// Produce side info: fill arrays with relevant transformation matrices
    void produceSideinfo();

    /// Reconstruct current volume
    void reconstructCurrent();

    /// Reconstruct volume with all projections
    void reconstruct();



};
//@}


