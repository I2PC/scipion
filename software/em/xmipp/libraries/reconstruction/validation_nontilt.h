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

#ifndef VALIDATION_NONTILT_H_
#define VALIDATION_NONTILT_H_
#define PI 3.14159265

#include <data/xmipp_program.h>
#include <math.h>

/**@defgroup Validation without tilt
   @ingroup ReconsLibrary */
//@{
class ProgValidationNonTilt: public XmippProgram
{


public:
    /** Filenames */
    FileName fnDir, fnSym, fnInit, fnParticles;

    MetaData mdPartial;

    size_t rank, Nprocessors;

    bool useSignificant;

    double significance_noise;

public:

    void readParams();

    void defineParams();

    void run();

public:

    ProgValidationNonTilt();

    void obtainSumU(const MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0);

    void obtainSumU_2(const MetaData & mdGallery, const MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0);

    void obtainSumW(const MetaData & tempMd, double & sum_W, std::vector<double> & sum_u, std::vector<double> & H, const double factor);

    /// Gather alignment
    virtual void gatherClusterability() {}

    /// Synchronize with other processors
    virtual void synchronize() {}

};
#endif /* VALIDATION_NONTILT_H_ */
