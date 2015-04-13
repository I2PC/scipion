/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#ifndef MPI_PROJECT_XR_H_
#define MPI_PROJECT_XR_H_

#include "parallel/xmipp_mpi.h"
#include "reconstruction/project_xray.h"

/**@defgroup ProgMPIXrayProject ProgMPIXrayProject
   @ingroup Programs */
//@{
/* Projection XR Program -------------------------------- */
/** Program class for the project program */
class ProgMPIXrayProject: public ProgXrayProject
{
    MpiNode *node;
public:

    ~ProgMPIXrayProject();
    void read(int argc, char** argv);

    void run();

protected:
    void defineParams();
};


/* Effectively project ===================================================== */
int PROJECT_mpi_XR_Effectively_project(
    ParametersProjectionTomography &prm,
    XrayProjPhantom &side,
    Projection &proj,
    XRayPSF &psf,
    MetaData &SF) ;

/** @} */
#endif /* MPI_PROJECT_XR_H_ */
