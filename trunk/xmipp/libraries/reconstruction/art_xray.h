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

#ifndef ART_XRAY_H_
#define ART_XRAY_H_

#include "project_XR.h"
#include "base_art_recons.h"

/**@defgroup ARTXray art_xray (ART for xrays)
   @ingroup ReconsLibrary */
//@{
/** ART+Xrays parameters.
    Here only those specific parameters for X-rays are found, the rest of
    parameters common with normal ART should be looked up in
    \ref BasicARTParameters */

class XrayARTRecons : public SinPartARTRecons
{
    FileName fnPSF;
    XRayPSF psf;

public:
    XrayARTRecons()
    {
    }
    virtual ~XrayARTRecons()
    {}

    static void defineParams(XmippProgram * program, const char* prefix=NULL, const char* comment=NULL);
    /** Read special parameters from command line. */
    void readParams(XmippProgram * proram);

    void preIterations(GridVolume &vol_basis0, int level = FULL, int rank = -1);

    void print(std::ostream &o)const;

    void singleStep(GridVolume &vol_in, GridVolume *vol_out,
                    Projection &theo_proj, Projection &read_proj,
                    int sym_no,
                    Projection &diff_proj, Projection &corr_proj, Projection &alig_proj,
                    double &mean_error, int numIMG, double lambda, int act_proj,
                    const FileName &fn_ctf, const MultidimArray<int> *maskPtr,
                    bool refine);

};


//@}

#endif /* ART_XRAY_H_ */
