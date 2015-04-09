///***************************************************************************
// * Authors:     Joaquin Oton (joton@cnb.csic.es)
// *
// *
// * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
// *
// * This program is free software; you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation; either version 2 of the License, or
// * (at your option) any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program; if not, write to the Free Software
// * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
// * 02111-1307  USA
// *
// *  All comments concerning this program package may be sent to the
// *  e-mail address 'xmipp@cnb.csic.es'
// ***************************************************************************/
//
//#include "art_xray.h"
//
//
//void XrayARTRecons::defineParams(XmippProgram * program, const char* prefix, const char* comment)
//{
//    char tempLine[256];
//
//    if(prefix == NULL)
//        sprintf(tempLine, "  [--xray <psf_param_file>]   : X-ray mode activation");
//    else
//        sprintf(tempLine,"%s   [--xray <psf_param_file>]   : X-ray mode activation", prefix);
//    if (comment != NULL)
//        sprintf(tempLine, "%s : %s", tempLine, comment);
//
//    program->addParamsLine(tempLine);
//
//
//    //    XRayPSF::defineParams(program);
//}
//
//void XrayARTRecons::readParams(XmippProgram * program)
//{
//    ARTReconsBase::readParams(program);
//    fnPSF = program->getParam("--xray");
//    psf.read(fnPSF);
//    //    psf.readParams(program);
//}
//
//void XrayARTRecons::preProcess(GridVolume &vol_basis0, int level, int rank)
//{
//    psf.calculateParams(artPrm.sampling*1.e-10); // sampling is read in angstrom
//
//    if (artPrm.basis.VolPSF == NULL)
//        artPrm.basis.VolPSF = new MultidimArray<double>;
//    psf.PSFGen().getImage(*artPrm.basis.VolPSF);
//
//    SinPartARTRecons::preProcess(vol_basis0);
//
//    //TODO: If Start volume is not loaded, then vol_basis (our mu in x-ray) must be estimated
//}
//
//void XrayARTRecons::print(std::ostream &o) const
//{
//    o << "X-rays information ------------------------------------------" << std::endl;
//    o << "Microscope parameters file: " << fnPSF.c_str() << std::endl;
//}
//
