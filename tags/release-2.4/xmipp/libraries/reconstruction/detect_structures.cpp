/***************************************************************************
 *
 * Authors:     Manuel Sanchez Pau 
 *              Carlos Oscar Sanchez Sorzano
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

#include "detect_structures.h"
#include <data/fft.h>
#include <data/volume.h>
#include <data/histogram.h>
#include <data/filters.h>
#include <data/morphology.h>


/* Detect structures ------------------------------------------------------- */
void Prog_Detect_Structures_Param::read(int argc, char **argv)
{
    fnIn=getParameter(argc,argv,"-i");
    fnOut=getParameter(argc,argv,"-oroot");
    filterType=getParameter(argc,argv,"-type");
    if (filterType!="wall" && filterType!="filament")
        REPORT_ERROR(1,(std::string)"Unrecognized filter type "+filterType);
    sigma0=textToInteger(getParameter(argc,argv,"-sigma0","1"));
    sigmaF=textToInteger(getParameter(argc,argv,"-sigmaF","-1"));
    sigmaStep=textToInteger(getParameter(argc,argv,"-sigmaStep","-1"));
    angStep=textToInteger(getParameter(argc,argv,"-angStep","5"));
    if (sigmaF<0) sigmaF=sigma0;
    if (sigmaStep<0) sigmaStep=1;
    removeBackground=checkParameter(argc,argv,"-removeBackground");
    removeMissingWedge=checkParameter(argc,argv,"-missing");
    if (removeMissingWedge)
    {
        int i=paremeterPosition(argc, argv, "-missing");
        if (i+4 >= argc)
            REPORT_ERROR(1, "Not enough parameters behind -missing");
        rot1  = textToFloat(argv[i+1]);
        tilt1 = textToFloat(argv[i+2]);
        rot2  = textToFloat(argv[i+3]);
        tilt2 = textToFloat(argv[i+4]);
    }
}

void Prog_Detect_Structures_Param::show() const
{
    std::cout
        << "Input volume:         " << fnIn               << std::endl
        << "Output rootname:      " << fnOut              << std::endl
        << "Filter type:          " << filterType         << std::endl
        << "Sigma0:               " << sigma0             << std::endl
        << "SigmaF:               " << sigmaF             << std::endl
        << "SigmaStep:            " << sigmaStep          << std::endl
        << "AngStep:              " << angStep            << std::endl
        << "Remove Background:    " << removeBackground   << std::endl
        << "Remove Missing Wedge: " << removeMissingWedge << std::endl
    ;
    if (removeMissingWedge)
        std::cout << "Plane 1: " << rot1 << " " << tilt1 << std::endl
                  << "Plane 2: " << rot2 << " " << tilt2 << std::endl;
}

void Prog_Detect_Structures_Param::usage() const
{
    std::cout << "Usage:\n"
        << "    -i <volume>         : Input volume\n"
        << "    -oroot <rootname>   : Output rootname\n"
        << "    -type <string>      : Filter type: wall or filament\n"
        << "   [-sigma0 <s=1>]      : Initial width\n"
        << "   [-sigmaF <s=-1>]     : Final width\n"
        << "   [-sigmaStep <s=-1>]  : Width step\n"
        << "   [-angStep <ang=5>]   : Angular step\n"
        << "   [-removeBackground]  : Remove background\n"
        << "   [-missing <rot1> <tilt1> <rot2> <tilt2>] : Remove missing wedge\n"
    ;
}

void Prog_Detect_Structures_Param::run()
{
    // Produce side info
    VolumeXmipp Vin(fnIn), Vout, Vaux, Vsigma;
    for (double sigma=sigma0; sigma<=sigmaF; sigma+=sigmaStep)
    {
        std::cout << "Filtering with sigma=" << sigma << std::endl;
        Vaux()=Vin();
        MissingWedge *MW=NULL;
        if (removeMissingWedge)
        {
            MW=new MissingWedge();
            MW->rotPos=rot1;
            MW->tiltPos=tilt1;
            MW->rotNeg=rot2;
            MW->tiltNeg=tilt2;
        }
        
        Steerable *filter=new Steerable(sigma,Vaux(),angStep,filterType,MW);
        
        // Compute energy percentage
        double totalEnergy=Vaux().sum2();
        Vaux()*=Vaux();
        Vaux()/=totalEnergy;

        if (XSIZE(Vout())==0)
        {
            Vout()=Vaux();
            Vsigma().resize(Vout());
            Vsigma().initConstant(sigma0);
        }
        else
            FOR_ALL_ELEMENTS_IN_MATRIX3D(Vout())
            {
                double vout=Vout(k,i,j);
                double vaux=Vaux(k,i,j);
                if (vout<vaux)
                {
                    Vout(k,i,j)=vaux;
                    Vsigma(k,i,j)=sigma;
                }
            }
        delete filter;
    }
    Vout.write(fnOut+"_energy.vol");

    if (removeBackground)
    {
        Matrix3D<double> Voutmask=Vout();
        EntropyOtsuSegmentation(Voutmask);
        dilate3D(Voutmask,Vout(),18,0,2);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(Voutmask)
            if (Vout(k,i,j)<0.5)
                Vsigma(k,i,j)=0;
    }
    Vsigma.write(fnOut+"_width.vol");
}
