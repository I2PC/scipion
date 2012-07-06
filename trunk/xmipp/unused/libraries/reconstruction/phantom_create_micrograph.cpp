/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "phantom_create_micrograph.h"
#include <data/projection.h>

#include <data/args.h>
#include <data/metadata.h>
#include <data/micrograph.h>

/* Empty constructor ------------------------------------------------------- */
ProgPhantomCreateMicrograph::ProgPhantomCreateMicrograph()
{
    microscope.command_line=false;
}

/* Read parameters --------------------------------------------------------- */
void ProgPhantomCreateMicrograph::read(int argc, char **argv)
{
    fn_vol=getParameter(argc,argv,"-vol");
    fn_root=getParameter(argc,argv,"-o","micrograph");
    Xdim=textToInteger(getParameter(argc,argv,"-dim","2048"));
    Nmicrographs=textToInteger(getParameter(argc,argv,"-N","1"));
    density=textToFloat(getParameter(argc,argv,"-density","10"));
    microscope.read(argc,argv);
}

/* Usage ------------------------------------------------------------------- */
void ProgPhantomCreateMicrograph::usage()
{
    std::cerr << "Usage: xmipp_phantom_create_micrograph\n"
              << "   -vol <Volume>            : Volume for the micrographs\n"
              << "  [-o <rootname=micrograph>]: Rootname for the micrographs\n"
              << "  [-dim <Xdim=2048>]        : Dimension of the micrographs\n"
              << "  [-N <N=1>]                : Number of micrographs\n"
              << "  [-density <d=10>]         : Density of particles (%)\n"
    ;
    microscope.usage();
}

/* Show -------------------------------------------------------------------- */
void ProgPhantomCreateMicrograph::show()
{
    std::cout
        << "Volume:       " << fn_vol       << std::endl
        << "Rootname:     " << fn_root      << std::endl
        << "Dimension:    " << Xdim         << std::endl
        << "Nmicrographs: " << Nmicrographs << std::endl
        << "Density:      " << density      << std::endl
    ;
    microscope.show();
}

/* Produce side information ------------------------------------------------ */
void ProgPhantomCreateMicrograph::produce_side_info()
{
    V.read(fn_vol);
    V().setXmippOrigin();
    microscope.Xdim=Xdim;
    microscope.Ydim=Xdim;
    microscope.produceSideInfo();
    Nproj=FLOOR(density/100.0*(double)(Xdim*Xdim)/(XSIZE(V())*XSIZE(V())));
}

/* Run -------------------------------------------------------------------.. */
void ProgPhantomCreateMicrograph::run()
{
    for (int n=0; n<Nmicrographs; n++)
    {
        // Rootname for this micrograph
        FileName dir_micrograph=fn_root+integerToString(n,3);
    
        // Create an empty micrograph ......................................
        system(((std::string)"mkdir "+dir_micrograph).c_str());
        FileName fn_micrograph=dir_micrograph;
        FileName fn_micrograph_clean=fn_micrograph+"_noiseless";

        Image<double> Md;
        Md().resize(Xdim,Xdim);
        Micrograph M;
        
        // Create the clean projections ....................................
        MetaData DF_out_clean;
        system(((std::string)"mkdir "+dir_micrograph+"/"+fn_micrograph_clean).c_str());
        FileName fn_proj_clean=fn_micrograph_clean+"/img_noiseless";
        int Xproj=XSIZE(V());
        std::cerr << "Creating micrograph " << n << "...\n";
        init_progress_bar(Nproj);
        for (int l=0; l<Nproj; l++)
        {
            // Create projection
            Projection P;
            double rot=rnd_unif(0,360);
            double tilt=rnd_unif(0,180);
            double psi=rnd_unif(0,360);
            projectVolume(V(), P, Xproj, Xproj, rot, tilt, psi);
            FileName fn_image=fn_proj_clean+"_"+integerToString(l,6)+".xmp";
            P.write(dir_micrograph+"/"+fn_image);
            DF_out_clean.addObject();
            DF_out_clean.setValue(MDL_IMAGE,fn_image);
            DF_out_clean.setValue(MDL_ANGLE_ROT,rot);
            DF_out_clean.setValue(MDL_ANGLE_TILT,tilt);
            DF_out_clean.setValue(MDL_ANGLE_PSI,psi);
            
            // Place this image in the micrograph
            int X=rnd_unif(Xproj,Xdim-Xproj);
            int Y=rnd_unif(Xproj,Xdim-Xproj);
            while (M.search_coord_near(X,Y,Xproj)!=-1)
            {
                X=rnd_unif(Xproj,Xdim-Xproj);
                Y=rnd_unif(Xproj,Xdim-Xproj);
            }
            FOR_ALL_ELEMENTS_IN_ARRAY2D(P())
                Md(Y+i,X+j)=P(i,j);
            M.add_coord(X,Y,0,0);
            
            progress_bar(l);
        }
        progress_bar(Nproj);
        Md.write(dir_micrograph+"/"+fn_micrograph+".xmp");
        DF_out_clean.write(dir_micrograph+"/"+fn_micrograph+".doc");
        M.write_coordinates(0,-1,dir_micrograph+"/"+fn_micrograph+".pos");

        // Create the micrograph with noise and CTF ........................
        microscope.apply(Md());
        
        // Write the projections
        system(((std::string)"mkdir "+dir_micrograph+"/"+fn_micrograph).c_str());
        FileName fn_proj=fn_micrograph+"/img";
        MetaData DF_out=DF_out_clean;
        DF_out.firstObject();
        STARTINGY(Md())=STARTINGX(Md())=0;
        Image<double> I;
        FileName fn_image;
        for (int l=0; l<Nproj; l++)
        {
            fn_image.compose(fn_proj, l, "xmp");
            Particle_coords particle=M.coord(l);
            I().resize(Xproj,Xproj);
            I().setXmippOrigin();
            FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
                I(i,j)=Md(particle.Y+i,particle.X+j);
            double rot; DF_out.getValue(MDL_ANGLE_ROT,rot); I.setRot(rot);
            double tilt; DF_out.getValue(MDL_ANGLE_TILT,tilt); I.setTilt(tilt);
            double psi; DF_out.getValue(MDL_ANGLE_PSI,psi); I.setPsi(psi);
            I.write(dir_micrograph+"/"+fn_image);
            DF_out.setValue(MDL_IMAGE,I.name());
            DF_out.setValue(MDL_CTF_MODEL,microscope.fn_ctf);
        }
        DF_out.write(dir_micrograph+"/"+fn_micrograph+".doc");
        Md.write(dir_micrograph+"/"+fn_micrograph+".xmp");
    }
}
