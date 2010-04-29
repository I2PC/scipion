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

#include "ctf_correct_idr.h"
#include <data/projection.h>

Prog_IDR_ART_Parameters::Prog_IDR_ART_Parameters()
{
    fn_vol=fn_ctfdat=fnRoot="";
    mu=1.8;
    MPIversion=false;
    numberOfProcessors=1;
    MPIrank=0;
}

void Prog_IDR_ART_Parameters::read(int argc, char **argv)
{
    fn_vol = getParameter(argc, argv, "-vol");
    fn_ctfdat = getParameter(argc, argv, "-ctfdat");
    fnRoot = getParameter(argc, argv, "-oroot");
    mu = textToFloat(getParameter(argc, argv, "-mu", "1.8"));
}

void Prog_IDR_ART_Parameters::produce_side_info()
{
    V.read(fn_vol);
    V().setXmippOrigin();
    ctfdat.read(fn_ctfdat);
}

void Prog_IDR_ART_Parameters::show()
{
    std::cout << "Input volume: " << fn_vol << std::endl
              << "CTFDat: " << fn_ctfdat << std::endl
              << "Oroot: " << fnRoot << std::endl
              << "Relaxation factor: " << mu << std::endl
              ;
}

void Prog_IDR_ART_Parameters::Usage()
{
    std::cerr << "Usage: IDR\n"
              << "   -vol <volume>        : Voxel volume with the current reconstruction\n"
              << "   -ctfdat <ctfdat>     : List of projections and CTFs\n"
              << "   -oroot <root>        : Rootname of the output files\n"
              << "  [-mu <mu=1.8>]        : Relaxation factor\n"
              ;
}

/* IDR correction ---------------------------------------------------------- */
//#define DEBUG
void Prog_IDR_ART_Parameters::IDR_correction()
{
    Projection Ireal, Inorm, Itheo, Itheo_CTF;

    if (MPIrank==0)
    {
        std::cerr << "Modifying input data ...\n";
        init_progress_bar(ctfdat.lineNo());
    }
    int istep = CEIL((double)ctfdat.lineNo() / 60.0);
    int imgs = 0;
    ctfdat.goFirstLine();
    while (!ctfdat.eof())
    {
        FileName fn_img, fn_ctf;
        ctfdat.getCurrentLine(fn_img,fn_ctf);
        if (fn_img!="" && (imgs%numberOfProcessors==MPIrank))
        {
            // Read current input image
            Ireal.read(fn_img);
            Ireal().selfTranslateBSpline(3,vectorR2(Ireal.Xoff(),Ireal.Yoff()));

            // Project the volume in the same direction
            project_Volume(V(), Itheo, YSIZE(Ireal()), XSIZE(Ireal()),
                           Ireal.rot(), Ireal.tilt(), Ireal.psi());

            // Copy to theo_CTF and resize
            Itheo_CTF() = Itheo();
            int Ydim = YSIZE(Itheo());
            int Xdim = XSIZE(Itheo());
            Itheo_CTF().setXmippOrigin();
            Itheo_CTF().window(FIRST_XMIPP_INDEX(2*Ydim), FIRST_XMIPP_INDEX(2*Xdim),
                               LAST_XMIPP_INDEX(2*Ydim), LAST_XMIPP_INDEX(2*Xdim));

            // Read CTF file
            FourierMask ctf;
            ctf.FilterBand = CTF;
            ctf.ctf.read(fn_ctf);
            ctf.ctf.enable_CTFnoise = false;
            ctf.ctf.Produce_Side_Info();
            ctf.generate_mask(Itheo_CTF());
            ctf.correct_phase();

            // Apply CTF
            ctf.apply_mask_Space(Itheo_CTF());
            Itheo_CTF().window(FIRST_XMIPP_INDEX(Ydim), FIRST_XMIPP_INDEX(Xdim),
                               LAST_XMIPP_INDEX(Ydim), LAST_XMIPP_INDEX(Xdim));

            // Center the all images
            Ireal().setXmippOrigin();
            Itheo().setXmippOrigin();
#ifdef DEBUG
            Itheo.write("PPPtheo.xmp");
            Itheo_CTF.write("PPPtheo_CTF.xmp");
            Ireal.write("PPPreal.xmp");
#endif

            // Apply IDR process
            FOR_ALL_ELEMENTS_IN_MATRIX2D(Ireal())
            IMGPIXEL(Itheo, i, j) = mu * IMGPIXEL(Ireal, i, j) +
                                    (IMGPIXEL(Itheo, i, j) - mu * IMGPIXEL(Itheo_CTF, i, j));

            // Save output image
            Itheo.write(fnRoot+integerToString(imgs,6)+".xmp");

#ifdef DEBUG
            ImageXmipp save;
            save() = Itheo() - mu * Itheo_CTF();
            save.write("PPPdiff.xmp");
            Itheo.write("PPPidr.xmp");
            std::cout << "Press any key to continue\n";
            char c;
            std::cin >> c;
#endif
        }

        if (imgs++ % istep == 0 && MPIrank==0) progress_bar(imgs);
        ctfdat.nextLine();
    }
    if (MPIrank==0) progress_bar(ctfdat.lineNo());
}
