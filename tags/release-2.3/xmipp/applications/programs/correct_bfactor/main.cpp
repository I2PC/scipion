/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/progs.h>
#include <data/args.h>
#include <reconstruction/correct_bfactor.h>

class Correct_Bfactor_parameters: public Prog_parameters
{
public:
    Prog_correct_bfactor_prm prmb;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        if (checkParameter(argc, argv, "-ref"))
        {
            if (checkParameter(argc, argv, "-allpoints"))
                prmb.mode = ALLPOINTS_REF;
            else
                prmb.mode = BFACTOR_REF;
            prmb.fn_ref= getParameter(argc, argv, "-ref");
        }
        else if (checkParameter(argc, argv, "-adhoc"))
        {
            prmb.mode = BFACTOR_ADHOC;
            prmb.adhocB = textToFloat(getParameter(argc, argv, "-adhoc"));
        }
        else if (checkParameter(argc, argv, "-auto"))
        {
            prmb.mode = BFACTOR_AUTO;
        }
        else
        {
            REPORT_ERROR(1,"Please provide -auto, -adhoc <B> or -ref <fn_ref> option");
        }
        prmb.sampling_rate = textToFloat(getParameter(argc, argv, "-sampling"));
        prmb.apply_maxres = textToFloat(getParameter(argc, argv, "-maxres"));
        prmb.fit_minres = textToFloat(getParameter(argc, argv, "-fit_minres","15"));
        prmb.fit_maxres = textToFloat(getParameter(argc, argv, "-fit_maxres","-1"));

        if (prmb.fit_maxres < 0.)
            prmb.fit_maxres = prmb.apply_maxres;
        prmb.fn_fsc = getParameter(argc, argv, "-fsc","");

    }

    void show()
    {
        Prog_parameters::show();
        std::cout << "Pixel size : "<<prmb.sampling_rate<<" Angstrom"<<std::endl;
        std::cout << "Maximum resolution: "<<prmb.apply_maxres<<" Angstrom"<<std::endl;
        if (prmb.mode == BFACTOR_REF || prmb.mode == BFACTOR_AUTO)
        {
            std::cerr<<"Fit within resolutions: "<<prmb.fit_minres<<" - "<<prmb.fit_maxres<<" Angstrom"<<std::endl;
        }
        if (prmb.mode == BFACTOR_REF)
        {
            std::cout << "Adjust B-factor according to reference "<<prmb.fn_ref<<std::endl;
        } 
        else if (prmb.mode == BFACTOR_ADHOC)
        {
            std::cout << "Apply ad-hoc B-factor of "<<prmb.adhocB<<" squared Angstroms"<<std::endl;
        }
        else
        {
            std::cout << "Use automated B-factor fit (Rosenthal and Henderson, 2003) "<<std::endl;
        }
        if (prmb.fn_fsc != "")
            std::cout << "Use signal-to-noise weighted based on "<<prmb.fn_fsc<<std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr 
            << " Use one of the three following modes: \n"
            << "  [-auto]                   : Use automated B-factor fit in flat Wilson region\n"
            << "  [-ref <fn_ref>]           : Fit B-factor according to the reference \n"
            << "  [-adhoc <B>]              : Use a user-provided (negative) B-factor\n"
            << " Specific parameters: \n"
            << "   -sampling <float>        : Pixel size (in Ang) \n"
            << "   -maxres <float>          : High-resolution limit for B-factor correction \n"
            << "  [-fit_minres <f=15>]      : Low-resolution  limit (in Ang) for fit in -auto or -ref \n"
            << "  [-fit_maxres <f=maxres>]  : High-resolution limit (in Ang) for fit in -auto or -ref \n"
            << "  [-allpoints]              : Do not fit B-factor, adjust power spectrum to reference \n"
            << "\n"
            << " Note: do not use the automated mode for maps with resolutions lower than 12-15 Angstroms!\n"
            ;
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Correct_Bfactor_parameters *eprm = (Correct_Bfactor_parameters *) prm;
    REPORT_ERROR(1,"B-factor correction for 2D images not implemented yet");

    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Correct_Bfactor_parameters *eprm = (Correct_Bfactor_parameters *) prm;
    FileName fn_guinier;
    fn_guinier = prm->fn_out + ".guinier";
    eprm->prmb.bfactor_correction(vol(), fn_guinier);
    
    return true;
}

int main(int argc, char **argv)
{
    Correct_Bfactor_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

