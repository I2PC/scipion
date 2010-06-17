/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>

class Similarity_parameters: public Prog_parameters
{
public:
    FileName    fn_ref, fn_msk;
    Image<double> ref, M;
    MultidimArray<int> mask;
    bool usemask, docc, doeu, domi, doco;

public:
    Similarity_parameters()
    {}
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        try
        {
            usemask = false;
            fn_ref = getParameter(argc, argv, "-ref");
            ref.read(fn_ref,true,-1,apply_geo);
            ref().setXmippOrigin();
            fn_msk = getParameter(argc, argv, "-mask", "");
            if (fn_msk != "")
            {
                usemask = true;
                M.read(fn_msk, true, -1, apply_geo);
                M().setXmippOrigin();
                typeCast(M(),mask);
            }
            doco=checkParameter(argc, argv, "-co");
            docc=checkParameter(argc, argv, "-cc");
            doeu=checkParameter(argc, argv, "-eu");
            domi=checkParameter(argc, argv, "-mi");
        }
        catch (Xmipp_error XE)
        {
            std::cout << XE;
            usage();
            exit(1);
        }
    }

    void show()
    {
        std::cout << "Reference file: " << fn_ref << std::endl;
        if (usemask) std::cout << "mask file: " << fn_msk << std::endl;
        Prog_parameters::show();
        std::cout << std::endl;
    }
    void usage()
    {
        std::cerr << "   -ref <input file>        : Filename for reference image/volume \n";
        Prog_parameters::usage();
        std::cerr << "  [-mask <input mask>]      : Restrict similarity calculation to region within the mask\n";
        std::cerr << "  [-co ]                    : Calculate correlation (i.e. signal product).\n";
        std::cerr << "  [-cc ]                    : Calculate cross-correlation coefficient \n";
        std::cerr << "  [-eu ]                    : Calculate euclidian distance \n";
        std::cerr << "  [-mi ]                    : Calculate mutual information\n";
    }
};


bool process_img(Image<double> &img, const Prog_parameters *prm)
{
    Similarity_parameters *eprm = (Similarity_parameters *) prm;

    double co, cc, eu, mi;
    if (!eprm->usemask)
    {
        if (eprm->doco) co = correlation(eprm->ref(), img());
        if (eprm->docc) cc = correlation_index(eprm->ref(), img());
        if (eprm->doeu) eu = euclidian_distance(eprm->ref(), img());
        if (eprm->domi) mi = mutual_information(eprm->ref(), img());
    }
    else
    {
        if (eprm->doco) co = correlation(eprm->ref(), img(), &eprm->mask);
        if (eprm->docc) cc = correlation_index(eprm->ref(), img(), &eprm->mask);
        if (eprm->doeu) eu = euclidian_distance(eprm->ref(), img(), &eprm->mask);
        if (eprm->domi) mi = mutual_information(eprm->ref(), img(), 0, 0, &eprm->mask);
    }

    std::cout << img.name() << ": ";
    if (eprm->doco) std::cout << " co= " << co;
    if (eprm->docc) std::cout << " cc= " << cc;
    if (eprm->doeu) std::cout << " eu= " << eu;
    if (eprm->domi) std::cout << " mi= " << mi;
    std::cout << std::endl;
    return true;
}

int main(int argc, char **argv)
{
    Similarity_parameters prm;
    prm.allow_time_bar = false;
    prm.each_image_produces_an_output = false;
    prm.apply_geo = true;
    SF_main(argc, argv, &prm, (void*)&process_img);
}
