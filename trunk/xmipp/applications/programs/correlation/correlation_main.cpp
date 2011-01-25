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

class ProgCorrelation: public XmippMetadataProgram
{
public:
    FileName    fn_ref, fn_msk;
    Image<double> ref, M;
    MultidimArray<int> mask;
    bool usemask, docc, doeu, domi, doco;

    void defineParams()
    {
        addParamsLine("   -ref <input file>        : Filename for reference image/volume ");
        XmippMetadataProgram::defineParams();
        addParamsLine("  [-mask <input mask=\"\">] : Restrict similarity calculation to region within the mask");
        addParamsLine("  [-co ]                    : Calculate correlation (i.e. signal product).");
        addParamsLine("  [-cc ]                    : Calculate cross-correlation coefficient ");
        addParamsLine("  [-eu ]                    : Calculate euclidian distance ");
        addParamsLine("  [-mi ]                    : Calculate mutual information");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        usemask = false;
        fn_ref = getParam("-ref");
        ref.read(fn_ref);
        ref().setXmippOrigin();
        fn_msk = getParam("-mask");
        if (fn_msk != "")
        {
            usemask = true;
            M.read(fn_msk);
            M().setXmippOrigin();
            typeCast(M(),mask);
        }
        doco = checkParam("-co");
        docc = checkParam("-cc");
        doeu = checkParam("-eu");
        domi = checkParam("-mi");
    }

    void show()
    {
        std::cout << "Reference file: " << fn_ref << std::endl;
        if (usemask)
            std::cout << "mask file: " << fn_msk << std::endl;
        XmippMetadataProgram::show();
        std::cout << std::endl;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, long int objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg,mdIn,objId);

        double co, cc, eu, mi;
        if (!usemask)
        {
            if (doco)
                co = correlation(ref(), img());
            if (docc)
                cc = correlation_index(ref(), img());
            if (doeu)
                eu = euclidian_distance(ref(), img());
            if (domi)
                mi = mutual_information(ref(), img());
        }
        else
        {
            if (doco)
                co = correlation(ref(), img(), &mask);
            if (docc)
                cc = correlation_index(ref(), img(), &mask);
            if (doeu)
                eu = euclidian_distance(ref(), img(), &mask);
            if (domi)
                mi = mutual_information(ref(), img(), 0, 0, &mask);
        }

        std::cout << img.name() << ": ";
        if (doco)
            std::cout << " co = " << co;
        if (docc)
            std::cout << " cc = " << cc;
        if (doeu)
            std::cout << " eu = " << eu;
        if (domi)
            std::cout << " mi = " << mi;
        std::cout << std::endl;
    }
}
;///end of class ProgCorrelation


/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  try
  {
      ProgCorrelation program;
      program.read(argc, argv);
      program.run();
  }
  catch (XmippError xe)
  {
      std::cerr << xe;
  }
}
