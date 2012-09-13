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
    Mask mask;
    bool usemask, docc, doeu, domi, doco;

    void defineParams()
    {
        addParamsLine("   --ref <input file>        : Filename for reference image/volume ");
        XmippMetadataProgram::defineParams();
        addParamsLine("  [--mask <input mask=\"\">] : Restrict similarity calculation to region within the mask");
        addParamsLine("  [--co ]                    : Calculate correlation (i.e. signal product).");
        addParamsLine("  [--cc ]                    : Calculate cross-correlation coefficient ");
        addParamsLine("  [--eu ]                    : Calculate euclidian distance ");
        addParamsLine("  [--mi ]                    : Calculate mutual information");
        mask.defineParams(this,INT_MASK,NULL,"Measure only within an area");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        fn_ref = getParam("--ref");
        ref.read(fn_ref);
        ref().setXmippOrigin();
        doco = checkParam("--co");
        docc = checkParam("--cc");
        doeu = checkParam("--eu");
        domi = checkParam("--mi");
        if (usemask = checkParam("--mask"))
            mask.readParams(this);
    }

    void show()
    {
    	if (verbose==0)
    		return;
        std::cout << "Reference file: " << fn_ref << std::endl;
        XmippMetadataProgram::show();
    }

    void preProcess() {
    	Image<double> Iref;
    	Iref.read(fn_ref);
    	mask.generate_mask(Iref());
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg,mdIn,objId);

        double co, cc, eu, mi;
        if (!usemask)
        {
            if (doco)
                co = correlation(ref(), img());
            if (docc)
                cc = correlationIndex(ref(), img());
            if (doeu)
                eu = euclidianDistance(ref(), img());
            if (domi)
                mi = mutualInformation(ref(), img());
        }
        else
        {
            if (doco)
                co = correlation(ref(), img(), &mask.get_binary_mask());
            if (docc)
                cc = correlationIndex(ref(), img(), &mask.get_binary_mask());
            if (doeu)
                eu = euclidianDistance(ref(), img(), &mask.get_binary_mask());
            if (domi)
                mi = mutualInformation(ref(), img(), 0, 0, &mask.get_binary_mask());
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
    ProgCorrelation program;
    program.read(argc, argv);
    return program.tryRun();
}
