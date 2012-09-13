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

#include <data/progs.h>
#include <data/args.h>
#include <data/geometry.h>

class ProgAverage: public XmippMetadataProgram
{
protected:
    Image<double>  sumI, sumI2;
    int         nI, nV;
    double      sumweight;
    bool        set_weight, weighted_avg, more_options, only_avg, keep_first_header, is_first;

    void defineParams()
    {
        XmippMetadataProgram::defineParams();
        addUsageLine("Calculate the average and standard deviation of a set of images or volumes.");
        addParamsLine( "  [-set_weight+]             : for 2D-images: set weight in header of average to nr. of particles");
        addParamsLine( "  [-weighted_avg+]           : for 2D-images: use header weights in weighted average calculation");
        addParamsLine( "  [-only_avg+]               : Skip stddev calculation; Output average will be called rootname.xmp");
        addParamsLine( "  [-keep_first_header+]      : Set header of output images equal to header of first image (only for 2D!) ");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        set_weight = checkParam("-set_weight");
        weighted_avg = checkParam("-weighted_avg");
        only_avg = checkParam("-only_avg");
        keep_first_header = checkParam("-keep_first_header");

        ///Some initializations
        sumweight = 0.;
        nI = nV = 0;
        is_first = true;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg,mdIn,objId);

        if (keep_first_header)
        {
            if (is_first)
            {
                sumI=img;
                sumI2=img;
                sumI().initZeros();
                sumI2().initZeros();
                is_first = false;
            }
        }
        else
        {
            sumI().resize(img());
            sumI2().resize(img());
        }
        if (weighted_avg)
        {
            img() *= img.weight();
            sumweight += img.weight();
        }
        FOR_ALL_ELEMENTS_IN_ARRAY3D(img())
        {
            A3D_ELEM(sumI(), k, i, j) += A3D_ELEM(img(), k, i, j);
        }
        if (!only_avg)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(img())
            {
                A3D_ELEM(sumI2(), k, i, j) += A3D_ELEM(img(), k, i, j) *
                                              A3D_ELEM(img(), k, i, j);
            }
        }
        ++nI;
    }

    void postProcess()
    {
        FileName fnt, fn_root = fn_in.withoutExtension();
        if (nI != 0)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(sumI())
            {
                A3D_ELEM(sumI(), k, i, j) /= nI;
            }
            if (!only_avg)
            {
                FOR_ALL_ELEMENTS_IN_ARRAY3D(sumI())
                {
                    A3D_ELEM(sumI2(), k, i, j) /= nI;
                    A3D_ELEM(sumI2(), k, i, j) -= A3D_ELEM(sumI(), k, i, j) *
                                                  A3D_ELEM(sumI(), k, i, j);
                    A3D_ELEM(sumI2(), k, i, j) = sqrt(ABS(A3D_ELEM(sumI2(), k, i, j)));
                }
            }
            if (weighted_avg)
            {
                sumI() /= sumweight;
                sumI.setWeight(sumweight);
            }
            else if (set_weight)
            {
                sumI.setWeight((double)nI);
                std::cerr << " Setting weight in the header of the average image to " << integerToString(nI) << std::endl;
            }
            if (only_avg)
            {
                sumI.write(fn_root + ".xmp");
            }
            else
            {
                sumI.write(fn_root + ".med.xmp");
                sumI2.write(fn_root + ".sig.xmp");
            }
        }
    }

}
;///end of class ProgAverage

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgAverage program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe << std::endl;
    }
    return 0;
}

