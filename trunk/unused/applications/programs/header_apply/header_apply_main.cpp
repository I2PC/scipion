/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

class ProgHeaderApply: public XmippMetadataProgram
{
protected:
    bool wrap;

    void defineParams()
    {
        addUsageLine("Apply the geometric transformation stored in the image header.");
        XmippMetadataProgram::defineParams();
        addParamsLine(" [--dont_wrap]              : By default, the image is wrapped");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        wrap = !checkParam("--dont_wrap");
    }

    void show()
    {
        XmippMetadataProgram::show();

        if (!wrap)
            std::cout << "Do not wrap"<<std::endl;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg,mdIn,objId);

        if (ZSIZE(img())!=1 || NSIZE(img())!=1)
            REPORT_ERROR(ERR_MULTIDIM_DIM, "This program is intended only for images");

        MultidimArray<double> Maux;
        Matrix2D<double> A;
        img.getTransformationMatrix(A);
        applyGeometry(BSPLINE3, Maux, img(), A, IS_INV, wrap);
        img()=Maux;
        //Reset in-plane transformations of the header
        img.setShifts(0,0);
        img.setPsi(0.);

        if (img.tilt() == 0)
            img.setRot(0.);
        img.setFlip(0.);

        img.write(fnImg);
    }
}
;///end of class ProgHeaderApply

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    try
    {
        ProgHeaderApply program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}
