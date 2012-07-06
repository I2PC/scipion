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
#include <data/metadata.h>
#include <data/geometry.h>

class Shift_Scale_parameters: public Prog_parameters
{
public:
    Matrix1D<double> shift;
    Matrix1D<double> scale;
    bool             wrap;
    bool             store_in_header;
    bool             center_mass;
    MetaData         DF_shifts;
    MDIterator      shiftsIter;
    MetaData         DF_scales;
    MDIterator      scalesIter;
    bool             Docfile;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        FileName fnShift = getParameter(argc, argv, "-shift","");
        FileName fnScale = getParameter(argc, argv, "-scale","");
        center_mass = checkParameter(argc, argv, "-center_mass");

        if (fnShift=="" && fnScale=="" && !center_mass)
            REPORT_ERROR(ERR_ARG_MISSING, "Shift_Scale:: Cannot find -shift or -scale");
        Docfile = (fnShift!="" && fnShift.isMetaData(false)) ||
                  (fnScale!="" && fnScale.isMetaData(false));
        if (Docfile)
        {
            if (fnShift!="" && fnShift.isMetaData(false))
            {
                DF_shifts.read(fnShift);
                shiftsIter = MDIterator(DF_shifts);
            }
            if (fnScale!="" && fnScale.isMetaData(false))
            {
                DF_scales.read(fnScale);
                scalesIter = MDIterator(DF_scales);
            }
        }
        else
        {
            int my_dim=0;
            if (checkParameter(argc, argv, "-shift"))
            {
                shift = getVectorParameter(argc, argv, "-shift", -1);
                my_dim = VEC_XSIZE(shift);
            }
            if (checkParameter(argc, argv, "-scale"))
            {
                scale = getVectorParameter(argc, argv, "-scale", -1);
                my_dim = VEC_XSIZE(scale);
            }

            if (!checkParameter(argc, argv, "-shift") && !center_mass)
                shift.initZeros(my_dim);
            if (!checkParameter(argc, argv, "-scale") && !center_mass)
            {
                scale.resize(my_dim);
                scale.initConstant(1);
            }
        }
        wrap = !checkParameter(argc, argv, "-dont_wrap");
        store_in_header = checkParameter(argc, argv, "-store_in_header");
    }

    void show()
    {
        Prog_parameters::show();
        if (wrap)
            std::cout << "Wrapping image/volume\n";
        else
            std::cout << "Not wrapping image/volume\n";
        if (store_in_header)
            std::cout << "Storing the shift in header\n";
        else
            std::cout << "Shifting image/volume\n";
        if (VEC_XSIZE(shift) > 1)
            std::cout << "Shift: " << shift.transpose() << std::endl;
        else if (DF_shifts.getFilename() != "")
            std::cout << "Shift docfile: " << DF_shifts.getFilename() << std::endl;
        if (VEC_XSIZE(scale) > 1)
            std::cout << "Scale: " << scale.transpose() << std::endl;
        else if (DF_scales.getFilename() != "")
            std::cout << "Scale: docfile: "       << DF_scales.getFilename() << std::endl;
        if (center_mass)
            std::cout << "Moving center of mass to origin\n";
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr
        << "   -shift \"[<x>,<y>[,<z>]]\" : Shift by (x,y,z) for volumes, (x,y) for images\n"
        << "   -scale \"[<x>,[<y>,<z>]]\" : Scale by (x,y,z)\n"
        << "   -shift <DocFile>         : Shifts are stored in a Docfile\n"
        << "   -scale <DocFile>         : Scales are stored in a Docfile (may be the same\n"
        << "                              Docfile used for shifts\n"
        << "  [-center_mass]            : Move the center of mass to the origin\n"
        << "  [-dont_wrap]              : By default, the image is wrapped\n"
        << "  [-store_in_header]        : Do not shift, but store the shift in the header\n"
        << "                              Not applicable to volumes\n";
    }
};

bool process_img(Image<double> &img, const Prog_parameters *prm)
{
    size_t id;
    Shift_Scale_parameters *eprm = (Shift_Scale_parameters *) prm;
    Matrix2D<double> A;
    int dim;
    if (ZSIZE(img())==1)
    {
        dim=2;
        A.resize(3, 3);
    }
    else
    {
        dim=3;
        A.resize(4, 4);
    }
    A.initIdentity();

    if (eprm->DF_shifts.getFilename() != "")
    {
        eprm->shift.resize(dim);
        id = eprm->shiftsIter.objId;
        eprm->DF_shifts.getValue(MDL_SHITF_X,XX(eprm->shift), id);
        eprm->DF_shifts.getValue(MDL_SHITF_Y,YY(eprm->shift), id);
        if (dim==3)
            eprm->DF_shifts.getValue(MDL_SHITF_Z,ZZ(eprm->shift), id);
        eprm->shiftsIter.moveNext();
    }
    else if (eprm->Docfile)
    {
        eprm->shift.resize(2);
        eprm->shift.initConstant(0.);
    }
    else if (eprm->center_mass)
    {
        img().centerOfMass(eprm->shift);
        eprm->shift = -eprm->shift;
    }

    if (dim==2)
    {
        A(0, 2) = XX(eprm->shift);
        A(1, 2) = YY(eprm->shift);
    }
    else
    {
        A(0, 3) = XX(eprm->shift);
        A(1, 3) = YY(eprm->shift);
        A(2, 3) = ZZ(eprm->shift);
    }
    if (eprm->DF_scales.getFilename() != "")
    {
        eprm->scale.resize(dim);
        id = eprm->scalesIter.objId;
        eprm->DF_scales.getValue(MDL_SCALE,XX(eprm->scale), id);
        eprm->DF_scales.getValue(MDL_SCALE,YY(eprm->scale), id);
        if (dim==3)
            eprm->DF_scales.getValue(MDL_SCALE,XX(eprm->scale), id);
        eprm->scalesIter.moveNext();
    }
    else if (eprm->Docfile || eprm->center_mass)
    {
        eprm->scale.resize(dim);
        eprm->scale.initConstant(1.);
    }
    A(0, 0) = XX(eprm->scale);
    A(1, 1) = YY(eprm->scale);
    if (dim==3)
        A(2,2) = ZZ(eprm->scale);
    if (!eprm->store_in_header)
        selfApplyGeometry(BSPLINE3,img(),A, IS_NOT_INV, eprm->wrap);
    else
    {
        if (dim==2)
            img.setShifts(XX(eprm->shift), YY(eprm->shift));
        else
            img.setShifts(XX(eprm->shift), YY(eprm->shift), ZZ(eprm->shift));
    }
    return true;
}

int main(int argc, char **argv)
{
    Shift_Scale_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img);
}
