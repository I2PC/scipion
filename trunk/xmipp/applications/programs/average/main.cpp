/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
#include <data/geometry.h>

class Statis_parameters: public Prog_parameters
{
public:
    ImageXmipp  sumI, sumI2;
    VolumeXmipp sumV, sumV2;
    int         nI, nV;
    double      sumweight;
    bool        set_weight, weighted_avg, more_options, only_avg, keep_first_header, is_first;
    
public:
    Statis_parameters()
    {
        nI = nV = 0;
	is_first = true;
    }
    void final_process();

    void read(int argc, char **argv)
    {
        more_options = checkParameter(argc, argv, "-more_options");
        Prog_parameters::read(argc, argv);
        set_weight = checkParameter(argc, argv, "-set_weight");
        weighted_avg = checkParameter(argc, argv, "-weighted_avg");
	only_avg = checkParameter(argc, argv, "-only_avg");
	keep_first_header = checkParameter(argc, argv, "-keep_first_header");
        sumweight = 0.;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-more_options]           : show additional options\n";
        if (more_options)
        {
            std::cerr << "  [-set_weight]             : for 2D-images: set weight in header of average to nr. of particles\n";
            std::cerr << "  [-weighted_avg]           : for 2D-images: use header weights in weighted average calculation\n";
	    std::cerr << "  [-only_avg]               : Skip stddev calculation; Output average will be called rootname.xmp\n";
	    std::cerr << "  [-keep_first_header]      : Set header of output images equal to header of first image (only for 2D!) \n";
        }
        std::cerr << std::endl;
        std::cerr << "Purpose: This program allows you to calculate the average and \n"
        << "         standard deviation of a set of images or volumes \n";

    }

};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Statis_parameters *eprm = (Statis_parameters *) prm;
    if (eprm->keep_first_header)
    {
	if (eprm->is_first)
	{
	    eprm->sumI=img;
	    eprm->sumI2=img;
	    eprm->sumI().initZeros();
	    eprm->sumI2().initZeros();
	    eprm->is_first = false;
	}
    }
    else
    {
	eprm->sumI().resize(img());
	eprm->sumI2().resize(img());
    }
    if (eprm->weighted_avg)
    {
        img() *= img.weight();
        eprm->sumweight += img.weight();
    }
    FOR_ALL_ELEMENTS_IN_MATRIX2D(img())
    {
        MAT_ELEM(eprm->sumI(), i, j) += MAT_ELEM(img(), i, j);
    }
    if (!eprm->only_avg)
    {
	FOR_ALL_ELEMENTS_IN_MATRIX2D(img())
	{
	    MAT_ELEM(eprm->sumI2(), i, j) += MAT_ELEM(img(), i, j) *
                MAT_ELEM(img(), i, j);
	}
    }
    eprm->nI++;
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Statis_parameters *eprm = (Statis_parameters *) prm;
    eprm->sumV().resize(vol());
    eprm->sumV2().resize(vol());
    FOR_ALL_ELEMENTS_IN_MATRIX3D(vol())
    {
        VOL_ELEM(eprm->sumV(), k, i, j) += VOL_ELEM(vol(), k ,i, j);
    }
    if (!eprm->only_avg)
    {
	FOR_ALL_ELEMENTS_IN_MATRIX3D(vol())
	{
	    VOL_ELEM(eprm->sumV2(), k, i, j) +=
                VOL_ELEM(vol(), k, i, j) * VOL_ELEM(vol(), k, i, j);
	}
    }
    eprm->nV++;
    return true;
}

void Statis_parameters::final_process()
{
    FileName fnt, fn_root = fn_in.without_extension();
    if (nI != 0)
    {
        FOR_ALL_ELEMENTS_IN_MATRIX2D(sumI())
        {
            MAT_ELEM(sumI(),  i, j) /= nI;
        }
	if (!only_avg)
	{
	    FOR_ALL_ELEMENTS_IN_MATRIX2D(sumI())
	    {
		MAT_ELEM(sumI2(), i, j) /= nI;
		MAT_ELEM(sumI2(), i, j) -= MAT_ELEM(sumI(), i, j) *
                    MAT_ELEM(sumI(), i, j);
		MAT_ELEM(sumI2(), i, j) = sqrt(ABS(MAT_ELEM(sumI2(), i, j)));
	    }
	}
        if (weighted_avg)
        {
            sumI() /= sumweight;
            sumI.set_weight(sumweight);
        }
        else if (set_weight)
        {
            sumI.set_weight((double)nI);
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
    if (nV != 0)
    {
        FOR_ALL_ELEMENTS_IN_MATRIX3D(sumV())
        {
            VOL_ELEM(sumV(),  k, i, j) /= nV;
        }
	if (!only_avg)
	{
	    FOR_ALL_ELEMENTS_IN_MATRIX3D(sumV())
	    {
		VOL_ELEM(sumV2(), k, i, j) /= nV;
		VOL_ELEM(sumV2(), k, i, j) -= VOL_ELEM(sumV(), k, i, j) *
                    VOL_ELEM(sumV(), k, i, j);
		VOL_ELEM(sumV2(), k, i, j) = sqrt(ABS(VOL_ELEM(sumV2(), k, i, j)));
	    }
	}
	if (only_avg)
	{
	    sumV.write(fn_root + ".vol");
	}
	else
	{
	    sumV.write(fn_root + ".med.vol");
	    sumV2.write(fn_root + ".sig.vol");
	}
    }
}

int main(int argc, char **argv)
{
    Statis_parameters prm;
    prm.each_image_produces_an_output = false;
    // Set default action for application of header transformation
    prm.apply_geo = true;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
    prm.final_process();
}
