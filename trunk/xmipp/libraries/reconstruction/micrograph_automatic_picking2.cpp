/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
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
#include <math.h>
#include "micrograph_automatic_picking2.h"
#include <data/filters.h>
#include <data/rotational_spectrum.h>
#include <reconstruction/denoise.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_filename.h>
#include <algorithm>


MultidimArray< std::complex<double> > arrayFourier(MultidimArray<double> inputArray)
{
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > fourierTransform;
    transformer.completeFourierTransform(inputArray,fourierTransform);
    return fourierTransform;
}

//Generate maps from the fourier transform of an array
void filterBankGenerator(MultidimArray<double> &inputMicrograph, const FileName &fnFilterBankStack)
{
	fnFilterBankStack.deleteFile();

	Image<double> Iaux;
	Iaux()=inputMicrograph;

	FourierFilter filter;
    filter.raised_w = 0.02;
    filter.FilterShape = RAISED_COSINE;
    filter.FilterBand = BANDPASS;

    MultidimArray< std::complex<double> > micrographFourier;
    FourierTransformer transformer;
    transformer.FourierTransform(Iaux(), micrographFourier, true);

    for (int i=0; i<7; i++)
    {
        filter.w1 =0.025*i;
        filter.w2 = (filter.w1)+0.025;

        transformer.setFourier(micrographFourier);
        filter.applyMaskFourierSpace(inputMicrograph,transformer.fFourier);
        transformer.inverseFourierTransform();

        Iaux.write(fnFilterBankStack,i+1,true,WRITE_APPEND);
    }
}

AutoParticlePicking2::AutoParticlePicking2(const FileName &fn, Micrograph *_m, bool _fast)
{
    __m = _m;
    fn_micrograph = fn;
    __numThreads = 1;
    piece_xsize = 512;
    particle_radius = 0;
    __piece_overlap = 0;
    __scan_overlap = 0;
    __learn_overlap = 0;
    __fast=_fast;
    __incore=false;
}

std::ostream & operator <<(std::ostream &_out, const Particle2 &_p)
{
    _out << _p.micrograph << " " << _p.x << " " << _p.y << " " << _p.cost << " ";
    FOR_ALL_ELEMENTS_IN_MATRIX1D(_p.vec)
    _out << VEC_ELEM(_p.vec,i) << " ";
    _out << std::endl;
    return _out;
}

void Particle2::read(std::istream &_in, int _vec_size)
{
    _in >> micrograph >> x >> y >> cost;
    if (cost>=0)
        status=1;
    else
        status=0;
    vec.resize(_vec_size);
    _in >> vec;
}

AutoParticlePicking2::~AutoParticlePicking2()
{}

void AutoParticlePicking2::learnParticles()
{
    int num_part = __m->ParticleNo();
    for (int i=0;i<num_part;i++)
    {
        int x = (__m->coord(i).X)/scaleRate;
        int y = (__m->coord(i).Y)/scaleRate;

    }

}

void  AutoParticlePicking2::buildPositiveVectors()
{}

void AutoParticlePicking2::saveAutoFeatureVectors(const FileName &fn, int Nvectors) const
{
    std::ofstream fh_auto;
    fh_auto.open(fn.c_str());
    if (!fh_auto)
        REPORT_ERROR(ERR_IO_NOWRITE,fn);
    fh_auto << Nvectors << " ";
    if (Nvectors == 0)
        fh_auto << "0\n";
    else
        fh_auto << VEC_XSIZE(__auto_candidates[0].vec) << std::endl;
    size_t Nauto = __auto_candidates.size();
    for (size_t i = 0; i < Nauto; i++)
    {
        const Particle2 &p = __auto_candidates[i];
        if (p.cost>0)
            fh_auto << p;
    }
    fh_auto.close();
}

void AutoParticlePicking2::loadAutoFeatureVectors(const FileName &fn)
{
    std::ifstream fh_auto;
    fh_auto.open(fn.c_str());
    if (!fh_auto)
        REPORT_ERROR(ERR_IO_NOTOPEN,fn);
    size_t Nauto;
    int vec_size;
    fh_auto >> Nauto >> vec_size;
    __auto_candidates.reserve(Nauto);
    Particle2 p;
    for (size_t i = 0; i < Nauto; i++)
    {
        p.read(fh_auto, vec_size);
        __auto_candidates.push_back(p);
    }
    fh_auto.close();
}

void AutoParticlePicking2::loadModels(const FileName &fn_root)
{
    FileName fn_training = fn_root + "_training.txt";
    if (!fn_training.exists())
        return;

    // Load training vectors
    String dummy;
    std::ifstream fh_training;
    fh_training.open(fn_training.c_str());
    if (!fh_training)
        REPORT_ERROR(ERR_IO_NOTOPEN, fn_training);
    fh_training
    >> dummy >> piece_xsize
    >> dummy >> particle_radius;
}

// ==========================================================================
// Section: Program interface ===============================================
// ==========================================================================
void ProgMicrographAutomaticPicking2::readParams()
{
    fn_micrograph = getParam("-i");
    fn_model = getParam("--model");
    size = getIntParam("--particleSize");
    mode = getParam("--mode");
    if (mode == "train")
        fn_train = getParam("--mode", 1);
    Nthreads = getIntParam("--thr");
    fn_root = getParam("--outputRoot");
    fast = checkParam("--fast");
    incore = checkParam("--in_core");
}

void ProgMicrographAutomaticPicking2::defineParams()
{
    addUsageLine("Automatic particle picking for micrographs");
    addUsageLine("+The algorithm is designed to learn the particles from the user, as well as from its own errors.");
    addUsageLine("+The algorithm is fully described in [[http://www.ncbi.nlm.nih.gov/pubmed/19555764][this paper]].");
    addParamsLine("  -i <micrograph>               : Micrograph image");
    addParamsLine("  --outputRoot <rootname>       : Output rootname");
    addParamsLine("  --mode <mode>                 : Operation mode");
    addParamsLine("         where <mode>");
    addParamsLine("                    try              : Try to autoselect within the training phase.");
    addParamsLine("                    train <posfile=\"\">  : posfile contains the coordinates of manually picked particles");
    addParamsLine("                                     : <rootname>_auto_feature_vectors.txt contains the particle structure created by this program when used in automatic selection mode");
    addParamsLine("                                     : <rootname>_false_positives.xmd contains the list of false positives among the automatically picked particles");
    addParamsLine("                    autoselect  : Autoselect");
    addParamsLine("  --model <model_rootname>      : Bayesian model of the particles to pick");
    addParamsLine("  --particleSize <size>         : Particle size in pixels");
    addParamsLine("  [--thr <p=1>]                 : Number of threads for automatic picking");
    addParamsLine("  [--fast]                      : Perform a fast preprocessing of the micrograph (Fourier filter instead of Wavelet filter)");
    addParamsLine("  [--in_core]                   : Read the micrograph in memory");
    addExampleLine("Automatically select particles during training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode try ");
    addExampleLine("Training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode train manual.pos");
    addExampleLine("Automatically select particles after training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode autoselect");
}

void ProgMicrographAutomaticPicking2::run()
{
    Micrograph m;
    m.open_micrograph(fn_micrograph);
    Image<double> I2;
    AutoParticlePicking2 *autoPicking = new AutoParticlePicking2(fn_micrograph,&m,fast);
    autoPicking->microImage.read(fn_micrograph);
    autoPicking->scaleRate=std::max(0.25,50.0/size);
    selfScaleToSizeFourier((m.Ydim)*autoPicking->scaleRate,(m.Xdim)*autoPicking->scaleRate,autoPicking->microImage(), 2);

    // Generating the filter bank
    FileName fnFilterBank=fn_micrograph.removeLastExtension()+"_filterbank.stk";
    filterBankGenerator(autoPicking->microImage(), fnFilterBank);

    autoPicking->loadModels(fn_model);
    FileName familyName=fn_model.removeDirectories();
    FileName fnAutoParticles=familyName+"@"+fn_root+"_auto.pos";
    FileName fnVectors=fn_root + "_auto_feature_vectors_"+familyName+".txt";
    //
    MetaData MD;
    // Insert all true positives
    if (fn_train!="")
    {
        MD.read(fn_train);
        int x, y;
        FOR_ALL_OBJECTS_IN_METADATA(MD)
        {
            MD.getValue(MDL_XINT, x, __iter.objId);
            MD.getValue(MDL_YINT, y, __iter.objId);
            m.add_coord(x, y, 0, 1);
        }
    }
    //autoPicking->learnParticles();

}


