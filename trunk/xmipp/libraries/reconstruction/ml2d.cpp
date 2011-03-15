/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2007)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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
#include "ml2d.h"

ML2DBaseProgram::ML2DBaseProgram()
{
    do_ML3D = false;
    refs_per_class = 1;
    defaultNiter = 100;
    defaultRoot = "ml2d";
    referenceExclusive = allowFastOption = true;
    allowThreads = allowRestart = allowIEM = true;
    defaultNiter = 100;
}

void ML2DBaseProgram::initSamplingStuff()
{
    // Set sampling stuff: flipping matrices, psi_step etc.
    Matrix2D<double> A(3, 3);
    psi_max = 90.;
    nr_psi = CEIL(psi_max / psi_step);
    psi_step = psi_max / nr_psi;
    nr_flip = nr_nomirror_flips = 4;
    // 0, 90, 180 & 270 degree flipping, as well as mirror
    A.initIdentity();
    F.push_back(A);
    A(0, 0) = 0.;
    A(1, 1) = 0.;
    A(1, 0) = 1.;
    A(0, 1) = -1;
    F.push_back(A);
    A(0, 0) = -1.;
    A(1, 1) = -1.;
    A(1, 0) = 0.;
    A(0, 1) = 0;
    F.push_back(A);
    A(0, 0) = 0.;
    A(1, 1) = 0.;
    A(1, 0) = -1.;
    A(0, 1) = 1;
    F.push_back(A);
    if (do_mirror)
    {
        nr_flip = 8;
        A.initIdentity();
        A(0, 0) = -1;
        F.push_back(A);
        A(0, 0) = 0.;
        A(1, 1) = 0.;
        A(1, 0) = 1.;
        A(0, 1) = 1;
        F.push_back(A);
        A(0, 0) = 1.;
        A(1, 1) = -1.;
        A(1, 0) = 0.;
        A(0, 1) = 0;
        F.push_back(A);
        A(0, 0) = 0.;
        A(1, 1) = 0.;
        A(1, 0) = -1.;
        A(0, 1) = -1;
        F.push_back(A);
    }
    // Set limit_rot
    limit_rot = (search_rot < 180.);
}

void ML2DBaseProgram::randomizeImagesOrder()
{
    //This static flag is for only randomize once
    static bool randomized = true;

    if (!randomized)
    {
        srand(seed);
        //-------Randomize the order of images
        std::random_shuffle(img_id.begin(), img_id.end());
        randomized = true;
    }
}//close function randomizeImagesOrder

void ML2DBaseProgram::setNumberOfLocalImages()
{
    nr_images_local = nr_images_global;
    //the following will be override in the MPI implementation.
    //  nr_images_local = divide_equally(nr_images_global, size, rank, myFirstImg,
    //                                   myLastImg);
}

// Check convergence
bool ML2DBaseProgram::checkConvergence()
{

#ifdef DEBUG
    std::cerr<<"entering checkConvergence"<<std::endl;
#endif

    if (iter == 0)
        return false;

    bool converged = true;
    double convv;
    MultidimArray<double> Maux;

    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    conv.clear();

    for (int refno = 0; refno < model.n_ref; refno++)
    {
        if (model.Iref[refno].weight() > 0.)
        {
            Maux = Iold[refno]() * Iold[refno]();
            convv = 1. / (Maux.computeAvg());
            Maux = Iold[refno]() - model.Iref[refno]();
            Maux = Maux * Maux;
            convv *= Maux.computeAvg();
            conv.push_back(convv);

            if (convv > eps)
                converged = false;
        }
        else
        {
            conv.push_back(-1.);
        }
    }

#ifdef DEBUG
    std::cerr<<"leaving checkConvergence"<<std::endl;

#endif

    return converged;
}//close function checkConvergence

//Standard run function of ML2D family
void ML2DBaseProgram::run()
{
    bool converged = false;

    produceSideInfo();
    //Do some initialization work
    produceSideInfo2();
    //Create threads to be ready for work
    createThreads();

    // Loop over all iterations
    for (iter = istart; !converged && iter <= Niter; iter++)
    {
        if (verbose)
            std::cout << "  Multi-reference refinement:  iteration " << iter << " of " << Niter << std::endl;

        for (int refno = 0;refno < model.n_ref; refno++)
            Iold[refno]() = model.Iref[refno]();

        //Perform an ML iteration
        iteration();

        // Check convergence
        converged = checkConvergence();

        // Write output files
        addPartialDocfileData(docfiledata, myFirstImg, myLastImg);
        writeOutputFiles(model, OUT_ITER);

    } // end loop iterations

    if (verbose)
    {
        std::cout << (converged ?
                      "--> Optimization converged!" :
                      "--> Optimization was stopped before convergence was reached!")
        << std::endl;
    }

    writeOutputFiles(model, OUT_FINAL);
    destroyThreads();
}

void ML2DBaseProgram::defineBasicParams(XmippProgram * prog)
{
    prog->addParamsLine("   -i <input_file>                : Metadata or stack with input images ");
    if (referenceExclusive)
    {
        prog->addParamsLine("   --nref <int=1>               : Number of references to generate automatically (recommended)");
        prog->addParamsLine("or --ref <reference_file=\"\">  : Image, stack or metadata with initial(s) references(s)");
    }
    else
        prog->addParamsLine("--ref <refence_file> <nref=1> : Volume, stack or metadata with initial volume references");

    prog->addParamsLine(formatString(" [ --oroot <rootname=%s> ]    : Output rootname", defaultRoot.c_str()));
    prog->addParamsLine(" [ --mirror ]                   : Also check mirror image of each reference ");

    if (allowFastOption)
    {
        prog->addParamsLine(" [ --fast ]                     : Use pre-centered images to pre-calculate significant orientations.");
        prog->addParamsLine(":++ If this flag is set part of the integration over all references, rotations and translations is skipped.");
        prog->addParamsLine(":++ The program will store all (=N_imgs*N_refs=) origin offsets that yield the maximum probability of observing");
        prog->addParamsLine(":++ each experimental image, given each of the references. In the first iterations a complete integration over");
        prog->addParamsLine(":++ all references, rotations and translations is performed for all images. In all subsequent iterations, for all");
        prog->addParamsLine(":++ combinations of experimental images, references and rotations, the probability of observing the image given");
        prog->addParamsLine(":++ the optimal origin offsets from the previous iteration is calculated. Then, if this probability is not");
        prog->addParamsLine(":++ considered \"significant\", we assume that none of the other translations will be significant, and we skip");
        prog->addParamsLine(":++ the integration over the translations. A combination of experimental image, reference and rotation is considered");
        prog->addParamsLine(":++ as \"significant\" if the probability at the corresponding optimal origin offsets is larger than C times the");
        prog->addParamsLine(":++ maximum of all these probabilities for that experimental image and reference (by default C=1e-12) This version");
        prog->addParamsLine(":++ may run up to ten times faster than the original, complete-search approach, while practically identical results may be obtained.");

    }
    if (allowThreads)
        prog->addParamsLine(" [ --thr <N=1> ]                : Use N parallel threads ");
    if (allowIEM)
        prog->addParamsLine(" [ --iem <blocks=1>]            : Number of blocks to be used with IEM");

}

void ML2DBaseProgram::defineAdditionalParams(XmippProgram * prog, const char * sectionLine)
{
    prog->addParamsLine(sectionLine);
    prog->addParamsLine(" [ --eps <float=5e-5> ]         : Stopping criterium");
    prog->addParamsLine(formatString(" [ --iter <int=%d> ]           : Maximum number of iterations to perform ", defaultNiter));
    prog->addParamsLine(" [ --psi_step <float=5.> ]       : In-plane rotation sampling interval [deg]");
    prog->addParamsLine(" [ --noise <float=1> ]          : Expected standard deviation for pixel noise ");
    prog->addParamsLine(" [ --offset <float=3.> ]         : Expected standard deviation for origin offset [pix]");
    prog->addParamsLine(" [ --frac <docfile=\"\"> ]      : Docfile with expected model fractions (default: even distr.)");
    prog->addParamsLine(" [ -C <double=1e-12> ]         : Significance criterion for fast approach ");
    prog->addParamsLine(" [ --zero_offsets ]             : Kick-start the fast algorithm from all-zero offsets ");
    if (allowRestart)
        prog->addParamsLine(" [ --restart <iter=1> ]  : restart a run with all parameters as in the logfile ");
    prog->addParamsLine(" [ --fix_sigma_noise ]           : Do not re-estimate the standard deviation in the pixel noise ");
    prog->addParamsLine(" [ --fix_sigma_offset ]          : Do not re-estimate the standard deviation in the origin offsets ");
    prog->addParamsLine(" [ --fix_fractions ]             : Do not re-estimate the model fractions ");
    prog->addParamsLine(" [ --student <df=6>]            : Use t-distributed instead of Gaussian model for the noise ");
    prog->addParamsLine("                                : df = Degrees of freedom for the t-distribution ");
    prog->addParamsLine(" [ --norm ]                     : Refined normalization parameters for each particle ");
    prog->addParamsLine(" [ --save_memA ]                : Save memory A(deprecated)");
    prog->addParamsLine(" [ --save_memB ]                : Save memory B(deprecated)");

}

void ML2DBaseProgram::defineHiddenParams(XmippProgram *prog)
{
  addParamsLine("==+++++ Hidden arguments ==");
  addParamsLine(" [--scratch <scratch=\"\">]");
  addParamsLine(" [--debug <int=0>]");
  addParamsLine(" [--no_sigma_trick]");
  addParamsLine(" [--trymindiff_factor <float=0.9>]");
  addParamsLine(" [--random_seed <int=-1>]");
  addParamsLine(" [--search_rot <float=999.>]");
  addParamsLine(" [--load <N=1>]");
}




///////////// ModelML2D Implementation ////////////
ModelML2D::ModelML2D()
{
    n_ref = -1;
    sumw_allrefs2 = 0;
    initData();

}//close default constructor

ModelML2D::ModelML2D(int n_ref)
{
    sumw_allrefs2 = 0;
    initData();
    setNRef(n_ref);
}//close constructor

void ModelML2D::initData()
{
    do_student = do_norm = false;
    do_student_sigma_trick = true;
    sumw_allrefs = sigma_noise = sigma_offset = LL = avePmax = 0;
    dim = 0;
}//close function initData

/** Before call this function model.n_ref should
 * be properly setted. */
void ModelML2D::setNRef(int n_ref)
{
    Image<double> Iempty;
    Iempty().initZeros(dim,dim);
    Iempty().setXmippOrigin();
    this->n_ref = n_ref;
    Iref.resize(n_ref, Iempty);
    alpha_k.resize(n_ref, 0.);
    mirror_fraction.resize(n_ref, 0.);
    scale.resize(n_ref, 1.);

}//close function setNRef

void ModelML2D::combineModel(ModelML2D model, int sign)
{
    if (n_ref != model.n_ref)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "Can not add models with different 'n_ref'");

    double sumw, sumw_mirror, sumwsc, sumweight;
    double wsum_sigma_offset = getWsumSigmaOffset() + sign
                               * model.getWsumSigmaOffset();
    double wsum_sigma_noise = getWsumSigmaNoise() + sign
                              * model.getWsumSigmaNoise();
    double local_sumw_allrefs = sumw_allrefs + sign * model.sumw_allrefs;
    double sumfracweight = getSumfracweight() + sign
                           * model.getSumfracweight();

    for (int refno = 0; refno < n_ref; refno++)
    {
        sumweight = Iref[refno].weight() + sign * model.Iref[refno].weight();
        if (sumweight > 0)
        {
            Iref[refno]() = (getWsumMref(refno) + sign * model.getWsumMref(
                                 refno)) / sumweight;
            Iref[refno].setWeight(sumweight);
        }
        else
        {
            //std::cerr << "sumweight: " << sumweight << std::endl;
            Iref[refno]().initZeros();
            Iref[refno].setWeight(0);
        }

        //Get all sums first, because function call will change
        //after updating model parameters.
        sumw = getSumw(refno) + sign * model.getSumw(refno);
        sumw_mirror = getSumwMirror(refno) + sign * model.getSumwMirror(
                          refno);
        sumwsc = getSumwsc(refno) + sign * model.getSumwsc(refno);

        //Update parameters
        //alpha_k[refno] = sumw / local_sumw_allrefs;
        //mirror_fraction[refno] = sumw_mirror / sumw;
        updateFractions(refno, sumw, sumw_mirror, local_sumw_allrefs);
        //scale[refno] = sumwsc / sumw;
        updateScale(refno, sumwsc, sumw);
    }

    sumw_allrefs = local_sumw_allrefs;
    sumw_allrefs2 += sign * model.sumw_allrefs2;

    updateSigmaNoise(wsum_sigma_noise);
    (wsum_sigma_offset);
    updateAvePmax(sumfracweight);
    LL += sign * model.LL;

}//close function combineModel

void ModelML2D::addModel(ModelML2D model)
{
    combineModel(model, 1);
}//close function addModel

void ModelML2D::substractModel(ModelML2D model)
{
    combineModel(model, -1);
}//close function substractModel

double ModelML2D::getSumw(int refno) const
{
    return alpha_k[refno] * sumw_allrefs;
}//close function sumw

double ModelML2D::getSumwMirror(int refno) const
{
    return getSumw(refno) * mirror_fraction[refno];
}//close function sumw_mirror

double ModelML2D::getSumwsc(int refno) const
{
    return scale[refno] * getSumw(refno);
}//close function get_sumwsc

MultidimArray<double> ModelML2D::getWsumMref(int refno) const
{
    return Iref[refno]() * Iref[refno].weight();
}//close function get_wsum_Mref

double ModelML2D::getWsumSigmaOffset() const
{
    return sigma_offset * sigma_offset * 2 * sumw_allrefs;
}//close function get_wsum_sigma_offset

double ModelML2D::getWsumSigmaNoise() const
{
    double sum = (do_student && do_student_sigma_trick) ? sumw_allrefs2
                 : sumw_allrefs;
    return sigma_noise * sigma_noise * dim * dim * sum;
}//close function get_wsum_sigma_noise

double ModelML2D::getSumfracweight() const
{
    return avePmax * sumw_allrefs;
}//close function get_sumfracweight

void ModelML2D::updateSigmaOffset(double wsum_sigma_offset)
{
    if (sumw_allrefs == 0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "'sumw_allrefs' couldn't be zero");
    sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs));
    if (wsum_sigma_offset < 0.)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "sqrt of negative 'wsum_sigma_offset'");
    if (sumw_allrefs < 0.)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "sqrt of negative 'wsum_sigma_offset'");
}//close function updateSigmaOffset

void ModelML2D::updateSigmaNoise(double wsum_sigma_noise)
{
    // The following converges faster according to McLachlan&Peel (2000)
    // Finite Mixture Models, Wiley p. 228!
    double sum = (do_student && do_student_sigma_trick) ? sumw_allrefs2
                 : sumw_allrefs;
    if (sum == 0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "'sumw_allrefs' couldn't be zero");

    double sigma_noise2 = wsum_sigma_noise / (sum * dim * dim);
    if (sigma_noise2 < 0.)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "sqrt of negative 'sigma_noise2'");
    sigma_noise = sqrt(sigma_noise2);
}//close function updateSigmaNoise

void ModelML2D::updateAvePmax(double sumfracweight)
{
    avePmax = sumfracweight / sumw_allrefs;
}//close function updateAvePmax

void ModelML2D::updateFractions(int refno, double sumw,
                                double sumw_mirror, double sumw_allrefs)
{
    if (sumw_allrefs == 0)
    {
        REPORT_ERROR(ERR_VALUE_INCORRECT, "updateFractions: sumw_allrefs == 0 ");
    }

    if (sumw > 0.)
    {
        alpha_k[refno] = sumw / sumw_allrefs;
        mirror_fraction[refno] = sumw_mirror / sumw;
    }
    else
    {
        alpha_k[refno] = 0.;
        mirror_fraction[refno] = 0.;
    }
}//close updateFractions

void ModelML2D::updateScale(int refno, double sumwsc, double sumw)
{
    if (do_norm)
        scale[refno] = (sumw > 0) ? sumwsc / sumw : 1;

}//close function updateScale

void ModelML2D::print() const
{
    std::cerr << "sumw_allrefs: " << sumw_allrefs << std::endl;
    std::cerr << "wsum_sigma_offset: " << getWsumSigmaOffset() << std::endl;
    std::cerr << "wsum_sigma_noise: " << getWsumSigmaNoise() << std::endl;
    std::cerr << "sigma_offset: " << sigma_offset << std::endl;
    std::cerr << "sigma_noise: " << sigma_noise << std::endl;
    std::cerr << "LL: " << LL << std::endl;

    for (int refno = 0; refno < n_ref; refno++)
    {
        std::cerr << "refno:       " << refno << std::endl;
        std::cerr << "sumw:        " << getSumw(refno) << std::endl;
        std::cerr << "sumw_mirror: " << getSumwMirror(refno) << std::endl;
        std::cerr << "alpha_k:        " << alpha_k[refno] << std::endl;
        std::cerr << "mirror_fraction: " << mirror_fraction[refno] << std::endl;

    }

}//close function print



