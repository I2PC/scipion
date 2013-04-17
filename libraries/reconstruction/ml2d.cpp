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

#include <algorithm>
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
    blocks = 1;
    factor_nref = 1;
    outRefsMd = "";
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
    static bool randomized = false;

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
//#define DEBUG_JM
#ifdef DEBUG_JM
    std::cerr<<">>>>>>>> entering checkConvergence"<<std::endl;
#endif

    if (iter == 0)
        return false;

    bool converged = true;
    double convv, avg2;
    MultidimArray<double> Maux;

    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    conv.clear();

    for (int refno = 0; refno < model.n_ref; refno++)
    {
        convv = -1;

        if (model.Iref[refno].weight() > 0.)
        {
            Maux = Iold[refno]() * Iold[refno]();
            avg2 = Maux.computeAvg();
            if (avg2 > 0.)
            {
              convv = 1. / avg2;
              Maux = Iold[refno]() - model.Iref[refno]();
              Maux = Maux * Maux;
              convv *= Maux.computeAvg();
              if (convv > eps)
                  converged = false;
            }
        }
        conv.push_back(convv);
    }

#ifdef DEBUG_JM
    std::cerr<<"<<<<<< leaving checkConvergence"<<std::endl;

#endif
#undef DEBUG_JM

    return converged;
}//close function checkConvergence

void ML2DBaseProgram::endIteration()
{
    // Write output files
    addPartialDocfileData(docfiledata, myFirstImg, myLastImg);
    writeOutputFiles(model, OUT_ITER);
}

//Standard run function of ML2D family
void ML2DBaseProgram::run()
{
    bool converged = false;
    CREATE_LOG(LOG_FN(fn_root));

    LOG(" starting produceSideInfo\n");
    produceSideInfo();

    //Do some initialization work
    LOG(" starting produceSideInfo2\n");
    produceSideInfo2();

    //Create threads to be ready for work
    LOG(" createThreads\n");
    createThreads();

    LOG(" starting iterations\n");
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

        // Do some task before ending iteration
        endIteration();

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
    CLOSE_LOG();
}

void ML2DBaseProgram::defineBasicParams(XmippProgram * prog)
{
    prog->addParamsLine("   -i <input_file>                : Metadata or stack with input images ");

    String orStr = "", lb = "", rb = "";
    if (!referenceExclusive)
    {
        orStr = "";
        lb = "[";
        rb = "]";
    }
    prog->addParamsLine(" [ --nref <int=1> ]                  : Number of references to generate automatically (recommended)");
    prog->addParamsLine(" [ --ref <reference_file=\"\">  ]    : Image, stack or metadata with initial(s) references(s)");

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
    prog->addParamsLine("==+++++ Hidden arguments ==");
    prog->addParamsLine(" [--scratch <scratch=\"\">]");
    prog->addParamsLine(" [--debug <int=0>]");
    prog->addParamsLine(" [--no_sigma_trick]");
    prog->addParamsLine(" [--trymindiff_factor <float=0.9>]");
    prog->addParamsLine(" [--random_seed <int=-1>]");
    prog->addParamsLine(" [--search_rot <float=999.>]");
    prog-> addParamsLine(" [--load <N=1>]");

    //fixme: only for debug
    prog->addParamsLine("[--no_iem] : bla bla bla");
}

FileName getIterExtraPath(const FileName &fn_root, int iter, bool makePath)
{
  FileName fn = formatString("%sextra/iter%03d/", fn_root.c_str(), iter);
  if (makePath)
    fn.makePath();
  fn += "iter_";
  return fn;
}



///////////// ModelML2D Implementation ////////////
ModelML2D::ModelML2D()
{
    initData();
}//close default constructor

ModelML2D::ModelML2D(int n_ref)
{
    initData();
    setNRef(n_ref);
}//close constructor

void ModelML2D::initData()
{
    n_ref = -1;
    sumw_allrefs2 = sumw_allrefs = 0.;
    do_student = do_norm = false;
    do_student_sigma_trick = true;
    sigma_noise = sigma_offset = LL = avePmax = 0.;
    wsum_sigma_noise = wsum_sigma_offset = 0.;
    sumfracweight = 0.;
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
    WsumMref.resize(n_ref, Iempty);
    alpha_k.resize(n_ref, 0.);
    mirror_fraction.resize(n_ref, 0.);
    scale.resize(n_ref, 1.);
    sumw_mirror.resize(n_ref, 0.);
    sumwsc.resize(n_ref, 0.);

}//close function setNRef



void ModelML2D::combineModel(const ModelML2D &model, int sign)
{
    if (n_ref != model.n_ref)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "Can not add models with different 'n_ref'");

#define COMBINE(var) var += sign * model.var
    COMBINE(wsum_sigma_offset);
    COMBINE(wsum_sigma_noise);
    COMBINE(sumw_allrefs);
    COMBINE(sumw_allrefs2);
    COMBINE(sumfracweight);
    COMBINE(LL);

    MultidimArray<double> tmp;
    double sumweight = 0., w1 = 0., w2 = 0.;

    for (int refno = 0; refno < n_ref; refno++)
    {
        w1 = WsumMref[refno].weight();
        w2 = sign * model.WsumMref[refno].weight();
        tmp = model.WsumMref[refno]();
        tmp *= sign;
        WsumMref[refno]() += tmp;
        sumweight = w1 + w2;

        if (sumweight > 0.)
        {
            Iref[refno].setWeight(sumweight);
            WsumMref[refno].setWeight(sumweight);
            COMBINE(sumw_mirror[refno]);
            COMBINE(sumwsc[refno]);
        }
        else
        {
            Iref[refno]().initZeros();
            WsumMref[refno]().initZeros();
            Iref[refno].setWeight(0);
            WsumMref[refno].setWeight(0.);
            sumw_mirror[refno] = 0.;
            sumwsc[refno] = 0.;
        }
    }
}//close function combineModel

void ModelML2D::addModel(const ModelML2D &model)
{
    combineModel(model, 1);
}//close function addModel

void ModelML2D::substractModel(const ModelML2D &model)
{
    combineModel(model, -1);
}//close function substractModel

void ModelML2D::update()
{
    if (sumw_allrefs < 0)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "updateFractions: sumw_allrefs should be greater than 0 ");

    //update sigma_noise
    double sum = (do_student && do_student_sigma_trick) ? sumw_allrefs2
                 : sumw_allrefs;
    double sigma_noise2 = wsum_sigma_noise / (sum * dim * dim);

    if (sigma_noise2 < 0.)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "sqrt of negative 'sigma_noise2'");
    sigma_noise = sqrt(sigma_noise2);

    //update sigma_offset
    if (wsum_sigma_offset < 0.)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "sqrt of negative 'wsum_sigma_offset'");
    sigma_offset = sqrt(wsum_sigma_offset / (2. * sumw_allrefs));

    //update avePmax
    avePmax = sumfracweight / sumw_allrefs;

    double weight, inv_weight;

//    std::cerr << "DEBUG_JM: sumw_allrefs: " << sumw_allrefs << std::endl;

    for (int refno = 0; refno < n_ref; ++refno)
    {
        weight = WsumMref[refno].weight();

//        std::cerr << "    DEBUG_JM: refno: " << refno << std::endl;

        if (weight > 0.)
        {
            //update weights
            Iref[refno].setWeight(weight);
            inv_weight = 1 / weight; //just to speed-up
            Iref[refno]() = WsumMref[refno]();
            Iref[refno]() *= inv_weight;
            //update fractions
            alpha_k[refno] = weight / sumw_allrefs;
            mirror_fraction[refno] = sumw_mirror[refno] * inv_weight;
            scale[refno] = sumwsc[refno] * inv_weight;
//            std::cerr << "    DEBUG_JM: weight:     " << weight << std::endl;
//            std::cerr << "    DEBUG_JM: sumw_mirror[refno]: " << sumw_mirror[refno] << std::endl;
        }
        else
        {
//          std::cerr << "    DEBUG_JM: all zeros... " << std::endl;
            //zero weights
            Iref[refno].setWeight(0.);
            Iref[refno]().initZeros(WsumMref[refno]());
            alpha_k[refno] = 0.;
            mirror_fraction[refno] = 0.;
            scale[refno] = 0.;
        }
    }
}

#define pp(s, x) s << std::setw(10) << x

void ModelML2D::print(int tabs) const
{
    String stabs = "";
    for (int t = 0; t < tabs; ++t)
        stabs += " ";

    std::cerr << "======================= Block ==================== " << std::endl;

    std::cerr << "sumw_allrefs: " << sumw_allrefs << std::endl;
    std::cerr << "wsum_sigma_offset: " << wsum_sigma_offset << std::endl;
    std::cerr << "wsum_sigma_noise: " << wsum_sigma_noise << std::endl;
    std::cerr << "sigma_offset: " << sigma_offset << std::endl;
    std::cerr << "sigma_noise: " << sigma_noise << std::endl;
    std::cerr << "LL: " << LL << std::endl;

    std::stringstream ss1, ss2, ss3, ss4, ss5, ss6, ss7;
    pp(ss1, "refno:");
    pp(ss3, "sumw_mirror:");
    pp(ss4, "alpha_k:");
    pp(ss5, "mirror_fraction");
    pp(ss6, "WsumMref.weight");
    pp(ss7, "Iref.weight");

    for (int refno = 0; refno < n_ref; refno++)
    {
        pp(ss1, refno);
        pp(ss3, sumw_mirror[refno]);
        pp(ss4, alpha_k[refno]);
        pp(ss5, mirror_fraction[refno]);
        pp(ss6, WsumMref[refno].weight());
        pp(ss7, Iref[refno].weight());
    }
    std::cerr << ss1.str() << std::endl<< ss3.str() << std::endl
    << ss4.str() << std::endl<< ss5.str() << std::endl << ss6.str() << std::endl << ss7.str()<< std::endl;

}//close function print





