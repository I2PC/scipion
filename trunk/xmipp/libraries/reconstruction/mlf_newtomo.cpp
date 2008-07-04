/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es
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
#include "mlf_newtomo.h"
#define DEBUG

// Read arguments ==========================================================
void Prog_mlf_tomo_prm::read(int argc, char **argv)
{

    // Read command line
    if (checkParameter(argc, argv, "-more_options"))
    {
        usage();
        extendedUsage();
    }
    SFi.read(getParameter(argc, argv, "-i"));
    SFw.read(getParameter(argc, argv, "-wedge"));
    fn_group = getParameter(argc, argv, "-groups","");
    nr_ref = textToInteger(getParameter(argc, argv, "-nref", "0"));
    fn_ref = getParameter(argc, argv, "-ref", "");
    fn_root = getParameter(argc, argv, "-o", "out");
    Niter = textToInteger(getParameter(argc, argv, "-iter", "100"));
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));
    fix_fractions = checkParameter(argc, argv, "-fix_fractions");
    fix_sigma_noise = checkParameter(argc, argv, "-fix_sigma_noise");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    highres = textToFloat(getParameter(argc, argv, "-highres", "-1"));
    lowres = textToFloat(getParameter(argc, argv, "-lowres", "-1"));

    debug = checkParameter(argc, argv, "-debug");

}

// Show ====================================================================
void Prog_mlf_tomo_prm::show()
{

    if (verb > 0)
    {
        // To screen
        std::cerr << "--> Maximum-likelihood sub-tomogram classification " << std::endl;
        std::cerr << "  Input subtomograms      : " << SFi.name() << " (" << SFi.ImgNo() << ")" << std::endl;
        if (fn_ref != "")
            std::cerr << "  References              : " << fn_ref << std::endl;
        else
            std::cerr << "  Number of references:   : " << nr_ref << std::endl;
        std::cerr << "  Output rootname         : " << fn_root << std::endl;
        std::cerr << "  Number of iterations    : " << Niter << std::endl;
        if (fix_fractions)
        {
            std::cerr << "  -> Do not update estimates of model fractions." << std::endl;
        }
        if (fix_sigma_noise)
        {
            std::cerr << "  -> Do not update sigma-estimate of noise." << std::endl;
        }
        if (lowres > 0. || highres > 0.)
        {
            std::cerr << "  -> Limit to resolutions between " << lowres << " and " << highres << std::endl;
        }
        std::cerr << " -----------------------------------------------------------------" << std::endl;
    }

}

// Usage ===================================================================
void Prog_mlf_tomo_prm::usage()
{
    std::cerr << "Usage:  mlf_tomo [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with the input sub-tomograms \n";
    std::cerr << "   -wedge <selfile>            : Selfile with the point-spread function of the missing wedges \n";
    std::cerr << "  [ -groups <selfile> ]         : Selfile with the group numbers for all subtomograms \n";
    std::cerr << "   -nref <int>                 : Number of references to generate automatically (recommended)\n";
    std::cerr << "   OR -ref <selfile/volume>         OR selfile with initial references/single reference image \n";
    std::cerr << " [ -o <rootname=\"out\"> ]       : Output rootname \n";
    std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_mlf_tomo_prm::extendedUsage()
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
    std::cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
    std::cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
    std::cerr << std::endl;
    exit(1);
}

// This routine is for SF-independent side-info calculations
void Prog_mlf_tomo_prm::produceSideInfo(double * dataRefs, double * dataSigma)
{

    headerXmipp head;
    Matrix3D<double> Maux, Faux_real, Faux_imag;
    Matrix3D<std::complex<double> > Faux, Faux2;
    int xdim, ydim, zdim, c, iaux, ifound;
    float xx, yy;
    double dum, avg;

    // Get SFr
    if (fn_ref != "")
    {
        if (Is_VolumeXmipp(fn_ref)) 
        {
            nr_ref = 1;
            SFr.insert(fn_ref);
        }
        else
        {
            SFr.read(fn_ref);
            nr_ref = SFr.ImgNo();
        }
    }
    else
    {
        generateInitialReferences();
    }

    // image sizes
    VolumeXmipp vol;
    SFr.go_beginning();
    vol.read(SFr.NextImg());
    Xdim = XSIZE(vol());
    Ydim = YSIZE(vol());
    Zdim = ZSIZE(vol());
    dim3 = Xdim * Ydim * Zdim;

    // Make FFTW objects and plans for forward and backward fftw objects
    int fNdim = 3;
    int * fN ;
    fN = new int[fNdim];
    fN[0] = XSIZE(vol());
    fN[1] = YSIZE(vol());
    fN[2] = ZSIZE(vol());
    forwfftw.myxmippFftw(fNdim, fN, false, NULL);
    forwfftw.Init("ES",FFTW_FORWARD,false);
    backfftw.myxmippFftw(fNdim, fN, false, NULL);
    backfftw.Init("ES",FFTW_BACKWARD,false);
    // Just for now (without resolution limits)
    // next line copied from fftw.cpp
    size = 2*int(double(forwfftw.GetSize())*(fN[fNdim-1]/2+1)/fN[fNdim-1]);

    // Get groups structure from SFg and store into SFi
    SelFile SFtmp;
    SelLine SL;
    nr_group = 1;
    if (fn_group != "")
    {
        SFg.read(fn_group);
        SFi.go_beginning();
        while (!SFi.eof())
        {
            SFg.search(SFi.NextImg());
            SL = SFg.current();
            SFtmp.insert(SL);
            nr_group = XMIPP_MAX(nr_group, SL.get_number());
        }
        SFi = SFtmp;
    }

 }


void Prog_mlf_tomo_prm::generateInitialReferences()
{

    SelFile SFtmp;
    VolumeXmipp Vave, Vtmp;
    double dummy;
    FileName fn_tmp;
    SelLine line;

    if (verb > 0)
    {
        std::cerr << "  Generating initial references by averaging over random subsets" << std::endl;
        init_progress_bar(nr_ref);
    }

    // Make random subsets and calculate average images
    FileName fnt;
    SFtmp = SFi.randomize();
    int Nsub = ROUND((double)SFtmp.ImgNo() / nr_ref);
    for (int refno = 0; refno < nr_ref; refno++)
    {
        SFtmp.go_beginning();
        if (Nsub*refno>0)
            SFtmp.jump_lines(Nsub*refno);
        if (refno == nr_ref - 1) Nsub = SFtmp.ImgNo() - refno * Nsub;
        for (int nn = 0; nn < Nsub; nn++)
        {
            fnt = SFtmp.NextImg();
            if (nn==0) 
                Vave.read(fnt);
            else
            {
                Vtmp.read(fnt);
                Vave() += Vtmp();
            }
        }
        Vave() /= Nsub;
        fn_tmp = fn_root + "_it";
        fn_tmp.compose(fn_tmp, 0, "");
        fn_tmp = fn_tmp + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".vol";
        Vave.write(fn_tmp);
        SFr.insert(fn_tmp, SelLine::ACTIVE);
        if (verb > 0) progress_bar(refno);
    }
    if (verb > 0) progress_bar(nr_ref);
    fn_ref = fn_root + "_it";
    fn_ref.compose(fn_ref, 0, "sel");
    SFr.write(fn_ref);

}

void Prog_mlf_tomo_prm::readAndFftwAllReferences(double * dataRefs)
{

#ifdef DEBUG
    std::cerr<<"start readAndFftwAllReferences"<<std::endl;
#endif

    // Read references in memory and calculate their FFTW
    VolumeXmipp vol;
    Matrix1D<double> data(size);

    alpha_k.clear();
    int refno = 0;
    SFr.go_beginning();
    while (!SFr.eof())
    {
        alpha_k.push_back(1./nr_ref);
        vol.read(SFr.NextImg());
        forwfftw.SetPoints(MULTIDIM_ARRAY(vol()));
        forwfftw.Transform();
        // TODO: GetPoints for certain resolution range...
        forwfftw.GetPoints(MULTIDIM_ARRAY(data));
        // Store all references in dataRefs structure
        for (int i = 0; i < size; i++)
        {
            dataRefs[refno*size + i] = data(i);
        }
        refno++;
    }
#ifdef DEBUG
    std::cerr<<"done readAndFftwAllReferences"<<std::endl;
#endif

}


// This routine is for (splitted) SF-dependent side-info calculations
void Prog_mlf_tomo_prm::calculateAllFFTWs(SelFile &SF)
{
#ifdef DEBUG
    std::cerr<<"start calculateAllFFTWs"<<std::endl;
#endif

    FileName fni, fno;
    VolumeXmipp Vin;
    int c, nn = SF.ImgNo(), imgno = 0;

    if (verb > 0)
    {
        std::cerr << "  Calculating FFTs for all input maps ... " << std::endl;
        init_progress_bar(nn);
        c = XMIPP_MAX(1, nn / 60);
    }

    SF.go_beginning();
    while (!SF.eof())
    {
        fni = SF.NextImg();
        fno = fni + ".fftw";
        if (!exists(fno))
        {
            Vin.read(fni);
            forwfftw.SetPoints(MULTIDIM_ARRAY(Vin()));
            forwfftw.Transform();
            forwfftw.write(fno);
        }
        if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
        imgno++;
    }
    if (verb > 0) progress_bar(nn);

#ifdef DEBUG
    std::cerr<<"done calculateAllFFTWs"<<std::endl;
#endif
}

// Estimate initial sigma2 for fourier-mode from power spectra of volumes
void Prog_mlf_tomo_prm::estimateInitialNoiseSpectra(double * dataSigma)
{
#ifdef DEBUG
    std::cerr<<"start estimateInitialNoiseSpectra"<<std::endl;
#endif

    FileName fn;
    Matrix1D<double> data;
    std::vector<Matrix1D<double> > sum, sum2;
    Matrix1D< int > radial_count;
    SelLine SL;

    if (verb > 0)
    {
        std::cerr << "  Estimating initial noise spectra ... " << std::endl;
    }

    // Initialize data, sum and sum2
    data.initZeros(size);
    for (int ig = 0; ig < nr_group; ig++)
    {
        sum.push_back(data);
        sum2.push_back(data);
    }

    SFi.go_beginning();
    while (!SFi.eof())
    {
        SL = SFi.current();
        int igg = SL.get_number() - 1;
        fn = SFi.NextImg();
        fn += ".fftw";
        backfftw.read(fn);
        backfftw.GetPoints(MULTIDIM_ARRAY(data), true);//TODO: add resolutionlimits!);
        for (int i=0; i<size; i++)
        {
            sum[igg](i) += data(i);
            sum2[igg](i) += data(i) * data(i);
        }
    }

    // Subtract squared amplitudes of the average subtomogram 
    // to prevent overestimated noise at low resolutions
    for (int ig = 0; ig < nr_group; ig++)
    {
        for (int i=0; i<size; i++)
        {
            sum2[ig](i) -= sum[ig](i) * sum[ig](i);
        }
        // Perform radial averaging in FFTW 
        // forwfftw.fftwRadialAverage(MULTIDIM_ARRAY(sum2[ig]), sigma_noise[ig], radial_count);
        // TODO: put radial average values in sum2[ig]
        // Now store in dataSigma structure (TODO: resolution limits!)
        for (int i=0; i<size; i++)
        {
            dataSigma[ig * size + i] = 2. * sum2[ig](i);
        }
        // TODO: Write sigma_noise[ig] to disc
    }

#ifdef DEBUG
    std::cerr<<"done estimateInitialNoiseSpectra"<<std::endl;
#endif
}


// Here perform the main probability-weighted integration over all
// rotations, translations and classes of the given image
void Prog_mlf_tomo_prm::expectationSingleImage(int igroup,
                                               double * dataImg,
                                               double * dataWedge,
                                               double * dataRefs,
                                               double * dataSigma,
                                               double * dataWsumRefs,
                                               double * dataWsumWedsPerRef,
                                               double * dataWsumWedsPerGroup,
                                               double * dataWsumDist,
                                               double * dataSumWRefs,
                                               int &opt_refno, 
                                               double &LL, 
                                               double &Pmax)

{

    double aux, weight, diff2, sumweight = 0., maxweight=0., mindiff2=99.e99;
    double * weights;
    weights = new double [nr_ref];

    for (int refno = 0; refno < nr_ref; refno++)
    {
        diff2 = 0.;
        for (int i = 0; i < size; i++)
        {
            if (dataWedge[i] > 0. && dataSigma[igroup * size + i] > 0.)
            {
                aux = dataImg[i] - dataRefs[refno * size + i];
                diff2 += (aux * aux)/dataSigma[igroup * size + i];
            }
        }
        weights[refno] = diff2;
        if (diff2 < mindiff2)
        {
            mindiff2 = diff2;
        }
    }

    // Now that we have mindiff2, calculate the actual weights
    for (int refno = 0; refno < nr_ref; refno++)
    {
        aux = weights[refno] - mindiff2;
        if (aux > 1000.) aux = 0.;
        else aux = exp(-aux) * alpha_k[refno];
        if (aux > maxweight)
        {
            maxweight = aux;
            opt_refno = refno;
        }
        weights[refno] = aux;
        sumweight += aux;
        std::cerr<<"refno= "<<refno<<" w= "<< aux<<std::endl;
    }


    // Store Pmax/sumP
    Pmax = maxweight / sumweight;
    std::cerr<<"Pmax= "<<Pmax<<std::endl;

    // Then, store the sum of all weights
    for (int refno = 0; refno < nr_ref; refno++)
    {
        weight = weights[refno] / sumweight;
        if (weight > SIGNIFICANT_WEIGHT_LOW)
        {
            dataSumWRefs[refno] += weight;
            for (int i = 0; i < size; i++)
            {
                int iiref = refno * size + i;
                int iig = igroup * size + i;
                if (dataWedge[i] > 0.)
                {
                    aux = dataImg[i] - dataRefs[iiref];
                    dataWsumDist[iig] += weight * aux * aux;
                    dataWsumWedsPerGroup[iig] += weight;
                    dataWsumRefs[iiref] += weight * dataImg[i];
                    dataWsumWedsPerRef[iiref] += weight;
                }
                if (i==320)
                    std::cerr<<" dataImg[320]= "<<dataImg[i]
                             <<" dataWsumDist[iig]= "<<dataWsumDist[iig]
                             <<" dataWsumWedsPerGroup[iig]= "<<dataWsumWedsPerGroup[iig]
                             <<" dataWsumRefs[iiref]= "<<dataWsumRefs[iiref]
                             <<" dataWsumWedsPerRef[iiref]= "<<dataWsumWedsPerRef[iiref]<<std::endl;
            }
        }
    }
    
    // Update the log-likelihood function value
    LL+= 1.; //TODO;

}



void Prog_mlf_tomo_prm::expectation(SelFile &mySFi, 
                                    SelFile &mySFw,
                                    double * dataRefs,
                                    double * dataSigma,
                                    double * dataWsumRefs,
                                    double * dataWsumWedsPerRef,
                                    double * dataWsumWedsPerGroup,
                                    double * dataWsumDist,
                                    double * dataSumWRefs,
                                    double * imgsPmax,
                                    int * imgsOptRefNos,
                                    double &LL)
{

#ifdef DEBUG
    std::cerr<<"start expectation"<<std::endl;
#endif
    int opt_refno, nn, imgno, nr_imgs, igroup;
    double Pmax;
    double *dataImg, *dataWedge;
    SelLine SL;
    FileName fn;

    // Check that selfiles are of equal length
    nr_imgs = mySFi.ImgNo();
    if (mySFi.ImgNo() != mySFw.ImgNo())
        REPORT_ERROR(1,"Selfiles of images and wedges have unequal lengths");

    // Reserve memory for dataImg and dataWedge
    try
    {
        dataImg = new double[size];
        dataWedge = new double[size];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in expectation");
    }
    // FOR NOW ALL WEGDES ARE 1:
    for (int i = 0; i < size; i++)
        dataWedge[i] = 1.;

    // Initialize weighted sums to zero
    for (int i = 0; i < nr_ref * size; i++)
    {
        dataWsumRefs[i] = 0.;
        dataWsumWedsPerRef[i] = 0.;
    }        
    for (int i = 0; i < nr_group * size; i++)
    {
        dataWsumWedsPerGroup[i] = 0.;
        dataWsumDist[i] = 0.;
    }
    for (int i = 0; i < nr_ref; i++)
    {
        dataSumWRefs[i] = 0.;
    }
    LL = 0.;

    // Loop over all images
    imgno = 0;
    nn = mySFi.ImgNo();
    if (verb > 0) init_progress_bar(nn);
    mySFi.go_beginning();
    mySFw.go_beginning();
    while ((!mySFi.eof()))
    {
        // Get the group
        SL = mySFi.current();
        igroup = SL.get_number() - 1;        

        // read tomogram from disc
        fn = mySFi.NextImg();
        backfftw.read(fn + ".fftw");
        backfftw.GetPoints(dataImg,true);
        //Again somehow limit resolution

/*
        VolumeXmipp vol(Xdim,Ydim,Zdim);
        double * data;
        data = new double[size];
        for (int i=0; i< size; i++)
        {
            data[i] = dataImg[i];
        }
        backfftw.SetPoints(dataImg);
        backfftw.Transform();
        backfftw.Normalize();
        backfftw.GetPoints(MULTIDIM_ARRAY(vol()));
        vol.write("img.vol"); exit(0);
        VolumeXmipp vol(Xdim,Ydim,Zdim);
        std::cerr<<"set dataimg in backfft"<<std::endl;
        backfftw.SetPoints(dataImg);
        backfftw.Transform();
        backfftw.Normalize();
        std::cerr<<"get the volume "<<std::endl;
        backfftw.GetPoints(MULTIDIM_ARRAY(vol()));
        vol.write("img.vol"); exit(0);
*/

        // read corresponding wedge as well
        //backfftw.read(mySFw.NextImg(), true);
        //Again someghow limit resolution
      
        // Perform expectation step 
        expectationSingleImage(igroup,
                               dataImg, 
                               dataWedge,
                               dataRefs,
                               dataSigma,
                               dataWsumRefs,
                               dataWsumWedsPerRef,
                               dataWsumWedsPerGroup,
                               dataWsumDist,
                               dataSumWRefs,
                               opt_refno, 
                               LL, 
                               Pmax);

        imgsOptRefNos[imgno] = opt_refno;
        imgsPmax[imgno] = Pmax;

#ifdef DEBUG
        std::cerr<<fn<<" belongs to class "<<opt_refno + 1 <<std::endl;
#endif

        imgno++;
        if (verb > 0) progress_bar(imgno);

    }
    if (verb > 0) progress_bar(nn);
#ifdef DEBUG
    std::cerr<<"done expectation"<<std::endl;
#endif
}

// Update all model parameters
void Prog_mlf_tomo_prm::maximization(double * dataRefs,
                                     double * dataSigma,
                                     double * dataWsumRefs,
                                     double * dataWsumWedsPerRef,
                                     double * dataWsumWedsPerGroup,
                                     double * dataWsumDist,
                                     double * dataSumWRefs,
                                     double * imgsPmax,
                                     int * imgsOptRefNos,
                                     double &sumw_allrefs )
{

#ifdef DEBUG
    std::cerr<<"start maximization"<<std::endl;
#endif

    // Update References
    sumw_allrefs = 0;
    for (int refno = 0;refno < nr_ref; refno++)
    {
        sumw_allrefs += dataSumWRefs[refno];
        for (int i = 0; i < size; i++)
        {
            int ii = refno * size + i;
            // Impute old reference for missing pixels
            dataRefs[ii] *= 1. - (dataWsumWedsPerRef[ii] / dataSumWRefs[refno]);
            // And sum the weighted sum for observed pixels
            dataRefs[ii] += dataWsumRefs[ii] / dataSumWRefs[refno];
        }
    }

    // Update fractions
    if (!fix_fractions)
    {
        for (int refno = 0; refno < nr_ref; refno++)
            alpha_k[refno] = dataSumWRefs[refno] / sumw_allrefs;
    }

    // Update sigma of the noise
    /////// TODO: radial averaging etc again
    for (int ig = 0; ig < nr_group; ig++)
    {
        for (int i = 0; i < size; i++)
        {
            int ii = ig * size + i;
            // Impute old sigma values for missing pixels
            dataSigma[ii] *= 1. - (dataWsumWedsPerGroup[ii] / sumw_allrefs);
            //  And sum the weighted sum for observedpixels
            dataSigma[ii] += dataWsumDist[ii] / sumw_allrefs;
        }
    }

#ifdef DEBUG
    std::cerr<<"done maximization"<<std::endl;
#endif
}


void Prog_mlf_tomo_prm::writeOutputFiles(int iter, 
                                         double * dataRefs,
                                         double &sumw_allrefs, 
                                         double &LL, 
                                         double &avePmax,
                                         std::vector<double> &conv)
{

#ifdef DEBUG
    std::cerr<<"start writeOutputFiles"<<std::endl;
#endif
    FileName fn_base, fn_tmp;
    Matrix1D<double>  fracline(2);
    DocFile           DFl;
    SelFile           SFo;
    std::string       comment;
    std::ofstream     fh;
    VolumeXmipp       vol;


    fn_base = fn_root;
    if (iter >= 0)
    {
        fn_base += "_it";
        fn_base.compose(fn_base, iter, "");
    }

    // Do backward FFTW to write out real-space maps again
    double * dataOut;
    try
    {
        dataOut = new double[size];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in maximization");
    }
    vol().resize(Xdim,Ydim,Zdim);
    for (int refno = 0; refno < nr_ref; refno++)
    {
        fn_tmp = fn_base + "_ref";
        fn_tmp.compose(fn_tmp, refno + 1, "");
        fn_tmp = fn_tmp + ".vol";
        
        for (int i = 0; i < size; i++)
        {
            dataOut[i] = dataRefs[refno*size + i];
        }
        // Somehow only setpoints within resolution limits
        backfftw.SetPoints(dataOut);
        backfftw.Transform();
        backfftw.Normalize();
        backfftw.GetPoints(MULTIDIM_ARRAY(vol()));
        vol.write(fn_tmp);
        // Fill selfile and docfile
        SFo.insert(fn_tmp, SelLine::ACTIVE);
        fracline(0) = alpha_k[refno];
        //fracline(1) = 1000 * conv[refno]; // Output 1000x the change for precision
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);
    }
    
    // Write out sel & log-file
    fn_tmp = fn_base + ".sel";
    SFo.write(fn_tmp);

    DFl.go_beginning();
    comment = "MLF_tomo: Number of images= " + floatToString(sumw_allrefs);
    comment += " LL= " + floatToString(LL, 15, 10) + " <Pmax/sumP>= " + floatToString(avePmax, 10, 5);
    comment = " -istart " + integerToString(iter + 1);
    //if (anneal > 1.) comment += " -anneal " + floatToString(anneal, 10, 7);
    DFl.insert_comment(comment);
    DFl.insert_comment("columns: model fraction (1); 1000x signal change (3)");
    fn_tmp = fn_base + ".log";
    DFl.write(fn_tmp);

#ifdef DEBUG
    std::cerr<<"done writeOutputFiles"<<std::endl;
#endif
}

