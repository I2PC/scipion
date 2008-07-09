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
    highres = textToFloat(getParameter(argc, argv, "-highres", "0.5"));
    lowres = textToFloat(getParameter(argc, argv, "-lowres", "0.0001"));
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
void Prog_mlf_tomo_prm::produceSideInfo()
{

    headerXmipp head;
    Matrix3D<double> Maux, Faux_real, Faux_imag;
    Matrix3D<std::complex<double> > Faux, Faux2;
    int xdim, ydim, zdim, c, iaux, ifound;
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
    dim3 = (double) (Xdim * Ydim * Zdim);

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
    // Get size right for the resolution limits
    fftw_hsize = forwfftw.GetSize()*(fN[fNdim-1]/2+1)/fN[fNdim-1];
    // Get an array of booleans whether Fourier components are within resolution range
    // TODO: Also exclude double x==0 axis!!
    try
    {
        is_in_range = new bool [fftw_hsize];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in main");
    }


    int zz, yy, xx;
    double sz,sy,sx,res;
    hsize = 0;
    for ( int z=0, ii=0; z<Zdim; z++ ) {
        if ( z > (Zdim - 1)/2 ) zz = Zdim-z;
        else zz = z;
        for ( int y=0; y<Ydim; y++ ) {
            if ( y > (Ydim - 1)/2 ) yy = Ydim-y;
            else yy = y;
            for ( int x=0; x<Xdim/2 + 1; x++,ii++ ) {
                if ( x > (Xdim - 1)/2 ) xx = Xdim-x;
                else xx = x;
                // Think about non-cubic volumes here: resolution in dig.freq.
                sx = (double)xx / (double)Xdim;
                sy = (double)yy / (double)Ydim;
                sz = (double)zz / (double)Zdim;
                res = sqrt(sx*sx+sy*sy+sz*sz);
                if (res <= highres && res >= lowres)
                {
                    is_in_range[ii] = true;
                    hsize++;
                }
                else
                    is_in_range[ii] = false;
            }
        }
    }
    size = 2 * hsize;
    std::cerr<<" lowres-highres hsize= "<<hsize<<" fftw hsize= "<<fftw_hsize<<std::endl;

    /*
    vol().initZeros();
    for(int i=0, ii=0;i<Zdim;i++)
        for(int j=0;j<Ydim;j++)
            for(int k=0;k<Xdim/2 +1;k++, ii++)
            {                                          //x,y
                vol(k,i,j)=dist_fftw[ii];
            }
    vol.write("dist.fft");
    std::cerr<<"written dist.fft"<<std::endl;
    */

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

#ifdef DEBUG
    std::cerr<<"start generateInitialReferences"<<std::endl;
#endif

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

#ifdef DEBUG
    std::cerr<<"done generateInitialReferences"<<std::endl;
#endif
}

void Prog_mlf_tomo_prm::readAndFftwAllReferences(double * dataRefs)
{

#ifdef DEBUG
    std::cerr<<"start readAndFftwAllReferences"<<std::endl;
#endif

    // Read references in memory and calculate their FFTW
    VolumeXmipp vol;

    alpha_k.clear();
    int refno = 0;
    SFr.go_beginning();
    while (!SFr.eof())
    {
        alpha_k.push_back(1./nr_ref);
        vol.read(SFr.NextImg());
        forwfftw.SetPoints(MULTIDIM_ARRAY(vol()));
        forwfftw.Transform();
        forwfftw.Normalize();
        // Get only points within resolution range
        for (int i = 0, ii=0; i< 2*fftw_hsize; i++)
        {
            if (is_in_range[i/2]) 
            {
                dataRefs[refno*size + ii] = forwfftw.fOut[i];
                ii++;
            }
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
    VolumeXmipp vol;
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
            vol.read(fni);
            forwfftw.SetPoints(MULTIDIM_ARRAY(vol()));
            forwfftw.Transform();
            forwfftw.Normalize();
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

    FileName                         fn, fnw;
    Matrix1D< int >                  radial_count, group_count;
    Matrix1D<double>                 radial_mean;
    SelLine                          SL;
    std::ofstream                    fh;
    std::complex<double>             * IMG, * SUM;
    double                           * sum2, * ave;

    // THIS IS NOT CORRECT YET: TAKE MISSING WEDGES INTO ACCOUNT!!!!

    if (verb > 0)
        std::cerr << "  Estimating initial noise spectra ... " << std::endl;

    // Initialize sigma_noise, sum and sum2
    SUM  = new std::complex<double> [hsize*nr_group];
    sum2 = new double [hsize*nr_group];
    ave = new double [fftw_hsize];
    for (int i=0; i<hsize*nr_group; i++)
    {
        SUM[i] = 0.;
        sum2[i] = 0.;
    }
    group_count.resize(nr_group);

    SFi.go_beginning();
    while (!SFi.eof())
    {
        SL = SFi.current();
        int igg = SL.get_number() - 1;
        group_count(igg)++;
        fn = SFi.NextImg();
        fn += ".fftw";
        backfftw.read(fn);
        IMG =  (std::complex<double> *) backfftw.fIn;
        for (int i=0, ii=0; i < fftw_hsize; i++)
        {
            if (is_in_range[i])
            {
                SUM[igg*hsize+ii] += IMG[i];
                sum2[igg*hsize+ii] += abs(IMG[i]) * abs(IMG[i]);
                ii++;
            }
        } 
    }

    // Subtract squared amplitudes of the average subtomogram 
    // to prevent overestimated noise at low resolutions
    for (int ig = 0; ig < nr_group; ig++)
    {
        for (int i=0; i<hsize; i++)
        {
            sum2[ig*hsize+i] /= (double)group_count(ig);
            SUM[ig*hsize+i] /= (double)group_count(ig);
            sum2[ig*hsize+i] -= abs(SUM[ig*hsize+i]) * abs(SUM[ig*hsize+i]);
        }
        for (int i=0, ii=0; i < fftw_hsize; i++)
        {
            if (is_in_range[i])
            {
                ave[i] = sum2[ig*hsize+ii];
                ii++;
            }
            else
            {
                ave[i] = 0;
            }
        }
        forwfftw.fftwRadialAverage(ave, radial_mean, radial_count, true, true);
        sigma_noise.push_back(radial_mean);

        // Now store in dataSigma structure
        for (int i = 0, ii=0; i< fftw_hsize; i++)
            if (is_in_range[i]) 
            {
                dataSigma[ig * hsize + ii] = 2. * ave[i]; // Store already TWO SIGMA^2
                ii++;
            }

        // Write sigma spectra to disc
        fn = fn_root + "_it";
        fn.compose(fn, istart - 1, "");
        if (nr_group > 1) 
        {
            fn.compose(fn+"_gr", ig + 1, "");
        }
        fn += ".noise";
        fh.open((fn).c_str(), std::ios::out);
        if (!fh) REPORT_ERROR(1, (std::string)"Error: Cannot write file: " + fn);
        for (int irr = 0; irr < XSIZE(sigma_noise[ig]); irr++)
        {
            fh << (double)irr/Xdim << " " << dVi(sigma_noise[ig], irr) << "\n";
        }
        fh.close();

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

    std::complex<double> *DATAREFS, *DATAWSUMREFS, *DATAIMG;
    DATAIMG      = (std::complex<double> *) dataImg;
    DATAREFS     = (std::complex<double> *) dataRefs;
    DATAWSUMREFS = (std::complex<double> *) dataWsumRefs;

    // TODO: Adapt dataSigma arrays to take into account that (x==0) is double in the fftw!
    for (int refno = 0; refno < nr_ref; refno++)
    {
        diff2 = 0.;
        for (int i = 0; i < hsize; i++)
        {
            int iiref = refno * hsize + i;
            int iig = igroup * hsize + i;
            if (dataWedge[i] > 0. && dataSigma[iig] > 0.)
            {
                aux = abs(DATAIMG[i] - DATAREFS[iiref]);
                diff2 += (aux * aux)/dataSigma[iig];
//#define DEBUG_EXPSINGLE 
#ifdef DEBUG_EXPSINGLE                
                std::cerr<<i<<" abs2= "<<aux<<" "<<DATAIMG[i]<<" "<<DATAREFS[iiref]<<" sigma2= "<<dataSigma[iig]<<" diff2= "<<diff2<<" mindiff2= "<<mindiff2<<std::endl;
#endif
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
    }


    // Store Pmax/sumP
    Pmax = maxweight / sumweight;
    std::cerr<<"Pmax= "<<Pmax<<" mindiff2= "<<mindiff2<<std::endl;

    // Then, store the sum of all weights
    for (int refno = 0; refno < nr_ref; refno++)
    {
        weight = weights[refno] / sumweight;
        std::cerr<<" refno= "<<refno<<" w= "<< weight<<" ";

        if (weight > SIGNIFICANT_WEIGHT_LOW)
        {
            dataSumWRefs[refno] += weight;
            for (int i = 0; i < hsize; i++)
            {
                int iiref = refno * hsize + i;
                int iig = igroup * hsize + i;
                if (dataWedge[i] > 0.)
                {
                    aux = abs(DATAIMG[i] - DATAREFS[iiref]);
                    dataWsumDist[iig] += weight * aux * aux;
                    dataWsumWedsPerGroup[iig] += weight;
                    DATAWSUMREFS[iiref] += weight * DATAIMG[i];
                    dataWsumWedsPerRef[iiref] += weight;
                }
#ifdef DEBUG_EXPSINGLE                
                if (i==7425)
                    std::cerr<<" dataImg[320]= "<<DATAIMG[i]
                             <<" dataRefs[320]= "<<DATAREFS[iiref]
                             <<" dataWsumDist[iig]= "<<dataWsumDist[iig]
                             <<" dataWsumWedsPerGroup[iig]= "<<dataWsumWedsPerGroup[iig]
                             <<" dataWsumRefs[iiref]= "<<DATAWSUMREFS[iiref]
                             <<" dataWsumWedsPerRef[iiref]= "<<dataWsumWedsPerRef[iiref]
                             <<std::endl;
#endif
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
                                    double &LL,
                                    double &avePmax)
{

#ifdef DEBUG
    std::cerr<<"start expectation"<<std::endl;
#endif
    int opt_refno, nn, imgno, nr_imgs, igroup, idum;
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
        dataWedge = new double[hsize];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in expectation");
    }
    // FOR NOW ALL WEGDES ARE 1:
    for (int i = 0; i < hsize; i++)
        dataWedge[i] = 1.;

    // Initialize weighted sums to zero
    for (int i = 0; i < nr_ref * size; i++)
    {
        dataWsumRefs[i] = 0.;
    }        
    for (int i = 0; i < nr_ref * hsize; i++)
    {
        dataWsumWedsPerRef[i] = 0.;
    }        
    for (int i = 0; i < nr_group * hsize; i++)
    {
        dataWsumWedsPerGroup[i] = 0.;
        dataWsumDist[i] = 0.;
    }
    for (int i = 0; i < nr_ref; i++)
    {
        dataSumWRefs[i] = 0.;
    }
    LL = 0.;
    avePmax = 0.;

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
        // Get only points within resolution range
        for (int i = 0, ii=0; i< 2*fftw_hsize; i++)
            if (is_in_range[i/2]) 
            {
                dataImg[ii] = backfftw.fIn[i];
                ii++;
            }

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
        avePmax += Pmax;

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
                                     double &sumw_allrefs,
                                     double &avePmax)
{

#ifdef DEBUG
    std::cerr<<"start maximization"<<std::endl;
#endif

    // Update References
    std::complex<double> *DATAREFS, *DATAWSUMREFS;
    DATAREFS     = (std::complex<double> *) dataRefs;
    DATAWSUMREFS = (std::complex<double> *) dataWsumRefs;
    sumw_allrefs = 0;
    for (int refno = 0;refno < nr_ref; refno++)
    {
        sumw_allrefs += dataSumWRefs[refno];
        if (dataSumWRefs[refno] > 0.)
        {
            for (int i = 0; i < hsize; i++)
            {
                int ii = refno * hsize + i;
                // Impute old reference for missing pixels
                DATAREFS[ii] *= 1. - (dataWsumWedsPerRef[ii] / dataSumWRefs[refno]);
                // And sum the weighted sum for observed pixels
                DATAREFS[ii] += DATAWSUMREFS[ii] / dataSumWRefs[refno];
            }
        }
        else
        {
            for (int i = 0; i < hsize; i++)
            {
                int ii = refno * hsize + i;
                DATAREFS[ii] = 0.;
            }
        }
    }

    // Average Pmax
    avePmax /= sumw_allrefs;

    // Update fractions
    if (!fix_fractions)
    {
        for (int refno = 0; refno < nr_ref; refno++)
            alpha_k[refno] = dataSumWRefs[refno] / sumw_allrefs;
    }


    // Update sigma of the noise
    Matrix1D<int> radial_count;
    double *ave;
    ave = new double [fftw_hsize];
    
    for (int ig = 0; ig < nr_group; ig++)
    {
        for (int i = 0; i < hsize; i++)
        {
            int ii = ig * hsize + i;
            // Return from two*SIGMA^2
            dataSigma[ii] /= 2;
            // Impute old sigma values for missing pixels
            dataSigma[ii] *= 1. - (dataWsumWedsPerGroup[ii] / sumw_allrefs);
            //  And sum the weighted sum for observedpixels
            dataSigma[ii] += dataWsumDist[ii] / sumw_allrefs;
        }

        // Set points within resolution range
        for (int i=0, ii=0; i<fftw_hsize; i++)
        {
            if (is_in_range[i])
            {
                ave[i] = dataSigma[ig * hsize + ii];
                ii++;
            }
            else
            {
                ave[i] = 0.;
            }
        }
        forwfftw.fftwRadialAverage(ave, sigma_noise[ig], radial_count, true, true);
        // Now store again in dataSigma structure
        for (int i = 0, ii=0; i< fftw_hsize; i++)
            if (is_in_range[i]) 
            {
                dataSigma[ig * hsize + ii] = 2. * ave[i]; // Store again TWO SIGMA^2
                ii++;
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
        dataOut = new double[2*fftw_hsize];
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
        // Set points within resolution range
        for (int i = 0, ii=0; i< 2*fftw_hsize; i++)
        {
            if (is_in_range[i/2]) 
            {
                dataOut[i] = dataRefs[refno*size + ii];
                ii++;
            }
            else
            {
                dataOut[i] = 0.;
            }
        }
        // Somehow only setpoints within resolution limits
        backfftw.SetPoints(dataOut);
        backfftw.Transform();
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
    //if (anneal > 1.) comment += " -anneal " + floatToString(anneal, 10, 7);
    DFl.insert_comment(comment);
    DFl.insert_comment("columns: model fraction (1); 1000x signal change (3)");
    fn_tmp = fn_base + ".log";
    DFl.write(fn_tmp);

    // Write out updated sigma2 vectors
    if (!fix_sigma_noise)
    {
        for (int ig = 0; ig < nr_group; ig++)
        {
	    fn_tmp = fn_base;
            if (nr_group > 1) 
            {
                fn_tmp.compose(fn_tmp+"_gr", ig + 1, "");
            }
            fn_tmp += ".noise";
            fh.open((fn_tmp).c_str(), std::ios::out);
            if (!fh) REPORT_ERROR(1, (std::string)"Error: Cannot write file: " + fn_tmp);
            for (int irr = 0; irr < XSIZE(sigma_noise[ig]); irr++)
            {
                fh << (double)irr/Xdim << " " << dVi(sigma_noise[ig], irr) << "\n";
            }
            fh.close();
	}
    }

#ifdef DEBUG
    std::cerr<<"done writeOutputFiles"<<std::endl;
#endif
}

