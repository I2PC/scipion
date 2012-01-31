/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.csic.es (2004)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "angular_projection_matching.h"

//#define DEBUG
//#define TIMING

// For blocking of threads
pthread_mutex_t update_refs_in_memory_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t debug_mutex = PTHREAD_MUTEX_INITIALIZER;


// Read arguments ==========================================================
void ProgAngularProjectionMatching::readParams()

{
    fn_exp  = getParam("-i");
    fn_out  = getParam("-o");
    fn_ref  = getParam("--ref");

    // Additional commands
    pad=XMIPP_MAX(1.,getDoubleParam("--pad"));
    Ri=getIntParam("--Ri");
    Ro=getIntParam("--Ro");
    search5d_shift  = getIntParam("--search5d_shift");
    search5d_step = getIntParam("--search5d_step");
    max_shift = getDoubleParam("--max_shift");
    avail_memory = getDoubleParam("--mem");
    if (checkParam("--ctf"))
        fn_ctf  = getParam("--ctf");
    phase_flipped = checkParam("--phase_flipped");
    threads = getIntParam("--thr");

    do_scale = checkParam("--scale");
    if (checkParam("--append"))
        do_overwrite = MD_APPEND;
    else
        do_overwrite = MD_OVERWRITE;

    if(do_scale)
    {
        scale_step = getDoubleParam("--scale",0);
        scale_nsteps = getDoubleParam("--scale",1);
    }

}

void ProgAngularProjectionMatching::defineParams()
{
    addUsageLine("Perform a discrete angular assignment using projection matching in real space.");
    addUsageLine("This program is relatively fast, using polar coordinates for the in-plane ");
    addUsageLine("angular searches and the 5-dimensional search of rotation angles and origin ");
    addUsageLine("offsets is broken in two: first the angles are search in a 3D-search; then, ");
    addUsageLine("for the optimal orientation the origin offsets are searched (2D).");
    addUsageLine(" ");
    addUsageLine("The output of the program consists of a document file with all assigned angles");
    addUsageLine("and rotations. This file also contains a column for the maximum cross-correlation ");
    addUsageLine("coefficient. Note that the program does not alter the image headers. ");
    addUsageLine("The recommended use of this program is within the python script of the ");
    addUsageLine("xmipp_protocol_projmatch.py");
    addSeeAlsoLine("angular_discrete_assign, angular_continuous_assign, angular_project_library");
    addExampleLine("Example of use: Sample at 2 pixel step size for 5D shift search",false);
    addExampleLine("xmipp_angular_projection_matching -i experimental.doc -o assigned_angles.doc --ref reference.stk --search5d_step 2");
    addParamsLine("   -i <doc_file>                : Docfile with input images");
    addParamsLine("   -o <output_filename>         : Output filename");
    addParamsLine("   -r <stackFile>               : Reference projections");
    addParamsLine("     alias --ref;");
    addParamsLine("  [--search5d_shift <s5dshift=0>]: Search range (in +/- pix) for 5D shift search");
    addParamsLine("  [--search5d_step <s5dstep=2>]  : Step size for 5D shift search (in pix)");
    addParamsLine("  [--Ri <ri=1>]               : Inner radius to limit rotational search");
    addParamsLine("  [--Ro <ro=-1>]              : Outer radius to limit rotational search");
    addParamsLine("                        : ro = -1 -> dim/2-1");
    addParamsLine("  [-s <step=1> <n_steps=3>]    : scale step factor (1 means 0.01 in/de-crements) and number of steps around 1.");
    addParamsLine("                               : with default values: 1 0.01 | 0.02 | 0.03");
    addParamsLine("    alias --scale;");
    addParamsLine("==+Extra parameters==");
    addParamsLine("  [--mem <mem=1>]             : Available memory for reference library (Gb)");
    addParamsLine("  [--max_shift <max_shift=-1>]   : Max. change in origin offset (+/- pixels; neg= no limit)");
    addParamsLine("  [--ctf <filename>]            : CTF to apply to the reference projections, either a");
    addParamsLine("                     : CTF parameter file or a 2D image with the CTF amplitudes");
    addParamsLine("  [--pad <pad=1>]             : Padding factor (for CTF correction only)");
    addParamsLine("  [--phase_flipped]            : Use this if the experimental images have been phase flipped");
    addParamsLine("  [--thr <threads=1>]           : Number of concurrent threads");
    addParamsLine("  [--append]                : Append (versus overwrite) data to the output file");
}

/* Show -------------------------------------------------------------------- */
void ProgAngularProjectionMatching::show()
{
    if (!verbose)
        return;

    std::cout
    << "  Input images            : "<< fn_exp << std::endl
    << "  Output file             : "<< fn_out << std::endl
    ;
    if (Ri>0)
        std::cout << "  Inner radius rot-search : " << Ri<< std::endl;
    if (Ro>0)
        std::cout << "  Outer radius rot-search : " << Ro << std::endl;
    if (max_nr_refs_in_memory<total_nr_refs)
    {
        std::cout << "  Number of references    : " << total_nr_refs << std::endl
        << "  Nr. refs in memory      : " << max_nr_refs_in_memory << " (using " << avail_memory <<" Gb)" << std::endl
        ;
    }
    else
    {
        std::cout << "  Number of references    : " << total_nr_refs << " (all stored in memory)" << std::endl;
    }
    std::cout << "  Max. allowed shift      : +/- " <<max_shift<<" pixels"<<std::endl;
    if (search5d_shift > 0)
    {
        std::cout << "  5D-search shift range   : "<<search5d_shift<<" pixels (sampled "<<nr_trans<<" times)"<<std::endl;
    }
    if (fn_ctf!="")
    {
        if (!fn_ctf.isMetaData())
        {
            std::cout << "  CTF image               :  " <<fn_ctf<<std::endl;
            if (pad > 1.)
                std::cout << "  Padding factor          : "<< pad << std::endl;
        }
        else
        {
            std::cout << "  CTF parameter file      :  " <<fn_ctf<<std::endl;
            if (phase_flipped)
                std::cout << "    + Assuming images have been phase flipped " << std::endl;
            else
                std::cout << "    + Assuming images have not been phase flipped " << std::endl;
        }
    }
    if (threads>1)
    {
        std::cout << "  -> Using "<<threads<<" parallel threads"<<std::endl;
    }
    std::cout << " ================================================================="<<std::endl;
}

/* Run --------------------------------------------------------------------- */
void ProgAngularProjectionMatching::run()
{
    produceSideInfo();

    show();

    processAllImages();

    writeOutputFiles();

    if (verbose)
        std::cout << "done!"<<std::endl;
}


// Side info stuff ===================================================================
void ProgAngularProjectionMatching::produceSideInfo()
{

    Image<double>    img,empty;
    Projection       proj;
    MetaData         DF;
    MetaData         SFr,emptySF;
    SymList          SL;
    FileName         fn_img;
    double           mean,stddev,psi=0.;
    MultidimArray<double> Maux;
    MultidimArray<double> dataline(3);
    int              nl;
    Polar<double>    P;
    Polar<std::complex <double> > fP;

    // Read Selfile and get dimensions
    DFexp.read(fn_exp);

    // Thread barrier
    barrier_init(&thread_barrier, threads);

    // Read one image to get dim
    DFexp.getValue(MDL_IMAGE,fn_img,DFexp.firstObject());
    img.read(fn_img);
    dim = XSIZE(img());

    // Check that the reference and the experimental images are of the same
    // size
    FileName fnt;
    fnt.compose(FIRST_IMAGE, fn_ref);
    Image<double> imgRef;
    imgRef.read(fnt,HEADER);

    if (!imgRef().sameShape(img()))
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "Check that the reference volume and the experimental images are of the same size");

    // Set padding dimension
    paddim=ROUND(pad*dim);

    // Set max_shift
    if (max_shift<0)
        max_shift = dim/2;

    // Set ring defaults
    if (Ri<1)
        Ri=1;
    if (Ro<0)
        Ro=(dim/2)-1;

    // Calculate necessary memory per image
    produceSplineCoefficients(BSPLINE3,Maux,img());
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    P.calculateFftwPlans(global_plans);
    fourierTransformRings(P,fP,global_plans,false);
    double memory_per_ref = 0.;
    for (int i = 0; i < fP.getRingNo(); i++)
    {
        memory_per_ref += (double) fP.getSampleNo(i) * 2 * sizeof(double);
    }
    memory_per_ref += dim * dim * sizeof(double);
    max_nr_imgs_in_memory = ROUND( 1024 * 1024 * 1024 * avail_memory / memory_per_ref);

    // Set up angular sampling
    mysampling.readSamplingFile(fn_ref.removeAllExtensions(),false);
    total_nr_refs = mysampling.no_redundant_sampling_points_angles.size();

    convert_refno_to_stack_position.resize(mysampling.numberSamplesAsymmetricUnit, -1);
    for (int i = 0; i < mysampling.no_redundant_sampling_points_index.size(); i++)
    	convert_refno_to_stack_position[mysampling.no_redundant_sampling_points_index[i]] = i;

    // Don't reserve more memory than necessary
    max_nr_refs_in_memory = XMIPP_MIN(max_nr_imgs_in_memory, total_nr_refs);

    // Initialize pointers for reference retrieval
    pointer_allrefs2refsinmem.resize(mysampling.numberSamplesAsymmetricUnit,-1);
    pointer_refsinmem2allrefs.resize(max_nr_refs_in_memory,-1);
    counter_refs_in_memory = 0;
    loop_forward_refs=true;

    // Initialize 5D search vectors
    search5d_xoff.clear();
    search5d_yoff.clear();
    // Make sure origin is included
    if(search5d_step == 0)
    {
        printf("***********************************************************\n");
        printf("*   ERROR: search step should be different from 0\n");
        printf("*   search step set to 1 \n");
        printf("***********************************************************\n");

        search5d_step = 1;
    }
    int myfinal=search5d_shift + search5d_shift%search5d_step;
    nr_trans = 0;//number translations in 5D
    for (int xoff = -myfinal; xoff <= myfinal; xoff+= search5d_step)
    {
        for (int yoff = -myfinal; yoff <= myfinal; yoff+= search5d_step)
        {
            // Only take a circle (not a square)
            if ( xoff*xoff + yoff*yoff <= search5d_shift*search5d_shift)
            {
                search5d_xoff.push_back(xoff);
                search5d_yoff.push_back(yoff);
                nr_trans++;
            }
        }
    }

    // Initialize all arrays
    try
    {
        fP_ref = new Polar<std::complex<double> >[max_nr_refs_in_memory];
        proj_ref = new MultidimArray<double>[max_nr_refs_in_memory];
        fP_img = new Polar<std::complex<double> >[nr_trans];
        fPm_img = new Polar<std::complex<double> >[nr_trans];

        stddev_ref = new double[max_nr_refs_in_memory];
        stddev_img = new double[nr_trans];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(ERR_MEM_BADREQUEST,"Error allocating memory in produceSideInfo");
    }

    // CTF stuff
    if (fn_ctf != "")
    {
        if (!fn_ctf.isMetaData())
        {
            Image<double> img;
            img.read(fn_ctf);
            Mctf=img();
            if (XSIZE(Mctf) != paddim)
            {
                std::cerr<<"image size= "<<dim<<" padding factor= "<<pad<<" padded image size= "<<paddim<<" Wiener filter size= "<<XSIZE(Mctf)<<std::endl;
                REPORT_ERROR(ERR_VALUE_INCORRECT,
                             "Incompatible padding factor for this CTF filter");
            }
        }
        else
        {
            CTFDescription ctf;
            MultidimArray<std::complex<double> >  ctfmask;
            ctf.read(fn_ctf);
            if (ABS(ctf.DeltafV - ctf.DeltafU) >1.)
            {
                REPORT_ERROR(ERR_VALUE_INCORRECT,
                             "ERROR!! Only non-astigmatic CTFs are allowed!");
            }
            ctf.enable_CTF = true;
            ctf.Produce_Side_Info();
            ctf.Generate_CTF(paddim, paddim, ctfmask);
            Mctf.resize(paddim,paddim);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mctf)
            {
                if (phase_flipped)
                    dAij(Mctf, i, j) = fabs(dAij(ctfmask, i, j).real());
                else
                    dAij(Mctf, i, j) = dAij(ctfmask, i, j).real();
            }
        }
    }

    //Store the id's of each experimental image from metadata
    DFexp.findObjects(ids);
}

int ProgAngularProjectionMatching::getCurrentReference(int refno,
        Polar_fftw_plans &local_plans)
{
    FileName                      fnt;
    Image<double>                 img;
    double                        mean,stddev;
    MultidimArray<double>         Maux;
    Polar<double>                 P;
    Polar<std::complex <double> > fP;
    FourierTransformer                     local_transformer;
    size_t _pointer;

    // Image was not stored yet: read it from disc and store
//    std::vector<size_t>::const_iterator found =
//        std::find((mysampling.no_redundant_sampling_points_index).begin(),
//                  (mysampling.no_redundant_sampling_points_index).end(),
//                  (size_t)refno
//                 );
//    //found = found - (mysampling.no_redundant_sampling_points_index).begin();
//    _pointer = found - (mysampling.no_redundant_sampling_points_index).begin();
//    if (found == (mysampling.no_redundant_sampling_points_index).end())
//        REPORT_ERROR(ERR_VALUE_INCORRECT, "Wrong reference number");
//    fnt.compose(_pointer + FIRST_IMAGE, fn_ref);
    fnt.compose(convert_refno_to_stack_position[refno] + FIRST_IMAGE, fn_ref);
    //!a delete _DATA_ALL
    img.read(fnt, _DATA_ALL);
    double rot_tmp,tilt_tmp,psi_tmp;
    img.getEulerAngles(rot_tmp,tilt_tmp,psi_tmp);
    img().setXmippOrigin();
#define DEBUG
#ifdef DEBUG

    {
        std::cerr << "refno: " << refno<<std::endl;
    	std::cerr << "index_found: " << convert_refno_to_stack_position[refno] << std::endl;
    	std::cerr << "reading image " << fnt << std::endl;
        std::cerr << "rot_tmp,tilt_tmp,psi_tmp: " << rot_tmp<< " "<< tilt_tmp<< " "<<psi_tmp<< std::endl;
        std::cerr << "XXXXno_redundant_sampling_points_indexXXXXXX" <<std::endl;
        for (std::vector<size_t>::iterator i =
                    mysampling.no_redundant_sampling_points_index.begin();
                i != mysampling.no_redundant_sampling_points_index.end();
                ++i)
            std::cerr << *i << " ";
        std::cerr << std::endl;
    }
#endif
    // Apply CTF
    if (fn_ctf!="")
    {
        MultidimArray<std::complex<double> > Faux;
        if (paddim > dim)
        {
            // pad real-space image
            int x0 = FIRST_XMIPP_INDEX(paddim);
            int xF = LAST_XMIPP_INDEX(paddim);
            img().selfWindow(x0, x0, xF,xF);
        }
        local_transformer.FourierTransform(img(),Faux);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Faux)
        {
            dAij(Faux,i,j) *= dAij(Mctf,i,j);
        }
        local_transformer.inverseFourierTransform(Faux,img());

        if (paddim > dim)
        {
            // de-pad real-space image
            int x0 = FIRST_XMIPP_INDEX(dim);
            int xF = LAST_XMIPP_INDEX(dim);
            img().selfWindow(x0, x0, xF,xF);
        }
    }

    // Calculate FTs of polar rings and its stddev
    produceSplineCoefficients(BSPLINE3,Maux,img());
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    P.computeAverageAndStddev(mean,stddev);
    P -= mean;
    fourierTransformRings(P,fP,local_plans,true);

    pthread_mutex_lock(  &update_refs_in_memory_mutex );

    int counter = counter_refs_in_memory % max_nr_refs_in_memory;
    pointer_allrefs2refsinmem[refno] = counter;
    if (pointer_refsinmem2allrefs[counter] != -1)
    {
        // This position was in use already
        // Images will be overwritten, so reset the
        // pointer_allrefs2refsinmem of the old images to -1
        pointer_allrefs2refsinmem[pointer_refsinmem2allrefs[counter]] = -1;
    }
    pointer_refsinmem2allrefs[counter] = refno;
    fP_ref[counter] = fP;
    stddev_ref[counter] = stddev;
    proj_ref[counter] = img();
#define DEBUG
#ifdef DEBUG

    std::cerr<<"counter= "<<counter<<"refno= "<<refno<<" stddev = "<<stddev;
    std::cerr<<" refsinmem2allrefs= "<<pointer_refsinmem2allrefs[counter];
    std::cerr<<" allrefs2refsinmem= "<<pointer_allrefs2refsinmem[pointer_refsinmem2allrefs[counter]] <<std::endl;
    std::cerr << "pointer_allrefs2refsinmem" <<std::endl;
    for (std::vector<int>::iterator i = pointer_allrefs2refsinmem.begin();
            i != pointer_allrefs2refsinmem.end();
            ++i)
        std::cerr << *i << " ";
    std::cerr <<std::endl;
    std::cerr << "pointer_refsinmem2allrefs" <<std::endl;
    for (std::vector<int>::iterator i = pointer_refsinmem2allrefs.begin();
            i != pointer_refsinmem2allrefs.end();
            ++i)
        std::cerr << *i << " ";
    std::cerr <<std::endl;
#endif

    counter_refs_in_memory++;
    pthread_mutex_unlock(  &update_refs_in_memory_mutex );

}

void * threadRotationallyAlignOneImage( void * data )
{
    structThreadRotationallyAlignOneImage * thread_data = (structThreadRotationallyAlignOneImage *) data;

    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage 1" <<std::endl;
    // Variables from above
    int thread_id = thread_data->thread_id;
    int thread_num = thread_data->thread_num;
    ProgAngularProjectionMatching *prm = thread_data->prm;
    MultidimArray<double> *img = thread_data->img;
    size_t *this_image = thread_data->this_image;
    int *opt_refno = thread_data->opt_refno;
    double *opt_psi = thread_data->opt_psi;
    bool *opt_flip = thread_data->opt_flip;
    double *maxcorr = thread_data->maxcorr;

    // Local variables
    MultidimArray<double>       Maux;
    MultidimArray<double>       ang, corr;
    int                         max_index, refno, myinit, myfinal, myincr;
    bool                        done_once=false;
    double                      mean, stddev;
    Polar<double>               P;
    Polar<std::complex <double> > fP,fPm;
    RotationalCorrelationAux    rotAux;
    Polar_fftw_plans            local_plans;
    size_t                         imgno = (*this_image) - FIRST_IMAGE;

#ifdef TIMING

    TimeStamp t0,t1,t2;
    time_config();
    annotate_time(&t0);
    annotate_time(&t2);
#endif

    *maxcorr = -99.e99;
    produceSplineCoefficients(BSPLINE3,Maux,*img);
    // Precalculate polar transform of each translation
    // This loop is also threaded
    myinit = thread_id;
    myincr = thread_num;
    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage 2 itrans LOOP" <<std::endl;
    for (int itrans = myinit; itrans < prm->nr_trans; itrans+=myincr)
    {
        //std::cerr << "DEBUG_JM:      itrans: " <<      itrans << std::endl;
        P.getPolarFromCartesianBSpline(Maux,prm->Ri,prm->Ro,3,
                                       (double)prm->search5d_xoff[itrans],
                                       (double)prm->search5d_yoff[itrans]);
        P.computeAverageAndStddev(mean,stddev);
        P -= mean; // for normalized cross-correlation coefficient
        if (itrans == myinit)
            P.calculateFftwPlans(local_plans);
        fourierTransformRings(P,prm->fP_img[itrans],local_plans,false);
        fourierTransformRings(P,prm->fPm_img[itrans],local_plans,true);
        prm->stddev_img[itrans] = stddev;
        done_once=true;
    }
    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage 3 after LOOP" <<std::endl;
    // If thread did not have to do any itrans, initialize fftw plans
    if (!done_once)
    {
        P.getPolarFromCartesianBSpline(Maux,prm->Ri,prm->Ro);
        P.calculateFftwPlans(local_plans);
    }
    // Prepare FFTW plan for rotational correlation
    corr.resize(P.getSampleNoOuterRing());
    rotAux.local_transformer.setReal(corr);
    rotAux.local_transformer.FourierTransform();

    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage 4 before thread WAIT" <<std::endl;
    // All threads have to wait until the itrans loop is done
    barrier_wait(&(prm->thread_barrier));

#ifdef TIMING

    float prepare_img = elapsed_time(t0);
    float get_refs = 0.;
    annotate_time(&t0);
#endif

    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage 5 after WAIT" <<std::endl;
    //pthread_mutex_lock(  &debug_mutex );
    // Switch the order of looping through the references every time.
    // That way, in case max_nr_refs_in_memory<total_nr_refs
    // the references read in memory for the previous image
    // will still be there when processing the next image
    // Every thread processes a part of the references
    if (prm->loop_forward_refs)
    {
        myinit = 0;
        myfinal = prm->mysampling.my_neighbors[imgno].size();
        myincr = +1;
    }
    else
    {
        myinit = prm->mysampling.my_neighbors[imgno].size() - 1;
        myfinal = -1;
        myincr = -1;
    }
    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage 6 i LOOP" <<std::endl;
    // Loop over all relevant "neighbours" (i.e. directions within the search range)
    //for (int i = myinit; i != myfinal; i+=myincr)
    for (size_t i = myinit; i != myfinal; i += myincr)
    {
        std::cerr << "DEBUG_JM:     i: " <<     i << std::endl;
        if (i%thread_num == thread_id)
        {

#ifdef DEBUG_THREADS
            pthread_mutex_lock(  &debug_mutex );
            std::cerr<<" thread_id= "<<thread_id<<" i= "<<i<<" "<<myinit<<" "<<myfinal<<" "<<myincr<<std::endl;
            pthread_mutex_unlock(  &debug_mutex );
#endif

#ifdef TIMING

            annotate_time(&t1);
#endif
            // Get pointer to the current reference image
#ifdef DEBUG

            if(prm->mysampling.my_neighbors[imgno][i]==58)
            {
                std::cerr << "XXXXpointer_allrefs2refsinmemXXXXXX" <<std::endl;
                for (std::vector<int>::iterator i = prm->
                                                    pointer_allrefs2refsinmem.begin();
                        i != prm->pointer_allrefs2refsinmem.end();
                        ++i)
                    std::cerr << *i << std::endl;
                std::cerr << "XXXXpointer_refsinmem2allrefsXXXXXX" <<std::endl;
                for (std::vector<int>
                        ::iterator i = prm->pointer_refsinmem2allrefs.begin();
                        i != prm->pointer_refsinmem2allrefs.end();
                        ++i)
                    std::cerr << *i << std::endl;
                std::cerr <<std::endl;

            }
#endif

            refno = prm->pointer_allrefs2refsinmem[prm->mysampling.my_neighbors[imgno][i]];
            if (refno == -1)
            {
                // Reference is not stored in memory (anymore): (re-)read from disc
                prm->getCurrentReference(prm->mysampling.my_neighbors[imgno][i],local_plans);
                refno = prm->pointer_allrefs2refsinmem[prm->mysampling.my_neighbors[imgno][i]];
            }


#ifdef TIMING
            get_refs += elapsed_time(t1);
#endif
            #define DEBUG
#ifdef DEBUG

            std::cerr << "imgno " << imgno <<std::endl;
            std::cerr<<"Got refno= "<<refno
            <<" pointer= "<<prm->mysampling.my_neighbors[imgno][i]<<std::endl;
#endif

            // Loop over all 5D-search translations
            for (int itrans = 0; itrans < prm->nr_trans; itrans++)
            {
#ifdef DEBUG

                std::cerr<< "prm->stddev_ref[refno], prm->stddev_img[itrans]: " <<
                prm->stddev_ref[refno] << " " <<
                prm->stddev_img[itrans];
#endif
                // A. Check straight image
                std::cerr << "threadRotationallyAlignOneImage.ONE" <<std::endl;
                std::cerr << "itrans, refno " << itrans << " " << refno
                <<std::endl;
		rotationalCorrelation(prm->fP_img[itrans],
				      prm->fP_ref[refno],
				      ang,rotAux);
                std::cerr << "threadRotationallyAlignOneImage.TWO" <<std::endl;
                corr /= prm->stddev_ref[refno] * prm->stddev_img[itrans]; // for normalized ccf
                for (int k = 0; k < XSIZE(corr); k++)
                {
                    if (DIRECT_A1D_ELEM(corr,k)> *maxcorr)
                    {
                        *maxcorr = DIRECT_A1D_ELEM(corr,k);
                        *opt_psi = DIRECT_A1D_ELEM(ang,k);
                        *opt_refno = prm->mysampling.my_neighbors[imgno][i];
                        *opt_flip = false;
                    }
                }
#ifdef DEBUG
                std::cerr<<"straight: corr "<<*maxcorr<<std::endl;
#endif
                // B. Check mirrored image
                rotationalCorrelation(prm->fPm_img[itrans],prm->fP_ref[refno],ang,rotAux);
                corr /= prm->stddev_ref[refno] * prm->stddev_img[itrans]; // for normalized ccf
                for (int k = 0; k < XSIZE(corr); k++)
                {
                    if (DIRECT_A1D_ELEM(corr,k)> *maxcorr)
                    {
                        *maxcorr = DIRECT_A1D_ELEM(corr,k);
                        *opt_psi = DIRECT_A1D_ELEM(ang,k);
                        *opt_refno = prm->mysampling.my_neighbors[imgno][i];
                        *opt_flip = true;
                    }
                }

#ifdef DEBUG
                std::cerr<<"mirror: corr "<<*maxcorr;
                if (*opt_flip)
                    std::cerr<<"**";
                std::cerr<<std::endl;
#endif
#undef DEBUG

            }
        }
    }

#ifdef TIMING
    float all_rot_align = elapsed_time(t0);
    float total_rot = elapsed_time(t2);
    std::cerr<<" rotal%% "<<total_rot
    <<" => prep: "<<prepare_img
    <<" all_refs: "<<all_rot_align
    <<" (of which "<<get_refs
    <<" to get "<< prm->mysampling.my_neighbors[imgno].size()
    <<" refs for imgno "<<imgno<<" )"
    <<std::endl;
#endif
    //pthread_mutex_unlock(  &debug_mutex );
    //std::cerr << "DEBUG_JM: threadRotationallyAlignOneImage END" <<std::endl;
}

void ProgAngularProjectionMatching::translationallyAlignOneImage(MultidimArray<double> &img,
        const int &opt_refno,
        const double &opt_psi,
        const bool &opt_flip,
        double &opt_xoff,
        double &opt_yoff,
        double &maxcorr)
{
    MultidimArray<double> Mtrans,Mimg,Mref;
    int refno;
    Mtrans.setXmippOrigin();
    Mimg.setXmippOrigin();
    Mref.setXmippOrigin();

#ifdef TIMING

    TimeStamp t0,t1,t2;
    time_config();
    annotate_time(&t0);
#endif

#ifdef DEBUG

    std::cerr<<"start trans: opt_refno= "<<opt_refno<<" pointer= "<<pointer_allrefs2refsinmem[opt_refno]<<" opt_psi= "<<opt_psi<<"opt_flip= "<<opt_flip<<std::endl;
#endif

    // Get pointer to the correct reference image in memory
    refno = pointer_allrefs2refsinmem[opt_refno];
    if (refno == -1)
    {
        // Reference is not stored in memory (anymore): (re-)read from disc
        getCurrentReference(opt_refno,global_plans);
        refno = pointer_allrefs2refsinmem[opt_refno];
    }

    // Rotate stored reference projection by phi degrees
    rotate(BSPLINE3,Mref,proj_ref[refno],opt_psi,DONT_WRAP);

#ifdef DEBUG

    std::cerr<<"rotated ref "<<std::endl;
#endif

    if (opt_flip)
    {
        // Flip experimental image
        Matrix2D<double> A(3,3);
        A.initIdentity();
        MAT_ELEM(A,0, 0) = -1;
        applyGeometry(LINEAR, Mimg, img, A, IS_INV, DONT_WRAP);
    }
    else
        Mimg = img;

    // Perform the actual search for the optimal shift
    if (max_shift>0)
    {
    	CorrelationAux aux;
        bestShift(Mref,Mimg,opt_xoff,opt_yoff,aux);
    }
    else
        opt_xoff = opt_yoff = 0.;
    if (opt_xoff * opt_xoff + opt_yoff * opt_yoff > max_shift * max_shift)
        opt_xoff = opt_yoff = 0.;

#ifdef DEBUG

    std::cerr<<"optimal shift "<<opt_xoff<<" "<<opt_yoff<<std::endl;
#endif

    // Calculate standard cross-correlation coefficient
    translate(LINEAR,Mtrans,Mimg,vectorR2(opt_xoff,opt_yoff),true);
    maxcorr = correlationIndex(Mref,Mtrans);

#ifdef DEBUG

    std::cerr<<"optimal shift corr "<<maxcorr<<std::endl;
#endif

    // Correct X-shift for mirrored images
    if (opt_flip)
        opt_xoff *= -1.;

#ifdef TIMING

    float total_trans = elapsed_time(t0);
    std::cerr<<" trans%% "<<total_trans <<std::endl;
#endif

}


void ProgAngularProjectionMatching::scaleAlignOneImage(MultidimArray<double> &img,
        const int &opt_refno,
        const double &opt_psi,
        const bool &opt_flip,
        const double &opt_xoff,
        const double &opt_yoff,
        double &opt_scale,
        double &maxcorr)
{
    MultidimArray<double> Mscale,Mtrans,Mref;
    int refno;

    Mscale.setXmippOrigin();
    Mtrans.setXmippOrigin();
    Mref.setXmippOrigin();

#ifdef TIMING

    TimeStamp t0,t1,t2;
    time_config();
    annotate_time(&t0);
#endif

    // ROTATE + SHIFT
    // Transformation matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();
    double ang, cosine, sine;
    ang = DEG2RAD(-opt_psi);
    cosine = cos(ang);
    sine = sin(ang);

    // Rotation
    MAT_ELEM(A,0, 0) = cosine;
    MAT_ELEM(A,0, 1) = -sine;
    MAT_ELEM(A,1, 0) = sine;
    MAT_ELEM(A,1, 1) = cosine;

    // Shift
    MAT_ELEM(A,0, 2) = -opt_xoff;
    MAT_ELEM(A,1, 2) = -opt_yoff;

    if (opt_flip)
    {
        MAT_ELEM(A,0, 0) *= -1.;
        MAT_ELEM(A,0, 1) *= -1.;
    }

    applyGeometry(LINEAR, Mtrans, img, A, IS_INV, DONT_WRAP);

    // SCALE
    // Get pointer to the correct reference image in memory
    refno = pointer_allrefs2refsinmem[opt_refno];
    if (refno == -1)
    {
        // Reference is not stored in memory (anymore): (re-)read from disc
        getCurrentReference(opt_refno,global_plans);
        refno = pointer_allrefs2refsinmem[opt_refno];
    }
    Mref = proj_ref[refno];

    // Means
    double ref_mean=Mref.computeAvg();
    Mref-=ref_mean;
    // Mtrans is already rotated and shifted
    double trans_mean=Mtrans.computeAvg();

    // Scale search
    double corr;
    opt_scale = 1;
    maxcorr = 0;

    // 1 (0.01 * scale_step * scale_nsteps)
    for(double scale = 1 - 0.01 * scale_step * scale_nsteps ;
            scale <= 1 + 0.01 * scale_step * (scale_nsteps + 1) ;
            scale += 0.01 * scale_step)
    {
        // apply current scale
        A.initIdentity();
        A *= scale;
        applyGeometry(LINEAR, Mscale, Mtrans, A, IS_INV, DONT_WRAP);

        // SUM (Mscale - trans_mean) * (Mref - ref_mean)
        corr=0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mscale)
        corr+=(DIRECT_MULTIDIM_ELEM(Mscale,n)-trans_mean)*DIRECT_MULTIDIM_ELEM(Mref,n);

        // best scale update
        if(corr > maxcorr)
        {
            opt_scale = scale;
            maxcorr = corr;
        }
    }

#ifdef TIMING

    float total_trans = elapsed_time(t0);
    std::cerr<<" trans%% "<<total_trans <<std::endl;
#endif

}

void ProgAngularProjectionMatching::processAllImages()
{
    //Init progress bar
    size_t total_number_of_images = DFexp.size();
    if (verbose)
    {
        progress_bar_step = XMIPP_MAX(1, total_number_of_images / 80);
        init_progress_bar(total_number_of_images);
    }
    processSomeImages(ids);
    if (verbose)
        progress_bar(total_number_of_images);
}

void ProgAngularProjectionMatching::processSomeImages(const std::vector<size_t> &imagesToProcess)
{
    Image<double> img;
    double opt_rot, opt_tilt, opt_psi, opt_xoff, opt_yoff, opt_scale, maxcorr=-99.e99;
    bool opt_flip;
    int opt_refno;

    size_t nr_images = imagesToProcess.size();
    size_t idNew, imgid;
    FileName fn;

    for (size_t imgno = 0; imgno < nr_images; imgno++)
    {
        imgid = imagesToProcess[imgno];
        //std::cerr << "DEBUG_JM: imgno: " << imgno << std::endl;
        //std::cerr << "DEBUG_JM: imgid: " << imgid << std::endl;
        //std::cerr << "DEBUG_JM: calling getCurrentImage" <<std::endl;
        getCurrentImage(imgid, img);
        //std::cerr << "DEBUG_JM:     after getCurrentImage" <<std::endl;
        //img.write("kk,spi");
        //exit(0);
        // Call threads to calculate the rotational alignment of each image in the selfile
        pthread_t * th_ids = (pthread_t *)malloc( threads * sizeof( pthread_t));

        structThreadRotationallyAlignOneImage * threads_d = (structThreadRotationallyAlignOneImage *)
                malloc ( threads * sizeof( structThreadRotationallyAlignOneImage ) );
        //std::cerr << "DEBUG_JM: creating threads..." <<std::endl;
        for( int c = 0 ; c < threads ; c++ )
        {
            threads_d[c].thread_id = c;
            threads_d[c].thread_num = threads;
            threads_d[c].prm = this;
            threads_d[c].img=&img();
            threads_d[c].this_image=&imgid;
            threads_d[c].opt_refno=&opt_refno;
            threads_d[c].opt_psi=&opt_psi;
            threads_d[c].opt_flip=&opt_flip;
            threads_d[c].maxcorr=&maxcorr;
            pthread_create( (th_ids+c), NULL, threadRotationallyAlignOneImage, (void *)(threads_d+c) );
        }
        //std::cerr << "DEBUG_JM: joining threads..." <<std::endl;
        // Wait for threads to finish and get optimal refno, psi, flip and maxcorr
        for( int c = 0 ; c < threads ; c++ )
        {
            if (*threads_d[c].maxcorr > maxcorr)
            {
                maxcorr = *threads_d[c].maxcorr;
                opt_refno = *threads_d[c].opt_refno;
                opt_psi = *threads_d[c].opt_psi;
                opt_flip = *threads_d[c].opt_flip;
            }
            pthread_join(*(th_ids+c),NULL);
        }
        //exit(1);      // add one because first number is number of elements in the array

        // Flip order to loop through references
        loop_forward_refs = !loop_forward_refs;

        std::cerr << "opt_refno: " << opt_refno <<std::endl;
        opt_rot  = XX(mysampling.no_redundant_sampling_points_angles[convert_refno_to_stack_position[opt_refno]]);
        opt_tilt = YY(mysampling.no_redundant_sampling_points_angles[convert_refno_to_stack_position[opt_refno]]);

        //std::cerr << "DEBUG_JM: calling translationallyAlignOneImage" <<std::endl;
        translationallyAlignOneImage(img(), opt_refno, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr);
        //std::cerr << "DEBUG_JM:     after translationallyAlignOneImage" <<std::endl;
        // Add previously applied translation to the newly found one
        opt_xoff += img.Xoff();
        opt_yoff += img.Yoff();

        opt_scale=1.0;

        if(do_scale)
        {
            //std::cerr << "DEBUG_JM: calling scaleAlignOneImage" <<std::endl;
            // Compute a better scale (scale_min -> scale_max)
            scaleAlignOneImage(img(), opt_refno, opt_psi, opt_flip, opt_xoff, opt_yoff, opt_scale, maxcorr);
            //std::cerr << "DEBUG_JM:     after scaleAlignOneImage" <<std::endl;
            //Add the previously applied scale to the newly found one
            opt_scale *= img.scale();
        }

        // Output
        DFexp.getValue(MDL_IMAGE, fn, imgid);

        idNew = DFo.addObject();
        DFo.setValue(MDL_IMAGE, fn,idNew);
        DFo.setValue(MDL_ANGLEROT, opt_rot,idNew);
        DFo.setValue(MDL_ANGLETILT,opt_tilt,idNew);
        DFo.setValue(MDL_ANGLEPSI, opt_psi,idNew);
        DFo.setValue(MDL_SHIFTX,   opt_xoff,idNew);
        DFo.setValue(MDL_SHIFTY,   opt_yoff,idNew);
        DFo.setValue(MDL_REF,      opt_refno + FIRST_IMAGE,idNew);
        DFo.setValue(MDL_FLIP,     opt_flip,idNew);
        DFo.setValue(MDL_SCALE,    opt_scale,idNew);
        DFo.setValue(MDL_MAXCC,    maxcorr,idNew);
        if (verbose && imgno % progress_bar_step == 0)
            progress_bar(imgno);

        //std::cerr << "DEBUG_JM: END OF ITERATION imgno: " << imgno <<std::endl;
    }
}

void ProgAngularProjectionMatching::getCurrentImage(size_t imgid, Image<double> &img)
{
    FileName fn_img;
    Matrix2D<double> A;

    // jump to line imgno+1 in DFexp, get data and filename
    DFexp.getValue(MDL_IMAGE,fn_img, imgid);

    // Read actual image
    img.read(fn_img);
    img().setXmippOrigin();

    // Store translation in header and apply it to the actual image
    double shiftX, shiftY;
    DFexp.getValue(MDL_SHIFTX,shiftX, imgid);
    DFexp.getValue(MDL_SHIFTY,shiftY, imgid);
    img.setShifts(shiftX,shiftY);
    img.setEulerAngles(0.,0.,0.);
    img.setFlip(0.);

    double scale;
    scale = 1.;
    if(DFexp.containsLabel(MDL_SCALE))
        DFexp.getValue(MDL_SCALE, scale, imgid);
    img.setScale(scale);

    img.getTransformationMatrix(A,true);
    if (!A.isIdentity())
        selfApplyGeometry(BSPLINE3, img(), A, IS_INV, WRAP);
}

void ProgAngularProjectionMatching::writeOutputFiles()
{
    DFo.write(fn_out,do_overwrite);
}
