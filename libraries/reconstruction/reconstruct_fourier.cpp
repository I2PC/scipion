/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Jose Roman Bilbao-Castro (jrbcast@ace.ual.es)
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

#include "reconstruct_fourier.h"

// Define params
void ProgRecFourier::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using direct Fourier interpolation with arbitrary geometry.");
    addUsageLine("Kaisser-windows are used for interpolation in Fourier space.");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  [-o <volume_file=\"rec_fourier.vol\">]  : Filename for output volume");
    addParamsLine("  [--iter <iterations=1>]      : Number of iterations for weight correction");
    addParamsLine("  [--sym <symfile=c1>]              : Enforce symmetry in projections");
    addParamsLine("  [--padding <proj=2.0> <vol=2.0>]  : Padding used for projections and volume");
    addParamsLine("  [--prepare_fsc <fscfile>]      : Filename root for FSC files");
    addParamsLine("  [--max_resolution <p=0.5>]     : Max resolution (Nyquist=0.5)");
    addParamsLine("  [--weight]                     : Use weights stored in the image metadata");
    addParamsLine("  [--thr <threads=1> <rows=1>]   : Number of concurrent threads and rows processed at time by a thread");
    addParamsLine("  [--blob <radius=1.9> <order=0> <alpha=15>] : Blob parameters");
    addParamsLine( "                                  : radius in pixels, order of Bessel function in blob and parameter alpha");
    addExampleLine("For reconstruct enforcing i3 symmetry and using stored weights:", false);
    addExampleLine("   xmipp_reconstruct_fourier  -i reconstruction.sel --sym i3 --weight");
}

// Read arguments ==========================================================
void ProgRecFourier::readParams()
{
    fn_sel = getParam("-i");
    fn_out = getParam("-o");
    fn_sym = getParam("--sym");
    if(checkParam("--prepare_fsc"))
        fn_fsc = getParam("--prepare_fsc");
    do_weights = checkParam("--weight");
    padding_factor_proj = getDoubleParam("--padding", 0);
    padding_factor_vol = getDoubleParam("--padding", 1);
    blob.radius   = getDoubleParam("--blob", 0);
    blob.order    = getIntParam("--blob", 1);
    blob.alpha    = getDoubleParam("--blob", 2);
    maxResolution = getDoubleParam("--max_resolution");
    numThreads = getIntParam("--thr");
    thrWidth = getIntParam("--thr", 1);
    NiterWeight = getIntParam("--iter");
}

// Show ====================================================================
void ProgRecFourier::show()
{
    if (verbose > 0)
    {
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Direct 3D reconstruction method using Kaiser windows as interpolators" << std::endl;
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Input selfile             : "  << fn_sel << std::endl;
        std::cout << " padding_factor_proj       : "  << padding_factor_proj << std::endl;
        std::cout << " padding_factor_vol        : "  << padding_factor_vol << std::endl;
        std::cout << " Output volume             : "  << fn_out << std::endl;
        if (fn_sym != "")
            std::cout << " Symmetry file for projections : "  << fn_sym << std::endl;
        if (fn_fsc != "")
            std::cout << " File root for FSC files: " << fn_fsc << std::endl;
        if (do_weights)
            std::cout << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cout << " Do NOT use weights" << std::endl;
        std::cout << "\n Interpolation Function"
        << "\n   blrad                 : "  << blob.radius
        << "\n   blord                 : "  << blob.order
        << "\n   blalpha               : "  << blob.alpha
        //<< "\n sampling_rate           : "  << sampling_rate
        << "\n max_resolution          : "  << maxResolution
        << "\n -----------------------------------------------------------------" << std::endl;
    }
}

// Main routine ------------------------------------------------------------
void ProgRecFourier::run()
{
    show();
    produceSideinfo();
    // Process all images in the selfile
    if (verbose)
        init_progress_bar(NiterWeight*SF.size());
    // Create threads stuff
    barrier_init( &barrier, numThreads+1 );
    pthread_mutex_init( &workLoadMutex, NULL );
    statusArray = NULL;
    th_ids = (pthread_t *)malloc( numThreads * sizeof( pthread_t));
    th_args = (ImageThreadParams *) malloc ( numThreads * sizeof( ImageThreadParams ) );

    // Create threads
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        th_args[nt].selFile = new MetaData(SF);
        pthread_create( (th_ids+nt) , NULL, processImageThread, (void *)(th_args+nt) );
    }

    processImages(0, SF.size() - 1, !fn_fsc.empty(), false);

    finishComputations(fn_out);

    threadOpCode = EXIT_THREAD;

    // Waiting for threads to finish
    barrier_wait( &barrier );
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
        pthread_join(*(th_ids+nt), NULL);
    barrier_destroy( &barrier );
}


void ProgRecFourier::produceSideinfo()
{
    // Translate the maximum resolution to digital frequency
    // maxResolution=sampling_rate/maxResolution;
    maxResolution2=maxResolution*maxResolution;

    // Read the input images
    SF.read(fn_sel);

    // Ask for memory for the output volume and its Fourier transform
    size_t objId = SF.firstObject();
    FileName fnImg;
    SF.getValue(MDL_IMAGE,fnImg,objId);
    Image<double> I;
    I.read(fnImg, HEADER);
    int Ydim=YSIZE(I());
    int Xdim=XSIZE(I());
    if (Ydim!=Xdim)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"This algorithm only works for squared images");
    imgSize=Xdim;
    volPadSizeX = volPadSizeY = volPadSizeZ=(int)(Xdim*padding_factor_vol);
    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);

    //use threads for volume inverse fourier transform, plan is created in setReal()
    transformerVol.setThreadsNumber(numThreads);
    transformerVol.setReal(Vout());

    Vout().clear(); // Free the memory so that it is available for FourierWeights
    transformerVol.getFourierAlias(VoutFourier);
    FourierWeights.initZeros(VoutFourier);

    // Ask for memory for the padded images
    size_t paddedImgSize=(size_t)(Xdim*padding_factor_proj);
    paddedImg.resize(paddedImgSize,paddedImgSize);
    paddedImg.setXmippOrigin();
    transformerImg.setReal(paddedImg);

    // Build a table of blob values
    blobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    fourierBlobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    Fourier_blob_table.resize(BLOB_TABLE_SIZE_SQRT);

    struct blobtype blobFourier,blobnormalized;
    blobFourier=blob;
    //Sjors 18aug10 blobFourier.radius/=(padding_factor_proj*Xdim);
    blobFourier.radius/=(padding_factor_vol*Xdim);
    blobnormalized=blob;
    blobnormalized.radius/=((double)padding_factor_proj/padding_factor_vol);
    double deltaSqrt     = (blob.radius*blob.radius) /(BLOB_TABLE_SIZE_SQRT-1);
    double deltaFourier  = (sqrt(3.)*Xdim/2.)/(BLOB_TABLE_SIZE_SQRT-1);

    // The interpolation kernel must integrate to 1
    double iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
    //Sjors 18aug10 double padXdim3 = padding_factor_proj * Xdim;
    double padXdim3 = padding_factor_vol * Xdim;
    padXdim3 = padXdim3 * padXdim3 * padXdim3;
    double blobTableSize = blob.radius*sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    //***
    //Following commented line seems to be the right thing but I do not understand it
    //double fourierBlobTableSize = (sqrt(3.)*Xdim*Xdim/2.)*blobFourier.radius *sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(blobTableSqrt)

    {
        //use a r*r sample instead of r
        //DIRECT_VEC_ELEM(blob_table,i)         = blob_val(delta*i, blob)  *iw0;
        VEC_ELEM(blobTableSqrt,i)    = blob_val(blobTableSize*sqrt((double)i), blob)  *iw0;
        //***
        //DIRECT_VEC_ELEM(fourierBlobTableSqrt,i) =
        //     blob_Fourier_val(fourierBlobTableSize*sqrt(i), blobFourier)*padXdim3  *iw0;
        VEC_ELEM(Fourier_blob_table,i) =
            blob_Fourier_val(deltaFourier*i, blobFourier)*padXdim3  *iw0;
        //#define DEBUG
#ifdef DEBUG

        std::cout << VEC_ELEM(Fourier_blob_table,i)
        << " " << VEC_ELEM(fourierBlobTableSqrt,i)
        << std::endl;
#endif
  #undef DEBUG

    }
    //iDelta        = 1/delta;
    iDeltaSqrt    = 1/deltaSqrt;
    iDeltaFourier = 1/deltaFourier;

    // Get symmetries
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    R_repository.push_back(Identity);
    if (fn_sym != "")
    {
        SymList SL;
        SL.readSymmetryFile(fn_sym);
        for (int isym = 0; isym < SL.symsNo(); isym++)
        {
            Matrix2D<double>  L(4, 4), R(4, 4);
            SL.getMatrices(isym, L, R);
            R.resize(3, 3);
            R_repository.push_back(R);
        }
    }
}

void ProgRecFourier::get_angles_for_image(const FileName &fn, double &rot,
        double &tilt, double &psi, double &xoff, double &yoff, bool &flip,
        double &weight, MetaData * docfile)
{
    std::vector<size_t> found;
    (*docfile).findObjects(found,MDValueEQ(MDL_IMAGE,(std::string)fn));

    if (found.size()==1)
    {
        (*docfile).getValue(MDL_ANGLE_ROT,rot,found[0]);
        (*docfile).getValue(MDL_ANGLE_TILT,tilt,found[0]);
        (*docfile).getValue(MDL_ANGLE_PSI,psi,found[0]);
        (*docfile).getValue(MDL_SHIFT_X,xoff,found[0]);
        (*docfile).getValue(MDL_SHIFT_Y,yoff,found[0]);
        flip=0;
        weight=0;
        (*docfile).getValue(MDL_FLIP,flip,found[0]);
        // COSS, ROB porque no coger weight?
        (*docfile).getValue(MDL_WEIGHT,weight,found[0]);
    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, (std::string)"Prog_RecFourier_prm: Cannot find " + fn + " in docfile " + fn_doc);
}

void * ProgRecFourier::processImageThread( void * threadArgs )
{

    ImageThreadParams * threadParams = (ImageThreadParams *) threadArgs;
    ProgRecFourier * parent = threadParams->parent;
    barrier_t * barrier = &(parent->barrier);

    int minSeparation;

    if ( (int)ceil(parent->blob.radius) > parent->thrWidth )
        minSeparation = (int)ceil(parent->blob.radius);
    else
        minSeparation = parent->thrWidth;

    minSeparation+=1;

    Matrix2D<double>  localA(3, 3), localAinv(3, 3);
    MultidimArray< std::complex<double> > localPaddedFourier;
    MultidimArray<double> localPaddedImg;
    FourierTransformer localTransformerImg;

    std::vector<size_t> objId;

    threadParams->selFile->findObjects(objId);
    ApplyGeoParams params;
    params.only_apply_shifts = true;
    MultidimArray<int> zWrapped(3*parent->volPadSizeZ),yWrapped(3*parent->volPadSizeY),xWrapped(3*parent->volPadSizeX),
    		zNegWrapped, yNegWrapped, xNegWrapped;
    zWrapped.initConstant(-1);
    yWrapped.initConstant(-1);
	xWrapped.initConstant(-1);
	zWrapped.setXmippOrigin();
    yWrapped.setXmippOrigin();
	xWrapped.setXmippOrigin();
    zNegWrapped=zWrapped;
    yNegWrapped=yWrapped;
    xNegWrapped=xWrapped;

    MultidimArray<double> x2precalculated(XSIZE(xWrapped)), y2precalculated(XSIZE(yWrapped)), z2precalculated(XSIZE(zWrapped));
    x2precalculated.initConstant(-1);
    y2precalculated.initConstant(-1);
    z2precalculated.initConstant(-1);
    x2precalculated.setXmippOrigin();
    y2precalculated.setXmippOrigin();
    z2precalculated.setXmippOrigin();
    do
    {
        barrier_wait( barrier );

        switch ( parent->threadOpCode )
        {
        case PRELOAD_IMAGE:
            {

                threadParams->read = 0;

                if ( threadParams->imageIndex >= 0 )
                {
                    // Read input image
                    double rot, tilt, psi, weight;
                    Projection proj;

                    //Read projection from selfile, read also angles and shifts if present
                    //but only apply shifts

                    proj.readApplyGeo(*(threadParams->selFile), objId[threadParams->imageIndex], params);
                    rot  = proj.rot();
                    tilt = proj.tilt();
                    psi  = proj.psi();
                    weight = proj.weight();

                    threadParams->weight = 1.;

                    if(parent->do_weights)
                        threadParams->weight = weight;
                    else if (!parent->do_weights)
                    {
                        weight=1.0;
                    }
                    else if (weight==0.0)
                    {
                        threadParams->read = 2;
                        break;
                    }

                    // Copy the projection to the center of the padded image
                    // and compute its Fourier transform
                    proj().setXmippOrigin();
                    size_t localPaddedImgSize=(size_t)(parent->imgSize*parent->padding_factor_proj);
                    if (threadParams->reprocessFlag)
                    	localPaddedFourier.initZeros(localPaddedImgSize,localPaddedImgSize/2+1);
                    else
                    {
                    	localPaddedImg.initZeros(localPaddedImgSize,localPaddedImgSize);
                    	localPaddedImg.setXmippOrigin();
                        FOR_ALL_ELEMENTS_IN_ARRAY2D(proj())
                    		A2D_ELEM(localPaddedImg,i,j)=weight*proj(i,j);
                    	CenterFFT(localPaddedImg,true);

                    	// Fourier transformer for the images
                    	localTransformerImg.setReal(localPaddedImg);
                    	localTransformerImg.FourierTransform();
                    	localTransformerImg.getFourierAlias(localPaddedFourier);
                    }

                    // Compute the coordinate axes associated to this image
                    Euler_angles2matrix(rot, tilt, psi, localA);
                    localAinv=localA.transpose();

                    threadParams->localweight = weight;
                    threadParams->localAInv = &localAinv;
                    threadParams->localPaddedFourier = &localPaddedFourier;
                    //#define DEBUG22
#ifdef DEBUG22

                    {//CORRECTO

                        if(threadParams->myThreadID%1==0)
                        {
                            proj.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                       integerToString(threadParams->imageIndex) + "proj.spi");

                            ImageXmipp save44;
                            save44()=localPaddedImg;
                            save44.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "local_padded_img.spi");

                            FourierImage save33;
                            save33()=localPaddedFourier;
                            save33.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "local_padded_fourier.spi");
                            FourierImage save22;
                            //save22()=*paddedFourier;
                            save22().alias(*(threadParams->localPaddedFourier));
                            save22.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "_padded_fourier.spi");
                        }

                    }
#endif
                    #undef DEBUG22

                    threadParams->read = 1;
                }
                break;
            }
        case EXIT_THREAD:
            return NULL;
        case PROCESS_WEIGHTS:
            {

                // Get a first approximation of the reconstruction
                double corr2D_3D=pow(parent->padding_factor_proj,2.)/
                                 (parent->imgSize* pow(parent->padding_factor_vol,3.));
                // Divide by Zdim because of the
                // the extra dimension added
                // and padding differences
                MultidimArray<double> &mFourierWeights=parent->preFourierWeights;
                for (int k=threadParams->myThreadID; k<=FINISHINGZ(mFourierWeights); k+=parent->numThreads)
                    for (int i=STARTINGY(mFourierWeights); i<=FINISHINGY(mFourierWeights); i++)
                        for (int j=STARTINGX(mFourierWeights); j<=FINISHINGX(mFourierWeights); j++)
                        {
                            double weight_kij=A3D_ELEM(mFourierWeights,k,i,j);
                            if (1.0/weight_kij>ACCURACY)
                                A3D_ELEM(parent->VoutFourier,k,i,j)*=corr2D_3D*A3D_ELEM(mFourierWeights,k,i,j);
                            else
                                A3D_ELEM(parent->VoutFourier,k,i,j)=0;
                        }
                break;
            }
        case PROCESS_IMAGE:
            {

                MultidimArray< std::complex<double> > *paddedFourier = threadParams->paddedFourier;
                bool reprocessFlag = threadParams->reprocessFlag;
                int * statusArray = parent->statusArray;

                int minAssignedRow;
                int maxAssignedRow;
                bool breakCase;
                bool assigned;
                do
                {
                    minAssignedRow = -1;
                    maxAssignedRow = -1;
                    breakCase = false;
                    assigned = false;

                    do
                    {
                        pthread_mutex_lock( &(parent->workLoadMutex) );

                        if ( parent->rowsProcessed == YSIZE(*paddedFourier) )
                        {
                            pthread_mutex_unlock( &(parent->workLoadMutex) );
                            breakCase = true;
                            break;
                        }

                        for (size_t w = 0 ; w < YSIZE(*paddedFourier) ; w++ )
                        {
                            if ( statusArray[w]==0 )
                            {
                                assigned = true;
                                minAssignedRow = w;
                                maxAssignedRow = w+minSeparation-1;

                                if ( maxAssignedRow > (int)(YSIZE(*paddedFourier)-1) )
                                    maxAssignedRow = YSIZE(*paddedFourier)-1;

                                for ( int in = (minAssignedRow - minSeparation) ; in < (int)minAssignedRow ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < (int)YSIZE(*paddedFourier) ))
                                    {
                                        if ( statusArray[in] > -1 )
                                            statusArray[in]++;
                                    }
                                }

                                for ( int in = minAssignedRow ; in <= (int)maxAssignedRow ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < (int)YSIZE(*paddedFourier) ))
                                    {
                                        if ( statusArray[in] == 0 )
                                        {
                                            statusArray[in] = -1;
                                            parent->rowsProcessed++;
                                        }
                                    }
                                }

                                for ( int in = maxAssignedRow+1 ; in <= (maxAssignedRow+minSeparation) ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < (int)YSIZE(*paddedFourier) ))
                                    {
                                        if ( statusArray[in] > -1 )
                                            statusArray[in]++;
                                    }
                                }

                                break;
                            }
                        }

                        pthread_mutex_unlock( &(parent->workLoadMutex) );
                    }
                    while ( !assigned );

                    if ( breakCase == true )
                    {
                        break;
                    }

                    Matrix2D<double> * A_SL = threadParams->symmetry;

                    // Loop over all Fourier coefficients in the padded image
                    Matrix1D<double> freq(3), gcurrent(3), real_position(3);
                    Matrix1D<int> corner1(3), corner2(3);

                    // Some alias and calculations moved from heavy loops
                    double blobRadiusSquared = parent->blob.radius * parent->blob.radius;
                    double iDeltaSqrt = parent->iDeltaSqrt;
                    Matrix1D<double> & blobTableSqrt = parent->blobTableSqrt;
                    int xsize_1 = XSIZE(parent->VoutFourier) - 1;
                    int zsize_1 = ZSIZE(parent->VoutFourier) - 1;
                    MultidimArray< std::complex<double> > &VoutFourier=parent->VoutFourier;
                    MultidimArray<double> &fourierWeights = parent->FourierWeights;
                    MultidimArray<double> &prefourierWeights = parent->preFourierWeights;
                    // Get i value for the thread
                    for (int i = minAssignedRow; i <= maxAssignedRow ; i ++ )
                    {
                        // Discarded rows can be between minAssignedRow and maxAssignedRow, check
                        if ( statusArray[i] == -1 )
                            for (int j=STARTINGX(*paddedFourier); j<=FINISHINGX(*paddedFourier); j++)
                            {
                                // Compute the frequency of this coefficient in the
                                // universal coordinate system
                                FFT_IDX2DIGFREQ(j,XSIZE(parent->paddedImg),XX(freq));
                                FFT_IDX2DIGFREQ(i,YSIZE(parent->paddedImg),YY(freq));
                                ZZ(freq)=0;
                                if (XX(freq)*XX(freq)+YY(freq)*YY(freq)>parent->maxResolution2)
                                    continue;
                                SPEED_UP_temps012;
                                M3x3_BY_V3x1(freq,*A_SL,freq);

                                // Look for the corresponding index in the volume Fourier transform
                                DIGFREQ2FFT_IDX_DOUBLE(XX(freq),parent->volPadSizeX,XX(real_position));
                                DIGFREQ2FFT_IDX_DOUBLE(YY(freq),parent->volPadSizeY,YY(real_position));
                                DIGFREQ2FFT_IDX_DOUBLE(ZZ(freq),parent->volPadSizeZ,ZZ(real_position));

                                // Put a box around that coefficient
                                XX(corner1)=CEIL (XX(real_position)-parent->blob.radius);
                                YY(corner1)=CEIL (YY(real_position)-parent->blob.radius);
                                ZZ(corner1)=CEIL (ZZ(real_position)-parent->blob.radius);
                                XX(corner2)=FLOOR(XX(real_position)+parent->blob.radius);
                                YY(corner2)=FLOOR(YY(real_position)+parent->blob.radius);
                                ZZ(corner2)=FLOOR(ZZ(real_position)+parent->blob.radius);

#ifdef DEBUG

                                std::cout << "Idx Img=(0," << i << "," << j << ") -> Freq Img=("
                                << freq.transpose() << ") ->\n    Idx Vol=("
                                << real_position.transpose() << ")\n"
                                << "   Corner1=" << corner1.transpose() << std::endl
                                << "   Corner2=" << corner2.transpose() << std::endl;
#endif
                                // Loop within the box
                                double *ptrIn;

                                ptrIn =(double *)&(A2D_ELEM(*paddedFourier, i,j));

                                // Some precalculations
                                for (int intz = ZZ(corner1); intz <= ZZ(corner2); ++intz)
                                {
									double z = intz - ZZ(real_position);
									A1D_ELEM(z2precalculated,intz)=z*z;
                                    if (A1D_ELEM(zWrapped,intz)<0)
                                    {
                                        int iz, izneg;
										fastIntWRAP(iz, intz, 0, zsize_1);
										A1D_ELEM(zWrapped,intz)=iz;
										int miz=-iz;
										fastIntWRAP(izneg, miz,0,zsize_1);
										A1D_ELEM(zNegWrapped,intz)=izneg;
                                    }
                                }
                                for (int inty = YY(corner1); inty <= YY(corner2); ++inty)
                                {
									double y = inty - YY(real_position);
									A1D_ELEM(y2precalculated,inty)=y*y;
                                    if (A1D_ELEM(yWrapped,inty)<0)
                                    {
										int iy, iyneg;
										fastIntWRAP(iy, inty, 0, zsize_1);
										A1D_ELEM(yWrapped,inty)=iy;
										int miy=-iy;
										fastIntWRAP(iyneg, miy,0,zsize_1);
										A1D_ELEM(yNegWrapped,inty)=iyneg;
                                    }
                                }
                                for (int intx = XX(corner1); intx <= XX(corner2); ++intx)
                                {
									double x = intx - XX(real_position);
									A1D_ELEM(x2precalculated,intx)=x*x;
                                    if (A1D_ELEM(xWrapped,intx)<0)
                                    {
                                    	int ix, ixneg;
										fastIntWRAP(ix, intx, 0, zsize_1);
										A1D_ELEM(xWrapped,intx)=ix;
										int mix=-ix;
										fastIntWRAP(ixneg, mix,0,zsize_1);
										A1D_ELEM(xNegWrapped,intx)=ixneg;
                                    }
                                }

                                // Actually compute
                                for (int intz = ZZ(corner1); intz <= ZZ(corner2); ++intz)
                                {
                                    double z2 = A1D_ELEM(z2precalculated,intz);
                                    int iz=A1D_ELEM(zWrapped,intz);
                                    int izneg=A1D_ELEM(zNegWrapped,intz);

                                    for (int inty = YY(corner1); inty <= YY(corner2); ++inty)
                                    {
                                        double y2z2 = A1D_ELEM(y2precalculated,inty) + z2;
                                        if (y2z2 > blobRadiusSquared)
                                            continue;
                                        int iy=A1D_ELEM(yWrapped,inty);
                                        int iyneg=A1D_ELEM(yNegWrapped,inty);

                                        for (int intx = XX(corner1); intx <= XX(corner2); ++intx)
                                        {
                                            // Compute distance to the center of the blob
                                            // Compute blob value at that distance
                                            double d2 = A1D_ELEM(x2precalculated,intx) + y2z2;

                                            if (d2 > blobRadiusSquared)
                                                continue;
                                            int aux = (int)(d2 * iDeltaSqrt + 0.5);//Same as ROUND but avoid comparison
                                            //double w = blobTableSqrt((int)aux);
                                            double w = VEC_ELEM(blobTableSqrt, aux);
                                            // Look for the location of this logical index
                                            // in the physical layout
#ifdef DEBUG

                                            std::cout << "   gcurrent=" << gcurrent.transpose()
                                            << " d=" << d << std::endl;
                                            std::cout << "   1: intx=" << intx
                                            << " inty=" << inty
                                            << " intz=" << intz << std::endl;
#endif

                                            int ix=A1D_ELEM(xWrapped,intx);
#ifdef DEBUG

                                            std::cout << "   2: ix=" << ix << " iy=" << iy
                                            << " iz=" << iz << std::endl;
#endif

                                            bool conjugate=false;
                                            int izp, iyp, ixp;
                                            if (ix > xsize_1)
                                            {
                                                izp = izneg;
                                                iyp = iyneg;
                                                ixp = A1D_ELEM(xNegWrapped,intx);
                                                conjugate=true;
                                            }
                                            else
                                            {
                                                izp=iz;
                                                iyp=iy;
                                                ixp=ix;
                                            }
#ifdef DEBUG
                                            std::cout << "   3: ix=" << ix << " iy=" << iy
                                            << " iz=" << iz << " conj="
                                            << conjugate << std::endl;
#endif

                                            // Add the weighted coefficient
                                            if (reprocessFlag)
                                                DIRECT_A3D_ELEM(fourierWeights, izp,iyp,ixp) += (w *  DIRECT_A3D_ELEM(prefourierWeights, izp,iyp,ixp));
                                            else
                                            {
                                                double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, izp,iyp,ixp));
                                                ptrOut[0] += w * ptrIn[0];

                                                if (conjugate)
                                                    ptrOut[1]-=w*ptrIn[1];
                                                else
                                                    ptrOut[1]+=w*ptrIn[1];
                                            }
                                        }
                                    }
                                }
                            }
                    }

                    pthread_mutex_lock( &(parent->workLoadMutex) );

                    for ( int w = (minAssignedRow - minSeparation) ; w < minAssignedRow ; w ++ )
                    {
                        if ( ( w >= 0 ) && ( w < (int)YSIZE(*paddedFourier) ))
                        {
                            if ( statusArray[w] > 0 )
                            {
                                statusArray[w]--;
                            }
                        }
                    }

                    for ( int w = maxAssignedRow+1 ; w <= (maxAssignedRow+minSeparation) ; w ++ )
                    {
                        if ( ( w >= 0 ) && ( w < (int)YSIZE(*paddedFourier) ))
                        {
                            if ( statusArray[w] > 0 )
                            {
                                statusArray[w]--;
                            }
                        }
                    }

                    pthread_mutex_unlock( &(parent->workLoadMutex) );

                }
                while (!breakCase);
                break;
            }
        default:
            break;
        }

        barrier_wait( barrier );
    }
    while ( 1 );
}

//#define DEBUG
void ProgRecFourier::processImages( int firstImageIndex, int lastImageIndex, bool saveFSC, bool reprocessFlag)
{
    MultidimArray< std::complex<double> > *paddedFourier;

    int repaint = (int)ceil((double)SF.size()/60);

    bool processed;
    int imgno = 0;
    int imgIndex = firstImageIndex;

    // This index tells when to save work for later FSC usage
    int FSCIndex = (firstImageIndex + lastImageIndex)/2;

    // Index of the image that has just been processed. Used for
    // FSC purposes
    int current_index;

    do
    {
        threadOpCode = PRELOAD_IMAGE;

        for ( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            if ( imgIndex <= lastImageIndex )
            {
                th_args[nt].imageIndex = imgIndex;
                th_args[nt].reprocessFlag = reprocessFlag;
                imgIndex++;
            }
            else
            {
                th_args[nt].imageIndex = -1;
            }
        }

        // Awaking sleeping threads
        barrier_wait( &barrier );
        // here each thread is reading a different image and compute fft
        // Threads are working now, wait for them to finish
        // processing current projection
        barrier_wait( &barrier );

        // each threads have read a different image and now
        // all the thread will work in a different part of a single image.
        threadOpCode = PROCESS_IMAGE;

        processed = false;

        for ( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            if ( th_args[nt].read == 2 )
                processed = true;
            else if ( th_args[nt].read == 1 )
            {
                processed = true;
                if (verbose && imgno++%repaint==0)
                    progress_bar(NiterWeight*imgno);

                double weight = th_args[nt].localweight;
                paddedFourier = th_args[nt].localPaddedFourier;
                current_index = th_args[nt].imageIndex;
                Matrix2D<double> *Ainv = th_args[nt].localAInv;

                //#define DEBUG22
#ifdef DEBUG22

                {
                    static int ii=0;
                    if(ii%1==0)
                    {
                        FourierImage save22;
                        //save22()=*paddedFourier;
                        save22().alias(*paddedFourier);
                        save22.write((std::string) integerToString(ii)  + "_padded_fourier.spi");
                    }
                    ii++;
                }
#endif
                #undef DEBUG22

                // Initialized just once
                if ( statusArray == NULL )
                {
                    statusArray = (int *) malloc ( sizeof(int) * paddedFourier->ydim );
                }

                // Determine how many rows of the fourier
                // transform are of interest for us. This is because
                // the user can avoid to explore at certain resolutions
                size_t conserveRows=(size_t)ceil((double)paddedFourier->ydim * maxResolution * 2.0);
                conserveRows=(size_t)ceil((double)conserveRows/2.0);

                // Loop over all symmetries
                for (size_t isym = 0; isym < R_repository.size(); isym++)
                {
                    rowsProcessed = 0;

                    // Compute the coordinate axes of the symmetrized projection
                    Matrix2D<double> A_SL=R_repository[isym]*(*Ainv);

                    // Poner lo necesario en la estructura de cada hilo.
                    for ( int th = 0 ; th < numThreads ; th ++ )
                    {
                        // Passing parameters to each thread
                        th_args[th].symmetry = &A_SL;
                        th_args[th].paddedFourier = paddedFourier;
                        th_args[th].weight = weight;
                        th_args[th].reprocessFlag = reprocessFlag;
                    }

                    // Init status array
                    for (size_t i = 0 ; i < paddedFourier->ydim ; i ++ )
                    {
                        if ( i >= conserveRows && i < (paddedFourier->ydim-conserveRows))
                        {
                            // -2 means "discarded"
                            statusArray[i] = -2;
                            rowsProcessed++;
                        }
                        else
                        {
                            statusArray[i] = 0;
                        }
                    }

                    // Awaking sleeping threads
                    barrier_wait( &barrier );
                    // Threads are working now, wait for them to finish
                    // processing current projection
                    barrier_wait( &barrier );

                    //#define DEBUG2
#ifdef DEBUG2

                    {
                        static int ii=0;
                        if(ii%1==0)
                        {
                            Image<double> save;
                            save().alias( FourierWeights );
                            save.write((std::string) integerToString(ii)  + "_1_Weights.vol");

                            Image< std::complex<double> > save2;
                            save2().alias( VoutFourier );
                            save2.write((std::string) integerToString(ii)  + "_1_Fourier.vol");
                        }
                        ii++;
                    }
#endif
                    #undef DEBUG2

                }

                if ( current_index == FSCIndex && saveFSC )
                {
                    // Save Current Fourier, Reconstruction and Weights
                    Image<double> save;
                    save().alias( FourierWeights );
                    save.write((std::string)fn_fsc + "_1_Weights.vol");

                    Image< std::complex<double> > save2;
                    save2().alias( VoutFourier );
                    save2.write((std::string) fn_fsc + "_1_Fourier.vol");

                    finishComputations(FileName((std::string) fn_fsc + "_1_recons.vol"));
                    Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                    transformerVol.setReal(Vout());
                    Vout().clear();
                    transformerVol.getFourierAlias(VoutFourier);
                    FourierWeights.initZeros(VoutFourier);
                    VoutFourier.initZeros();
                }
            }
        }
    }
    while ( processed );

    if( saveFSC )
    {
        // Save Current Fourier, Reconstruction and Weights
        Image<double> auxVolume;
        auxVolume().alias( FourierWeights );
        auxVolume.write((std::string)fn_fsc + "_2_Weights.vol");

        Image< std::complex<double> > auxFourierVolume;
        auxFourierVolume().alias( VoutFourier );
        auxFourierVolume.write((std::string) fn_fsc + "_2_Fourier.vol");

        finishComputations(FileName((std::string) fn_fsc + "_2_recons.vol"));

        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
        transformerVol.setReal(Vout());
        Vout().clear();
        transformerVol.getFourierAlias(VoutFourier);
        FourierWeights.initZeros(VoutFourier);
        VoutFourier.initZeros();

        auxVolume.sumWithFile(fn_fsc + "_1_Weights.vol");
        auxVolume.sumWithFile(fn_fsc + "_2_Weights.vol");
        auxFourierVolume.sumWithFile(fn_fsc + "_1_Fourier.vol");
        auxFourierVolume.sumWithFile(fn_fsc + "_2_Fourier.vol");
        remove((fn_fsc + "_1_Weights.vol").c_str());
        remove((fn_fsc + "_2_Weights.vol").c_str());
        remove((fn_fsc + "_1_Fourier.vol").c_str());
        remove((fn_fsc + "_2_Fourier.vol").c_str());

        /*
        //Save SUM
                                    //this is an image but not an xmipp image
                                    auxFourierVolume.write((std::string)fn_fsc + "_all_Fourier.vol",
                                            false,VDOUBLE);
                                    auxVolume.write((std::string)fn_fsc + "_all_Weight.vol",
                                            false,VDOUBLE);
        //
        */
    }
}

void ProgRecFourier::correctWeight()
{

    preFourierWeights.initZeros(FourierWeights);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(preFourierWeights)
    A3D_ELEM(preFourierWeights,k,i,j) = 1;
    for (int i=0;i<NiterWeight;i++)
    {
        FourierWeights.initZeros(preFourierWeights);
        processImages(0, SF.size() - 1, !fn_fsc.empty(), true);
        forceWeightSymmetry(FourierWeights);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FourierWeights)
        {
            double weight_kij=A3D_ELEM(FourierWeights,k,i,j);
            double preweight_kij=A3D_ELEM(preFourierWeights,k,i,j);
            A3D_ELEM(preFourierWeights,k,i,j) = preweight_kij/weight_kij;
        }
    }
}

void ProgRecFourier::finishComputations( const FileName &out_name )
{
    //#define DEBUG_VOL
#ifdef DEBUG_VOL
    {
        VolumeXmipp save;
        save().alias( FourierWeights );
        save.write((std::string) fn_out + "Weights.vol");

        FourierVolume save2;
        save2().alias( VoutFourier );
        save2.write((std::string) fn_out + "FourierVol.vol");
    }
#endif

    // Method for correcting the weights
    correctWeight();
    // Enforce symmetry in the Fourier values as well as the weights
    // Sjors 19aug10 enforceHermitianSymmetry first checks ndim...
    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);
    transformerVol.setReal(Vout());
    transformerVol.enforceHermitianSymmetry();
	//forceWeightSymmetry(preFourierWeights);

    // Tell threads what to do
    //#define DEBUG_VOL1
#ifdef DEBUG_VOL1

    {
        Image<double> save;
        save().alias( FourierWeights );
        save.write((std::string) fn_out + "hermiticWeights.vol");

        Image< std::complex<double> > save2;
        save2().alias( VoutFourier );
        save2.write((std::string) fn_out + "hermiticFourierVol.vol");
    }
#endif
    threadOpCode = PROCESS_WEIGHTS;
    // Awake threads
    barrier_wait( &barrier );
    // Threads are working now, wait for them to finish
    barrier_wait( &barrier );

    transformerVol.inverseFourierTransform();
    CenterFFT(Vout(),false);

    // Correct by the Fourier transform of the blob
    Vout().setXmippOrigin();
    Vout().selfWindow(FIRST_XMIPP_INDEX(imgSize),FIRST_XMIPP_INDEX(imgSize),
                      FIRST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize),
                      LAST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize));
    double pad_relation= ((double)padding_factor_proj/padding_factor_vol);
    pad_relation = (pad_relation * pad_relation * pad_relation);

    MultidimArray<double> &mVout=Vout();
    double ipad_relation=1.0/pad_relation;
    double meanFactor2=0;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
    {
        double radius=sqrt((double)(k*k+i*i+j*j));
        double aux=radius*iDeltaFourier;
        double factor = Fourier_blob_table(ROUND(aux));
        double factor2=(pow(Sinc(radius/(2*(imgSize))),2));
        A3D_ELEM(mVout,k,i,j) /= (ipad_relation*factor2*factor);
        meanFactor2+=factor2;
    }
    meanFactor2/=MULTIDIM_SIZE(mVout);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
    A3D_ELEM(mVout,k,i,j) *= meanFactor2;
    Vout.write(out_name);
}

void ProgRecFourier::setIO(const FileName &fn_in, const FileName &fn_out)
{
    this->fn_sel = fn_in;
    this->fn_out = fn_out;
}

void ProgRecFourier::forceWeightSymmetry(MultidimArray<double> &FourierWeights)
{
    int yHalf=YSIZE(FourierWeights)/2;
    if (YSIZE(FourierWeights)%2==0)
        yHalf--;
    int zHalf=ZSIZE(FourierWeights)/2;
    if (ZSIZE(FourierWeights)%2==0)
        zHalf--;
    int zsize=(int)ZSIZE(FourierWeights);
    int zsize_1=zsize-1;
    int ysize_1=(int)YSIZE(FourierWeights)-1;
    for (int k=0; k<zsize; k++)
    {
        int ksym=intWRAP(-k,0,zsize_1);
        for (int i=1; i<=yHalf; i++)
        {
            int isym=intWRAP(-i,0,ysize_1);
            double mean=0.5*(
                            DIRECT_A3D_ELEM(FourierWeights,k,i,0)+
                            DIRECT_A3D_ELEM(FourierWeights,ksym,isym,0));
            DIRECT_A3D_ELEM(FourierWeights,k,i,0)=
                DIRECT_A3D_ELEM(FourierWeights,ksym,isym,0)=mean;
        }
    }
    for (int k=1; k<=zHalf; k++)
    {
        int ksym=intWRAP(-k,0,zsize_1);
        double mean=0.5*(
                        DIRECT_A3D_ELEM(FourierWeights,k,0,0)+
                        DIRECT_A3D_ELEM(FourierWeights,ksym,0,0));
        DIRECT_A3D_ELEM(FourierWeights,k,0,0)=
            DIRECT_A3D_ELEM(FourierWeights,ksym,0,0)=mean;
    }
}
