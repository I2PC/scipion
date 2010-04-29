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
#include <data/args.h>
#include <data/fft.h>
#include <sys/time.h>
// Read arguments ==========================================================
void Prog_RecFourier_prm::read(int argc, char **argv)
{
    fn_sel = getParameter(argc, argv, "-i");
    fn_doc = getParameter(argc, argv, "-doc","");
    fn_out = getParameter(argc, argv, "-o", "rec_fourier.vol");
    fn_sym = getParameter(argc, argv, "-sym", "");
    fn_fsc = getParameter(argc, argv, "-prepare_fsc","");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    do_weights = checkParameter(argc, argv, "-weight");
    padding_factor_proj = textToFloat(getParameter(argc, argv, "-pad_proj","2"));
    padding_factor_vol = textToFloat(getParameter(argc, argv, "-pad_vol","2"));
    fn_control    = getParameter(argc, argv, "-control", "");
    blob.radius   = textToFloat(getParameter(argc, argv, "-r","1.9"));
    blob.order    = textToFloat(getParameter(argc, argv, "-m","0"));
    blob.alpha    = textToFloat(getParameter(argc, argv, "-a","15"));
    //sampling_rate = textToFloat(getParameter(argc, argv, "-sampling_rate", "1"));
    maxResolution = textToFloat(getParameter(argc, argv,
                                             "-max_resolution",".5"));
    numThreads = textToInteger(getParameter(argc, argv, "-thr", "1"));
    thrWidth = textToInteger(getParameter(argc,argv, "-thr_width", "-1"));


}

// Show ====================================================================
void Prog_RecFourier_prm::show()
{
    if (verb > 0)
    {
        std::cerr << " =====================================================================" << std::endl;
        std::cerr << " Direct 3D reconstruction method using Kaiser windows as interpolators" << std::endl;
        std::cerr << " =====================================================================" << std::endl;
        std::cerr << " Input selfile             : "  << fn_sel << std::endl;
        std::cerr << " padding_factor_proj       : "  << padding_factor_proj << std::endl;
        std::cerr << " padding_factor_vol        : "  << padding_factor_vol << std::endl;
        if (fn_doc != "")
            std::cerr << " Input docfile         : "  << fn_doc << std::endl;
        std::cerr << " Output volume             : "  << fn_out << std::endl;
        if (fn_sym != "")
            std::cerr << " Symmetry file for projections : "  << fn_sym << std::endl;
        if (fn_fsc != "")
            std::cerr << " File root for FSC files: " << fn_fsc << std::endl;
        if (do_weights)
            std::cerr << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cerr << " Do NOT use weights" << std::endl;
        std::cerr << "\n Interpolation Function" 
        << "\n   blrad                 : "  << blob.radius
        << "\n   blord                 : "  << blob.order
        << "\n   blalpha               : "  << blob.alpha
        //<< "\n sampling_rate           : "  << sampling_rate
        << "\n max_resolution          : "  << maxResolution
        << "\n -----------------------------------------------------------------" << std::endl;
    }
}

// Usage ====================================================================
void Prog_RecFourier_prm::usage()
{

    // To screen
    std::cerr << "  Usage:\n";
    std::cerr << "  reconstruct_fourier  <options>\n";
    std::cerr << "   -i <input selfile>          : Selection file with input images \n";
    std::cerr << " [ -o <\"rec_fourier.vol\">]     : Filename for output volume \n";
    std::cerr << " [ -sym <symfile> ]            : Enforce symmetry in projections\n";
    std::cerr << " [ -pad_proj <p=2.0> ]         : Projection padding factor \n";
    std::cerr << " [ -pad_vol  <p=2.0> ]         : Volume padding factor \n";
    std::cerr << " [ -prepare_fsc <fscfile> ]    : Filename root for FSC files \n";
    std::cerr << " [ -doc <docfile> ]            : Ignore headers and get angles from this docfile \n";
    std::cerr << " [ -max_resolution <p=0.5> ]   : Max resolution (Nyquist=0.5) \n";
    std::cerr << " [ -weight ]                   : Use weights stored in the image headers or doc file\n";
    std::cerr << " [ -thr <threads=1> ]          : Number of concurrent threads\n";
    std::cerr << " [ -thr_width <width=blob_radius> : Number of image rows processed at a time by a thread\n";
    std::cerr << " -----------------------------------------------------------------\n";
    std::cerr << "Interpolation Function\n ";
    std::cerr << " [-r blrad=1.9]                : blob radius in pixels\n";
    std::cerr << " [-m blord=0]                  : order of Bessel function in blob\n";
    std::cerr << " [-a blalpha=15]               : blob parameter alpha\n";
    std::cerr << " -----------------------------------------------------------------"<< std::endl;
}

void Prog_RecFourier_prm::produce_Side_info()
{
    // Translate the maximum resolution to digital frequency
    // maxResolution=sampling_rate/maxResolution;
    maxResolution2=maxResolution*maxResolution;

    // Read the input images
    SF.read(fn_sel);
    
    // Read docfile and get column numbers
    if (fn_doc != "")
    {
        DF.read(fn_doc);
        if (SF.size() != DF.size())
            REPORT_ERROR(1, "docfile and corresponding selfile have unequal (active) entries");
    }


    // Ask for memory for the output volume and its Fourier transform
    SF.firstObject();
    FileName fnImg; SF.getValue(MDL_IMAGE,fnImg);
    ImageXmipp I(fnImg);
    int Ydim=YSIZE(I());
    int Xdim=XSIZE(I());
    if (Ydim!=Xdim)
        REPORT_ERROR(1,"This algorithm only works for squared images");
    imgSize=Xdim;
    volPadSizeX = volPadSizeY = volPadSizeZ=Xdim*padding_factor_vol;
    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);

    //use threads for volume inverse fourier transform, plan is created in setReal()
    transformerVol.setThreadsNumber(numThreads);
    transformerVol.setReal(Vout());

    Vout().clear(); // Free the memory so that it is available for FourierWeights
    transformerVol.getFourierAlias(VoutFourier);
    FourierWeights.initZeros(VoutFourier);
   
    // Ask for memory for the padded images
    paddedImg.resize(Xdim*padding_factor_proj,Xdim*padding_factor_proj);
    paddedImg.setXmippOrigin();
    transformerImg.setReal(paddedImg);

    // Build a table of blob values
    blobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    fourierBlobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    Fourier_blob_table.resize(BLOB_TABLE_SIZE_SQRT);

    struct blobtype blobFourier,blobnormalized;
    blobFourier=blob;
    blobFourier.radius/=(padding_factor_proj*Xdim);
    blobnormalized=blob;
    blobnormalized.radius/=((double)padding_factor_proj/padding_factor_vol);
    double deltaSqrt     = (blob.radius*blob.radius) /(BLOB_TABLE_SIZE_SQRT-1);
    double deltaFourier  = (sqrt(3.)*Xdim/2.)/(BLOB_TABLE_SIZE_SQRT-1);

    // The interpolation kernel must integrate to 1
    double iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
    double padXdim3 = padding_factor_proj * Xdim;
    padXdim3 = padXdim3 * padXdim3 * padXdim3;
    double blobTableSize = blob.radius*sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    //***
    //Following commented line seems to be the right thing but I do not understand it
    //double fourierBlobTableSize = (sqrt(3.)*Xdim*Xdim/2.)*blobFourier.radius *sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(blobTableSqrt)
    {
    	//use a r*r sample instead of r
    	//DIRECT_VEC_ELEM(blob_table,i)         = blob_val(delta*i, blob)  *iw0;
        DIRECT_VEC_ELEM(blobTableSqrt,i)    = blob_val(blobTableSize*sqrt((double)i), blob)  *iw0;
        //***
        //DIRECT_VEC_ELEM(fourierBlobTableSqrt,i) =
        //     blob_Fourier_val(fourierBlobTableSize*sqrt(i), blobFourier)*padXdim3  *iw0;
        DIRECT_VEC_ELEM(Fourier_blob_table,i) =
             blob_Fourier_val(deltaFourier*i, blobFourier)*padXdim3  *iw0;
		//#define DEBUG
		#ifdef DEBUG
			std::cout << DIRECT_VEC_ELEM(Fourier_blob_table,i)
					  << " " << DIRECT_VEC_ELEM(fourierBlobTableSqrt,i)
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
        SL.read_sym_file(fn_sym);
        for (int isym = 0; isym < SL.SymsNo(); isym++)
        {
            Matrix2D<double>  L(4, 4), R(4, 4);
            SL.get_matrices(isym, L, R);
            R.resize(3, 3);
            R_repository.push_back(R);
        }
    }
}

void Prog_RecFourier_prm::get_angles_for_image(const FileName &fn, double &rot,
                                               double &tilt, double &psi, double &xoff, double &yoff, double &flip,
                                               double &weight, MetaData * docfile)
{

    std::vector<long int> found=(*docfile).findObjects(MDL_IMAGE,(std::string)fn);
    if (found.size()==1)
    {
    	(*docfile).goToObject(found[0]);
        (*docfile).getValue(MDL_ANGLEROT,rot);
        (*docfile).getValue(MDL_ANGLETILT,tilt);
        (*docfile).getValue(MDL_ANGLEPSI,psi);
        (*docfile).getValue(MDL_SHIFTX,xoff);
        (*docfile).getValue(MDL_SHIFTY,yoff);
        flip=0;
        weight=0;
        bool iflip;
        (*docfile).getValue(MDL_FLIP,iflip);
        flip=iflip;
        // COSS (*docfile).getValue(MDL_WEIGHT,weight);
    }
    else
        REPORT_ERROR(1, (std::string)"Prog_RecFourier_prm: Cannot find " + fn + " in docfile " + fn_doc);
}

void * Prog_RecFourier_prm::processImageThread( void * threadArgs )
{
    ImageThreadParams * threadParams = (ImageThreadParams *) threadArgs;
    Prog_RecFourier_prm * parent = threadParams->parent;
    barrier_t * barrier = &(parent->barrier);
    
    int minSeparation;
    
    if ( (int)ceil(parent->blob.radius) > parent->thrWidth )
        minSeparation = (int)ceil(parent->blob.radius);
    else
        minSeparation = parent->thrWidth;

    minSeparation+=1;

    Matrix2D<double>  localA(3, 3), localAinv(3, 3);
    Matrix2D< std::complex<double> > localPaddedFourier;
    Matrix2D<double> localPaddedImg;
    XmippFftw localTransformerImg;

    MetaData * docFile = threadParams->docFile;

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
                    FileName fn_img;
                    parent->SF.getValue(MDL_IMAGE, fn_img, threadParams->imageIndex );

                    // Read input image
                    double rot, tilt, psi, xoff,yoff,flip,weight;
                    Projection proj;

                    if (parent->fn_doc == "")
                    {
                        proj.read(fn_img, true); //true means apply shifts 
                        rot  = proj.rot();
                        tilt = proj.tilt();
                        psi  = proj.psi();
                        weight = proj.weight();
                    }
                    else
                    {
                        proj.read(fn_img, false); // do not apply shifts since they are not in
                        // the header
                        parent->get_angles_for_image(fn_img, rot, tilt, psi, xoff, yoff, flip, weight, docFile);

                        proj.set_angles(rot, tilt, psi); 
                        proj.set_Xoff(xoff);
                        proj.set_Yoff(yoff);
                        proj.set_flip(flip);
                        proj.set_weight(weight);
                        localA = proj.get_transformation_matrix(true);
                        if (!localA.isIdentity())
                            proj().selfApplyGeometryBSpline(localA, 3, IS_INV, WRAP);
                    }

                    threadParams->weight = 1.;
                    
                    if(parent->do_weights)
                        threadParams->weight = weight;
                    
                    if (!parent->do_weights)
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
                    localPaddedImg.resize(parent->imgSize*parent->padding_factor_proj,
                                          parent->imgSize*parent->padding_factor_proj);
                    localPaddedImg.setXmippOrigin();
                    //added ROB
                    FOR_ALL_ELEMENTS_IN_MATRIX2D(localPaddedImg)
                         localPaddedImg(i,j)=0.;
                    FOR_ALL_ELEMENTS_IN_MATRIX2D(proj())
                         localPaddedImg(i,j)=weight*proj(i,j);
                    CenterFFT(localPaddedImg,true);

                    // Fourier transformer for the images
                    localTransformerImg.setReal(localPaddedImg);
                    localTransformerImg.FourierTransform();
                    localTransformerImg.getFourierAlias(localPaddedFourier);

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

                for (size_t k=threadParams->myThreadID; k<=FINISHINGZ(parent->FourierWeights); k+=parent->numThreads)
                    for (size_t i=STARTINGY(parent->FourierWeights); i<=FINISHINGY(parent->FourierWeights); i++)
                        for (size_t j=STARTINGX(parent->FourierWeights); j<=FINISHINGX(parent->FourierWeights); j++)
                        {
                            if (parent->FourierWeights(k,i,j)>ACCURACY)
                            {
                                parent->FourierWeights(k,i,j)=1/parent->FourierWeights(k,i,j);
                                VOL_ELEM(parent->VoutFourier,k,i,j)*=corr2D_3D*VOL_ELEM(parent->FourierWeights,k,i,j);
                            }
                            else
                                VOL_ELEM(parent->VoutFourier,k,i,j)=0;

                        }

                break;
            }
        case PROCESS_IMAGE:
            {   
                Matrix2D< std::complex<double> > *paddedFourier = threadParams->paddedFourier;
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

                        if ( parent->rowsProcessed == paddedFourier->ydim )
                        {
                            pthread_mutex_unlock( &(parent->workLoadMutex) );
                            breakCase = true;
                            break;
                        }

                        for (int w = 0 ; w < paddedFourier->ydim ; w++ )
                        {
                            if ( statusArray[w]==0 )
                            {
                                assigned = true;
                                minAssignedRow = w;
                                maxAssignedRow = w+minSeparation-1;

                                if ( maxAssignedRow > (paddedFourier->ydim-1) )
                                    maxAssignedRow = paddedFourier->ydim-1;

                                for ( int in = (minAssignedRow - minSeparation) ; in < minAssignedRow ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < paddedFourier->ydim ))
                                    {
                                        if ( statusArray[in] > -1 )
                                            statusArray[in]++;
                                    }
                                }                                    

                                for ( int in = minAssignedRow ; in <= maxAssignedRow ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < paddedFourier->ydim ))
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
                                    if ( ( in >= 0 ) && ( in < paddedFourier->ydim ))
                                    {
                                        if ( statusArray[in] > -1 )
                                            statusArray[in]++;
                                    }
                                }

                                break;
                            }
                        }

                        pthread_mutex_unlock( &(parent->workLoadMutex) );
                    }while ( !assigned );

                    if ( breakCase == true )
                    {
                        break;
                    }

                    Matrix2D<double> * A_SL = threadParams->symmetry;

                    // Loop over all Fourier coefficients in the padded image
                    Matrix1D<double> freq(3), gcurrent(3), real_position(3);
                    Matrix1D<int> corner1(3), corner2(3);

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
                                SPEED_UP_temps;
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

                                double *ptrIn =(double *)&((*paddedFourier)(i,j));
                                double blobRadiusSquared = parent->blob.radius * parent->blob.radius;
                                for (int intz = ZZ(corner1); intz <= ZZ(corner2); intz++)
                                {
                                    for (int inty = YY(corner1); inty <= YY(corner2); inty++)
                                    {
                                        for (int intx = XX(corner1); intx <= XX(corner2); intx++)
                                        {
                                            // Compute distance to the center of the blob
                                            // Compute blob value at that distance
                                            VECTOR_R3(gcurrent, intx, inty, intz);
                                            V3_MINUS_V3(gcurrent, real_position, gcurrent);
                                            double d = (XX(gcurrent) * XX(gcurrent) +
                                                        YY(gcurrent) * YY(gcurrent) +
                                                        ZZ(gcurrent) * ZZ(gcurrent));
                                            if (d > blobRadiusSquared) continue;
                                            // COSS: *** AVOID THE SQUARE ROOTS
                                            double w = parent->blobTableSqrt(ROUND(d*parent->iDeltaSqrt));
                                            //double w = parent->blob_table(ROUND(gcurrent.module()*parent->iDelta));
                                            //if(w<MINIMUMWEIGHT)
                                            //   continue;
                                            // Look for the location of this logical index
                                            // in the physical layout
#ifdef DEBUG
                                            std::cout << "   gcurrent=" << gcurrent.transpose()
                                            << " d=" << d << std::endl;
                                            std::cout << "   1: intx=" << intx
                                            << " inty=" << inty
                                            << " intz=" << intz << std::endl;
#endif
                                            int iz=intWRAP(intz,0,ZSIZE(parent->VoutFourier)-1);
                                            int iy=intWRAP(inty,0,ZSIZE(parent->VoutFourier)-1);
                                            int ix=intWRAP(intx,0,ZSIZE(parent->VoutFourier)-1);
#ifdef DEBUG
                                            std::cout << "   2: ix=" << ix << " iy=" << iy
                                            << " iz=" << iz << std::endl;
#endif
                                            bool conjugate=false;
                                            if (ix>=XSIZE(parent->VoutFourier))
                                            {
                                                iz=intWRAP(-iz,0,ZSIZE(parent->VoutFourier)-1);
                                                iy=intWRAP(-iy,0,ZSIZE(parent->VoutFourier)-1);
                                                ix=intWRAP(-ix,0,ZSIZE(parent->VoutFourier)-1);
                                                conjugate=true;
                                            }
#ifdef DEBUG
                                            std::cout << "   3: ix=" << ix << " iy=" << iy
                                            << " iz=" << iz << " conj="
                                            << conjugate << std::endl;
#endif

                                            // Add the weighted coefficient
                                            double *ptrOut=(double *)&((parent->VoutFourier)(iz,iy,ix));
                                            ptrOut[0]+=w*ptrIn[0];
                                            if (conjugate) ptrOut[1]-=w*ptrIn[1];
                                            else           ptrOut[1]+=w*ptrIn[1];

                                            (parent->FourierWeights)(iz,iy,ix)+=w*threadParams->weight ;
                                        }
                                    }
                                }
                            }
                    }   

                    pthread_mutex_lock( &(parent->workLoadMutex) );

                    for ( int w = (minAssignedRow - minSeparation) ; w < minAssignedRow ; w ++ )
                    {
                        if ( ( w >= 0 ) && ( w < paddedFourier->ydim ))
                        {
                            if ( statusArray[w] > 0 )
                            {
                                statusArray[w]--;
                            }
                        }
                    }    

                    for ( int w = maxAssignedRow+1 ; w <= (maxAssignedRow+minSeparation) ; w ++ )
                    {
                        if ( ( w >= 0 ) && ( w < paddedFourier->ydim ))
                        {
                            if ( statusArray[w] > 0 )
                            {
                                statusArray[w]--;
                            }
                        }
                    }

                    pthread_mutex_unlock( &(parent->workLoadMutex) );

                }while (!breakCase);  
                break;
            }
        default: break;
        }

        barrier_wait( barrier );

    }while ( 1 );
}

//#define DEBUG
void Prog_RecFourier_prm::processImages( int firstImageIndex, int lastImageIndex, bool saveFSC )
{
    Matrix2D< std::complex<double> > *paddedFourier;

    int repaint = ceil((double)SF.size()/60);

    bool processed;
    int imgno = 0;
    int imgIndex = firstImageIndex;
    struct timeval start_time, end_time;
    long int total_usecs;
    double total_time;
    
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
                if (verb && imgno++%repaint==0) progress_bar(imgno);
            
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
                int conserveRows= ceil((double)paddedFourier->ydim * maxResolution * 2.0);
                conserveRows= ceil((double)conserveRows/2.0);

                // Loop over all symmetries
                for (int isym = 0; isym < R_repository.size(); isym++)
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
                    }

                    // Init status array
                    for (int i = 0 ; i < paddedFourier->ydim ; i ++ )
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
                                VolumeT save;
                                save().alias( FourierWeights );
                                save.write((std::string) integerToString(ii)  + "_1_Weights.vol");

                                FourierVolume save2;
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
                    VolumeT<double> save;
                    save().alias( FourierWeights );
                    save.write((std::string)fn_fsc + "_1_Weights.vol",
                            false,VDOUBLE);
                    
                    FourierVolume save2;
                    save2().alias( VoutFourier );
                    save2.write((std::string) fn_fsc + "_1_Fourier.vol",
                            false,VDOUBLE);
                        
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
    }while ( processed );
                      
    if( saveFSC ) 
    {
        // Save Current Fourier, Reconstruction and Weights
        VolumeT<double> auxVolume;
        auxVolume().alias( FourierWeights );
        auxVolume.write((std::string)fn_fsc + "_2_Weights.vol",
                false,VDOUBLE);

        FourierVolume auxFourierVolume;
        auxFourierVolume().alias( VoutFourier );
        auxFourierVolume.write((std::string) fn_fsc + "_2_Fourier.vol",
                false,VDOUBLE);

        finishComputations(FileName((std::string) fn_fsc + "_2_recons.vol"));

        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
        transformerVol.setReal(Vout());
        Vout().clear();
        transformerVol.getFourierAlias(VoutFourier);
        FourierWeights.initZeros(VoutFourier);
        VoutFourier.initZeros();   

        auxVolume.sumWithFile((std::string) fn_fsc + "_1_Weights.vol",false, VDOUBLE);
        auxVolume.sumWithFile((std::string) fn_fsc + "_2_Weights.vol",false, VDOUBLE);

        auxFourierVolume.sumWithFile((std::string) fn_fsc + "_1_Fourier.vol",false, VDOUBLE);
        auxFourierVolume.sumWithFile((std::string) fn_fsc + "_2_Fourier.vol",false, VDOUBLE);
        remove(((std::string) fn_fsc + "_1_Weights.vol").c_str());
        remove(((std::string) fn_fsc + "_2_Weights.vol").c_str());
	remove(((std::string) fn_fsc + "_1_Fourier.vol").c_str());
	remove(((std::string) fn_fsc + "_2_Fourier.vol").c_str());

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

// Main routine ------------------------------------------------------------
void Prog_RecFourier_prm::run()
{
    // Process all images in the selfile
    if (verb) init_progress_bar(SF.size());

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
        th_args[nt].docFile = new MetaData(DF);
        pthread_create( (th_ids+nt) , NULL, processImageThread, (void *)(th_args+nt) );
    }

    if( fn_fsc != "" )
        processImages(0,SF.size()-1,true);
    else
        processImages(0,SF.size()-1,false);

    finishComputations(fn_out);

    threadOpCode = EXIT_THREAD;
    barrier_wait( &barrier );

    // Waiting for threads to finish
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
    {
        pthread_join(*(th_ids+nt), NULL);
    }

    barrier_destroy( &barrier );

}

void Prog_RecFourier_prm::finishComputations( FileName out_name )
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

    // Enforce symmetry in the Fourier values as well as the weights
    transformerVol.enforceHermitianSymmetry();
    int yHalf=YSIZE(FourierWeights)/2;
    if (YSIZE(FourierWeights)%2==0) yHalf--;
    int zHalf=ZSIZE(FourierWeights)/2;
    if (ZSIZE(FourierWeights)%2==0) zHalf--;
    for (int k=0; k<ZSIZE(FourierWeights); k++)
    {
        int ksym=intWRAP(-k,0,ZSIZE(FourierWeights)-1);
        for (int i=1; i<=yHalf; i++)
        {
            int isym=intWRAP(-i,0,YSIZE(FourierWeights)-1);
            double mean=0.5*(
                            DIRECT_VOL_ELEM(FourierWeights,k,i,0)+
                            DIRECT_VOL_ELEM(FourierWeights,ksym,isym,0));
            DIRECT_VOL_ELEM(FourierWeights,k,i,0)=
            DIRECT_VOL_ELEM(FourierWeights,ksym,isym,0)=mean;
        }
    }
    for (int k=1; k<=zHalf; k++)
    {
        int ksym=intWRAP(-k,0,ZSIZE(FourierWeights)-1);
        double mean=0.5*(
                        DIRECT_VOL_ELEM(FourierWeights,k,0,0)+
                        DIRECT_VOL_ELEM(FourierWeights,ksym,0,0));
        DIRECT_VOL_ELEM(FourierWeights,k,0,0)=
        DIRECT_VOL_ELEM(FourierWeights,ksym,0,0)=mean;
    }

    // Tell threads what to do
//#define DEBUG_VOL1
#ifdef DEBUG_VOL1
	{
         VolumeXmipp save;
         save().alias( FourierWeights );
         save.write((std::string) fn_out + "hermiticWeights.vol");

         FourierVolume save2;
         save2().alias( VoutFourier );
         save2.write((std::string) fn_out + "hermiticFourierVol.vol");
	}
#endif
    threadOpCode = PROCESS_WEIGHTS;
    // Awake threads
    barrier_wait( &barrier );
    // Threads are working now, wait for them to finish
    barrier_wait( &barrier );

    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);
    transformerVol.setReal(Vout());
    
    transformerVol.inverseFourierTransform();
    CenterFFT(Vout(),false);

    // Correct by the Fourier transform of the blob
    Vout().setXmippOrigin();
    Vout().window(FIRST_XMIPP_INDEX(imgSize),FIRST_XMIPP_INDEX(imgSize),
                  FIRST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize),
                  LAST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize));
    double pad_relation= ((double)padding_factor_proj/padding_factor_vol);
    pad_relation = (pad_relation * pad_relation * pad_relation);
    //THIS IS WRONG CHANGE int by size_t
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vout())
    {
        // COSS: *** Avoid the square root
        double factor = Fourier_blob_table(ROUND(sqrt((double)(k*k+i*i+j*j))
                                                 *iDeltaFourier));
        Vout(k,i,j) /= (factor/pad_relation);
    }

    Vout.write(out_name);

}
