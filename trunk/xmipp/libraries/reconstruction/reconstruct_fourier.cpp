/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "reconstruct_fourier.h"
#include <data/args.h>
#include <data/fft.h>

// Read arguments ==========================================================
void Prog_RecFourier_prm::read(int argc, char **argv)
{
    fn_sel = getParameter(argc, argv, "-i");
    fn_doc = getParameter(argc, argv, "-doc","");
    fn_out = getParameter(argc, argv, "-o", "rec_fourier.vol");
    fn_sym = getParameter(argc, argv, "-sym", "");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    do_weights = checkParameter(argc, argv, "-weight");
    padding_factor_proj = textToFloat(getParameter(argc, argv, "-pad_proj","2"));
    padding_factor_vol = textToFloat(getParameter(argc, argv, "-pad_vol","2"));
    fn_control    = getParameter(argc, argv, "-control", "");
    blob.radius   = textToFloat(getParameter(argc, argv, "-r","1.9"));
    blob.order    = textToFloat(getParameter(argc, argv, "-m","0"));
    blob.alpha    = textToFloat(getParameter(argc, argv, "-a","15"));
    sampling_rate = textToFloat(getParameter(argc, argv, "-sampling_rate", "1"));
    maxResolution = textToFloat(getParameter(argc, argv,
        "-max_resolution","2"));
    NiterWeight = textToInteger(getParameter(argc, argv, "-n","20"));
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
        if (do_weights)
            std::cerr << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cerr << " Do NOT use weights" << std::endl;
        std::cerr << " Iterations weight         : " << NiterWeight << std::endl;
        std::cerr << "\n Interpolation Function" 
                  << "\n   blrad                 : "  << blob.radius
                  << "\n   blord                 : "  << blob.order
                  << "\n   blalpha               : "  << blob.alpha
                  << "\n sampling_rate           : "  << sampling_rate
                  << "\n max_resolution          : "  << maxResolution
                  << "\n -----------------------------------------------------------------" << std::endl;
    }
}

// Usage ====================================================================
void Prog_RecFourier_prm::usage()
{

    // To screen
    std::cerr << "  Usage:\n";
    std::cerr << "  reconstruct_fourier_interpolation  <options>\n";
    std::cerr << "   -i <input selfile>          : selection file with input images \n";
    std::cerr << " [ -pad_proj <p=2.0> ]         : projection padding factor \n";
    std::cerr << " [ -pad_vol  <p=2.0> ]         : volume padding factor \n";
    std::cerr << " [ -o <name=\"rec_fourier.vol\">  : filename for output volume \n";
    std::cerr << " [ -doc <docfile>              : Ignore headers and get angles from this docfile \n";
    std::cerr << " [ -sym     <symfile> ]        : Enforce symmetry in projections\n";
    std::cerr << " [ -sym_vol <symfile> ]        : Enforce symmetry in volume \n";
    std::cerr << " [ -n <iter=20>]               : Iterations for computing the weight\n";
    std::cerr << " -----------------------------------------------------------------" << std::endl;
    std::cerr << " [ -do_weights ]               : Use weights stored in the image headers or doc file" << std::endl;
    std::cerr << "\n Interpolation Function"
              << "\n   [-r blrad=1.9]        blob radius in pixels"
              << "\n   [-m blord=0]          order of Bessel function in blob"
              << "\n   [-a blalpha=15]       blob parameter alpha"
              << "\n   [-sampling_rate =1>]            : Sampling rate (Angstroms/pixel)\n"
              << "\n   [-max_resolution 2>]            : Max resolution in "
              << "\n\t\tAngstroms, 2*sampling_rate is the maximum resolution)\n"
              << " -----------------------------------------------------------------" << std::endl;
}

void Prog_RecFourier_prm::produce_Side_info()
{
    // Translate the maximum resolution to digital frequency
    maxResolution=sampling_rate/maxResolution;
    maxResolution2=maxResolution*maxResolution;

    // Read docfile and get column numbers
    if (fn_doc != "")
    {
        DF.read(fn_doc);
        col_rot    = DF.getColNumberFromHeader("rot")  - 1;
        col_tilt   = DF.getColNumberFromHeader("tilt") - 1;
        col_psi    = DF.getColNumberFromHeader("psi")  - 1;
        col_xoff   = DF.getColNumberFromHeader("Xoff") - 1;
        col_yoff   = DF.getColNumberFromHeader("Yoff") - 1;
        col_flip   = DF.getColNumberFromHeader("Flip") - 1;
        if (do_weights)
            col_weight = DF.getColNumberFromHeader("Weight") - 1;
    }

    // Read the input images
    SF.read(fn_sel);
    SF.go_beginning();

    // Ask for memory for the output volume and its Fourier transform
    int Ydim, Xdim;
    SF.ImgSize(Ydim, Xdim);
    if (Ydim!=Xdim)
        REPORT_ERROR(1,"This algorithm only works for squared images");
    imgSize=Xdim;
    Vout().initZeros(Xdim*padding_factor_vol,Xdim*padding_factor_vol,
        Xdim*padding_factor_vol);
    transformerVol.setReal(Vout());
    transformerVol.getFourierAlias(VoutFourier);
    FourierWeights.initZeros(VoutFourier);

    // Ask for memory for the padded images
    paddedImg.resize(Xdim*padding_factor_proj,Xdim*padding_factor_proj);
    paddedImg.setXmippOrigin();
    transformerImg.setReal(paddedImg);

    // Build a table of blob values
    blob_table.resize(BLOB_TABLE_SIZE);
    Fourier_blob_table.resize(BLOB_TABLE_SIZE);

    struct blobtype blobFourier,blobnormalized;
    blobFourier=blob;
    blobFourier.radius/=(padding_factor_proj*Xdim);
    blobnormalized=blob;
    blobnormalized.radius/=((double)padding_factor_proj/padding_factor_vol);
    double delta = blob.radius/(BLOB_TABLE_SIZE-1);
    double deltaFourier = (sqrt(3.)*Xdim/2.)/(BLOB_TABLE_SIZE-1);
    
    // The interpolation kernel must integrate to 1
    double iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
    double padXdim3 = padding_factor_proj * Xdim;
    padXdim3 = padXdim3 * padXdim3 * padXdim3;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(blob_table)
    {
        DIRECT_VEC_ELEM(blob_table,i) = blob_val(delta*i, blob)  *iw0;
        DIRECT_VEC_ELEM(Fourier_blob_table,i) =
            blob_Fourier_val(deltaFourier*i, blobFourier)*padXdim3  *iw0;
    }
    iDelta=1/delta;
    iDeltaFourier=1/deltaFourier;

    // Kernel for the weight correction
    int L=CEIL(blob.radius);
    kernel.resize(2*L+1,2*L+1,2*L+1);
    kernel.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(kernel)
    {
        double r=sqrt(k*k+i*i+j*j);
        if (r>=blob.radius) kernel(k,i,j)=0;
        else kernel(k,i,j)=blob_table(ROUND(r*iDelta));
    }

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
    double &weight)
{
    if (fn_doc == "")
    {
        headerXmipp head;
        head.read(fn);
        rot    = head.Phi();
        tilt   = head.Theta();
        psi    = head.Psi();
        xoff   = head.fXoff();
        yoff   = head.fYoff();
        flip   = head.Flip();
        weight = head.Weight();
    } 
    else
    {
        if (DF.search_comment(fn))
        {
            rot    = DF(col_rot);
            tilt   = DF(col_tilt);
            psi    = DF(col_psi);
            xoff   = DF(col_xoff);
            yoff   = DF(col_yoff);
            if (col_flip < 0)
                flip   = 0.;
            else
                flip   = DF(col_flip);
            if (col_weight < 0)
                weight = 0.;
            else
                weight = DF(col_weight);
        }
        else
        {
            REPORT_ERROR(1, (std::string)"Prog_RecFourier_prm: Cannot find " + fn + " in docfile " + fn_doc);
        }
    }
}

//#define DEBUG
void Prog_RecFourier_prm::processImage(const FileName &fn_img)
{
    // Read input image
    double rot, tilt, psi, xoff,yoff,flip,weight;
    Projection proj;
    if (fn_doc == "")
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
        get_angles_for_image(fn_img, rot, tilt, psi, xoff, yoff, flip, weight);
        proj.set_rot(rot);
        proj.set_tilt(tilt);
        proj.set_psi(psi);
        proj.set_Xoff(xoff);
        proj.set_Yoff(yoff);
        proj.set_flip(flip);
        proj.set_weight(weight);
        Matrix2D<double> A;
        A = proj.get_transformation_matrix(true);
        if (!A.isIdentity())
            proj().selfApplyGeometryBSpline(A, 3, IS_INV, WRAP);
    }
    if (!do_weights) weight=1.0;
    else if (weight==0.0) return;

    // Copy the projection to the center of the padded image
    // and compute its Fourier transform
    proj().setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(proj())
        paddedImg(i,j)=weight*proj(i,j);
    CenterFFT(paddedImg,true);
    transformerImg.FourierTransform();
    Matrix2D< std::complex<double> > paddedFourier;
    transformerImg.getFourierAlias(paddedFourier);

    // Compute the coordinate axes associated to this image
    Matrix2D<double>  A(3, 3), Ainv;
    Euler_angles2matrix(rot, tilt, psi, A);
    Ainv=A.transpose();

    // Loop over all symmetries
    for (int isym = 0; isym < R_repository.size(); isym++)
    {
        // Compute the coordinate axes of the symmetrized projection
        Matrix2D<double> A_SL=Ainv * R_repository[isym];
        
        // Loop over all Fourier coefficients in the padded image
        Matrix1D<double> freq(3), gcurrent(3), real_position(3);
        Matrix1D<int> corner1(3), corner2(3);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(paddedFourier)
        {
            // Compute the frequency of this coefficient in the
            // universal coordinate system
            FFT_IDX2DIGFREQ(j,XSIZE(paddedImg),XX(freq));
            FFT_IDX2DIGFREQ(i,YSIZE(paddedImg),YY(freq));
            ZZ(freq)=0;
            if (XX(freq)*XX(freq)+YY(freq)*YY(freq)>maxResolution2)
                continue;
            SPEED_UP_temps;
            M3x3_BY_V3x1(freq,A_SL,freq);

            // Look for the corresponding index in the volume Fourier transform
            DIGFREQ2FFT_IDX_DOUBLE(XX(freq),XSIZE(VOLMATRIX(Vout)),XX(real_position));
            DIGFREQ2FFT_IDX_DOUBLE(YY(freq),YSIZE(VOLMATRIX(Vout)),YY(real_position));
            DIGFREQ2FFT_IDX_DOUBLE(ZZ(freq),ZSIZE(VOLMATRIX(Vout)),ZZ(real_position));
            
            // Put a box around that coefficient
            XX(corner1)=CEIL (XX(real_position)-blob.radius);
            YY(corner1)=CEIL (YY(real_position)-blob.radius);
            ZZ(corner1)=CEIL (ZZ(real_position)-blob.radius);
            XX(corner2)=FLOOR(XX(real_position)+blob.radius);
            YY(corner2)=FLOOR(YY(real_position)+blob.radius);
            ZZ(corner2)=FLOOR(ZZ(real_position)+blob.radius);

            #ifdef DEBUG
                std::cout << "Idx Img=(0," << i << "," << j << ") -> Freq Img=("
                          << freq.transpose() << ") ->\n    Idx Vol=("
                          << real_position.transpose() << ")\n"
                          << "   Corner1=" << corner1.transpose() << std::endl
                          << "   Corner2=" << corner2.transpose() << std::endl;
            #endif

            // Loop within the box
            double *ptrIn =(double *)&(paddedFourier(i,j));
            for (int intz = ZZ(corner1); intz <= ZZ(corner2); intz++)
                for (int inty = YY(corner1); inty <= YY(corner2); inty++)
                    for (int intx = XX(corner1); intx <= XX(corner2); intx++)
                    {
                        // Compute distance to the center of the blob
                        // Compute blob value at that distance
                        VECTOR_R3(gcurrent, intx, inty, intz);
                        V3_MINUS_V3(gcurrent, real_position, gcurrent);
                        double d = sqrt(XX(gcurrent) * XX(gcurrent) +
                                 YY(gcurrent) * YY(gcurrent) +
                                 ZZ(gcurrent) * ZZ(gcurrent));
                        if (d > blob.radius) continue;
                        double w = blob_table(ROUND(gcurrent.module()*iDelta));
                        
                        // Look for the location of this logical index
                        // in the physical layout
                        #ifdef DEBUG
                            std::cout << "   gcurrent=" << gcurrent.transpose()
                                      << " d=" << d << std::endl;
                            std::cout << "   1: intx=" << intx
                                      << " inty=" << inty
                                      << " intz=" << intz << std::endl;
                        #endif
                        int iz=intWRAP(intz,0,ZSIZE(VoutFourier)-1);
                        int iy=intWRAP(inty,0,ZSIZE(VoutFourier)-1);
                        int ix=intWRAP(intx,0,ZSIZE(VoutFourier)-1);
                        #ifdef DEBUG
                            std::cout << "   2: ix=" << ix << " iy=" << iy
                                      << " iz=" << iz << std::endl;
                        #endif
                        bool conjugate=false;
                        if (ix>=XSIZE(VoutFourier))
                        {
                            iz=intWRAP(-iz,0,ZSIZE(VoutFourier)-1);
                            iy=intWRAP(-iy,0,ZSIZE(VoutFourier)-1);
                            ix=intWRAP(-ix,0,ZSIZE(VoutFourier)-1);
                            conjugate=true;
                        }
                        #ifdef DEBUG
                            std::cout << "   3: ix=" << ix << " iy=" << iy
                                      << " iz=" << iz << " conj="
                                      << conjugate << std::endl;
                        #endif
                        
                        // Add the weighted coefficient
                        double *ptrOut=(double *)&(VoutFourier(iz,iy,ix));
                        ptrOut[0]+=w*ptrIn[0];
                        if (conjugate) ptrOut[1]-=w*ptrIn[1];
                        else           ptrOut[1]+=w*ptrIn[1];
                        FourierWeights(iz,iy,ix)+=w;
                    }
        }
    }
}
#undef DEBUG

// Correct weight ----------------------------------------------------------
void Prog_RecFourier_prm::correctWeight()
{
    FourierWeights.printStats(); std::cout << std::endl;
    FourierWeightsConvolved.initZeros(FourierWeights);

    // Compute the convolution of the kernel with the weights
    FOR_ALL_ELEMENTS_IN_MATRIX3D(kernel)
    {
        double wkij=kernel(k,i,j);
        if (wkij==0) continue;
        for (int kp=STARTINGZ(FourierWeights); kp<=FINISHINGZ(FourierWeights); kp++)
        {
            int iz=intWRAP(kp-k,0,ZSIZE(VoutFourier)-1);
            for (int ip=STARTINGY(FourierWeights); ip<=FINISHINGY(FourierWeights); ip++)
            {
                int iy=intWRAP(ip-i,0,ZSIZE(VoutFourier)-1);
                for (int jp=STARTINGX(FourierWeights); jp<=FINISHINGX(FourierWeights); jp++)
                {
                    int ix=intWRAP(jp-j,0,ZSIZE(VoutFourier)-1);
                    if (ix>=XSIZE(VoutFourier))
                    {
                        iz=intWRAP(-iz,0,ZSIZE(VoutFourier)-1);
                        iy=intWRAP(-iy,0,ZSIZE(VoutFourier)-1);
                        ix=intWRAP(-ix,0,ZSIZE(VoutFourier)-1);
                    }
                    FourierWeightsConvolved(iz,iy,ix)+=wkij*
                        VOL_ELEM(FourierWeights,kp,ip,jp);
                }
            }
        }
    }

/*
    VolumeXmipp save;
    save()=FourierWeights; 
    FOR_ALL_ELEMENTS_IN_MATRIX3D(save())
        save(k,i,j)=log10(save(k,i,j)+1);
    save.write("PPPweights.vol");
    save()=FourierWeightsConvolved;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(save())
        save(k,i,j)=log10(save(k,i,j)+1);
    save.write("PPPweightsConvolved.vol");
    std::cout << "kernel=\n" << kernel << std::endl;
    std::cout << "Press any key\n";
    char c; std::cin >> c;
*/
    
    // Update the weights with the convolved values
    FOR_ALL_ELEMENTS_IN_MATRIX3D(FourierWeights)
        if (FourierWeightsConvolved(k,i,j)>XMIPP_EQUAL_ACCURACY)
            FourierWeights(k,i,j)/=FourierWeightsConvolved(k,i,j);
}

// Main routine ------------------------------------------------------------
void Prog_RecFourier_prm::run()
{
    // Process all images in the selfile
    if (verb) init_progress_bar(SF.ImgNo());
    int imgno = 0;
    int repaint = ceil(SF.ImgNo()/60);
    while (!SF.eof())
    {
        exit_if_not_exists(fn_control);

        FileName fn_img = SF.NextImg();
        if (fn_img=="") break;
        processImage(fn_img);

        if (verb && imgno++%repaint==0) progress_bar(imgno);
    }
    if (verb > 0) progress_bar(SF.ImgNo());
}

void Prog_RecFourier_prm::finishComputations()
{
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

    // Correct weights
    std::cerr << "Correcting the weight ...\n";
    FOR_ALL_ELEMENTS_IN_MATRIX3D(FourierWeights)
        if (FourierWeights(k,i,j)>XMIPP_EQUAL_ACCURACY)
            FourierWeights(k,i,j)=1/FourierWeights(k,i,j);
/*
    verb=false;
    if (verb) init_progress_bar(NiterWeight);
    for (int n=0; n<NiterWeight; n++)
    {
        correctWeight();
        if (verb) progress_bar(n);
    }
    if (verb > 0) progress_bar(NiterWeight);
*/
    // Get a first approximation of the reconstruction
    double corr2D_3D=pow(padding_factor_proj,2.)/
                        (imgSize* pow(padding_factor_vol,3.)); 
                                   // Divide by Zdim because of the
                                   // the extra dimension added
                                   // and padding differences
    FOR_ALL_ELEMENTS_IN_MATRIX3D(FourierWeights)
        if (VOL_ELEM(FourierWeights,k,i,j)!=0)
            VOL_ELEM(VoutFourier,k,i,j)*=corr2D_3D*VOL_ELEM(FourierWeights,k,i,j);
//    std::cout << "Aqui2\n" << FourierWeights << std::endl;
//    std::cout << "Aqui3\n" << VoutFourier << std::endl;
    transformerVol.inverseFourierTransform();
    CenterFFT(Vout(),false);

    // Correct by the Fourier transform of the blob
    Vout().setXmippOrigin();
    Vout().window(FIRST_XMIPP_INDEX(imgSize),FIRST_XMIPP_INDEX(imgSize),
        FIRST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize),
        LAST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize));
    double pad_relation= ((double)padding_factor_proj/padding_factor_vol);
    pad_relation = (pad_relation * pad_relation * pad_relation);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vout())
    {
        double factor = Fourier_blob_table(ROUND(sqrt(k*k+i*i+j*j)
                                                *iDeltaFourier));
        Vout(k,i,j) /= (factor/pad_relation);
    }
    Vout.write(fn_out);
}
