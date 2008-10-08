/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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

#include "reconstruct_fourier.h"

// Read arguments ==========================================================
void Prog_RecFourier_prm::read(int argc, char **argv)
{

    fn_sel = getParameter(argc, argv, "-i");
    fn_doc = getParameter(argc, argv, "-doc","");
    fn_out = getParameter(argc, argv, "-o", "rec_fourier.vol");
    fn_sym = getParameter(argc, argv, "-sym", "");
    fn_sym_vol = getParameter(argc, argv, "-sym_vol", "");
    do_resolution = checkParameter(argc, argv, "-do_resolution");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    do_weights = checkParameter(argc, argv, "-weight");
    padding_factor_proj = textToFloat(getParameter(argc, argv, "-pad_proj","2"));
    //padding_factor_vol  = textToFloat(getParameter(argc, argv, "-pad_vol" ,"2"));
    padding_factor_vol=padding_factor_proj;    
    // For improved killing control
    fn_control = getParameter(argc, argv, "-control", "");
    blob.radius        = textToFloat(getParameter(argc, argv,  "-r","2.0"));
    blob.order         = textToFloat(getParameter(argc, argv,    "-m","0"));
    blob.alpha         = textToFloat(getParameter(argc, argv,   "-a","15"));
    sampling_rate      = textToFloat(getParameter(argc, argv, "-sampling_rate", "1"));
    maxResolution      = textToFloat(getParameter(argc, argv,
                                                    "-max_resolution","2"));
    maxResolution_normalize =  1/(maxResolution/sampling_rate);
    if(maxResolution_normalize > 0.5000001)
    {
        std::cerr << "this resolution cannot be obtained with this sampling rate";
        exit(0);
    }    
}

// Show ====================================================================
void Prog_RecFourier_prm::show()
{

    if (verb > 0)
    {
        // To screen
        std::cerr << " =================================================================" << std::endl;
        std::cerr << " Direct 3D reconstruction method using Kaiser windows as interpolators" << std::endl;
        std::cerr << " =================================================================" << std::endl;
        std::cerr << " Input selfile             : "  << fn_sel << std::endl;
        std::cerr << " padding_factor_proj       : "  << padding_factor_proj << std::endl;
        std::cerr << " padding_factor_vol        : "  << padding_factor_vol << std::endl;
        if (fn_doc != "")
            std::cerr << " Input docfile         : "  << fn_doc << std::endl;
        std::cerr << " Output volume             : "  << fn_out << std::endl;
        if (fn_sym != "")
            std::cerr << " Symmetry file for projections : "  << fn_sym << std::endl;
        if (fn_sym_vol != "")
            std::cerr << " Symmetry file for volume      : "  << fn_sym_vol << std::endl;
        if (do_resolution)
            std::cerr << " Compute resolution" << std::endl;
        else
            std::cerr << " Do NOT compute resolution" << std::endl;
        if (do_weights)
            std::cerr << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cerr << " Do NOT use weights" << std::endl;
        std::cerr << "\n Interpolation Function" 
                  << "\n   blrad                 : "  << blob.radius
                  << "\n   blord                 : "  << blob.order
                  << "\n   blalpha               : "  << blob.alpha
                  << "\n sampling_rate           : "  << sampling_rate
                  << "\n max_resolution          : "  << maxResolution
       
                  << "\n -----------------------------------------------------------------" << std::endl;
    }
// NOTE about blob parameters: following Matej et al.IEEE trans medical imaging
// 25(7)845 año 2006 the above parameters ar optimal for padding from 0 to 100%
//* do not multiply by padding factor
}

// Usage ====================================================================
void Prog_RecFourier_prm::usage()
{

    // To screen
    std::cerr << "  Usage:\n";
    std::cerr << "  reconstruct_fourier_interpolation  <options>\n";
    std::cerr << "   -i <input selfile>          : selection file with input images \n";
    std::cerr << "   -pad_proj <2.0>             : projection padding factor \n";
    std::cerr << " [ -o <name=\"rec_fourier.vol\">       : filename for output volume \n";
    std::cerr << " [ -doc <docfile>              : Ignore headers and get angles from this docfile \n";
    std::cerr << " [ -sym     <symfile> ]        : Enforce symmetry in projections\n";
    std::cerr << " [ -sym_vol <symfile> ]        : Enforce symmetry in volume \n";
    std::cerr << " [ -do_resolution]             : compute resolution while you reconstruct \n";
    std::cerr << " -----------------------------------------------------------------" << std::endl;
    if (do_weights)
        std::cerr << " --> Use weights stored in the image headers or doc file" << std::endl;
    std::cerr << "\n Interpolation Function"
              << "\n   [-r blrad=2.0]        blob radius in pixels"
              << "\n   [-m blord=0]          order of Bessel function in blob"
              << "\n   [-a blalpha=15]       blob parameter alpha"
              << "\n   [-sampling_rate =1>]            : Sampling rate (Angstroms/pixel)\n"
              << "\n   [-max_resolution 2>]            : Max resolution in "
              << "\n\t\tAngstroms, 2*sampling_rate is the maximum resolution)\n"

              << " -----------------------------------------------------------------" << std::endl;
}

void Prog_RecFourier_prm::produce_Side_info()
{
    double dum, weight;
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

    //remove images with weight=0
    if (do_weights)
    {
        SelFile SF_aux;
        SF_aux.read(fn_sel);
        SF_aux.go_beginning();
        while (!SF_aux.eof())
        {
            get_angles_for_image(SF_aux.get_current_file(), dum, dum, dum, dum, dum, dum, weight);
            if (weight != 0)
            {
                SF.insert(SF_aux.current());
            }
            SF_aux.NextImg();
        }
        if (SF.ImgNo() == 0)
        {
            std::cerr << "there is no input file with weight!=0" << std::endl;
            exit(1);
        }
    }
    else
        SF.read(fn_sel);

    SF.ImgSize(dim, dim);
    if (fn_sym != "")     SL.read_sym_file(fn_sym);
    if (fn_sym_vol != "") SL_vol.read_sym_file(fn_sym_vol);
    //precompute blob matrix
    blob_table = new double [BLOB_TABLE_SIZE];
    //convert index to pixel
    //factor 3. I do not know why
    double ww = blob.radius/((double)BLOB_TABLE_SIZE-1);
    //ww /= ((double)dim*padding_factor_proj);
    //acces to the blob table should take into account 
    //pad factor but not sampling
    double w;
    for (int i=0; i<BLOB_TABLE_SIZE; i++)
    {
      w = ww*(double)i;
      blob_table[i] =  blob_val(w, blob);
      //#define DEBUG
      #ifdef DEBUG
      std::cout.setf(std::ios::scientific); 
      std::cout.width(12); 
      std::cout << i << " " 
                << w <<  "\t" 
                << "\t"<< blob_table[i] << " "
                <<  blob.radius
                << std::endl;
      #endif
      #undef DEBUG
    }
    fill_L_R_repository(SL);
}

void Prog_RecFourier_prm::fill_L_R_repository(SymList & mySL)
{
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    R_repository.push_back(Identity);
    //L_repository.push_back(Identity);
    for (int isym = 0; isym < mySL.SymsNo(); isym++)
    {
        mySL.get_matrices(isym, L, R);
        R.resize(3, 3);
        //L.resize(3, 3);
        R_repository.push_back(R);
        //L_repository.push_back(L);
    }
//#define DEBUG3
#ifdef  DEBUG3
    for (int isym = 0; isym < R_repository.size(); isym++)
        {
        std::cout << R_repository[isym];
        //std::cout << L_repository[isym];
        }
#endif
#undef DEBUG3
}


void Prog_RecFourier_prm::get_angles_for_image(FileName fn, double &rot, double &tilt, double &psi,
                                        double &xoff, double &yoff, double &flip, double &weight)
{
    if (fn_doc == "")
    {
        headerXmipp      head;
        head.read(fn);
        rot    = head.Phi();
        tilt   = head.Theta();
        psi    = head.Psi();
        xoff   = head.fXoff();
        yoff   = head.fYoff();
        flip   = head.Flip();
        weight = head.Weight();
    //#define ANGLES
    #ifdef ANGLES
    std::cerr << "\nFrom header\n";
    #endif
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
            #ifdef ANGLES
            std::cerr << "\nFrom doc " << DF << "col_rot " << col_rot << " " << DF(col_rot) << "\n";
            #endif
        }
        else
        {
            REPORT_ERROR(1, (std::string)"Prog_RecFourier_prm: Cannot find " + fn + " in docfile " + fn_doc);
        }
    }
    //#define ANGLES
    #ifdef ANGLES
    std::cerr << "\nrot, tilt,psi " << rot << " " << tilt << " " << psi << "\n";
    #endif
    #undef ANGLES
    
}

void Prog_RecFourier_prm::ProcessOneImage(FileName &fn_img, 
                                          xmippFftw &fftPaddedImg,//esta no hace falta
                                          int paddim_proj,
                                          int paddim_vol,
                                          Projection &proj,
                                          std::complex<double> * FOURIERVOL,
                                          double * FourierVolWeight,
                                          std::complex<double> * FOURIERPROJ)
{
    double           rot, tilt, psi, xoff,yoff,flip,weight;
    //APPLY MIRROR WHEN
    if (fn_doc == "")
    {
        proj.read(fn_img, true);//true means apply shifts 
        rot  = proj.rot();
        tilt = proj.tilt();
        psi  = proj.psi();
        weight = proj.weight();
    }
    else
    {
        proj.read(fn_img, false);//do not apply shifts since they are not in the header
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
        //apply  shifts only
        if (!A.isIdentity())
            proj().selfApplyGeometryBSpline(A, 3, IS_INV, WRAP);
    }
    //pad
    //copy image to fft structure
    //I do not think is important to center the image
    //so padding will not be made at image center
    //note the transpose problem
    int Xdim=dim;//NOT padded image
    int Ydim=dim;
    int center_shift = (paddim_proj-Xdim)/2+Xdim%2;// add one if odd
    //#define DEBUG_PADD
    #ifdef DEBUG_PADD
    {
        ImageXmipp test(proj);
        for(int i=0;i<Ydim;i++)
           for(int j=0;j<Xdim;j++)
               {                                          //x,y
               test(j,i)=proj(j,i);
               }
       test.write("image.xmp");
    }
    ////std::cerr << "center_shift: " << center_shift << std::endl;
    std::cerr << "paddim_proj: "  << paddim_proj << std::endl;
    #endif
    for(int i=0;i<Ydim;i++)
        for(int j=0;j<Xdim;j++)
            {                                          //x,y
            fftPaddedImg.fIn[(center_shift+j) + 
                             (center_shift+i) * paddim_proj]=(double)proj(i,j) / 
                             fftPaddedImg.fTotalSize;//x,y
////            fftPaddedImg.fIn[(j)                  + 
////                             (i)  * paddim_proj]=(double)proj(j,i) / 
////                             fftPaddedImg.fTotalSize;//x,y;
            }
    //#define DEBUG_PADD
    #ifdef DEBUG_PADD
    {
        ImageXmipp test(paddim_proj,paddim_proj);
        for(int i=0;i<paddim_proj;i++)
           for(int j=0;j<paddim_proj;j++)
               {                                          //x,y
               test(j,i)=fftPaddedImg.fIn[j +i * paddim_proj];
               }
       test.write("padded_image.xmp");
    }
    #endif
    //#ifdef MYNEW
    fftPaddedImg.CenterFourierTransformInRealSpace(false);//ok MOVE TRANSFORM TO THE CENTER
    //#else
    //fftPaddedImg.CenterRealDataBeforeTransform();
    //#endif
    //#define DEBUG_PADD_FFT
    #ifdef DEBUG_PADD_FFT
    {
        ImageXmipp test(paddim_proj,paddim_proj);
        int ysize = YSIZE(test());
        int xsize = XSIZE(test());
        int ii=0;
        for(int i=0;i<xsize;i++)
            for(int j=0;j<ysize;j++)
                {                    //x,y
                DIRECT_MAT_ELEM(test(),j,i) = fftPaddedImg.fIn[ii];
                ii++;
                }
        test.write("shift_padded_image.xmp");
    }
    #endif
    fftPaddedImg.Transform();

    //#ifdef MYNEW
    fftPaddedImg.CenterRealImageInFourierSpace(true) ;//Change phases 
    //#endif

    //#define DEBUG_PADD_FFT
    #ifdef DEBUG_PADD_FFT
    {
        ImageXmipp test(paddim_proj,paddim_proj);
        std::complex<double> * IMG;
        IMG = (std::complex<double> *)fftPaddedImg.fOut;
        int ysize = (int)(((double)YSIZE(test()) / 2.) +1.);
        int xsize = XSIZE(test());
        int ii=0;
        for(int i=0;i<xsize;i++)
            for(int j=0;j<ysize;j++)
                {                    //x,y
                DIRECT_MAT_ELEM(test(),i,j) = log(1+abs(IMG[ii]));
                //std::cerr<< i << " " << j << " " << IMG[ii] << "\n";  
                ii++;
                }
        test.write("shift_padded_image.fft");
    }
    #endif
   
    
        
    Matrix2D<double>  A(3, 3);
    Matrix2D<double>  A_inv(3, 3);
    Matrix2D<double>  A_SL(3, 3);
    //APPLY MIRROR WHEN: if flip activated then selfApplyGeometryBSpline will apply mirror
    //compute euler matrix
    Euler_angles2matrix(rot, tilt, psi, A);
    //A_inv = A.inv();
    for (int isym = 0; isym < R_repository.size(); isym++)
        {
        A_SL=  A * R_repository[isym] ;
        
        //std::cout << "A*R" << A_SL << A_SL.det() << std::endl;
        ////std::cout << "R"   << R_repository[isym];
        // L is only needed for h cases  L_repository[isym];
        if(do_weights==false)
             weight=1.0;
        placeInFourerVolumeOneImage(A_SL,//3*3 euler matrix
                                    FOURIERVOL,//vector
                                    fftPaddedImg,//fft object
                                    FourierVolWeight,//vector
                                    paddim_proj,
                                    paddim_vol,
                                    FOURIERPROJ,
                                    weight
                                    );
        }
}

//place one projection in the 3D volume (Fourier space )
void Prog_RecFourier_prm::placeInFourerVolumeOneImage(
                                     Matrix2D<double> &A,//3*3 euler matrix
                                     std::complex<double> * FOURIERVOL,//vector
                                     xmippFftw &fftPaddedImg,//fft object
                                     double * FourierVolWeight,//vector
                                     int paddim_proj,
                                     int paddim_vol,
                                     std::complex<double> * FOURIERPROJ,
                                     double proj_weight)
{
    //dim, no padded proj dimension, vol dimension
    int dim2 = dim / 2;
    //remember size transposition
    
    int Xdim,Xsize;
    int Ydim=1,Ysize=1;
    int Zdim=1,Zsize=1;
    Zdim=Zsize=paddim_vol;
    Ydim=Ysize=paddim_vol;
    Xsize=paddim_vol;
    Xdim = (int)  paddim_vol/2  +1;
    int XCenter=Xsize/2;
    int YCenter=Ysize/2;
    int ZCenter=Zsize/2;


    
    int Udim,Usize;
    int Vdim,Vsize;
    Vdim=Vsize=paddim_proj;
    Usize=paddim_proj;
    Udim = (int)  paddim_proj/2 +1;
    int UCenter=Usize/2;
    int VCenter=Vsize/2;
    int uu, vv;
    double su2, sv2,s;
    double paddedMaxResolution2 = maxResolution_normalize;
    paddedMaxResolution2 *= paddedMaxResolution2;
    
    //create a 3D vector with dimensions
    Matrix1D<double> _3Dposition(3),_2Dposition(3);
    double xp,yp,zp;
    //note that this only works with cubic vols
    //#define DEBUG
    #ifdef DEBUG
    std::cerr << "A" << A(0,0) << " " << A(1,0) << A << std::endl ;
    #endif
    #undef DEBUG

    //#define CHIMERA
    #ifdef CHIMERA
    static std::ofstream filestr;
    static int my_open=0;
    if(filestr.is_open()==false)
        filestr.open ("fft_volume.bild");
    filestr << ".color yellow\n" ;
    int my_counter2=0;
    #endif
    //#define CHIMERA2
    #ifdef CHIMERA2
    static std::ofstream filestr2;
    static int my_open2=0;
    if(filestr2.is_open()==false)
        filestr2.open ("fft_volume_rectified.bild");
    filestr2 << ".color magenta\n" ;
    int my_counter22=0;
    #endif
    //#define CHIMERA3
    #ifdef CHIMERA3
    static std::ofstream filestr3;
    static int my_open3=0;
    if(filestr3.is_open()==false)
        filestr3.open ("arround_volume_point.bild");
    filestr3 << ".color magenta\n" ;
    int my_counter23=0;
    #endif

    double padding_ratio=padding_factor_proj/padding_factor_vol;
    //matrix convert from padded vol to padded projection
    A = A / padding_ratio;
//std::cerr << "A Vdim Udim" << A << Vdim << " " << Udim << std::endl;    
    for ( int v=0,i=0; v<Vdim; v++ ) {
	    //if ( v > (Vsize - 1)/2 ) vv = v-Vsize;
	    //else vv = v;
        vv = v -VCenter;
        sv2 = (double)vv/Vsize;
	    sv2 *= sv2;
	    for ( int u=0; u<Udim; u++, i++ ) {
		    //if ( u > (Usize - 1)/2 ) uu = u-Usize;
		    //else uu = u;
            uu = u -UCenter;
            su2 = (double)uu/Usize;
		    su2 *= su2;
		    s = (su2 + sv2);
            if (s >= paddedMaxResolution2) 
				continue;
            xp = (double)vv * A(1,0) + (double)uu * A(0,0); 
            yp = (double)vv * A(1,1) + (double)uu * A(0,1);
            zp = (double)vv * A(1,2) + (double)uu * A(0,2);
            #ifdef CHIMERA
            {
                if(my_counter2%2==0)
                {
                    filestr << ".sphere " <<  xp<< " " <<
                                            yp<< " " <<
                                            zp<< " " /*<<
                                            u << " " <<
                                            v << " " <<
                                            FOURIERPROJ[i]*/ 
                                            << " 0.5\n" ;

                //    filestr << ".sphere " <<  -xp<< " " <<
                //                            -yp<< " " <<
                //                            -zp<< " " /*<<
                //                            u << " " <<
                //                            v << " " <<
                //                            FOURIERPROJ[i]*/ 
                //                            << " 0.25\n" ;
                }
                my_counter2++;
            }
            #endif
            bool conjugate_flag;
            conjugate_flag=false;
            /*as we are inside a sphere you cannot be outside 
              a circle of diameter dim */
            if(xp > 0)//x is short dim  
            {
///
///continue;
///
               conjugate_flag=true;
               xp = -xp;
               yp = -yp;
               zp = -zp;
            }

            #ifdef CHIMERA2
            {
                if(my_counter22%2==0)
                {
                    /*
                    if(conjugate_flag==true)
                        filestr2 << ".color green" << "\n";
                    else
                        filestr2 << ".color red" << "\n";
                    */
                    filestr2 << ".sphere " <<  xp << " " <<
                                               yp << " " <<
                                               zp << " " 
                                                  << " 0.5\n" ;
                    if(xp==0)
                    filestr2 << ".sphere " <<  xp << " " <<
                                               -yp << " " <<
                                               -zp << " " 
                                                  << " 0.75\n" ;
                                                  

                //    std::cout << ".sphere " <<  -xp<< " " <<
                //                                -yp<< " " <<
                //                                -zp<< " " /*<<
                //                                 u << " " <<
                //                                 v << " " <<
                //                            FOURIERPROJ[i]*/ 
                //                            << " 0.25\n" ;
                }
                my_counter22++;
            }
            #endif
            //LOOP AROUND POINT
            int xx, yy, zz;
            int xxCentered,yyCentered,zzCentered;
            double r2,r;
            double BRadius=blob.radius;
            double BRadius2=BRadius*BRadius;
            double d_zp, d_yp, d_xp;
            int xx0 = int(xp-BRadius);
            int yy0 = int(yp-BRadius);
            int zz0 = int(zp-BRadius);
            int xxF = int(xp+BRadius);
            int yyF = int(yp+BRadius);
            int zzF = int(zp+BRadius);
            bool my_conjugate_flag;
/*
SEGMENTATION FAULT if symmetry one.sel
i3 does no symmetrice properlly z OK            
*/
            for (zz=zz0; zz <= zzF; zz++)
                for (yy=yy0; yy <= yyF; yy++)
                   for (xx=xx0; xx <= xxF; xx++)
//zz=int(zp);
//yy=int(yp);
//xx=int(xp);
                   {
                       my_conjugate_flag = conjugate_flag;
                       xxCentered = xx+XCenter;
                       //border of x axis

                       if(xxCentered>=Xdim)
                       {
///
//continue;
///                       
                           xxCentered = XCenter-xx;
                           yyCentered = YCenter-yy;
                           zzCentered = ZCenter-zz;
                           my_conjugate_flag = !my_conjugate_flag;
                       } 
                       else
                       {
                           yyCentered = yy+YCenter;
                           zzCentered = zz+ZCenter;
                       }
                       if (  xxCentered < 0     || 
                             yyCentered < 0     || 
                             zzCentered < 0     || 
                             /*xxCentered >= Xdim || */
                             yyCentered >= Ydim  || 
                             zzCentered >= Zdim )//Xdim half 
                           continue;
//OUT SPHERE


                       r2=    (xp-(double)xx)*(xp-(double)xx)+
                              (yp-(double)yy)*(yp-(double)yy)+
                              (zp-(double)zz)*(zp-(double)zz);//distance from projection grid point to volume grid point projection
                       if (  r2 >= BRadius2)
                           continue;
                       #ifdef CHIMERA3
                       {
                               if(u==63 && v== 80)
                               {
                                    if(my_conjugate_flag==true)
                                     {
                                        filestr3 << ".color green" << "\n";
                                        filestr3 << ".sphere " <<  xx<< " " <<
                                                            yy<< " " <<
                                                            zz<< " " 
                                                            << " 0.5\n" ;
                                     }
                                    else
                                     {   
                                        filestr3 << ".color red" << "\n";
                                        filestr3 << ".sphere " <<  -xx<< " " <<
                                                            -yy<< " " <<
                                                            -zz<< " " 
                                                            << " 0.5\n" ;
                                     }
                                }                            

                       }
                       #endif
                       r = sqrt(r2);
                       double weight = blob_table[(int)(r*BLOB_TABLE_SIZE/BRadius)];
                       weight *= proj_weight;
                       /*
                       if(weight < MINIMUMWEIGHT)//max is 1
                          continue;
                       */
                       std::complex<double> aux_complex;
                       if(my_conjugate_flag)
                           aux_complex = conj(FOURIERPROJ[v*Udim+u]);//Udimshort
                       else
                           aux_complex =      FOURIERPROJ[v*Udim+u];
                       FOURIERVOL[xxCentered+
                                  yyCentered*Xdim+
                                  zzCentered*Xdim*Ydim] += aux_complex*weight;// aux_complex*weight;//may be
                                                           //transpose
                       
                       FourierVolWeight[xxCentered+
                                        yyCentered*Xdim+
                                        zzCentered*Xdim*Ydim] += weight;
                       if(xx==0)//the line xx=0 mut be antiymmetric
                       {
                           if(my_conjugate_flag)
                           {
                              yyCentered = +yy+YCenter;
                              zzCentered = +zz+ZCenter;
                           }
                           else
                           {
                              yyCentered = -yy+YCenter;
                              zzCentered = -zz+ZCenter;
                           }
                           aux_complex = conj(aux_complex);
                           if (  yyCentered >= Ydim  || 
                                 zzCentered >= Zdim ) 
                               continue;

                           FOURIERVOL[xxCentered+
                                      yyCentered*Xdim+
                                      zzCentered*Xdim*Ydim] += aux_complex*weight;// aux_complex*weight;//may be
                                                               //transpose

                           FourierVolWeight[xxCentered+
                                            yyCentered*Xdim+
                                            zzCentered*Xdim*Ydim] += weight;

                       }
                       //LINE x=0 and warp arround this line
                       
                       //#define DEBUG
                       #ifdef DEBUG
                       {
         
                           if(xxCentered+
                              yyCentered*Xdim+
                              zzCentered*Xdim*Ydim==66949)
                               std::cerr << "FOURIERVOL FourierVolWeight "
                                         << FOURIERVOL[xxCentered+
                                                       yyCentered*Xdim+
                                                       zzCentered*Xdim*Ydim]
                                         << " "               
                                         << FourierVolWeight[xxCentered+
                                                             yyCentered*Xdim+
                                                             zzCentered*Xdim*Ydim]              
                                         << "\n";
                       }
                       #endif
                       #undef DEBUG
                       #ifdef DEBUG_CURIOSITY
                       my_counter++;
                       #endif
                   }//xx,yy,zz
        }//u
    }//v
  
    #ifdef CHIMERA
    //filestr.close();
    #endif
    #undef CHIMERA   
//We need to divide by the interpolation kernel but that may wait;
//A normalization 2D - 3D by volume size is needed
//symmetrization
//#endif
}
// Main Loop ======================================
void Prog_RecFourier_prm::MainLoop(VolumeXmipp &vol)
{

    int               c, nn, imgno;
    double            rot, tilt, psi, newrot, newtilt, newpsi, xoff, yoff, flip, weight;
    Projection        proj;
    Matrix2D<double>  L(4, 4), R(4, 4), A;
    Mask_Params       mask_prm;
    FileName          fn_img;            

    //vol().resize(dim, dim, dim);
    //vol().initZeros();
    //alloc memory for volume in Fourier space
    //the allocation is a bit wier but can be reused for fftw
    //I do not want to create a fourier object yet
    //because I only need it at the end of the program
    //and requires to alloc the doble amount of memory
    double * FourierVol;
    double * FourierVolWeight;
    int paddim_vol=ROUND(dim*padding_factor_vol);
    if(paddim_vol%2!=0)
        REPORT_ERROR(1, "\n\nvolume dim * padding_factor must be even.\n\n");
    padding_factor_vol=(double)paddim_vol/(double)dim;
    int sizeout;
    {
        int ndim    = 3;
        int * myfN;
        try
        {
            myfN = new int [ndim];
        }
        catch (std::bad_alloc&)
        {
          std::cout << "Error allocating memory." << std::endl;
          exit(1);
        }    
        myfN[0]=paddim_vol;
        myfN[1]=paddim_vol;
        myfN[2]=paddim_vol;
        // Set padding dimension
        int fTotalSize=myfN[0]*myfN[1]*myfN[2];
        sizeout = int(double(fTotalSize)*(int)(myfN[ndim-1]/2+1)/myfN[ndim-1]);
        try
        {
            FourierVol = new double[2*sizeout];
        }
        catch (std::bad_alloc&)
        {
          std::cout << "Error allocating memory." << std::endl;
          exit(1);
        }
        for (int i=0; i< 2*sizeout; i++) FourierVol[i]=0.;
        
        //we will need a volume for weight with same size in x and y and half in z
        try
        {                                 //the 2 factor is not needed
            FourierVolWeight = new double[sizeout];
        }
        catch (std::bad_alloc&)
        {
            std::cout << "Error allocating memory." << std::endl;
            exit(1);
        } 
    }
    //
    // init volumes with zeros
    std::complex<double> * FOURIERVOL;
    FOURIERVOL = (std::complex<double> *)FourierVol;
    for (int i=0; i<sizeout;i++)
    {
        FOURIERVOL[i]=(std::complex<double>)0.0;
        FourierVolWeight[i]=0.;
    }


    SF.go_beginning();
    imgno = 0;
    //create Fourier object of padded size
    //I want to reuse it so I will create it by hand
    int paddim_proj=ROUND(dim*padding_factor_proj);
    if(paddim_proj%2!=0)
        REPORT_ERROR(1, "\n\nprojection dim * padding_factor must be even.\n\n");
    padding_factor_proj=(double)paddim_proj/(double)dim;
    ImageXmipp paddedImg(paddim_proj,paddim_proj);
    paddedImg().initZeros();
    xmippFftw fftPaddedImg(paddedImg(), false);
    fftPaddedImg.Init("ES",FFTW_FORWARD,false);
    
    // pointer for output matrix
    std::complex<double> * FOURIERPROJ;
    FOURIERPROJ = (std::complex<double> *)fftPaddedImg.fOut;
    int i_counter=0;
    int screen_size=XMIPP_MIN(60,SF.ImgNo());
    int step_counter=(int)(SF.ImgNo()/screen_size);
    if (verb)
    {
        if (step_counter<=0) step_counter=1;
        init_progress_bar(screen_size);
    }

    while (!SF.eof())
    {
	// Check whether to kill job
        exit_if_not_exists(fn_control);
        fn_img = SF.NextImg();
        if (SF.eof()) break;

        ProcessOneImage(fn_img,
                        fftPaddedImg,
                        paddim_proj,
                        paddim_vol,
                        proj,
                        FOURIERVOL,
                        FourierVolWeight,
                        FOURIERPROJ);

        imgno++;
        if (verb && imgno%step_counter==0)
        {
            progress_bar(i_counter++);
        }
    }
    if (verb > 0) progress_bar(nn);
    //normalize and save volume
    int Xdim,Xsize;
    int Ydim=1,Ysize=1;
    int Zdim=1,Zsize=1;
    Zdim=Zsize=paddim_vol;
    Ydim=Ysize=paddim_vol;
    Xsize=paddim_vol;
    Xdim = (int)  paddim_vol/2  +1;
    //#define DEBUG_VOL
    #ifdef DEBUG_VOL
    {
        double aux_double;
        VolumeXmipp test(paddim_vol,
                         paddim_vol,
                         paddim_vol);
        for(int i=0;i<Zdim;i++)
           for(int j=0;j<Ydim;j++)
               for(int k=0;k<Xdim;k++)
               {
               aux_double  = abs(FOURIERVOL[k + j*Xdim +i*Ydim*Xdim ]) ; //x,y
               if (aux_double > 0)
                   test(i,j,k)=log(1+aux_double);
               else
                   test(i,j,k)=0.;
               }
       test.write("out_before.vol");
    }
    #endif
    #undef VOL
    /*Let us do some magic now. 
    what we have implement (afer dividing by FourierVolWeight[i] is a 
    Nadaraya-Watson kernel regression as described in 
    http://en.wikipedia.org/wiki/Kernel_regression
    
    I do not know why but there is a clear dumping efect in the reconstructed
    volume as we move from the center.
    
    I am going to divide by the FT¯1 of the kernel function
    since this seems to fix this problem.
    This division makes sense if you use a gridding technique
    but  simply I can not justify it in a kernel regresion approach.
    
    By the way the difference between kernel regression
    and gridding is that in the second case you do not divide by
    the sum of the weight and you DO multiply
    by a term related with the voronoi region of each projection
    Fourier transform term (that is, take into account the Fourier
    space density).
    
    spider in its command bp 3f has implement this same incongruence
    
    */
    
    double dZdim = (double) Zdim;
    for ( int z=0, i=0; z<Zdim; z++ ) 
		for ( int y=0; y<Ydim; y++ ) 
			for ( int x=0; x<Xdim; x++, i++ ) 
            {
                if( FourierVolWeight[i] > 0)//MINIMUMWEIGHT)
                     FOURIERVOL[i] /= FourierVolWeight[i] * dZdim ;
            }
    //#define DEBUG_VOL
    #ifdef DEBUG_VOL
    {
        double aux_double;
        VolumeXmipp test(paddim_vol,
                         paddim_vol,
                         paddim_vol);
        for(int i=0;i<Zdim;i++)
           for(int j=0;j<Ydim;j++)
               for(int k=0;k<Xdim;k++)
               {                                          
               aux_double  = abs(FOURIERVOL[k + j*Xdim +i*Ydim*Xdim ]) ; //x,y
               if (aux_double > 0)
                   test(i,j,k)=log(1+aux_double);
               else
                   test(i,j,k)=0.;
               }
       test.write("out_after.vol");
        for(int i=0;i<Zdim;i++)
           for(int j=0;j<Ydim;j++)
               for(int k=0;k<Xdim;k++)
               {                                          
               aux_double  = FourierVolWeight[k + j*Xdim +i*Ydim*Xdim ] ; //x,y
               if (aux_double > 0)
                   test(i,j,k)=log(1+aux_double);
               else
                   test(i,j,k)=0.;
               }
       test.write("weight.vol");
    }
    #endif
    #undef VOL

    // release memory for weights and create xmippfourier object for volume
    delete [] FourierVolWeight;
    {
        int ndim    = 3;
        int * myfN;
        try
        {
            myfN = new int [ndim];
        }
        catch (std::bad_alloc&)
        {
          std::cout << "Error allocating memory." << std::endl;
          exit(1);
        }    
        myfN[0]=paddim_vol;
        myfN[1]=paddim_vol;
        myfN[2]=paddim_vol;
        // Set padding dimension
        int fTotalSize=myfN[0]*myfN[1]*myfN[2];
        bool inplace = false;
        xmippFftw Volfft(ndim, myfN, inplace,FourierVol);
        Volfft.Init("ES",FFTW_BACKWARD,false);
        Volfft.CenterRealImageInFourierSpace(false) ;//Change phases
        Volfft.Transform();
        Volfft.CenterFourierTransformInRealSpace(true);//ok MOVE TRANSFORM TO THE CENTER
        
        /* create table with 3D blob fourier transform */
        //precompute blob fourier transform values close to origin
        fourier_blob_table = new double [BLOB_TABLE_SIZE];
        /*it should be divided by 2 not by 4
          but 4 seems to work. I do not know why */
        double ww = 1/(4*(double)(BLOB_TABLE_SIZE-1));
        
        //acces to the blob table should take into account 
        //pad factor but not sampling
        double w;
        double w0= blob_Fourier_val(0., blob);
        for (int i=0; i<BLOB_TABLE_SIZE; i++)
        {
          w = ww*(double)i;
          fourier_blob_table[i] =  blob_Fourier_val(w, blob)/w0;
          //#define DEBUG
          #ifdef DEBUG
          std::cout.setf(std::ios::scientific); 
          std::cout.width(12); 
          std::cout /*<< i << " " 
                    */<< w*paddim_vol <<  "\t" 
                    << "\t"<< fourier_blob_table[i] << " "
                    <<  blob.radius
                    << std::endl;
          #endif
          #undef DEBUG
        }
        //copy volume to original volume
        int xdim=dim;
        int ydim=dim;
        int zdim=dim;
        int center_shiftx = (paddim_vol-xdim)/2+xdim%2;// add one if odd
        int center_shifty = (paddim_vol-ydim)/2+ydim%2;// add one if odd
        int center_shiftz = (paddim_vol-zdim)/2+zdim%2;// add one if odd
        vol().resize(dim,dim,dim);
        for(int i=0;i<zdim;i++)
           for(int j=0;j<ydim;j++)
               for(int k=0;k<xdim;k++)
                {           
                vol(i,j,k)=Volfft.fOut[(center_shiftz + k) + 
                                       (center_shifty + j) * paddim_vol +
                                       (center_shiftx + i) * paddim_vol * paddim_vol ];
                                               //x,y
                }
        vol().setXmippOrigin();
        for (int k = STARTINGZ(vol()); k <= FINISHINGZ(vol()); k++)
        {
            for (int i = STARTINGY(vol()); i <= FINISHINGY(vol()); i++)
            {
                for (int j = STARTINGX(vol()); j <= FINISHINGX(vol()); j++)
                {
                    double r      = sqrt(k*k+i*i+j*j);
                    if(r>paddim_vol/2)
                        vol(i,j,k)=0.;
                    else
                    {
                        double factor = fourier_blob_table[(int)(r*BLOB_TABLE_SIZE/(paddim_vol/2.))];
                        if (factor > 0.001)
                        {
                            vol(i,j,k)   /=  factor;
                        }
                    }
                }
            }
        }
        /*        
        CENTER XMIPP VOLUME   
        APPLY FILTER
        */     
#ifdef NEVERDEFINED
        vol().resize(paddim_vol,paddim_vol,paddim_vol);
        for(int i=0;i<paddim_vol;i++)
           for(int j=0;j<paddim_vol;j++)
               for(int k=0;k<paddim_vol;k++)
                {           
                vol(i,j,k)=Volfft.fOut[(  k) + 
                                       (  j) * paddim_vol +
                                       (  i) * paddim_vol * paddim_vol ];
                                               //x,y
                }
#endif
                
    }    
    //symmetrize in reaL SPACE IF Needed 
    if (fn_sym_vol == "")
    {
        vol.write(fn_out);
    }
    else 
    {
        if (verb > 0) 
            std::cout << std::endl << "Symmetrizing volume (using Bsplines)" << std::endl;
        VolumeXmipp     V_out;
        symmetrize_Bspline(SL_vol, vol, V_out, 3, false, true);
        V_out.write(fn_out);
    }
        /*
       for (int i=0;i<10;i++)
           std::cerr << Volfft.fOut[i] << " " << std::endl;
       std::cerr << "paddim_vol " << paddim_vol << std::endl;
       std::cerr << " myfN[0] myfN[1]  myfN[2]" << myfN[0] << " " 
                                                << myfN[1] << " " 
                                                << myfN[2] <<std::endl;     
       */
    // free memory
    //free(mat_g);
    //free(mat_f);
}
