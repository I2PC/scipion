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
//#define DEBUG

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
    fn_group = getParameter(argc, argv, "-groups","");
    fn_doc = getParameter(argc,argv,"-doc","");
    nr_ref = textToInteger(getParameter(argc, argv, "-nref", "0"));
    fn_prior = getParameter(argc, argv, "-prior", "");
    fn_mask = getParameter(argc, argv, "-mask", "");
    fn_ref = getParameter(argc, argv, "-ref", "");
    //fn_sym = getParameter(argc, argv, "-sym", "");
    fn_root = getParameter(argc, argv, "-o", "out");
    Niter = textToInteger(getParameter(argc, argv, "-iter", "25"));
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));
    fix_fractions = checkParameter(argc, argv, "-fix_fractions");
    fix_sigma_noise = checkParameter(argc, argv, "-fix_sigma_noise");
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    highres = textToFloat(getParameter(argc, argv, "-highres", "0.45"));
    lowres = textToFloat(getParameter(argc, argv, "-lowres", "0.02"));
    debug = checkParameter(argc, argv, "-debug");
    reg0 = textToFloat(getParameter(argc, argv, "-reg0", "0.99"));
    regF = textToFloat(getParameter(argc, argv, "-regF", "0"));
    reg_steps = textToInteger(getParameter(argc, argv, "-steps", "5"));
    eps =  textToFloat(getParameter(argc, argv, "-eps", "0"));
    dont_recalc_fftw = checkParameter(argc, argv, "-dont_recalc_fftw");
    use_tom_conventions = checkParameter(argc, argv, "-tom_conventions");
    do_impute = !checkParameter(argc, argv, "-dont_impute");
    do_ravg_sigma = checkParameter(argc, argv, "-ravg");
    do_som = checkParameter(argc, argv, "-som");
    som_xdim = textToInteger(getParameter(argc, argv, "-xdim", "5"));
    som_ydim = textToInteger(getParameter(argc, argv, "-ydim", "4"));
    do_norm = checkParameter(argc, argv, "-norm");
    do_scale = checkParameter(argc, argv, "-scale");
   
    ang= textToFloat(getParameter(argc, argv, "-ang", "0"));
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
        std::cerr << "  Prior map               : " << fn_prior << std::endl;
        std::cerr << "  Output rootname         : " << fn_root << std::endl;
        std::cerr << "  Low resolution limit    : " <<lowres << " pix^-1 "<< std::endl;
        std::cerr << "  High resolution limit   : " <<highres << " pix^-1 "<< std::endl;
        std::cerr << "  Number of iterations    : " << Niter << std::endl;
        std::cerr << "  Initial regularisation  : " << reg0 << std::endl;
        std::cerr << "  Final regularisation    : " << regF << std::endl;
        std::cerr << "  Convergence criterium   : " << eps <<std::endl;
        std::cerr << "  Number of reg. steps    : " << reg_steps << std::endl;
        if (use_tom_conventions)
        {
            std::cerr << "  -> Using TOM Toolbox's geometrical conventions "<<std::endl;
        }
        if (fix_fractions)
        {
            std::cerr << "  -> Do not update estimates of model fractions." << std::endl;
        }
        if (fix_sigma_noise)
        {
            std::cerr << "  -> Do not update sigma-estimate of noise." << std::endl;
        }
        if (!do_impute)
        {
            std::cerr << "  -> Do not use imputation EM-agorithm." << std::endl;
        }
        if (do_ravg_sigma)
        {
            std::cerr << "  -> Radial average the sigma estimates." << std::endl;
        }
        if (do_norm)
        {
            std::cerr << "  -> Internally normalize the subtomograms." << std::endl;
        }
        if (do_scale)
        {
            std::cerr << "  -> Internally rescale brightness of the subtomograms." << std::endl;
        }

        std::cerr << " -----------------------------------------------------------------" << std::endl;
    }

}

// Usage ===================================================================
void Prog_mlf_tomo_prm::usage()
{
    std::cerr << "Usage:  mlf_tomo [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with the input sub-tomograms \n";
    std::cerr << "   -nref <int>                 : Number of references to generate automatically (recommended)\n";
    std::cerr << "   OR -ref <selfile/volume>         OR selfile with initial references/single reference image \n";
    std::cerr << "   OR -som -xdim <int> -ydim <int>  OR self-organizing map of given dimensions \n";
    std::cerr << " [ -prior <volume> ]           : Name for prior map (only relevant for -nref option) \n";
    std::cerr << " [ -o <rootname=\"out\"> ]       : Output rootname \n";
    std::cerr << " [ -doc <docfile>]             : Docfile with angles, offsets and wedge information \n";
    std::cerr << " [ -lowres <=0.02> ]           : Low-resolution limit (in 1/pixel) \n";
    std::cerr << " [ -highres <=0.5> ]           : High-resolution limit (in 1/pixel) \n";
    std::cerr << " [ -tom_conventions ]          : Use TOM Toolbox's geometrical conventions \n";
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
#ifdef DEBUG
    std::cerr<<"start produceSideInfo"<<std::endl;
#endif

    headerXmipp head;
    Matrix3D<double> Maux, Faux_real, Faux_imag;
    Matrix3D<std::complex<double> > Faux, Faux2;
    int xdim, ydim, zdim, c, iaux, ifound;
    double dum, avg;

    // image sizes
    VolumeXmipp vol;
    SFi.go_beginning();
    vol.read(SFi.NextImg());
    Xdim = XSIZE(vol());
    Ydim = YSIZE(vol());
    Zdim = ZSIZE(vol());
    dim3 = (double) (Xdim * Ydim * Zdim);

    // Get number of references if selfile/volume given
    if (fn_ref != "")
    {
        // A. User-provided reference(s)
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
    else if (do_som)
    {
        nr_ref = som_xdim * som_ydim;
    }
    else if (nr_ref <= 0)
        REPORT_ERROR(1,"Please provide either -ref, -nref or -som");

    // Initialize regularization matrix
    initializeRegularizationMatrix(false);

    // Read SymList
    //if (fn_sym!="") 
    //    SL.read_sym_file(fn_sym);
    //else
    //    SL.read_sym_file("c1");

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
    try
    {
        is_in_range = new bool [fftw_hsize];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in main");
    }


    int zz, yy;
    double sz,sy,sx,res;
    hsize = 0;
    for ( int z=0, ii=0; z<Zdim; z++ ) {
        if ( z > (Zdim - 1)/2 ) zz = z-Zdim;
        else zz = z;
        for ( int y=0; y<Ydim; y++ ) {
            if ( y > (Ydim - 1)/2 ) yy = y-Ydim;
            else yy = y;
            for ( int xx=0; xx<Xdim/2 + 1; xx++,ii++ ) {
                // Think about non-cubic volumes here: resolution in dig.freq.
                sx = (double)xx / (double)Xdim;
                sy = (double)yy / (double)Ydim;
                sz = (double)zz / (double)Zdim;
                res = sqrt(sx*sx+sy*sy+sz*sz);

// TODO: Also exclude double x==0 axis!!

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
#ifdef DEBUG
    std::cerr<<" lowres-highres hsize= "<<hsize<<" fftw hsize= "<<fftw_hsize<<std::endl;
#endif

    // Get groups structure from SFg and store into SFi
    SelFile SFtmp;
    SelLine SL;
    int igroup;
    nr_group = 1;
    nr_imgs = (double)SFi.ImgNo();
    nr_imgs_per_group.resize(1);
    nr_imgs_per_group(0) = nr_imgs;
    if (fn_group != "")
    {
        SFg.read(fn_group);
        SFi.go_beginning();
        while (!SFi.eof())
        {
            FileName fn_img=SFi.NextImg();
            if (SFi.eof()) break;
            SFg.search(fn_img);
            SL = SFg.current();
            SFtmp.insert(SL);
            igroup = SL.get_number();
            if (igroup > nr_group)
            {
                nr_group = igroup;
                nr_imgs_per_group.resize(nr_group);
            }
            nr_imgs_per_group(igroup-1) += 1.;
        }
        SFi = SFtmp;
    }

    // Initialize optimal scales to one.
    if (do_scale)
    {
	imgs_scale.clear();
        for (int i = 0; i < SFi.ImgNo(); i++)
	    imgs_scale.push_back(1.);
        for (int i = 0; i < nr_ref; i++)
	    refs_avgscale.push_back(1.);
        average_scale = 1.;
    }

    // Initialize reg
    reg = reg0;

#ifdef DEBUG
    std::cerr<<"done produceSideInfo"<<std::endl;
#endif

}

void Prog_mlf_tomo_prm::initializeRegularizationMatrix(bool is_som)
{
    som_reg_matrix.resize(nr_ref,nr_ref);
    som_reg_matrix.initZeros();
    if (is_som)
    {
        for (int iref1=0; iref1<nr_ref; iref1++)
        {
            int x1 = iref1%som_xdim;
            int y1 = (iref1 - x1)/som_xdim;
            int count = 0;
            for (int iref2=0; iref2<nr_ref; iref2++)
            {
                int x2 = iref2%som_xdim;
                int y2 = (iref2 - x2)/som_xdim;
                int r2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
                if (r2 == 1)
                {
                    dMij(som_reg_matrix,iref1,iref2) = 1.;
                    count++;
                }
            }
            for (int iref2=0; iref2<nr_ref; iref2++)
            {
                dMij(som_reg_matrix,iref1,iref2) /= (double)count;
            }
        }
    }
    else
    {
        for (int iref1=0; iref1<nr_ref; iref1++)
            for (int iref2=0; iref2<nr_ref; iref2++)
                dMij(som_reg_matrix,iref1,iref2) = 1./(double)(nr_ref);
    }

//#define DEBUG_REG_MATRIX
#ifdef DEBUG_REG_MATRIX
    std::cerr<<"Regularization matrix= "<<som_reg_matrix<<std::endl;
#endif

}


void Prog_mlf_tomo_prm::produceSideInfo2()
{
#ifdef DEBUG
    std::cerr<<"start produceSideInfo2"<<std::endl;
#endif

    DocFile DF;
    FileName fn_vol;
    int iaux;
    bool found;

    // Store tomogram angles, offset vectors and missing wedge parameters
    if (fn_doc!="") {
        DF.read(fn_doc);
        SFi.go_beginning();
        while (!SFi.eof()) 
        {
            fn_vol=SFi.NextImg();
            if (SFi.eof()) break;
            DF.go_beginning();
            if (DF.search_comment(fn_vol)) 
            {
                img_ang1.push_back(DF(0));
                img_ang2.push_back(DF(1));
                img_ang3.push_back(DF(2));
                img_xoff.push_back(DF(3));
                img_yoff.push_back(DF(4));
                img_zoff.push_back(DF(5));
                img_th0.push_back( DF(6));
                img_thF.push_back( DF(7));
            } 
            else 
            {
                std::cerr << "ERROR% "<<fn_vol
                          <<" not found in document file"
                          <<std::endl;
                exit(0);
            }
        }
    } 
    else 
    {
        SFi.go_beginning();
        while (!SFi.eof()) 
        {
            SFi.NextImg();
            if (SFi.eof()) break;
            img_ang1.push_back( 0.);
            img_ang2.push_back( 0.);
            img_ang3.push_back( 0.);
            img_xoff.push_back( 0.);
            img_yoff.push_back( 0.);
            img_zoff.push_back( 0.);
            img_th0.push_back(  0.);
            img_thF.push_back(  0.);
        }
    }

#ifdef DEBUG
    std::cerr<<"done produceSideInfo2"<<std::endl;
#endif
}


/// Get Transformation matrix 
Matrix2D< double > Prog_mlf_tomo_prm::getTransformationMatrix(double ang1, 
                                                              double ang2, 
                                                              double ang3, 
                                                              double xoff, /* = 0. */
                                                              double yoff, /* = 0. */
                                                              double zoff) /* = 0. */
{
    Matrix2D<double>  A(4,4);
    if (use_tom_conventions)
    {
        // Use TOM Toolbox conventions
        double cosphi = COSD(ang1);
        double cospsi = COSD(ang2);
        double costheta = COSD(ang3);
        double sinphi = SIND(ang1);
        double sinpsi = SIND(ang2);
        double sintheta = SIND(ang3);

        /* The following piece of code is fromTOM: Sptrans/tom_rotatec.c
          rm00=cospsi*cosphi-costheta*sinpsi*sinphi;
          rm10=sinpsi*cosphi+costheta*cospsi*sinphi;
          rm20=sintheta*sinphi;
          rm01=-cospsi*sinphi-costheta*sinpsi*cosphi;
          rm11=-sinpsi*sinphi+costheta*cospsi*cosphi;
          rm21=sintheta*cosphi;
          rm02=sintheta*sinpsi;
          rm12=-sintheta*cospsi;
          rm22=costheta;
        */

        A(0, 0) =  cospsi*cosphi-costheta*sinpsi*sinphi;
        A(1, 0) =  sinpsi*cosphi+costheta*cospsi*sinphi;
        A(2, 0) =  sintheta*sinphi;
        A(0, 1) = -cospsi*sinphi-costheta*sinpsi*cosphi;
        A(1, 1) = -sinpsi*sinphi+costheta*cospsi*cosphi; 
        A(2, 1) =  sintheta*cosphi; 
        A(0, 2) =  sintheta*sinpsi;
        A(1, 2) = -sintheta*cospsi;
        A(2, 2) =  costheta;

        A(0,3) = xoff;
        A(1,3) = yoff;
        A(2,3) = zoff;

//#define DEBUG_TOM_CONVENTIONS
#ifdef  DEBUG_TOM_CONVENTIONS
/*
XMIPP  =      TOM
rot    =      90. + psi
tilt   =      -theta 
psi    =      -90 + phi
xoff   =      -x
yoff   =      -y
zoff   =      -z
*/
            Matrix2D<double> B=A;
            double r,t,p;
            B.resize(3,3);
            std::cerr<<"------"<<std::endl;
            std::cerr<< "angles TOM= "<<ang1<<" "<<ang2<<" "<<ang3<<std::endl;
            std::cerr<< "angles XMIPP= "<<90. + ang2<<" "<<-ang3<<" "<<-90. + ang1<<std::endl;
            std::cerr<<"TOM matrix = "<<B<<std::endl;
            Euler_angles2matrix(90. + ang2, -ang3, -90. + ang1, B);
            std::cerr<<"XMIPP matrix = "<<B<<std::endl;
            Euler_matrix2angles(B,r,t,p);
            std::cerr<<"angle XMIPP again= "<<r<<" "<<t<<" "<<p<<std::endl;
            B=B.transpose();
            std::cerr<<"XMIPP transposed matrix = "<<B<<std::endl;
            Euler_matrix2angles(B,r,t,p);
            std::cerr<<"angles XMIPP transposed matrix= "<<r<<" "<<t<<" "<<p<<std::endl;
          
#endif

    }
    else
    {
        // Use XMIPP conventions
        A = Euler_rotation3DMatrix(ang1, ang2, ang3);
        A = A.inv();
        A(0,3) = -xoff;
        A(1,3) = -yoff;
        A(2,3) = -zoff;
    }

    return A;
}
 
/// Get binary missing wedge (or pyramid) 
void Prog_mlf_tomo_prm::getMissingWedge(bool * measured,
                                        Matrix2D<double> A,
                                        const double theta0_alongy, 
                                        const double thetaF_alongy,
                                        const double theta0_alongx /*= 0.*/,
                                        const double thetaF_alongx /*= 0.*/)
{

    if (theta0_alongy==0. && thetaF_alongy==0. && theta0_alongx==0. && thetaF_alongx==0.)
    {
        for (int i=0; i < hsize; i++)
            measured[i] = true;
        return;
    }

    double xp, yp, zp;
    double tg0_y, tgF_y, tg0_x, tgF_x, limx0, limxF, limy0, limyF;

    tg0_y = -tan(PI * (-90. - thetaF_alongy) / 180.);
    tgF_y = -tan(PI * (90. - theta0_alongy) / 180.);
    tg0_x = -tan(PI * (-90. - thetaF_alongx) / 180.);
    tgF_x = -tan(PI * (90. - theta0_alongx) / 180.);

//#define DEBUG_WEDGE
#ifdef DEBUG_WEDGE
    std::cerr<<"tg0_y= "<<tg0_y<<std::endl;
    std::cerr<<"tgF_y= "<<tgF_y<<std::endl;
    std::cerr<<"tg0_x= "<<tg0_x<<std::endl;
    std::cerr<<"tgF_x= "<<tgF_x<<std::endl;
#endif

    int zz, yy;
    for ( int z=0, i=0, ii=0; z<Zdim; z++ ) {
        if ( z > (Zdim - 1)/2 ) 
            zz = z-Zdim;
        else 
            zz = z;
        for ( int y=0; y<Ydim; y++ ) {
            if ( y > (Ydim - 1)/2 ) yy = y-Ydim;
            else yy = y;
            for ( int xx=0; xx<Xdim/2 + 1; xx++,ii++ ) {
                // Only concerned with those components within resolution limits
                if (is_in_range[ii])
                {
                    // Rotate the wedge
                    xp = dMij(A, 0, 0) * xx + dMij(A, 0, 1) * yy + dMij(A, 0, 2) * zz;
                    yp = dMij(A, 1, 0) * xx + dMij(A, 1, 1) * yy + dMij(A, 1, 2) * zz;
                    zp = dMij(A, 2, 0) * xx + dMij(A, 2, 1) * yy + dMij(A, 2, 2) * zz;
                    // Calculate the limits
                    limx0 = tg0_y * zp;
                    limxF = tgF_y * zp;
                    limy0 = tg0_x * zp;
                    limyF = tgF_x * zp;
                    
                    if (zp >= 0)
                    {
                        if ((xp <= limx0 || xp >= limxF) && (yp <= limy0 || yp >= limyF))
                            measured[i] = true; 
                        else
                            measured[i] = false;
                    }
                    else
                    {
                         if ((xp <= limxF || xp >= limx0) && (yp <= limyF || yp >= limy0))
                            measured[i] = true; 
                        else
                            measured[i] = false;
                    }
//#define DEBUG_WEDGE2
#ifdef DEBUG_WEDGE2
                    std::cerr<<" xx,yy,zz= "<<xx<<" "<<yy<<" "<<zz<<std::endl;
                    std::cerr<<" xp,yp,zp= "<<xp<<" "<<yp<<" "<<zp<<std::endl;
                    std::cerr<<"limx0, limxF= "<<limx0<<" "<<limxF<<std::endl;
                    if (measured[i])
                        std::cerr<<"true"<<std::endl;
                    else 
                        std::cerr<<"false"<<std::endl;
#endif
                    i++;
                } 
            }
        }
    }

#ifdef DEBUG_WEDGE
    VolumeXmipp test(Xdim,Ydim,Zdim), rot(Xdim,Ydim,Zdim);
    test().initZeros();

    for(int i=0,ii=0,iii=0;i<Zdim;i++)
        for(int j=0;j<Ydim;j++)
            for(int k=0;k<Xdim/2 + 1; k++,ii++)
                if (is_in_range[ii])
                {
                    if (measured[iii])
                    {
                        test(i,j,k)=1.;
                    }
                    iii++;
                }
    test.write("wedge.ftt");
    Matrix1D<double> off(3);
    off.initConstant(Xdim/2);
    test().selfTranslate(off);
    test.write("Fwedge.vol");
    test().setXmippOrigin();
    Matrix3D<int> ress(Xdim,Ydim,Zdim);
    ress.setXmippOrigin();
    A=A.inv();
    BinaryWedgeMask(test(),theta0_alongy,thetaF_alongy, A);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(test())
    {
        double res = (double)(i*i+j*j+k*k);
        res = sqrt(res)/Xdim;
        if (!(res <= highres && res >= lowres))
            VOL_ELEM(test(),k,i,j) = 0;
    }
    test.write("Mwedge.vol");
    exit(0);
#endif


}


void Prog_mlf_tomo_prm::generateInitialReferences(double * dataRefs, double * dataWsumWedsPerRef)
{

#ifdef DEBUG
    std::cerr<<"start generateInitialReferences"<<std::endl;
#endif

    if (fn_ref != "")
    {
        // A. User-provided reference(s)

        // Read Xmipp volumes and calculate their FTTws        
        std::complex<double> *DATAREFS, *REF;
        DATAREFS  = (std::complex<double> *) dataRefs;
        VolumeXmipp vol;
        int refno = 0;
        SFr.go_beginning();
        while (!SFr.eof())
        {
            FileName fn_img=SFr.NextImg();
            if (SFr.eof()) break;
            vol.read(fn_img);
            if (XSIZE(vol()) != Xdim ||
                XSIZE(vol()) != Xdim ||
                XSIZE(vol()) != Xdim )
                REPORT_ERROR(1,"Wrong dimensions for reference volume!");
            forwfftw.SetPoints(MULTIDIM_ARRAY(vol()));
            forwfftw.Transform();
            forwfftw.Normalize();
            // Get only points within resolution range
            REF = (std::complex<double> *) forwfftw.fOut;
            for (int i = 0, ii=0; i< fftw_hsize; i++)
            {
                if (is_in_range[i]) 
                {
                    DATAREFS[refno*hsize + ii] = REF[i];
                    ii++;
                }
            }
            refno++;
        }

    }
    else
    {

        // B. Random subsets of the data
        double               *dataImg, *dataPrior;
        std::complex<double> *DATAIMG, *DATAREFS, *DATAPRIOR, *IMG;
        bool                 *dataMeasured;
        FileName             fn;
        std::vector<double>  count_refs;
        Matrix2D<double>     A_img(4,4);
        double               th0,thF;
        int                  refno, c, nn = SFi.ImgNo(), imgno = 0;

        // Reserve memory for dataImg and dataWedge
        try
        {
            dataImg      = new double[size];
            dataPrior    = new double[size];
            dataMeasured = new bool[hsize];
        }
        catch (std::bad_alloc&)
        {
            REPORT_ERROR(1,"Error allocating memory in expectation");
        }
        DATAREFS  = (std::complex<double> *) dataRefs;
        DATAIMG   = (std::complex<double> *) dataImg;
        DATAPRIOR = (std::complex<double> *) dataPrior;

        // Read and transform the prior map
        if (fn_prior != "")
        {
            VolumeXmipp prior;
            prior.read(fn_prior);
            if (XSIZE(prior()) != Xdim ||
                XSIZE(prior()) != Xdim ||
                XSIZE(prior()) != Xdim )
                REPORT_ERROR(1,"Wrong dimensions for prior volume!");

            forwfftw.SetPoints(MULTIDIM_ARRAY(prior()));
            forwfftw.Transform();
            forwfftw.Normalize();
            IMG =  (std::complex<double> *) forwfftw.fOut;
            // Get only points within resolution range
            for (int i = 0, ii=0; i < fftw_hsize; i++)
                if (is_in_range[i]) 
                {
                    DATAPRIOR[ii] = IMG[i];
                    ii++;
                }
        }
        else
        {
            // If no prior is given: assume all-zero prior
            for (int i=0; i<hsize; i++)
                DATAPRIOR[i] = 0.;
        }

        // Initialize sums to zero
        for (int i = 0; i < nr_ref * hsize; i++)
        {
            DATAREFS[i] = 0.;
            dataWsumWedsPerRef[i] = 0.;
        }
        count_refs.resize(nr_ref);

        if (verb > 0)
        {
            std::cerr << "  Generating initial references by averaging over random subsets" << std::endl;
            init_progress_bar(nn);
            c = XMIPP_MAX(1, nn / 60);
        }
        SFi.go_beginning();
        randomize_random_generator();
        while (!SFi.eof())
        {
            // read tomogram from disc
            fn = SFi.NextImg();
            if (SFi.eof()) break;
            backfftw.read(fn + ".fftw");
            IMG =  (std::complex<double> *) backfftw.fIn;
            // Get only points within resolution range
            for (int i = 0, ii=0; i < fftw_hsize; i++)
                if (is_in_range[i]) 
                {
                    DATAIMG[ii] = IMG[i];
                    ii++;
                }

            // Get the missing wedge
            getMissingWedge(dataMeasured,
                            getTransformationMatrix(img_ang1[imgno],img_ang2[imgno],img_ang3[imgno]),
                            img_th0[imgno],
                            img_thF[imgno]);

            // Randomly choose a reference
            refno = FLOOR(rnd_unif(0,nr_ref));
            count_refs[refno] += 1.;

            // Store sum and sum of weights
            for (int i = 0; i < hsize; i++)
            {
                int ii = refno*hsize + i;
                DATAREFS[ii] += DATAIMG[i];
                if (dataMeasured[i])
                    dataWsumWedsPerRef[ii] += 1.;
            }

            if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
            imgno++;
        }
        if (verb > 0) progress_bar(nn);

#ifdef DEBUG
        std::cerr<<" Distribution of initial random sets= "<<std::endl;
        for (int r=0; r<nr_ref; r++)
            std::cerr<<" -> reference "<<r+1<<" has "<<count_refs[r]<<" particles."<<std::endl;
#endif

        // Calculate dataRefs with imputation of the prior
        for (int refno = 0; refno < nr_ref; refno++)
        {
            for (int i=0; i < hsize; i++)
            {
                int ii = refno * hsize + i;
                // Observed values
                DATAREFS[ii] /= count_refs[refno];
                // Impute missing values with the prior
                DATAREFS[ii] += DATAPRIOR[i] * (1. - dataWsumWedsPerRef[ii] / count_refs[refno]);
            }
        }
    }

#ifdef DEBUG
    VolumeXmipp test(Xdim,Ydim,Zdim);
    for (int refno=0; refno < nr_ref; refno++)
    {
        test().initZeros();
        for(int i=0,ii=0,iii=0;i<Zdim;i++)
            for(int j=0;j<Ydim;j++)
                for(int k=0;k<Xdim/2 + 1; k++,ii++)
                    if (is_in_range[ii])
                    {
                        test(i,j,k)=dataWsumWedsPerRef[refno*hsize + iii];
                        iii++;
                    }
        CenterFFT(test(),true);
        FileName fnt;
        fnt.compose("sumwedges_ref",refno,"vol");
        test.write(fnt);
    }
#endif

    // Initialize model fractions
    alpha_k.clear();
    for (int refno = 0; refno < nr_ref; refno++)
    {
        alpha_k.push_back(1./nr_ref);
    }

#ifdef DEBUG
    std::cerr<<"done generateInitialReferences"<<std::endl;
#endif
}

// This routine is for (splitted) SF-dependent side-info calculations
void Prog_mlf_tomo_prm::calculateAllFFTWs()
{
#ifdef DEBUG
    std::cerr<<"start calculateAllFFTWs"<<std::endl;
#endif

    FileName fni, fno;
    VolumeXmipp vol, ave(Xdim,Ydim,Zdim), mask(Xdim,Ydim,Zdim);
    Matrix2D<double> A_img(4,4);
    int c, nn = SFi.ImgNo(), imgno = 0;

    ave().setXmippOrigin();
    if (fn_mask != "")
    {
        mask.read(fn_mask);
        mask().setXmippOrigin();
    }
    else
    {
        mask().setXmippOrigin();
        RaisedCosineMask(mask(),(double)Xdim-2,(double)Xdim);
    }

    if (verb > 0)
    {
        std::cerr << "  Calculating FFTs for all input maps ... " << std::endl;
        init_progress_bar(nn);
        c = XMIPP_MAX(1, nn / 60);
    }
    SFi.go_beginning();
    while (!SFi.eof())
    {
        fni = SFi.NextImg();
        if (SFi.eof()) break;
        if ( ! (exists(fni) && dont_recalc_fftw) )
        {
            vol.read(fni);
            vol().setXmippOrigin();
            if (XSIZE(vol()) != Xdim ||
                XSIZE(vol()) != Xdim ||
                XSIZE(vol()) != Xdim )
                REPORT_ERROR(1,"Wrong dimensions for input volume!");

            // Get the transformation matrix
            A_img = getTransformationMatrix(img_ang1[imgno], img_ang2[imgno], img_ang3[imgno],
                                            img_xoff[imgno], img_yoff[imgno], img_zoff[imgno]);

            // IS_INV to be consistent with Xmipp conventions
            // if some day want to switch to IS_NOT_INV, do A=A.inv() in getMissingWedge!!!!
            vol().selfApplyGeometryBSpline(A_img,3,IS_INV,DONT_WRAP,0.);

            // Mask if necessary
            vol() *= mask();

            // Normalize if necessary
            if (do_norm)
            {
                double avg = 0., stddev = 0., c = 0.;
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(vol())
                {
                    if (dVkij(mask(),k,i,j) > 0.)
                    {
                        avg += dVkij(vol(),k,i,j);
                        stddev += dVkij(vol(),k,i,j) * dVkij(vol(),k,i,j);
                        c+= 1.;
                    }
                }
                avg /= c;
                stddev = stddev / c - avg * avg;
                stddev *= c / (c - 1.);
                vol() -= avg;
                vol() /= stddev;
            }
            // Store average of the (masked and/or normalized) maps
            ave() += vol();

            forwfftw.SetPoints(MULTIDIM_ARRAY(vol()));
            forwfftw.Transform();
            forwfftw.Normalize();
            fno = fni + ".fftw";
            forwfftw.write(fno);
            if (verb > 0) if (imgno % c == 0) progress_bar(imgno);
            imgno++;

        }
    }
    if (verb > 0) progress_bar(nn);

    if (!dont_recalc_fftw)
    {
        // Write out normal average map of all subtomograms
        ave() /= (double)imgno;
        fno = fn_root + "_average.vol"; 
        std::cerr<<" Writing normal average of all subtomograms as "<<fno<<std::endl; 
        ave.write(fno);
    }

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
        if (SFi.eof()) break;
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
    }

#ifdef DEBUG
    std::cerr<<"done estimateInitialNoiseSpectra"<<std::endl;
#endif
}


// Here perform the main probability-weighted integration over all
// rotations, translations and classes of the given image
void Prog_mlf_tomo_prm::expectationSingleImage(int igroup,
                                               double * dataImg,
                                               bool   * dataMeasured,
                                               double * dataRefs,
                                               double * dataSigma,
                                               double * dataWsumRefs,
                                               double * dataWsumWedsPerRef,
                                               double * dataWsumWedsPerGroup,
                                               double * dataWsumDist,
                                               double * dataSumWRefs,
                                               double * dataWsumScale,
                                               double * dataWsumScale2,
                                               int    & opt_refno, 
                                               double & LL, 
                                               double & Pmax,
                                               double & opt_scale )

{

    double aux, weight, diff2, sumweight = 0., maxweight = 0., mindiff2 = 99.e99;
    double ref_scale, wsum_sc = 0., wsum_sc2 = 0.;
    double * weights;
    weights = new double [nr_ref];

    std::complex<double> *DATAREFS, *DATAWSUMREFS, *DATAIMG;
    DATAIMG      = (std::complex<double> *) dataImg;
    DATAREFS     = (std::complex<double> *) dataRefs;
    DATAWSUMREFS = (std::complex<double> *) dataWsumRefs;

    if (!do_scale) 
    {
        opt_scale = 1.;
        ref_scale = 1.;
    }
    // TODO: Adapt dataSigma arrays to take into account that (x==0) is double in the fftw!
    for (int refno = 0; refno < nr_ref; refno++)
    {
        diff2 = 0.;
        if (do_scale) ref_scale = opt_scale / refs_avgscale[refno];
        for (int i = 0; i < hsize; i++)
        {
            int iiref = refno * hsize + i;
            int iig = igroup * hsize + i;
            if (dataMeasured[i])
            {
                aux = abs(DATAIMG[i] - ref_scale * DATAREFS[iiref]);
                diff2 += (aux * aux)/dataSigma[iig];
//#define DEBUG_ALOT_EXPSINGLE 
#ifdef DEBUG_ALOT_EXPSINGLE                
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

        if (do_scale)
        {
            for (int i = 0; i < hsize; i++)
            {
                int iiref = refno * size + 2*i;
                if (dataMeasured[i])
                {
                    wsum_sc  += weights[refno] * dataImg[2*i] * dataRefs[iiref];
                    wsum_sc  += weights[refno] * dataImg[2*i+1] * dataRefs[iiref+1];
                    wsum_sc2 += weights[refno] * dataRefs[iiref] * dataRefs[iiref];
                    wsum_sc2 += weights[refno] * dataRefs[iiref+1] * dataRefs[iiref+1];
                }
            }
        }
    }

    // Update the internal scale
    if (do_scale)
    {
        opt_scale = wsum_sc / wsum_sc2;
    }

    // Store Pmax/sumP
    Pmax = maxweight / sumweight;
//#define DEBUG_EXPSINGLE 
#ifdef DEBUG_EXPSINGLE                
    std::cerr<<"Pmax= "<<Pmax<<" mindiff2= "<<mindiff2<<std::endl;
#endif

    
    // Then, store the sum of all weights
    for (int refno = 0; refno < nr_ref; refno++)
    {
        weight = weights[refno] / sumweight;
#ifdef DEBUG_EXPSINGLE                
        std::cerr<<" refno= "<<refno<<" w= "<< weight<<" ";
#endif
        if (weight > SIGNIFICANT_WEIGHT_LOW)
        {
            dataSumWRefs[refno] += weight;
            dataWsumScale[refno] += weight * opt_scale;
            dataWsumScale2[refno] += weight * opt_scale * opt_scale;
            for (int i = 0; i < hsize; i++)
            {
                int iiref = refno * hsize + i;
                int iig = igroup * hsize + i;
                if (dataMeasured[i])
                {
                    aux = abs(DATAIMG[i] - opt_scale * DATAREFS[iiref]);
                    dataWsumDist[iig] += weight * aux * aux;
                    dataWsumWedsPerGroup[iig] += weight;
                    DATAWSUMREFS[iiref] += weight * opt_scale * DATAIMG[i];
                    dataWsumWedsPerRef[iiref] += weight;
                }
            }
        }
    }
    
    // Precalculate normalization constant for LL update
    double logsigma2 = 0.;
    for (int i = 0; i < hsize; i++)
        if (dataMeasured[i])
            logsigma2 += 2 * log( sqrt(PI * dataSigma[igroup * hsize + i]));

    // Update the log-likelihood function value
    // 1st term: log(refw_i)
    // 2nd term: for subtracting mindiff2
    // 3rd term: for missing normalization constant
    LL += log(sumweight) - mindiff2 - logsigma2;

}



void Prog_mlf_tomo_prm::expectation(double * dataRefs,
                                    double * dataSigma,
                                    double * dataWsumRefs,
                                    double * dataWsumWedsPerRef,
                                    double * dataWsumWedsPerGroup,
                                    double * dataWsumDist,
                                    double * dataSumWRefs,
                                    double * dataWsumScale,
                                    double * dataWsumScale2,
                                    double & LL,
                                    double & avePmax,
                                    DocFile & DFo)
{

#ifdef DEBUG
    std::cerr<<"start expectation"<<std::endl;
#endif
    int              opt_refno, nn, imgno, igroup, idum;
    double           *dataImg;
    bool             *dataMeasured;
    SelLine          SL;
    FileName         fn;
    Matrix1D<double> dataline(11);
    Matrix2D<double> A_img(4,4);
    double           Pmax, th0,thF;
    double           opt_scale = 1.;

    // Reserve memory for dataImg and dataWedge
    try
    {
        dataImg      = new double[size];
        dataMeasured = new bool[hsize];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in expectation");
    }

    // Initialize docfile header
    DFo.clear();
    if (use_tom_conventions)
        DFo.append_comment("Headerinfo columns: phi (1), psi (2), theta (3), Xoff (4), Yoff (5), Zoff (6), WedTh0 (7), WedThF (8), Ref (9), Pmax/sumP (10), brightness (11)");
    else
        DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Zoff (6), WedTh0 (7), WedThF (8), Ref (9), Pmax/sumP (10), brightness (11)");

    // Initialize weighted sums to zero
    for (int i = 0; i < nr_ref * size; i++)
        dataWsumRefs[i] = 0.;
    for (int i = 0; i < nr_ref * hsize; i++)
        dataWsumWedsPerRef[i] = 0.;
    for (int i = 0; i < nr_group * hsize; i++)
        dataWsumWedsPerGroup[i] = 0.;
    for (int i = 0; i < nr_group * hsize; i++)
        dataWsumDist[i] = 0.;
    for (int i = 0; i < nr_ref; i++)
    {
        dataSumWRefs[i] = 0.;
        dataWsumScale[i] = 0.;
        dataWsumScale2[i] = 0.;
    }
    LL = 0.;
    avePmax = 0.;

//#define DEBUG_TOM
#ifdef DEBUG_TOM
    use_tom_conventions=true;
    std::cerr<<"tom angles = 31, 152, 47";
    A_img=getTransformationMatrix(31,152,47);
    std::cerr<<"TOM matrix= "<<A_img;
    std::cerr<<"inverse TOM matrix= "<<A_img.inv();
    std::cerr<<"tom angles = -152, -31, -47";
    A_img=getTransformationMatrix(-152,-31,-47);
    std::cerr<<"TOM matrix= "<<A_img;
    std::cerr<<"inverse TOM matrix= "<<A_img.inv();

    use_tom_conventions=false;
    std::cerr<<"xmipp angles = 31, 47, 152";
    A_img=getTransformationMatrix(31,47,152);
    std::cerr<<"XMIPP matrix= "<<A_img;
    std::cerr<<"inverse XMIPP matrix= "<<A_img.inv();
    std::cerr<<"xmipp angles = -152, -47, -31";
    A_img=getTransformationMatrix(-152,-47,-31);
    std::cerr<<"XMIPP matrix= "<<A_img;
    std::cerr<<"inverse XMIPP matrix= "<<A_img.inv();
    exit(0);

#endif 

    // Loop over all images
    imgno = 0;
    //nn = SFi.ImgNo();
    //if (verb > 0) init_progress_bar(nn);
    SFi.go_beginning();
    while ((!SFi.eof()))
    {
        // Get the group
        SL = SFi.current();
        igroup = SL.get_number() - 1;        
        fn = SFi.NextImg();
        if (SFi.eof()) break;

        // Get the geometrical information
        A_img=getTransformationMatrix(img_ang1[imgno],img_ang2[imgno],img_ang3[imgno]);
        th0=img_th0[imgno];
        thF=img_thF[imgno];

        // read tomogram from disc
        backfftw.read(fn + ".fftw");
        // Get only points within resolution range
        for (int i = 0, ii=0; i< 2*fftw_hsize; i++)
            if (is_in_range[i/2]) 
            {
                dataImg[ii] = backfftw.fIn[i];
                ii++;
            }

        // Get optimal internal scale factor for this subtomogram
	if (do_scale)
            opt_scale = imgs_scale[imgno];

        // get missing wedge
        getMissingWedge(dataMeasured,A_img,th0,thF);

//#define EXP_DEBUG_WEDGE 
#ifdef EXP_DEBUG_WEDGE
    VolumeXmipp img(Xdim,Ydim,Zdim), wed(Xdim,Ydim,Zdim);
    img().initZeros();
    wed().initZeros();

    std::complex<double> * iaux;
    iaux = (std::complex<double> *) dataImg;
    for(int i=0,ii=0,iii=0;i<Zdim;i++)
        for(int j=0;j<Ydim;j++)
            for(int k=0;k<Xdim/2 + 1; k++,ii++)
                if (is_in_range[ii])
                {
                    img(i,j,k)=abs(iaux[iii]);
                    if (dataMeasured[iii])
                    {
                        wed(i,j,k)=1.;
                    }
                    iii++;
                }
    CenterFFT(wed(),true);
    CenterFFT(img(),true);
    wed.write("wed.ftt");
    img.write("img.ftt");
    std::cerr<<" Euler angles= "<<img_ang1[imgno]<<" "<<img_ang2[imgno]<<" "<<img_ang3[imgno]<<std::endl;
    std::cerr<<" Wedge angles= "<<th0<<" "<<thF<<std::endl;
    std::cerr<<"EXP_DEBUG_WEDGE: wrote wed.fft and img.fft for subtomogram "<<fn<<std::endl;
    std::cerr<<"Press any key to continue.. "<<std::endl;
    char a;
    std::cin >>a;
#endif

        // Perform expectation step 
        expectationSingleImage(igroup,
                               dataImg, 
                               dataMeasured,
                               dataRefs,
                               dataSigma,
                               dataWsumRefs,
                               dataWsumWedsPerRef,
                               dataWsumWedsPerGroup,
                               dataWsumDist,
                               dataSumWRefs,
                               dataWsumScale,
                               dataWsumScale2,
                               opt_refno, 
                               LL, 
                               Pmax,
                               opt_scale);

	// Store optimal scale in memory
	if (do_scale)
	    imgs_scale[imgno] = opt_scale;

        // Output to docfile                            XMIPP or  TOM
        dataline(0)=img_ang1[imgno];                 // rot   or  phi
        dataline(1)=img_ang2[imgno];                 // tilt  or  psi
        dataline(2)=img_ang3[imgno];                 // psi   or  theta
        dataline(3)=img_xoff[imgno];                 // Xoff
        dataline(4)=img_yoff[imgno];                 // Yoff
        dataline(5)=img_zoff[imgno];                 // Zoff
        dataline(6)=img_th0[imgno];                  // Wedge th0
        dataline(7)=img_thF[imgno];                  // Wedge thF
        dataline(8)=(double)(opt_refno+1);           // Ref
        dataline(9)=Pmax;                            // P_max/P_tot
        dataline(10)=opt_scale;                      // internal scale
        DFo.append_comment(fn);
        DFo.append_data_line(dataline);
        avePmax += Pmax;

        imgno++;
        //if (verb > 0) progress_bar(imgno);

    }
    //if (verb > 0) progress_bar(nn);
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
                                     double * dataWsumScale,
                                     double * dataWsumScale2,
                                     double & sumw_allrefs,
                                     double & avePmax)
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
        if (do_impute)
        {
            if (dataSumWRefs[refno] > 0.)
            {
                for (int i = 0; i < hsize; i++)
                {
                    int ii = refno * hsize + i;
                    // Impute old reference for missing pixels
                    //if (i==0) std::cerr<<abs(DATAREFS[ii])<<" "<<dataWsumWedsPerRef[ii]<<" "<<dataSumWRefs[refno]<<" "<<abs(DATAWSUMREFS[ii])/ dataSumWRefs[refno]<<" ";
                    DATAREFS[ii] *= 1. - (dataWsumWedsPerRef[ii] / dataSumWRefs[refno]);
                    // And sum the weighted sum for observed pixels
                    
                    //DATAREFS[ii] += DATAWSUMREFS[ii] / dataSumWRefs[refno];
                    DATAREFS[ii] += DATAWSUMREFS[ii] / dataWsumScale2[refno];
                    //if (i==0) std::cerr<<" new= "<<abs(DATAREFS[ii]) <<std::endl;
                }
            }
            else
            {
                // Do nothing (i.e. impute the old DATAREFS)
                /*
                  for (int i = 0; i < hsize; i++)
                  {
                  int ii = refno * hsize + i;
                  DATAREFS[ii] = 0.;
                  }
                */
            }
        }
        else // No imputation: divide by number of times a pixel has been observed 
        { 
            for (int i = 0; i < hsize; i++)
            {
                int ii = refno * hsize + i;
                if (dataWsumWedsPerRef[ii] > OBSERVED_THRESHOLD)
                {
                    DATAREFS[ii] = DATAWSUMREFS[ii] / dataWsumWedsPerRef[ii];
                }
                else
                {
                    DATAREFS[ii] = 0.;
                }
            }
        }
    }

    // Adjust average scale
    if (do_scale) 
    {
        for (int refno = 0; refno < nr_ref; refno++)
        {
            average_scale += dataWsumScale[refno];
            if (dataWsumScale[refno]>0.)
            {
                refs_avgscale[refno] =  dataWsumScale[refno] / dataSumWRefs[refno];
                for (int i = 0; i < hsize; i++)
                {
                    int ii = refno * hsize + i;
                    DATAREFS[ii] *= refs_avgscale[refno];
                }
            }
            else
                refs_avgscale[refno] = 1.;
        }
        average_scale /= sumw_allrefs;
    }

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
    
    // TODO replace sumw_allrefs by number of images per group!!!
    for (int ig = 0; ig < nr_group; ig++)
    {
        for (int i = 0; i < hsize; i++)
        {
            int ii = ig * hsize + i;
            if (do_impute)
            {
                // Return from two*SIGMA^2
                dataSigma[ii] /= 2;
                // Impute old sigma values for missing pixels
                dataSigma[ii] *= 1. - (dataWsumWedsPerGroup[ii] / nr_imgs_per_group(ig));
                //  And sum the weighted sum for observedpixels
                dataSigma[ii] += dataWsumDist[ii] / nr_imgs_per_group(ig);
            }
            else
            {
                if (dataWsumWedsPerGroup[ii] > OBSERVED_THRESHOLD)
                    dataSigma[ii] = dataWsumDist[ii] / dataWsumWedsPerGroup[ii];
                else
                    dataSigma[ii] = 0.;
            }
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

        // If I am not imputing, then just taking a radial average here will go WRONG!!
        if (do_ravg_sigma) 
            forwfftw.fftwRadialAverage(ave, sigma_noise[ig], radial_count, true, true);

        // Now store again in dataSigma structure
        for (int i = 0, ii=0; i< fftw_hsize; i++)
            if (is_in_range[i]) 
            {
                dataSigma[ig * hsize + ii] = 2. * ave[i]; // Store again TWO SIGMA^2
                ii++;
            }

    }

    // Apply regularisation
    regularise(dataRefs,dataSigma);

    // Average Pmax
    avePmax /= sumw_allrefs;

#ifdef DEBUG
    std::cerr<<"done maximization"<<std::endl;
#endif
}


void Prog_mlf_tomo_prm::regularise(double * dataRefs,
                                   double * dataSigma)
{
#ifdef DEBUG
    std::cerr<<"start regularise"<<std::endl;
#endif

    // Do a kerdenSOM like regularization (see Pascual et al., JSB, 2001)
    if (reg > 0.)
    {
        double * regRef, *regSigma;
        // ??? Divide by nr_ref*nr_ref??
        //double reg_norm = reg * nr_imgs / (nr_ref * nr_ref);
        double reg_norm = reg * nr_imgs / nr_ref;
        std::complex<double> * REGREF, * DATAREFS;
        try
        {
            regRef = new double [nr_ref * size];
            regSigma = new double [nr_ref * hsize];
        }
        catch (std::bad_alloc&)
        {
            REPORT_ERROR(1,"Error allocating memory in regularise");
        }
        REGREF = (std::complex<double> *)regRef;
        DATAREFS = (std::complex<double> *)dataRefs;
        // Initialize
        for (int i=0; i< nr_ref * size; i++)
            regRef[i] = 0.; 
        for (int i=0; i< nr_ref * hsize; i++)
            regSigma[i] = 0.; 

        // First, calculate averages of updated references and the squared distance between them
        for (int refno = 0; refno < nr_ref; refno++)
        {
            for (int refno2 = 0; refno2 < nr_ref; refno2++)
            {
                for (int i=0; i< size; i++)
                    regRef[refno*size + i] += dMij(som_reg_matrix,refno,refno2) * dataRefs[refno2*size + i];
                        
                for (int i=0; i< hsize; i++)
                {
                    double aux = abs(DATAREFS[refno*hsize + i] - DATAREFS[refno2*hsize + i]);
                    regSigma[refno*hsize + i] += dMij(som_reg_matrix,refno,refno2) * aux*aux;
                }
            }
        }

        // Then, update the references
        for (int refno = 0; refno < nr_ref; refno++)
        {
            double sumw = alpha_k[refno] * nr_imgs;
//#define DEBUG_REGULARISE
#ifdef DEBUG_REGULARISE
            std::cerr<<"refno= "<<refno<<" sumw = "<<sumw<<" reg_norm= "<<reg_norm<<" data11= "<<DATAREFS[refno * hsize + 11]<<" reg11= "<<REGREF[refno * hsize + 11];
#endif
            for (int i = 0; i < hsize; i++)
            {
                int ii = refno * hsize + i;
                DATAREFS[ii] = DATAREFS[ii] * sumw + reg_norm * REGREF[ii];
                DATAREFS[ii] /= sumw + reg_norm;
            }
#ifdef DEBUG_REGULARISE
            std::cerr<<" reged= "<<DATAREFS[refno * hsize + 11]<<std::endl;
#endif
        }
        // and the noise spectra
        for (int ig = 0; ig < nr_group; ig++)
            for (int i = 0; i < hsize; i++)
            {
                //int ii = ig * hsize + i;
                //double sumw = nr_imgs_per_group(ig);
                //dataSigma[ii] = dataSigma[ii] * sumw + reg_norm * regSigma[i];
                //dataSigma[ii] /= sumw;
                int ii = ig * hsize + i;
                dataSigma[ii] += (reg_norm / nr_imgs) * regSigma[ii];
            }
    }

#ifdef DEBUG
    std::cerr<<"done regularise"<<std::endl;
#endif

}

    // Convergence check
bool Prog_mlf_tomo_prm::checkConvergence(double * dataRefs,
                                         double * oldDataRefs)
{

    double * dataOut;
    double diff2;
    VolumeXmipp vol;
    bool converged = true;
    vol().resize(Xdim,Ydim,Zdim);

    try
    {
        dataOut  = new double[2*fftw_hsize];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in maximization");
    }

    // Do convergence check in real space...
    // Calculate difference volumes for all references
    for (int refno = 0; refno < nr_ref; refno++)
    {
        // Set points within resolution range
        for (int i = 0, ii=0; i< 2*fftw_hsize; i++)
        {
            if (is_in_range[i/2]) 
            {
                dataOut[i] = oldDataRefs[refno*size + ii] - dataRefs[refno*size + ii];
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
        diff2 = vol().sum2() / dim3;
#define DEBUG_CONVERGENCE
#ifdef DEBUG_CONVERGENCE
        FileName fnt;
        fnt.compose("diff",refno,"vol");
        vol.write(fnt);
        //std::cerr<<"diff volume= "<<fnt<<" diff2= "<<diff2<<" eps= "<<eps<<std::endl;
#endif
        if (diff2 > eps)
        {
            converged = false;
            break;
        }
    }

    // Set oldDataRefs again to dataRefs
    for (int refno = 0; refno < nr_ref; refno++)
        for (int i = 0; i < size; i++)
        {
            int ii = refno * size + i;
            oldDataRefs[ii] = dataRefs[ii];
        }

    return converged;

}        


void Prog_mlf_tomo_prm::writeOutputFiles(int step,
                                         int iter, 
                                         double  * dataRefs,
                                         double  * dataWsumWedsPerRef,
                                         DocFile & DFo,
                                         double  & sumw_allrefs, 
                                         double  & LL, 
                                         double  & avePmax)
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
    if (step >= 0)
    {
        fn_base += "_step";
        fn_base.compose(fn_base, step, "");
    }
    if (iter >= 0)
    {
        fn_base += "_it";
        fn_base.compose(fn_base, iter, "");
    }

    // Do backward FFTW to write out real-space maps again
    double * dataOut;
    try
    {
        dataOut  = new double[2*fftw_hsize];
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
        DFl.insert_comment(fn_tmp);
        DFl.insert_data_line(fracline);

        if (iter != 0)
        {
            // Also write out sum of wedges
            vol().initZeros();
            for(int i=0,ii=0,iii=0;i<Zdim;i++)
                for(int j=0;j<Ydim;j++)
                    for(int k=0;k<Xdim/2 + 1; k++,ii++)
                        if (is_in_range[ii])
                        {
                            vol(i,j,k)=dataWsumWedsPerRef[refno*hsize + iii];
                            iii++;
                        }
            CenterFFT(vol(),true);
            fn_tmp = fn_base + "_ref";
            fn_tmp.compose(fn_tmp, refno + 1, "");
            fn_tmp = fn_tmp + "_sumwedge.vol";
            vol.write(fn_tmp);
        }
    }
    
    // Write out sel & log-file
    fn_tmp = fn_base + ".sel";
    SFo.write(fn_tmp);

    DFl.go_beginning();
    comment = "MLF_tomo: nr. imgs= " + floatToString(sumw_allrefs,8,8);
    comment += " LL= " + floatToString(LL, 13, 10) + " <Pmax>= " + floatToString(avePmax, 8, 5) + " reg= " + floatToString(reg, 5, 3);
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

    if (iter != 0)
    {
        // Write out docfile with optimal transformation & references
        fn_tmp=fn_base + ".doc";
        DFo.write(fn_tmp);
        
        // Write out class selfiles
        for (int refno = 0;refno < nr_ref; refno++)
        {
            DFo.go_beginning();
            SFo.clear();
            for (int n = 0; n < DFo.dataLineNo(); n++)
            {
                DFo.next();
                fn_tmp = ((DFo.get_current_line()).get_text()).erase(0, 3);
                DFo.adjust_to_data_line();
                if ((refno + 1) == (int)DFo(8)) SFo.insert(fn_tmp, SelLine::ACTIVE);
            }
            fn_tmp = fn_base + "_class";
            fn_tmp.compose(fn_tmp, refno + 1, "sel");
            SFo.write(fn_tmp);
        }
    }

#ifdef DEBUG
    std::cerr<<"done writeOutputFiles"<<std::endl;
#endif
}

