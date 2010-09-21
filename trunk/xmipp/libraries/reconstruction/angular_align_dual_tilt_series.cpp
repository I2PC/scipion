/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *
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

#include "angular_align_dual_tilt_series.h"
#include <data/args.h>
#include <data/filters.h>
#include <data/fftw.h>

double wrapperDualAligment(double *p, void *prm)
{
    Prog_align_dual *eprm=(Prog_align_dual *)prm;
    Matrix1D<double> aux(6);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(aux)
    aux(i)=p[i+1];
    return eprm->evaluateAlignment(aux);
}

/// Read parameters
void Prog_align_dual::read(int argc, char **argv)
{
    fnRef=getParameter(argc,argv,"-ref");
    fnDual=getParameter(argc,argv,"-dual");
    fnOut=getParameter(argc,argv,"-o","");
    scaleFactor=textToFloat(getParameter(argc,argv,"-scale","0.25"));
    verbose=false;
}

/// Show parameters
void Prog_align_dual::show()
{
    std::cout << "Reference: " << fnRef       << std::endl
    << "Dual:      " << fnDual      << std::endl
    << "Output:    " << fnOut       << std::endl
    << "Scale:     " << scaleFactor << std::endl
    ;
}

/// Usage
void Prog_align_dual::usage()
{
    std::cout << "Usage:\n"
    << "  -ref  <selfile>  : Reference tilt series\n"
    << "  -dual <selfile>  : Dual tilt series\n"
    << " [-o <rootname>]   : Rootname for the aligned tilt series\n"
    << " [-scale <s=0.25>] : Scale for performing the common line comparisons\n"
    ;
}

/// Read dual series
void Prog_align_dual::readDual()
{
    imgDual.clear();
    tiltDual.initZeros(SFDual.size());
    int i=0;
    FOR_ALL_OBJECTS_IN_METADATA(SFDual)
    {
        Image<double> I;
        FileName fnImg;
        SFDual.getValue(MDL_IMAGE,fnImg);
        I.read(fnImg);
        tiltDual(i++)=I.tilt();
        selfScaleToSize(BSPLINE3,I(),ROUND(YSIZE(I())*scaleFactor),
                        ROUND(XSIZE(I())*scaleFactor));
        I().setXmippOrigin();
        Xdim=XSIZE(I());
        Ydim=YSIZE(I());
        imgDual.push_back(I());
    }
}

/// Produce side info
void Prog_align_dual::produceSideInfo()
{
    if (fnOut=="")
        fnOut=fnDual.withoutExtension();

    std::cout << "Reading data...\n";
    SFRef.read(fnRef);
    SFDual.read(fnDual);

    // Read Reference series
    tiltRef.initZeros(SFRef.size());
    int i=0;
    FOR_ALL_OBJECTS_IN_METADATA(SFRef)
    {
        Image<double> I;
        FileName fnImg;
        SFRef.getValue(MDL_IMAGE,fnImg);
        I.read(fnImg);
        tiltRef(i++)=I.tilt();
        selfScaleToSize(BSPLINE3,I(),ROUND(YSIZE(I())*scaleFactor),
                        ROUND(XSIZE(I())*scaleFactor));
        I().setXmippOrigin();
        imgRef.push_back(I());
    }

    // Read Dual series
    readDual();

    normali.resize(3);
    normalj.resize(3);
    commonline.resize(3);
    commonlinei.resize(3);
    commonlinej.resize(3);
    commonlinejE.resize(3);
    A.initIdentity(3);
    int minDim=XMIPP_MIN(Xdim,Ydim);
    profilei.resize(minDim);
    profilei.setXmippOrigin();
    profilej=profilei;
}

/// Find parameters (shift+rotation) at 0 degrees
//#define DEBUG
void Prog_align_dual::findParametersAt0degrees(bool rotateDual)
{
    // Look for the images at 0 degrees
    int idxRef0=-1;
    double minAbsTilt=1000;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(tiltRef)
    if (ABS(tiltRef(i))<minAbsTilt)
    {
        minAbsTilt=ABS(tiltRef(i));
        idxRef0=i;
    }

    // Look for the images at 0 degrees
    int idxDual0=-1;
    minAbsTilt=1000;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(tiltDual)
    if (ABS(tiltDual(i))<minAbsTilt)
    {
        minAbsTilt=ABS(tiltDual(i));
        idxDual0=i;
    }

    FileName fnRef, fnDual;
    SFRef.getValue(MDL_IMAGE,fnRef,idxRef0);
    SFDual.getValue(MDL_IMAGE,fnDual,idxDual0);
    Image<double> Iref;
    Iref.read(fnRef);
    Image<double> Idual;
    Idual.read(fnDual);
    Iref().setXmippOrigin();
    Idual().setXmippOrigin();
    if (rotateDual)
        selfRotate(BSPLINE3,Idual(),180);
    Matrix2D<double> M0,M0inv;
#ifdef DEBUG

    Image<double> save;
    save()=Iref();
    save.write("PPPref0.xmp");
    save()=Idual();
    save.write("PPPdual0.xmp");
#endif

    alignImages(Iref(), Idual(), M0);
    M0inv=M0.inv();
#ifdef DEBUG

    save()=Idual();
    save.write("PPPdual0Aligned.xmp");
    std::cout << "M0\n" << M0 << std::endl;
    std::cout << "M0inv\n" << M0inv << std::endl;
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif

    alignment.resize(6);
    alignment(0)=-RAD2DEG(atan2(M0inv(1,0),M0inv(1,1)));  // rot
    alignment(1)=0;                                // tilt
    alignment(2)=0;                                // psi
    alignment(3)=-M0inv(0,2)*scaleFactor;          // X
    alignment(4)=-M0inv(1,2)*scaleFactor;          // Y
    alignment(5)=0;                                // Z

    /*
            alignment(0)=90;
            alignment(3)=-25*scaleFactor;          // X
            alignment(4)= 15*scaleFactor;          // Y
            alignment(5)=5;
    */
    if (rotateDual)
        std::cout << "First estimate (0) of (rot,tilt,psi,x,y,z)=\n"
        << alignment.transpose() << std::endl;
    else
        std::cout << "First estimate (180) of (rot,tilt,psi,x,y,z)=\n"
        << alignment.transpose() << std::endl;
    /*
            verbose=true;
    */
}
#undef DEBUG

/// Distance between a pair of common lines
//#define DEBUG
double Prog_align_dual::distanceBetweenCommonLines(
    int refi, int dualj, const Matrix2D<double> &E,
    double X, double Y, double Z)
{
    SPEED_UP_temps;

    // Compute the direction of the common line in the
    // universal coordinate system
    Euler_direction(0, tiltRef(refi), 0, normali);
    Euler_direction(0, tiltDual(dualj), 0, normalj);
    M3x3_BY_V3x1(normalj, Et, normalj);

    vectorProduct(normali,normalj,commonline);
    if (commonline.module()<XMIPP_EQUAL_ACCURACY)
        return -2;
    commonline.selfNormalize();

    // Compute the direction of the common line in each
    // of the projection coordinate systems
    Uproject_to_plane(commonline, 0, tiltRef(refi), 0, commonlinei);
    M3x3_BY_V3x1(commonlinejE, E, commonline);
    Uproject_to_plane(commonlinejE, 0, tiltDual(dualj), 0, commonlinej);
    commonlinei.selfNormalize();
    commonlinej.selfNormalize();

    // Compute the angle of the common line in i and j images
    double angi=-RAD2DEG(atan2(YY(commonlinei),XX(commonlinei)));
    double angj=-RAD2DEG(atan2(YY(commonlinej),XX(commonlinej)));

    I=imgRef[refi];
#ifdef DEBUG

    Image<double> save;
    save()=I;
    save.write("PPPref.xmp");
#endif

    selfRotate(LINEAR,I, angi, DONT_WRAP);
    profilei.initZeros();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    A1D_ELEM(profilei,j)+=A2D_ELEM(I,i,j);
#ifdef DEBUG

    save()=I;
    save.write("PPPrefaligned.xmp");
#endif

    I=imgDual[dualj];
#ifdef DEBUG

    save()=I;
    save.write("PPPdual.xmp");
#endif

    // Modify Idual according to z
    shiftProjectionInZ(I, dualj, Z/scaleFactor);

    // Modify Idual according to x,y,angj
    // First translation, then rotation
    double cj=COSD(angj);
    double sj=SIND(angj);
    MAT_ELEM(A,0,0)=MAT_ELEM(A,1,1)=cj;
    MAT_ELEM(A,1,0)=sj;
    MAT_ELEM(A,0,1)=-sj;
    MAT_ELEM(A,0,2)=cj*X-sj*Y;
    MAT_ELEM(A,1,2)=sj*X+cj*Y;
    selfApplyGeometry(LINEAR, I, A, IS_NOT_INV, DONT_WRAP);
    profilej.initZeros();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    A1D_ELEM(profilej,j)+=A2D_ELEM(I,i,j);
#ifdef DEBUG

    save()=I;
    save.write("PPPdualAligned.xmp");
#endif

    if (verbose)
    {
        std::cout << "refi=" << refi << " dualj=" << dualj << std::endl
        << "   tiltRef=" << tiltRef(refi) << " tiltDual=" << tiltDual(dualj) << std::endl
        << "   normali=" << normali.transpose() << std::endl
        << "   normalj=" << normalj.transpose() << std::endl
        << "   commonline=" << commonline.transpose() << std::endl
        << "   commonlinei=" << commonlinei.transpose() << std::endl
        << "   commonlinejE=" << commonlinejE.transpose() << std::endl
        << "   commonlinej=" << commonlinej.transpose() << std::endl
        << "   angi=" << angi << " angj=" << angj << std::endl
        << "   A\n" << A << std::endl
        << "   corr=" << correlation_index(profilei,profilej) << std::endl
        ;

        profilei.write("PPPi.txt");
        profilej.write("PPPj.txt");

        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
    }
    return 1-correlation_index(profilei,profilej);
}
#undef DEBUG

/// Evaluate alignment
double Prog_align_dual::evaluateAlignment(const Matrix1D<double> &_alignment)
{
    Euler_angles2matrix(_alignment(0),_alignment(1),_alignment(2),E);
    Et=E.transpose();
    double X=_alignment(3);
    double Y=_alignment(4);
    double Z=_alignment(5);

    int Nref=imgRef.size();
    int Ndual=imgDual.size();
    int Ncomparisons=0;
    double avgDistance=0;
    for (int nref=0; nref<Nref; nref++)
        for (int ndual=0; ndual<Ndual; ndual++)
        {
            double d=distanceBetweenCommonLines(nref,ndual,E,X,Y,Z);
            if (d>0)
            {
                avgDistance+=d;
                Ncomparisons++;
            }
        }
    if (Ncomparisons>0)
        avgDistance/=Ncomparisons;
    return avgDistance;
}

//#define DEBUG
double Prog_align_dual::optimizeAlignment()
{
    Matrix1D<double> steps;
    steps.resize(alignment);
    steps.initConstant(1);
    double fitness;
    int iter;
    powellOptimizer(alignment,1,6,&wrapperDualAligment,this,0.01,fitness,
                    iter,steps,true);
#ifdef DEBUG

    verbose=true;
    evaluateAlignment(alignment);
    verbose=false;
#endif

    return fitness;
}
#undef DEBUG

/// Align dual
void Prog_align_dual::alignDual()
{
/* There might be a problem, the code seemingly working in TomoJ is

        for(int i=0;i<ts.length;i++){ // ts is an array of tilt series
            double[] ali=ts[i].getGlobalOrientation(); // I get the
                alignment computed between the tilt series with common lines
            if(ali!=null){ // is there an alignment?
                double cxtmp=ali[3]; //this is the shift of the center of
                                        volume original is (0,0,0)
                double cytmp=ali[4];
                //Now where will be this center of tilt axis in the volume?
                DoubleMatrix2D et= TiltSeries.eulerAngles2Matrix(ali[0],
                    ali[1], ali[2]).viewDice(); // I compute the euler angles of the rotation
                                                    between the tilt series
                cx[i]=cxtmp*et.getQuick(0,0)+cytmp*et.getQuick(0,1)+width/2;
                // I compute the new position of the center of the volume + add width/2 for
                    correct indexation of arrays (in Xmipp should not be there)
 
                cy[i]=cxtmp*et.getQuick(1,0)+cytmp*et.getQuick(1,1)+height/2;


            }else{ //there is no alignment so I define it as classical
                cx[i]=width/2;
                cy[i]=height/2;
            }
        }
*/

    Euler_angles2matrix(alignment(0),alignment(1),alignment(2),E);
    int Ndual=imgDual.size();
    Matrix1D<double> shift3D=vectorR3(alignment(3),alignment(4),alignment(5));
    Matrix1D<double> shift2D=vectorR2(alignment(3),alignment(4));
    shift3D/=scaleFactor;
    shift2D/=scaleFactor;
    MetaData SFout;
    int n=0;
    //std::cout << "Aligning dual" << std::endl;
    SFDual.write(std::cout);
    FOR_ALL_OBJECTS_IN_METADATA(SFDual)
    {
        Image<double> Idual;
        FileName fnImg;
        SFDual.getValue(MDL_IMAGE,fnImg);
        Idual.read(fnImg);
        Idual().setXmippOrigin();
        if (rotatedDual)
            selfRotate(BSPLINE3,Idual(),180);
        selfTranslate(BSPLINE3,Idual(),shift2D);
        shiftProjectionInZ(Idual(), n, ZZ(shift3D));
        Matrix2D<double> Edual;
        Euler_angles2matrix(0, tiltDual(n), 0, Edual);
        Edual=Edual*E;
        double rot, tilt, psi;
        Euler_matrix2angles(Edual, rot, tilt, psi);
        Idual.setRot((float)rot);
        Idual.setTilt((float)tilt);
        Idual.setPsi((float)psi);
        FileName fn=fnOut+"_aligned"+integerToString(n,4)+".xmp";
        Idual.write(fn);
        SFout.addObject();
        SFout.setValue(MDL_IMAGE,fn);
        n++;
    }
    SFout.write(fnOut.removeDirectories()+"_aligned.sel");
}

/// Run
void Prog_align_dual::run()
{
    produceSideInfo();

    // Check if there is a rotation of 180 degrees
    findParametersAt0degrees(false);
    double objective0=evaluateAlignment(alignment);
    Matrix1D<double> alignmentBackup=alignment;

    int Ndual=imgDual.size();
    for (int ndual=0; ndual<Ndual; ndual++)
        selfRotate(BSPLINE3,imgDual[ndual],180);
    findParametersAt0degrees(true);
    double objective180=evaluateAlignment(alignment);
    std::cout << "objective0=" << objective0 << " objective180=" << objective180 << std::endl;
    if (objective0<objective180)
    {
        readDual();
        rotatedDual=false;
        alignment=alignmentBackup;
    }
    else
        rotatedDual=true;

    optimizeAlignment();
    alignDual();
}

void Prog_align_dual::shiftProjectionInZ(MultidimArray<double> &I, int dualj, double Z) const
{
    if (Z==0)
        return;
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > Ifft;
    Matrix1D<double> w(3);
    Matrix1D<int> idx(3);
    Matrix2D<double> Edual, Edualt;
    Euler_angles2matrix(0, tiltDual(dualj), 0, Edual);
    Edualt=Edual.transpose();
    transformer.FourierTransform(I,Ifft,false);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(Ifft)
    {
        SPEED_UP_temps;
        VECTOR_R2(idx,j,i);
        FFT_idx2digfreq(I,idx,w);
        ZZ(w)=0;
        M3x3_BY_V3x1(w,Edualt,w);
        double arg=ZZ(w)*Z;
        Ifft(i,j)*=std::complex<double>(cos(arg),sin(arg));
    }
    transformer.inverseFourierTransform();
}
