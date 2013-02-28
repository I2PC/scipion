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

#include "tomo_align_dual_tilt_series.h"
#include <data/args.h>
#include <data/filters.h>
#include <data/xmipp_fftw.h>

double wrapperDualAligment(double *p, void *prm)
{
    ProgAlignDualTiltSeries *eprm=(ProgAlignDualTiltSeries *)prm;
    Matrix1D<double> aux(6);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(aux)
    aux(i)=p[i+1];
    return eprm->evaluateAlignment(aux);
}

/// Read parameters
void ProgAlignDualTiltSeries::readParams()
{
    fnRef=getParam("--ref");
    fnDual=getParam("--dual");
    fnOut=getParam("-o");
    scaleFactor=getDoubleParam("--scale");
}

/// Show parameters
void ProgAlignDualTiltSeries::show()
{
    std::cout
    << "Reference: " << fnRef       << std::endl
    << "Dual:      " << fnDual      << std::endl
    << "Output:    " << fnOut       << std::endl
    << "Scale:     " << scaleFactor << std::endl
    ;
}

/// Usage
void ProgAlignDualTiltSeries::defineParams()
{
    addUsageLine("Align two dual tilt series that have been previously internally aligned");
    addUsageLine("with xmipp_tomo_align_tilt_series");
    addSeeAlsoLine("tomo_align_tilt_series");
    addParamsLine("  --ref  <selfile>  : Reference tilt series");
    addParamsLine("                    : The selfile must contain the list of micrographs");
    addParamsLine("                    : and its tilt angles");
    addParamsLine("  --dual <selfile>  : Dual tilt series");
    addParamsLine(" [-o <rootname=\"\">]  : Rootname for the aligned tilt series");
    addParamsLine("                    : If not given, the dual tilt series+_aligned");
    addParamsLine(" [--scale <s=0.25>] : Scale for performing the common line comparisons");
}

/// Read dual series
void ProgAlignDualTiltSeries::readDual()
{
    imgDual.clear();
    tiltDual.initZeros(SFDual.size());
    int i=0;
    Image<double> I;
    FileName fnImg;
    double minAbsTilt=1000;
    FOR_ALL_OBJECTS_IN_METADATA(SFDual)
    {
        SFDual.getValue(MDL_IMAGE,fnImg, __iter.objId);
        SFDual.getValue(MDL_ANGLE_TILT,tiltDual(i), __iter.objId);
        if (fabs(tiltDual(i))<minAbsTilt)
        {
            minAbsTilt=fabs(tiltDual(i));
            fnDual0=fnImg;
        }
        I.read(fnImg);
        selfScaleToSize(BSPLINE3,I(),ROUND(YSIZE(I())*scaleFactor),
                        ROUND(XSIZE(I())*scaleFactor));
        I().setXmippOrigin();
        Xdim=XSIZE(I());
        Ydim=YSIZE(I());
        imgDual.push_back(I());
        i++;
    }
}

/// Produce side info
void ProgAlignDualTiltSeries::produceSideInfo()
{
	debugging=false;
    if (fnOut=="")
        fnOut=fnDual.withoutExtension()+"_aligned";

    std::cout << "Reading data...\n";
    SFRef.read(fnRef);
    SFDual.read(fnDual);

    // Read Reference series
    tiltRef.initZeros(SFRef.size());
    int i=0;
    Image<double> I;
    FileName fnImg;
    double minAbsTilt=1000;
    FOR_ALL_OBJECTS_IN_METADATA(SFRef)
    {
        SFRef.getValue(MDL_IMAGE,fnImg,__iter.objId);
        SFRef.getValue(MDL_ANGLE_TILT,tiltRef(i),__iter.objId);
        if (fabs(tiltRef(i))<minAbsTilt)
        {
            minAbsTilt=fabs(tiltRef(i));
            fnRef0=fnImg;
        }
        I.read(fnImg);
        selfScaleToSize(BSPLINE3,I(),ROUND(YSIZE(I())*scaleFactor),
                        ROUND(XSIZE(I())*scaleFactor));
        I().setXmippOrigin();
        imgRef.push_back(I());
        i++;
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
void ProgAlignDualTiltSeries::findParametersAt0degrees(bool rotateDual)
{
    Image<double> Iref, Idual;
    Iref.read(fnRef0);
    Idual.read(fnDual0);
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
        std::cout << "First estimate (180) of (rot,tilt,psi,x,y,z)=\n"
        << alignment.transpose() << std::endl;
    else
        std::cout << "First estimate (0) of (rot,tilt,psi,x,y,z)=\n"
        << alignment.transpose() << std::endl;
    /*
            debugging=true;
    */
}
#undef DEBUG

/// Distance between a pair of common lines
//#define DEBUG
double ProgAlignDualTiltSeries::distanceBetweenCommonLines(
    int refi, int dualj, const Matrix2D<double> &E,
    double X, double Y, double Z)
{
    SPEED_UP_temps012;

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
    debugging=true;
    Image<double> save;
    save()=I;
    save.write("PPPref.xmp");
#endif

    selfRotate(LINEAR,I, -angi, DONT_WRAP);
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

    if (debugging)
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
        << "   corr=" << correlationIndex(profilei,profilej) << std::endl
        ;

        profilei.write("PPPi.txt");
        profilej.write("PPPj.txt");

        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
    }
    return 1-correlationIndex(profilei,profilej);
}
#undef DEBUG

/// Evaluate alignment
double ProgAlignDualTiltSeries::evaluateAlignment(const Matrix1D<double> &_alignment)
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
double ProgAlignDualTiltSeries::optimizeAlignment()
{
    Matrix1D<double> steps;
    steps.resize(alignment);
    steps.initConstant(1);
    double fitness;
    int iter;
    powellOptimizer(alignment,1,6,&wrapperDualAligment,this,0.01,fitness,
                    iter,steps,true);
#ifdef DEBUG

    debugging=true;
    evaluateAlignment(alignment);
    debugging=false;
#endif

    return fitness;
}
#undef DEBUG

/// Align dual
void ProgAlignDualTiltSeries::alignDual()
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
    Matrix1D<double> shift3D=vectorR3(alignment(3),alignment(4),alignment(5));
    Matrix1D<double> shift2D=vectorR2(alignment(3),alignment(4));
    shift3D/=scaleFactor;
    shift2D/=scaleFactor;
    MetaData SFout;
    int n=0;
    //std::cout << "Aligning dual" << std::endl;
    Image<double> Idual;
    Matrix2D<double> Edual;
    FileName fnImg, fn;
    FOR_ALL_OBJECTS_IN_METADATA(SFDual)
    {
        SFDual.getValue(MDL_IMAGE,fnImg,__iter.objId);
        Idual.read(fnImg);
        Idual().setXmippOrigin();
        if (rotatedDual)
            selfRotate(BSPLINE3,Idual(),180,DONT_WRAP);
        selfTranslate(BSPLINE3,Idual(),shift2D,DONT_WRAP);
        shiftProjectionInZ(Idual(), n, ZZ(shift3D));
        Euler_angles2matrix(0, tiltDual(n), 0, Edual);
        Edual=Edual*E;
        double rot, tilt, psi;
        Euler_matrix2angles(Edual, rot, tilt, psi);
        fn.compose(n+1,fnOut+".stk");
        Idual.write(fn);
        size_t id = SFout.addObject();
        SFout.setValue(MDL_IMAGE, fn, id);
        SFout.setValue(MDL_ANGLE_ROT, rot, id);
        SFout.setValue(MDL_ANGLE_TILT, tilt, id);
        SFout.setValue(MDL_ANGLE_PSI, psi, id);
        n++;
    }
    SFout.write(fnOut.removeDirectories()+".sel");
}

/// Run
void ProgAlignDualTiltSeries::run()
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

void ProgAlignDualTiltSeries::shiftProjectionInZ(MultidimArray<double> &I, int dualj, double Z) const
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
        SPEED_UP_temps012;
        VECTOR_R2(idx,j,i);
        FFT_idx2digfreq(I,idx,w);
        ZZ(w)=0;
        M3x3_BY_V3x1(w,Edualt,w);
        double arg=ZZ(w)*Z;
        A2D_ELEM(Ifft,i,j)*=std::complex<double>(cos(arg),sin(arg));
    }
    transformer.inverseFourierTransform();
}
