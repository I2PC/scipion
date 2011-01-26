/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano coss@cnb.csic.es
 *              Slavica Jonic        Slavica.Jonic@impmc.jussieu.fr
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

#include "reconstruct_art_pseudo.h"
#include <fstream>
#include <data/histogram.h>

#define FORWARD   1
#define BACKWARD -1

/** Projection of a pseudoatom volume */
void project_Pseudo(const std::vector< Matrix1D<double> > &atomPosition,
                    std::vector<double> &atomWeight, double sigma,
                    MultidimArray<double> &proj, MultidimArray<double> &norm_proj,
                    Matrix2D<double> &Euler, double shiftX, double shiftY,
                    const std::vector<double> &lambda,
                    const std::vector< Matrix2D<double> > &NMA,
                    int direction, const Matrix1D<double> &gaussianProjectionTable,
                    const Matrix1D<double> &gaussianProjectionTable2)
{
    // Project all pseudo atoms ............................................
    int nmax=atomPosition.size();
    Matrix1D<double> actprj(3);
    double sigma4=4*sigma;
    Matrix1D<double> actualAtomPosition;
    int lambdaSize=lambda.size();
    for (int n=0; n<nmax; n++)
    {
        actualAtomPosition=atomPosition[n];
        for (int mode=0; mode<lambdaSize; mode++)
        {
            const Matrix2D<double> &NMAmode=NMA[mode];
            double lambdam=lambda[mode];
            FOR_ALL_ELEMENTS_IN_MATRIX1D(actualAtomPosition)
            VEC_ELEM(actualAtomPosition,i)+=
                lambdam*MAT_ELEM(NMAmode,n,i);
        }

        Uproject_to_plane(actualAtomPosition, Euler, actprj);
        XX(actprj)+=shiftX;
        YY(actprj)+=shiftY;
#ifdef DEBUG

        bool condition = true;
        if (condition)
        {
            std::cout << "Projecting point " << atomPositions[n].transpose() << std::endl;
            std::cout << "Vol there = " << atomWeight[n] << std::endl;
            std::cout << " Center of the basis proj (2D) " << XX(actprj) << "," << YY(actprj) << std::endl;
        }
#endif

        // Search for integer corners for this basis
        int XX_corner1 = CEIL(XMIPP_MAX(STARTINGX(proj), XX(actprj) - sigma4));
        int YY_corner1 = CEIL(XMIPP_MAX(STARTINGY(proj), YY(actprj) - sigma4));
        int XX_corner2 = FLOOR(XMIPP_MIN(FINISHINGX(proj), XX(actprj) + sigma4));
        int YY_corner2 = FLOOR(XMIPP_MIN(FINISHINGY(proj), YY(actprj) + sigma4));

#ifdef DEBUG

        if (condition)
        {
            std::cout << "Clipped and rounded Corner 1 " << XX_corner1
            << " " << YY_corner1 << " " << std::endl;
            std::cout << "Clipped and rounded Corner 2 " << XX_corner2
            << " " << YY_corner2 << " " << std::endl;
        }
#endif

        // Check if the basis falls outside the projection plane
        if (XX_corner1 <= XX_corner2 && YY_corner1 <= YY_corner2)
        {
            double vol_corr=0;

            // Effectively project this basis
            for (int y = YY_corner1; y <= YY_corner2; y++)
            {
                double y_diff2=y-YY(actprj);
                y_diff2=y_diff2*y_diff2;
                for (int x = XX_corner1; x <= XX_corner2; x++)
                {
                    double x_diff2=x-XX(actprj);
                    x_diff2=x_diff2*x_diff2;
                    double r=sqrt(x_diff2+y_diff2);
                    int idx=ROUND(r*1000);
                    double a=gaussianProjectionTable(idx);
                    double a2=gaussianProjectionTable(idx);
#ifdef DEBUG

                    if (condition)
                        std::cout << "=" << a << " , " << a2;
#endif

                    if (direction==FORWARD)
                    {
                        A2D_ELEM(proj, y, x) += atomWeight[n] * a;
                        A2D_ELEM(norm_proj, y, x) += a2;
#ifdef DEBUG
                        if (condition)
                        {
                            std::cout << " proj= " << A2D_ELEM(proj, y, x)
                            << " norm_proj=" << A2D_ELEM(norm_proj, y, x) << std::endl;
                            std::cout.flush();
                        }
#endif
                    }
                    else
                    {
                        vol_corr += A2D_ELEM(norm_proj, y, x) * a;
#ifdef DEBUG
                        if (condition)
                        {
                            std::cout << " corr_img= " << A2D_ELEM(norm_proj, y, x)
                            << " correction=" << vol_corr << std::endl;
                            std::cout.flush();
                        }
#endif
                    }
                }
            }

            if (direction==BACKWARD)
            {
                atomWeight[n] += vol_corr;
#ifdef DEBUG
                if (condition)
                {
                    std::cout << "\nFinal value at ( " << n << ") = "
                    << atomWeight[n] << std::endl;
                }
#endif
            }
        } // If not collapsed
    }
}

void ProgARTPseudo::readParams()
{
    fnDoc=getParam("-i");
    fnPseudo=getParam("-pseudo");
    fnRoot=getParam("-o");
    lambdaART=getDoubleParam("-l");
    Nit=getIntParam("-n");
    sigma=getDoubleParam("-sigma");
    fnNMA=getParam("-nma");
    sampling=getDoubleParam("-sampling");
}

void ProgARTPseudo::show() const
{
    std::cout
    << "Input images:    " << fnDoc     << std::endl
    << "Pseudoatoms:     " << fnPseudo  << std::endl
    << "Sigma:           " << sigma     << std::endl
    << "Sampling rate:   " << sampling  << std::endl
    << "NMA:             " << fnNMA     << std::endl
    << "Output rootname: " << fnRoot    << std::endl
    << "Lambda ART:      " << lambdaART << std::endl
    << "N. Iterations:   " << Nit       << std::endl
    ;
}

void ProgARTPseudo::defineParams()
{
    addUsageLine("Perform 3D reconstruction using ART on pseudoatoms");
    addParamsLine("   -i <selfile>         : Input images");
    addParamsLine("   -pseudo <pseudofile> : Pseudo atomic structure (PDB format)");
    addParamsLine("   -o <rootname>        : Output rootname");
    addParamsLine("  [-sigma <s=-1>]       : Pseudoatom sigma. By default, from pseudo file");
    addParamsLine("  [-sampling <Ts=1>]    : Sampling rate in Angstroms/pixel");
    addParamsLine("  [-l <lambda=0.1>]     : Relaxation factor");
    addParamsLine("  [-n <N=1>]            : Number of iterations");
    addParamsLine("  [-nma <selfile=\"\">] : Selfile with NMA");
}

void ProgARTPseudo::produceSideInfo()
{
    DF.read(fnDoc);
    std::ifstream fhPseudo;
    fhPseudo.open(fnPseudo.c_str());
    if (!fhPseudo)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnPseudo);
    while (!fhPseudo.eof())
    {
        std::string line;
        getline(fhPseudo, line);
        if (line.length() == 0)
            continue;
        if (line.substr(7,13)=="fixedGaussian" && sigma<0)
        {
            std::vector < std::string> results;
            splitString(line," ",results);
            sigma=textToFloat(results[2]);
            sigma/=sampling;
        }
        else if (line.substr(0,4)=="ATOM")
        {
            Matrix1D<double> v(3);
            v(0)=textToFloat(line.substr(30,8));
            v(1)=textToFloat(line.substr(38,8));
            v(2)=textToFloat(line.substr(46,8));
            v/=sampling;
            atomPosition.push_back(v);
            atomWeight.push_back(0);
        }
    }
    fhPseudo.close();

    double sigma4=4*sigma;
    gaussianProjectionTable.resize(CEIL(sigma4*sqrt(2)*1000));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(gaussianProjectionTable)
    gaussianProjectionTable(i)=gaussian1D(i/1000.0,sigma);
    gaussianProjectionTable*=gaussian1D(0,sigma);
    gaussianProjectionTable2=gaussianProjectionTable;
    gaussianProjectionTable2*=gaussianProjectionTable;

    // NMA
    if (fnNMA!="")
    {
        MetaData DFNMA(fnNMA);
        FOR_ALL_OBJECTS_IN_METADATA(DFNMA)
        {
            Matrix2D<double> mode;
            mode.initZeros(atomPosition.size(),3);
            FileName fnMode;
            DFNMA.getValue( MDL_IMAGE, fnMode,__iter.objId);
            mode.read(fnMode);
            NMA.push_back(mode);
        }
    }
}

void ProgARTPseudo::run()
{
    produceSideInfo();
    Image<double> Iexp;
    for (int it=0; it<Nit; it++)
    {
        double itError=0;
        FOR_ALL_OBJECTS_IN_METADATA(DF)
        {
            FileName fnExp;
            DF.getValue( MDL_IMAGE, fnExp,__iter.objId);
            double rot;
            DF.getValue( MDL_ANGLEROT, rot,__iter.objId);
            double tilt;
            DF.getValue( MDL_ANGLETILT, tilt,__iter.objId);
            double psi;
            DF.getValue( MDL_ANGLEPSI, psi,__iter.objId);
            double shiftX;
            DF.getValue( MDL_SHIFTX, shiftX,__iter.objId);
            double shiftY;
            DF.getValue( MDL_SHIFTY, shiftY,__iter.objId);
            std::vector<double> lambda;
            if (NMA.size()>0)
                DF.getValue( MDL_NMA, lambda,__iter.objId);

            Iexp.read(fnExp);
            Iexp().setXmippOrigin();
            itError+=ART_single_step(Iexp(),rot,tilt,psi,shiftX,shiftY,lambda);
        }
        if (DF.size()>0)
            itError/=DF.size();
        std::cerr << "Error at iteration " << it << " = " << itError << std::endl;
    }
    writePseudo();
}

void ProgARTPseudo::writePseudo()
{
    // Convert from pseudoatoms to volume
    Image<double> V;
    size_t objId = DF.firstObject();
    FileName fnExp;
    DF.getValue( MDL_IMAGE, fnExp, objId);
    Image<double> I;
    I.read(fnExp,false);
    V().resize(XSIZE(I()),XSIZE(I()),XSIZE(I()));
    V().setXmippOrigin();

    int nmax=atomPosition.size();
    double sigma4=4*sigma;
    for (int n=0; n<nmax; n++)
    {
        int XX_corner1 = CEIL(XMIPP_MAX(STARTINGX(V()), XX(atomPosition[n]) - sigma4));
        int YY_corner1 = CEIL(XMIPP_MAX(STARTINGY(V()), YY(atomPosition[n]) - sigma4));
        int ZZ_corner1 = CEIL(XMIPP_MAX(STARTINGY(V()), ZZ(atomPosition[n]) - sigma4));
        int XX_corner2 = FLOOR(XMIPP_MIN(FINISHINGX(V()), XX(atomPosition[n]) + sigma4));
        int YY_corner2 = FLOOR(XMIPP_MIN(FINISHINGY(V()), YY(atomPosition[n]) + sigma4));
        int ZZ_corner2 = FLOOR(XMIPP_MIN(FINISHINGY(V()), ZZ(atomPosition[n]) + sigma4));
        if (XX_corner1 <= XX_corner2 && YY_corner1 <= YY_corner2 &&
            ZZ_corner1 <= ZZ_corner2)
        {
            for (int z = ZZ_corner1; z <= ZZ_corner2; z++)
                for (int y = YY_corner1; y <= YY_corner2; y++)
                    for (int x = XX_corner1; x <= XX_corner2; x++)
                        V(z,y,x)+=
                            gaussian1D(z-ZZ(atomPosition[n]),sigma)*
                            gaussian1D(y-YY(atomPosition[n]),sigma)*
                            gaussian1D(x-XX(atomPosition[n]),sigma);
        }
    }
    V.write(fnRoot+".vol");

    // Histogram of the intensities
    MultidimArray<double> intensities(atomWeight);
    Histogram1D hist;
    compute_hist(intensities, hist, 100);
    hist.write(fnRoot+"_intensities.hist");
}

double ProgARTPseudo::ART_single_step(const MultidimArray<double> &Iexp,
                                      double rot, double tilt, double psi, double shiftX, double shiftY,
                                      const std::vector<double> &lambda)
{
    MultidimArray<double> Itheo, Icorr, Idiff;
    Itheo.initZeros(Iexp);
    Icorr.initZeros(Iexp);
    Matrix2D<double> Euler;
    Euler_angles2matrix(rot, tilt, psi, Euler);
    project_Pseudo(atomPosition, atomWeight, sigma, Itheo, Icorr,
                   Euler, shiftX, shiftY, lambda, NMA, FORWARD,
                   gaussianProjectionTable, gaussianProjectionTable2);
    Idiff.initZeros(Iexp);

    double mean_error=0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iexp)
    {
        // Compute difference image and error
        DIRECT_A2D_ELEM(Idiff, i, j) = DIRECT_A2D_ELEM(Iexp, i, j) - DIRECT_A2D_ELEM(Itheo, i, j);
        mean_error += DIRECT_A2D_ELEM(Idiff, i, j) * DIRECT_A2D_ELEM(Idiff, i, j);

        // Compute the correction image
        DIRECT_A2D_ELEM(Icorr, i, j) = XMIPP_MAX(DIRECT_A2D_ELEM(Icorr, i, j), 1);
        DIRECT_A2D_ELEM(Icorr, i, j) =
            lambdaART * DIRECT_A2D_ELEM(Idiff, i, j) / DIRECT_A2D_ELEM(Icorr, i, j);
    }
    mean_error /= YXSIZE(Iexp);

    project_Pseudo(atomPosition, atomWeight, sigma, Itheo, Icorr,
                   Euler, shiftX, shiftY, lambda, NMA, BACKWARD,
                   gaussianProjectionTable, gaussianProjectionTable2);
    return mean_error;
}
