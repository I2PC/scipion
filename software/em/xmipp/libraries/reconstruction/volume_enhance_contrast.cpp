/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2010)
 *             Fernando Fuentes
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

#include "volume_enhance_contrast.h"

//#define DEBUG

// Define paramsUsage -------------------------------------------------------------------s
void ProgVolumeEnhanceContrast::defineParams()
{
	addUsageLine("Enhance volume contrast by applying a non-linear transformation to gray levels");
	addUsageLine("+The contrast enhancement is based on an evolution of an article of");
	addUsageLine("+[[http://www.ncbi.nlm.nih.gov/pubmed/16579377][Sean Matz]] and is further ");
	addUsageLine("+explained [[http://biocomp.cnb.csic.es/~coss/Articulos/Fuentes2010.pdf][here]]");
	addUsageLine("+The algorithm detects the molecule borders and applies a nonlinear transformation ");
	addUsageLine("+of the gray values");
    addSeeAlsoLine("volume_correct_bfactor");
    addParamsLine("   -i <file>                  : Input volume");
    addParamsLine("  [-o <file=\"\">]                 : Output volume");
    addParamsLine("  [--removeBackground]        : Remove the noise of the background");
    addParamsLine("  [--alpha+ <a=0.01>]         : Confidence interval for background identification");
    addParamsLine("  [--lowerIntensity+ <intensity=0>]     : Only process if the gray value is higher than this value");
    addParamsLine("  [--saveMask+ <filename=\"\">]  : Filename for the background mask");
    addExampleLine("xmipp_volume_enhance_contrast -i volume.vol -o volumeEnhanced.vol");
}

// Read from command line --------------------------------------------------
void ProgVolumeEnhanceContrast::readParams()
{
	fnIn=getParam("-i");
	fnOut=getParam("-o");
	if (fnOut=="")
		fnOut=fnIn;
    alpha = getDoubleParam("--alpha");
    lowerIntensity = getDoubleParam("--lowerIntensity");
    removeBg = checkParam("--removeBackground");
    fnMask = getParam("--saveMask");
}

void ProgVolumeEnhanceContrast::run()
{
  Image<double> img;
  img.read(fnIn);
  enhance(img());
  img.write(fnOut);
}

// Show --------------------------------------------------------------------
void ProgVolumeEnhanceContrast::show()
{
	if (verbose==0)
		return;
    std::cout
    << "Input: " << fnIn << std::endl
    << "Output: " << fnOut << std::endl
    << "Alpha: " << alpha << std::endl
    << "lowerIntensity: " << lowerIntensity << std::endl
    << "removeBackground: " << removeBg << std::endl
    << "saveMask: " << fnMask << std::endl;
}

// Enhance volume ----------------------------------------------------------
void ProgVolumeEnhanceContrast::enhance(MultidimArray<double> &vol)
{
    // 1.-Scale volume between 0 and 255----------------------------------
    double minVal=0., maxVal=0.;
    vol.computeDoubleMinMax(minVal,maxVal);
    vol.rangeAdjust(0,255);

    double vol_avg = vol.computeAvg();

    // Padd volume with two new rows and cols
    if( FINISHINGZ(vol) - STARTINGZ(vol) > 0 )
    {
        vol.selfWindow( STARTINGZ(vol)-2,STARTINGY(vol)-2,STARTINGX(vol)-2,
                    FINISHINGZ(vol)+2, FINISHINGY(vol)+2,FINISHINGX(vol)+2, vol_avg );
    }
    else
    {
        vol.selfWindow( STARTINGZ(vol),STARTINGY(vol)-2,STARTINGX(vol)-2,
                    FINISHINGZ(vol), FINISHINGY(vol)+2,FINISHINGX(vol)+2, vol_avg );
    }

    Image<double> save;
#ifdef DEBUG

    save()=vol;
    save.write("PPP1_init.vol");
#endif

    // 2.-Elimination of the background-----------------------------------
    MultidimArray<double> mask;
    double bg_mean;
    detectBackground(vol,mask,0.01,bg_mean);
#ifdef DEBUG

    save()=mask;
    save.write("PPPmask.vol");
#endif

    // We change 0<->1
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        if (A3D_ELEM(mask,k,i,j)==1)
        {
            A3D_ELEM(mask,k,i,j)=0;
        }
        else
        {
            A3D_ELEM(mask,k,i,j)=1;
        }
    } // Now if 0:background and if 1:mol
#ifdef DEBUG
    save()=mask;
    save.write("PPPmask_c.vol");
#endif

    if (fnMask!="")
    {
        save()=mask;
        save.write(fnMask);
    }
    if (removeBg==true)
    {
        // We put all the background with the same value
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol)
        if (DIRECT_A3D_ELEM(mask,k,i,j)==0)
            DIRECT_A3D_ELEM(vol,k,i,j)=bg_mean;
    }
#ifdef DEBUG
    save()=vol;
    save.write("PPP2_no_bg.vol");
#endif

    // 3.-BSplines + diff = edges-------------------------------------------
    MultidimArray<double> vol_edge;
    computeEdges(vol,vol_edge);

#ifdef DEBUG
    save()=vol_edge;
    save.write("PPP3_edge.vol");
#endif
    // 4.-MEAN EDGE GRAY VALUE

    // 4.1.-Variable Neighbourhood
#define COMPUTE_STATISTICS(V,N,k,i,j,avg,stddev,cubeSize) \
     {\
   int k0,kF,i0,iF,j0,jF; \
     k0=XMIPP_MAX(k-N,STARTINGZ(V)); \
     kF=XMIPP_MIN(k+N,FINISHINGZ(V)); \
     i0=XMIPP_MAX(i-N,STARTINGY(V)); \
     iF=XMIPP_MIN(i+N,FINISHINGY(V)); \
     j0=XMIPP_MAX(j-N,STARTINGX(V)); \
         jF=XMIPP_MIN(j+N,FINISHINGX(V)); \
         double sum=0, sum2=0; \
         for (int kk=k0; kk<=kF; ++kk) \
    for (int ii=i0; ii<=iF; ++ii) \
     for (int jj=j0; jj<=jF; ++jj) \
     { \
      double v=A3D_ELEM(V,kk,ii,jj); \
      sum+=v; \
      sum2+=v*v; \
     } \
   cubeSize=(kF-k0+1)*(iF-i0+1)*(jF-j0+1); \
   avg=sum/cubeSize; \
   stddev=sum2/cubeSize-avg*avg; \
   stddev=sqrt(XMIPP_MAX(stddev,0)); \
  }

    // We use gaussian because we have a T-Student with more than 100
    float z=icdf_gauss(1-(0.01/2));  // degrees of freedom
    // We create a volume with the neighbor sizes
    MultidimArray<int> vol_tam;
    vol_tam.resize(vol);
    vol_tam.initConstant(4);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
    {
        if (A3D_ELEM(mask,k,i,j)==1)
        {
            // only compute if the pixel is not background
            int cubeDim=1; // 3x3x3
            int newCubeDim=1; // 3x3x3
            int max_tam=4; // 9x9x9
            double avgSmall, stddevSmall, NSmall;
            COMPUTE_STATISTICS(vol,cubeDim,k,i,j,avgSmall,stddevSmall,NSmall);
            bool found=false;
            while (!found)
            {
                // Big Cube
                newCubeDim=newCubeDim+1;
                double avgBig, stddevBig, NBig;
                COMPUTE_STATISTICS(vol,newCubeDim,k,i,j,avgBig,stddevBig,NBig);
                // Computing the degrees of freedom
                //double num, den;
                //num=(((stddevSmall*stddevSmall)/NSmall)+((stddevBig*stddevBig)/NBig))*
                // (((stddevSmall*stddevSmall)/NSmall)+((stddevBig*stddevBig)/NBig));
                //den=((stddevSmall*stddevSmall)/((NSmall-1)*(NSmall*NSmall)))+
                // ((stddevBig*stddevBig)/((NBig-1)*(NBig*NBig)));
                //double df;
                //df=num/den; // Degrees of freedom
                // Confidence Intervals --> Comparison of two cubes
                double K, d, A, B;
                K=sqrt(((stddevSmall*stddevSmall)/NSmall)+((stddevBig*stddevBig)/NBig));
                d=avgSmall-avgBig;
                A=d-(K*z);
                B=d+(K*z);
                // If the interval [A,B] contains the cero there are the same
                if ((A<0) && (0<B))
                {  // equal continue or not
                    if (newCubeDim>=max_tam)
                    { // not continue
                        found=true;
                        A3D_ELEM(vol_tam,k,i,j)=newCubeDim;
                    }
                    else
                    { // continue searching
                        avgSmall=avgBig;
                        stddevSmall=stddevBig;
                        NSmall=NBig;
                    }
                }
                else // Not equal -> we stop searching
                {
                    found=true;  // Size found
                    A3D_ELEM(vol_tam,k,i,j)=newCubeDim-1; // it's the latest
                }

            } // end of while
        } // end of if
    } // end of FOR_ALL_ELEMENTS
#ifdef DEBUG
    save()=vol_tam;
    save.write("PPP4_vol_tam.vol");
#endif

    // 4.2.- Compute the mean edge gray value
    MultidimArray<double> VolxEdge;
    VolxEdge=vol;
    VolxEdge*=vol_edge;
    // Macro for compute de Mean Edge Gray Value
#define COMPUTE_MEGV(VxE,V_E,N,k,i,j,pixel) \
     {\
   int k0,kF,i0,iF,j0,jF; \
     k0=XMIPP_MAX(k-N,STARTINGZ(V_E)); \
     kF=XMIPP_MIN(k+N,FINISHINGZ(V_E)); \
     i0=XMIPP_MAX(i-N,STARTINGY(V_E)); \
     iF=XMIPP_MIN(i+N,FINISHINGY(V_E)); \
     j0=XMIPP_MAX(j-N,STARTINGX(V_E)); \
         jF=XMIPP_MIN(j+N,FINISHINGX(V_E)); \
         double sum1=0, sum2=0; \
         for (int kk=k0; kk<=kF; kk++) \
    for (int ii=i0; ii<=iF; ii++) \
     for (int jj=j0; jj<=jF; jj++) \
     { \
      double v1=A3D_ELEM(VxE,kk,ii,jj); \
      sum1=sum1+v1; \
      double v2=A3D_ELEM(V_E,kk,ii,jj); \
      sum2=sum2+v2; \
     } \
   pixel=sum1/sum2; \
  }

    MultidimArray<double> Vol_E;
    Vol_E=vol;
    int tam;
    double pixel;
    // We compute the MEGV for all the volume
    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
    { // Only compute if no background
        if (A3D_ELEM(mask,k,i,j)==1)
        {
            tam=A3D_ELEM(vol_tam,k,i,j);
            COMPUTE_MEGV(VolxEdge,vol_edge,tam,k,i,j,pixel);
            A3D_ELEM(Vol_E,k,i,j)=pixel;
        }
    }
#ifdef DEBUG
    save()=Vol_E;
    save.write("PPP5_vol_megv.vol");
#endif

    //5.-Nonlinear function

    MultidimArray<double> vol_f=vol;
    double a, b, x, E;
    double threshold=lowerIntensity*255;

    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
    {
        x=A3D_ELEM(vol,k,i,j);
        if (A3D_ELEM(mask,k,i,j)==1 && x>threshold)
        {
            E=A3D_ELEM(Vol_E,k,i,j);
            // We search the interval
            if ((x>=0)  && (x<=3))
            {
                a=0;
                b=3;
            }
            else if ((x>3)   && (x<=8))
            {
                a=3;
                b=8;
            }
            else if ((x>8)   && (x<=16))
            {
                a=8;
                b=16;
            }
            else if ((x>16)  && (x<=30))
            {
                a=16;
                b=30;
            }
            else if ((x>30)  && (x<=49))
            {
                a=30;
                b=49;
            }
            else if ((x>49)  && (x<=75))
            {
                a=49;
                b=75;
            }
            else if ((x>75)  && (x<=107))
            {
                a=75;
                b=107;
            }
            else if ((x>107) && (x<=147))
            {
                a=107;
                b=147;
            }
            else if ((x>147) && (x<=196))
            {
                a=147;
                b=196;
            }
            else
            {
                a=196;
                b=255;
            }
            // nonlinear function
            double x_new;
            if (x<=E)
            {
                x_new=E-sqrt(((E-a)*(E-a))-((x-a)*(x-a)));
            }
            if (x>E)
            {
                x_new=E+sqrt(((b-E)*(b-E))-((x-b)*(x-b)));
            }
            A3D_ELEM(vol_f,k,i,j)=x_new;
        }
    }
#ifdef DEBUG
    save()=vol_f;
    save.write("PPP6_vol_f.vol");
#endif

    // Unpadd volume with two new rows and cols
    if( FINISHINGZ(vol_f) - STARTINGZ(vol_f) > 0 )
    {
        vol_f.selfWindow( STARTINGZ(vol_f)+2,STARTINGY(vol_f)+2,STARTINGX(vol_f)+2,
                      FINISHINGZ(vol_f)-2, FINISHINGY(vol_f)-2,FINISHINGX(vol_f)-2 );
    }
    else
    {
        vol_f.selfWindow( STARTINGZ(vol_f),STARTINGY(vol_f)+2,STARTINGX(vol_f)+2,
                      FINISHINGZ(vol_f), FINISHINGY(vol_f)-2,FINISHINGX(vol_f)-2 );
    }

    //6.-Re-Scale (now we put the original values different from [0-255])
    double x_s,x_ns;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_f)
    {
        x_s=A3D_ELEM(vol_f,k,i,j);
        x_ns=((maxVal-minVal)/255)*(x_s+((255*minVal)/(maxVal-minVal)));
        A3D_ELEM(vol_f,k,i,j)=x_ns;
    }
    // OUTPUT
    vol=vol_f;
} // End of enhance

