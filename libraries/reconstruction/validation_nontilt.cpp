/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#include "validation_nontilt.h"
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time
#include <algorithm>


void ProgValidationNonTilt::readParams()
{

    fnIn = getParam("-i");
    fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--volume");
    alpha0 = getDoubleParam("--alpha0");
    //alphaF = getDoubleParam("--alphaF");
    //Niter = getIntParam("--iter");
    //keepIntermediateVolumes = checkParam("--keepIntermediateVolumes");
    angularSampling=getDoubleParam("--angularSampling");
    sampling_rate = getDoubleParam("--sampling_rate");
    //maxShift=getDoubleParam("--maxShift");
    //tilt0=getDoubleParam("--minTilt");
    //tiltF=getDoubleParam("--maxTilt");
    //useImed=checkParam("--useImed");
    //strict=checkParam("--strictDirection");
    //angDistance=getDoubleParam("--angDistance");
    //Nvolumes=getIntParam("--numberOfVolumes");
    Nvolumes = 1;

}

void ProgValidationNonTilt::defineParams()
{
    //usage
    addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    //addParamsLine("  [--numberOfVolumes <N=1>]    : Number of volumes to reconstruct");
    addParamsLine("  [--volume <md_file=\"\">]    : Volume to validate");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    //addParamsLine("  [--iter <N=10>]              : Number of iterations");
    addParamsLine("  [--alpha0 <N=0.05>]          : Significance");
    //addParamsLine("  [--alphaF <N=0.005>]         : Final significance");
    //addParamsLine("  [--keepIntermediateVolumes]  : Keep the volume of each iteration");
    addParamsLine("  [--angularSampling <a=5>]    : Angular sampling in degrees for generating the projection gallery");
    addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate in A/px");
    //addParamsLine("  [--maxShift <s=-1>]          : Maximum shift allowed (+-this amount)");
    //addParamsLine("  [--minTilt <t=0>]            : Minimum tilt angle");
    //addParamsLine("  [--maxTilt <t=90>]           : Maximum tilt angle");
    //addParamsLine("  [--useImed]                  : Use Imed for weighting");
    //addParamsLine("  [--strictDirection]          : Images not significant for a direction are also discarded");
    //addParamsLine("  [--angDistance <a=10>]       : Angular distance");

}

void ProgValidationNonTilt::run()
{
    //Clustering Tendency and Cluster Validity
    //Stephen D. Scott

    randomize_random_generator();
    char buffer[400];
    sprintf(buffer, "xmipp_reconstruct_significant -i %s  --initvolumes %s --odir %s --sym  %s --iter 1 --alpha0 %f --angularSampling %f",fnIn.c_str(), fnInit.c_str(),fnDir.c_str(),fnSym.c_str(),alpha0,angularSampling);
    system(buffer);

    MetaData md,mdOut,mdOut2,tempMd2,mdWeight;
    FileName fnMd,fnOut,fnFSC,fnOut2;
    fnMd = fnDir+"/angles_iter01_00.xmd";
    fnOut = fnDir+"/kk.xmd";
    fnOut2 = fnDir+"/kk2.xmd";
    fnFSC = fnDir+"/fsc.xmd";
    size_t nSamplesRandom = 20;

    md.read(fnMd);
    size_t maxNImg;
    size_t sz = md.size();
    md.getValue(MDL_IMAGE_IDX,maxNImg,sz);

    String expression;

    //double H,H0,d;
    double W,sumW;
    double tempW;
    MDRow row;

    init_progress_bar(maxNImg);
    for (size_t i=0; i<=maxNImg;i++)
    {
    	MetaData tempMd;
        expression = formatString("imageIndex == %lu",1);
        tempMd.importObjects(md, MDExpression(expression));

        if (tempMd.size()==0)
            continue;

        std::vector<double> sum_u(nSamplesRandom);
        std::vector<double> sum_w(nSamplesRandom);
        std::vector<double> H0(nSamplesRandom);
        std::vector<double> H(nSamplesRandom);

        obtainSumU(tempMd,sum_u,H0);
        obtainSumW(tempMd,sum_w,sum_u,H);

        std::sort(H0.begin(),H0.end());
        std::sort(H.begin(),H.end());

        double P = 0;
        for(int i=0; i<sum_u.size();i++)
        {
            if (H0.at(i)>H.at(i) )
                P += 1;
        }

        P = (P/nSamplesRandom);
        std::cout << "P : " << P << std::endl;

        row.setValue(MDL_IMAGE_IDX,i);
        row.setValue(MDL_ZSCORE,P);
        mdOut.addRow(row);

        sum_u.clear();
		sum_w.clear();
		H0.clear();
		H.clear();
		tempMd.clear();
        progress_bar(i+1);

    }

    mdOut.write(fnOut);

    FileName fnVolume=fnDir+"/volume_projMatch.vol";
    String args=formatString("-i %s -o %s --sym %s --weight -v 0",fnMd.c_str(),fnVolume.c_str(),fnSym.c_str());
    String cmd=(String)"xmipp_reconstruct_fourier "+args;

    if (system(cmd.c_str())==-1)
        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");


    // Size of the images
    size_t Xdim, Ydim, Zdim,Ndim;
    getImageSize(fnVolume,Xdim,Ydim,Zdim,Ndim);
    args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-Xdim/2);
    cmd=(String)"xmipp_transform_mask "+args;
    if (system(cmd.c_str())==-1)
        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

    FileName fnVolumeSig=fnDir+"/volume_iter01_00.vol";
    sprintf(buffer, "xmipp_resolution_fsc --ref %s -i %s -s %f -o %s",fnVolume.c_str(),fnVolumeSig.c_str(),sampling_rate,fnFSC.c_str());
    system(buffer);

    /*



    size_t sz = md.size();
    md.getValue(MDL_IMAGE_IDX,maxNImg,sz);
    std::cout << maxNImg << "\n";

    std::vector<double> clustering(maxNImg+1);
    std::vector<double> clusteringRandom(maxNImg+1);

    size_t nSamplesRandom = 30;
    double sumWRan;

    init_progress_bar(maxNImg);

    for (size_t i=0; i<=maxNImg;i++)
{


        //For each projection we have different number of orientations with weights within
        //the significant value. We have for each projection obtain the null hypotesis with
        // points uniformily distributed in the sphere
        //null Hypothesis



        // sumWRanArray[n] = sumWRan;
        //sumWRan = 1e3;
        // for (int var = 0; var < nSamplesRandom; var++)
        //if (sumWRanArray[var] < sumWRan)
        //  sumWRan += sumWRanArray[var];
        //sumWRan += sumWRanArray[var];

        //clusteringRandom.at(i)=sumWRan/nSamplesRandom;
        clusteringRandom.at(i)=sumWRan/nSamplesRandom;
        sumWRan = 0;

        //END null Hypothesis
        double a;
        sumW=0;
        FOR_ALL_OBJECTS_IN_METADATA(tempMd)
        {
            tempMd2.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
            tempMd2.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
            tempMd2.getValue(MDL_WEIGHT,w,__iter.objId);
            x = sin(tilt*PI/180)*cos(rot*PI/180);
            y = sin(tilt*PI/180)*sin(rot*PI/180);
            z = std::abs(cos(tilt*PI/180));

            tempW = 20;
            FOR_ALL_OBJECTS_IN_METADATA(tempMd)
            {
                tempMd2.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
                tempMd2.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
                tempMd2.getValue(MDL_WEIGHT,w2,__iter.objId);
                xx = sin(tilt*PI/180)*cos(rot*PI/180);
                yy = sin(tilt*PI/180)*sin(rot*PI/180);
                zz = std::abs(cos(tilt*PI/180));
                a = std::abs(std::acos(x*xx+y*yy+z*zz));

                if ( (a<tempW) && (a != 0))
                {
                    //W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));
                    W = a;
                    tempW = a;
                }
            }
            sumW += W;
        }

        clustering.at(i)=sumW;
        progress_bar(i+1);
}

    std::vector<double> H(maxNImg+1);
    std::vector<double> H0(maxNImg+1);

    for (size_t i=0; i<=maxNImg;i++)
{
        H[i] = clustering.at(i)/(clusteringRandom.at(i)+clustering.at(i));
        row.setValue(MDL_IMAGE_IDX,i);
        row.setValue(MDL_ZSCORE,clustering.at(i));
        row.setValue(MDL_ZSCORE_SHAPE1,clusteringRandom.at(i));
        row.setValue(MDL_ZSCORE_SHAPE2,H[i]);
        mdOut.addRow(row);
}

    for (size_t i=0; i<=maxNImg;i++)
{

        std::srand(std::time(0));
        std::random_shuffle(clusteringRandom.begin(), clusteringRandom.end());

        H0[i] = clusteringRandom.at(0)/(clusteringRandom.at(0)+clusteringRandom.at(i));
        row.setValue(MDL_IMAGE_IDX,i);
        row.setValue(MDL_ZSCORE,clusteringRandom.at(0));
        row.setValue(MDL_ZSCORE_SHAPE1,clusteringRandom.at(i));
        row.setValue(MDL_ZSCORE_SHAPE2,H0[i]);
        mdOut2.addRow(row);
}

    std::sort(H.begin(), H.end());
    std::sort(H0.begin(), H0.end());

    double p5Rand =  H0.at(std::floor(maxNImg*0.05));
    double pValue;
    for (size_t i=0; i<=maxNImg;i++)
        if (H[i] <= p5Rand)
        {
            pValue = double(i)/maxNImg;
        }

    mdOut.write(fnOut);
    mdOut2.write(fnOut2);

    std :: cout << " p5Rand : " << p5Rand << std::endl;
    std :: cout << " pValue : " << pValue << std::endl;


       size_t p50,p15,p85;
       sz = mdOut.size();
       p50 = std::floor(sz/2);
       p15 = std::floor(sz/2);
       p85 = std::floor(sz/2);

       //Reconstruct the "projection matching" volume and after that compare with significant volume
       FileName fnAngles=fnDir+"/images_significant_iter01_00.xmd";
       if (!fnAngles.exists())
           REPORT_ERROR(ERR_UNCLASSIFIED,"Angles file does not exist");





    */
}

void ProgValidationNonTilt::obtainSumU(MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0)
{
    //MetaData md,tempMd,mdWeight;
    //FileName fnMd;
    //randomize_random_generator();
    double xRan,yRan,zRan,norm;
    double sumWRan;
    double xRanArray[tempMd.size()];
    double yRanArray[tempMd.size()];
    double zRanArray[tempMd.size()];
    std::vector<double> weightV;
    double a;

    for (size_t n=0; n<sum_u.size(); n++)
    {
        sumWRan = 0;
        for (size_t nS=0; nS<tempMd.size(); nS++)
        {
            xRan=(std::rand()-RAND_MAX/2);
            yRan=(std::rand()-RAND_MAX/2);
            zRan=(std::rand()-RAND_MAX/2);
            norm = std::sqrt(xRan*xRan+yRan*yRan+zRan*zRan);
            xRan = (xRan/norm);
            yRan = (yRan/norm);
            zRan = std::abs(zRan/norm);

            xRanArray[nS] = xRan;
            yRanArray[nS] = yRan;
            zRanArray[nS] = zRan;
            tempMd.getColumnValues(MDL_WEIGHT, weightV);

            std::random_shuffle(weightV.begin(), weightV.end());
        }

        sumWRan = 0;
        double WRan, tempWRan, tempW1, tempW2;
        for (size_t nS1=0; nS1<tempMd.size(); nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<tempMd.size(); nS2++)
            {
                a = std::abs(std::acos(xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2]));
                if ( (a<tempWRan) && (a != 0))
                {
                    tempWRan = a;
                    tempW2 = weightV[nS2];
                    tempW1 = weightV[nS1];
                    WRan = a*std::exp(std::abs(tempW1-tempW2))*std::exp(-(tempW1+tempW2));
                    //WRan = a;
                }
            }
            sumWRan += WRan;
        }
        sum_u.at(n)=sumWRan;
    }

    size_t idx = 0;
    while (idx < sum_u.size())
    {
        std::random_shuffle(sum_u.begin(), sum_u.end());

        if(sum_u.at(0) != sum_u.at(1))
        {
            H0[idx] = sum_u.at(0)/(sum_u.at(0)+sum_u.at(1));
            idx += 1;
        }
    }

}

void ProgValidationNonTilt::obtainSumW(MetaData & tempMd,std::vector<double> & sum_W,std::vector<double> & sum_u,std::vector<double> & H)
{
    double a;
    double rot,tilt,w;
    double x,y,z;
    double xx,yy,zz;
    double w2;
    double tempW;
    double W;
    double sumW;

    sumW = 0;
    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        tempMd.getValue(MDL_WEIGHT,w,__iter.objId);
        x = sin(tilt*PI/180)*cos(rot*PI/180);
        y = sin(tilt*PI/180)*sin(rot*PI/180);
        z = std::abs(cos(tilt*PI/180));

        tempW = 1e3;
        FOR_ALL_OBJECTS_IN_METADATA(tempMd)
        {
            tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
            tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
            tempMd.getValue(MDL_WEIGHT,w2,__iter.objId);
            xx = sin(tilt*PI/180)*cos(rot*PI/180);
            yy = sin(tilt*PI/180)*sin(rot*PI/180);
            zz = std::abs(cos(tilt*PI/180));
            a = std::abs(std::acos(x*xx+y*yy+z*zz));

            if ( (a<tempW) && (a != 0))
            {
                W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));
                //W = a;
                tempW = a;
            }
        }
        sumW +=  W;
    }

    size_t idx = 0;
    for (size_t n=0; n<sum_u.size(); n++)
    {
        H[n] = sumW/(sumW+sum_u.at(n));
    }
}

