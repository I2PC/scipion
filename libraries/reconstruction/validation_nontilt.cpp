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
	char buffer[400];
	//sprintf(buffer, "xmipp_reconstruct_significant -i %s  --initvolumes %s --odir %s --sym  %s --iter 1 --alpha0 %f --angularSampling %f",fnIn.c_str(), fnInit.c_str(),fnDir.c_str(),fnSym.c_str(),alpha0,angularSampling);
	//system(buffer);

	MetaData md,tempMd,mdOut,tempMd2;
	FileName fnMd,fnOut,fnFSC;
	fnMd = fnDir+"/angles_iter01_00.xmd";
	fnOut = fnDir+"/kk.xmd";
	fnFSC = fnDir+"/fsc.xmd";

	md.read(fnMd);
	size_t nImg,maxNImg;
	String expression;
	double rot,tilt,w;
	double x,y,z;
	double R,R0;
	double x_,y_,z_;
	double x2,y2,z2,w2;
	double k,d;
	bool test;
	MDRow row;

	size_t sz = md.size();
	md.getValue(MDL_IMAGE_IDX,maxNImg,sz);
	std::cout << maxNImg << "\n";
	for (size_t i=0; i<=maxNImg;i++)
	{
		expression = formatString("imageIndex == %lu",i);
		tempMd.importObjects(md, MDExpression(expression));
		if (tempMd.size()==0)
			continue;

		R = 0;
		x2=0;
		y2=0;
		z2=0;
		w2=0;

		tempMd2.sort(tempMd,MDL_WEIGHT,false, -1, 0);
	    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
	    {
	    	tempMd2.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
	    	tempMd2.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
	    	tempMd2.getValue(MDL_WEIGHT,w,__iter.objId);
	    	x = sin(tilt*PI/180)*cos(rot*PI/180);
	    	y = sin(tilt*PI/180)*sin(rot*PI/180);
	    	z = (cos(tilt*PI/180));

	    	if (__iter.objId == 1)
	    	{
	    		// __iter.objId =1 is the first element with the largest weight
	    		x_ = x;
	    		y_ = y;
	    		z_ = z;
	    	}

	    	if (z < 0)
	    	{
	    		x = -x;
	    		y = -y;
	    		z = -z;
	    	}

	    	d = std::sqrt( (x-x_)*(x-x_)+(y-y_)*(y-y_)+(z-z_)*(z-z_));

	    	if ( d>1.85 ) //we have mirror the projection is at 180º
	    	{
	    		x = -x;
	    		y = -y;
	    	}

	    	x2 += x*w;
	    	y2 += y*w;
	    	z2 += z*w;
	    	w2 += w;
	    }

	    R = std::sqrt(x2*x2+y2*y2+z2*z2);
	    R0 = std::sqrt(3.782*tempMd2.size());

	    x_=x2/R;
	    y_=y2/R;
	    z_=z2/R;
	    k=(w2-1)/(w2-R);

	    row.setValue(MDL_IMAGE_IDX,i);
	    row.setValue(MDL_ZSCORE,k);
	    row.setValue(MDL_KSTEST,R);
	    row.setValue(MDL_FLIP,R>R0);

	    mdOut.addRow(row);

	}

	mdOut.write(fnOut);
	size_t p50,p15,p85;
	sz = mdOut.size();
	p50 = std::floor(sz/2);
	p15 = std::floor(sz/2);
	p85 = std::floor(sz/2);

	//Reconstruct the "projection matching" volume and after that compare with significant volume
	FileName fnAngles=fnDir+"/images_significant_iter01_00.xmd";
	if (!fnAngles.exists())
		REPORT_ERROR(ERR_UNCLASSIFIED,"Angles file does not exist");

	FileName fnVolume=fnDir+"/volume_projMatch.vol";
	String args=formatString("-i %s -o %s --sym %s --weight -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str());
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

}


