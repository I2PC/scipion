/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

#include "angular_distance.h"

#include <data/args.h>
#include <data/histogram.h>

// Read arguments ==========================================================
void ProgAngularDistance::readParams()
{
    fn_ang1 = getParam("--ang1");
    fn_ang2 = getParam("--ang2");
    fn_out = getParam("--oroot");
    fn_sym = getParam("--sym");
    check_mirrors = checkParam("--check_mirrors");
    object_rotation = checkParam("--object_rotation");
    compute_weights = checkParam("--compute_weights");
    if (compute_weights)
    {
    	minSigma = getDoubleParam("--compute_weights");
    	idLabel = getParam("--compute_weights",1);
    	minSigmaD = getDoubleParam("--compute_weights",2);
    }
    set = getIntParam("--set");
}

// Show ====================================================================
void ProgAngularDistance::show()
{
    std::cout
    << "Angular docfile 1: " << fn_ang1       << std::endl
    << "Angular docfile 2: " << fn_ang2       << std::endl
    << "Angular output   : " << fn_out    << std::endl
    << "Symmetry file    : " << fn_sym        << std::endl
    << "Check mirrors    : " << check_mirrors << std::endl
    << "Object rotation  : " << object_rotation<<std::endl
    << "Compute weights  : " << compute_weights << std::endl
    << "Min sigma        : " << minSigma << std::endl
    << "Min sigmaD       : " << minSigmaD << std::endl
    << "IdLabel          : " << idLabel << std::endl
    << "Set              : " << set << std::endl
    ;
}

// usage ===================================================================
void ProgAngularDistance::defineParams()
{
    addUsageLine("Computes the angular distance between two angle files. The angular distance ");
    addUsageLine("is defined as the average angular distance between the 3 vectors of the ");
    addUsageLine("coordinate system defined by the Euler angles (taking into account any ");
    addUsageLine("possible symmetry). ");
    addParamsLine("   --ang1 <Metadata1>        : Angular document file 1");
    addParamsLine("   --ang2 <Metadata2>        : Angular document file 2");
    addParamsLine("  [--oroot <rootname=\"\">]  : Output rootname");
    addParamsLine("                             : rootname.xmd Angular comparison;");
    addParamsLine("                             : rootname_vec_diff_hist.txt Histogram of the differences in vector directions;");
    addParamsLine("                             : rootname_shift_diff_hist.txt Histogram of the differences in shifts;");
    addParamsLine("                             :+ rootname_rot_diff_hist.txt (verbose>=2) Histogram of the differences in rot;");
    addParamsLine("                             :+ rootname_tilt_diff_hist.txt (verbose>=2) Histogram of the differences in tilt;");
    addParamsLine("                             :+ rootname_psi_diff_hist.txt (verbose>=2) Histogram of the differences in psi;");
    addParamsLine("                             :+ rootname_X_diff_hist.txt (verbose>=2) Histogram of the differences in shiftX;");
    addParamsLine("                             :+ rootname_Y_diff_hist.txt (verbose>=2) Histogram of the differences in shiftY;");
    addParamsLine("  [--sym <symmetry=\"\">]    : Symmetry file if any");
    addParamsLine("                             :+The definition of the symmetry is described at [[transform_symmetrize_v3][transform_symmetrize]]");
    addParamsLine("  [--check_mirrors]          : Check if mirrored projections give better results");
    addParamsLine("  [--object_rotation]        : Use object rotations rather than projection directions");
    addParamsLine("                             : fit (Spider, APMQ)");
    addParamsLine("  [--compute_weights <minSigma=1> <idLabel=particleId> <minSigmaD=-1>] : If this flag is given, images in ang2 are given a weight according to their ");
    addParamsLine("                             : distance to the same image in ang1. The ang2 file is rewritten.");
    addParamsLine("                             : Ang1 and ang2 are supposed to have a label called itemId");
    addParamsLine("                             : The output is written to oroot+_weights.xmd");
    addParamsLine("                             : Min sigma is the minimum angular standard deviation, by default, 1 degree");
    addParamsLine("                             : Min sigmaD is the minimum shift standard deviation, set to -1 for not using shifts for weighting");
    addParamsLine("  [--set <set=1>]            : Set of distance to compute (angular_diff0 and jumper_weight0, angular_diff and jumper_weight,");
    addParamsLine("                             : or angular_diff2 and jumper_weight2)");
}

// Produce side information ================================================
void ProgAngularDistance::produce_side_info()
{
    if (fn_sym != "")
        SL.readSymmetryFile(fn_sym);

    // Check that both docfiles are of the same length
    if (fn_ang1!="" && fn_ang2!="")
    {
        DF1.read(fn_ang1);
        DF2.read(fn_ang2);
        if (DF1.size() != DF2.size() && !compute_weights)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "Angular_distance: Input Docfiles with different number of entries");
    }
}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgAngularDistance::run()
{
    produce_side_info();
    if (compute_weights)
    {
    	computeWeights();
    	return;
    }

    MetaData DF_out;
    double angular_distance=0;
    double shift_distance=0;

    MultidimArray<double> rot_diff, tilt_diff, psi_diff, vec_diff,
    X_diff, Y_diff, shift_diff;
    rot_diff.resize(DF1.size());
    tilt_diff.resize(rot_diff);
    psi_diff.resize(rot_diff);
    vec_diff.resize(rot_diff);
    X_diff.resize(rot_diff);
    Y_diff.resize(rot_diff);
    shift_diff.resize(rot_diff);

    // Build output comment
    /////DF_out.setComment("image rot1 rot2 diff_rot tilt1 tilt2 diff_tilt psi1 psi2 diff_psi ang_dist X1 X2 Xdiff Y1 Y2 Ydiff ShiftDiff");

    int i = 0;
    size_t id;
    FileName fnImg;
    std::vector<double> output;
    output.resize(17,0);
    bool fillOutput=fn_out!="";
    MDRow row;
    FOR_ALL_OBJECTS_IN_METADATA2(DF1, DF2)
    {
        // Read input data
        double rot1,  tilt1,  psi1;
        double rot2,  tilt2,  psi2;
        double rot2p, tilt2p, psi2p;
        double distp;
        double X1, X2, Y1, Y2;
        DF1.getValue(MDL_IMAGE,fnImg,__iter.objId);

        DF1.getValue(MDL_ANGLE_ROT,rot1,__iter.objId);
        DF1.getValue(MDL_ANGLE_TILT,tilt1,__iter.objId);
        DF1.getValue(MDL_ANGLE_PSI,psi1,__iter.objId);
        DF1.getValue(MDL_SHIFT_X,X1,__iter.objId);
        DF1.getValue(MDL_SHIFT_Y,Y1,__iter.objId);

        DF2.getValue(MDL_ANGLE_ROT,rot2,__iter2.objId);
        DF2.getValue(MDL_ANGLE_TILT,tilt2,__iter2.objId);
        DF2.getValue(MDL_ANGLE_PSI,psi2,__iter2.objId);
        DF2.getValue(MDL_SHIFT_X,X2,__iter2.objId);
        DF2.getValue(MDL_SHIFT_Y,Y2,__iter2.objId);

        // Bring both angles to a normalized set
        rot1 = realWRAP(rot1, -180, 180);
        tilt1 = realWRAP(tilt1, -180, 180);
        psi1 = realWRAP(psi1, -180, 180);

        rot2 = realWRAP(rot2, -180, 180);
        tilt2 = realWRAP(tilt2, -180, 180);
        psi2 = realWRAP(psi2, -180, 180);

        // Apply rotations to find the minimum distance angles
        rot2p = rot2;
        tilt2p = tilt2;
        psi2p = psi2;
        distp = SL.computeDistance(rot1, tilt1, psi1,
        		                   rot2p, tilt2p, psi2p, false,
                                   check_mirrors, object_rotation);
        angular_distance += distp;

        // Compute angular difference
        rot_diff(i) = rot1 - rot2p;
        tilt_diff(i) = tilt1 - tilt2p;
        psi_diff(i) = psi1 - psi2p;
        vec_diff(i) = distp;
        X_diff(i) = X1 - X2;
        Y_diff(i) = Y1 - Y2;
        shift_diff(i) = sqrt(X_diff(i)*X_diff(i)+Y_diff(i)*Y_diff(i));
        shift_distance += shift_diff(i);

        // Fill the output result
        if (fillOutput)
        {
            //output[0]=rot1;
            row.setValue(MDL_ANGLE_ROT, rot1);
            //output[1]=rot2p;
            row.setValue(MDL_ANGLE_ROT2, rot2p);
            //output[2]=rot_diff(i);
            row.setValue(MDL_ANGLE_ROT_DIFF, rot_diff(i));
            //output[3]=tilt1;
            row.setValue(MDL_ANGLE_TILT, tilt1);
            //output[4]=tilt2p;
            row.setValue(MDL_ANGLE_TILT2, tilt2p);
            //output[5]=tilt_diff(i);
            row.setValue(MDL_ANGLE_TILT_DIFF,tilt_diff(i) );
            //output[6]=psi1;
            row.setValue(MDL_ANGLE_PSI, psi1);
            //output[7]=psi2p;
            row.setValue(MDL_ANGLE_PSI2, psi2);
            //output[8]=psi_diff(i);
            row.setValue(MDL_ANGLE_PSI_DIFF, psi_diff(i));
            //output[9]=distp;
            if (set==1)
            	row.setValue(MDL_ANGLE_DIFF, distp);
            else
            	row.setValue(MDL_ANGLE_DIFF2, distp);
            //output[10]=X1;
            row.setValue(MDL_SHIFT_X,X1);
            //output[11]=X2;
            row.setValue(MDL_SHIFT_X2, X2);
            //output[12]=X_diff(i);
            row.setValue(MDL_SHIFT_X_DIFF, X_diff(i));
            //output[13]=Y1;
            row.setValue(MDL_SHIFT_Y, Y1);
            //output[14]=Y2;
            row.setValue(MDL_SHIFT_Y2, Y2);
            //output[15]=Y_diff(i);
            row.setValue(MDL_SHIFT_Y_DIFF, Y_diff(i));
            //output[16]=shift_diff(i);
            row.setValue(MDL_SHIFT_DIFF,shift_diff(i));

            id = DF_out.addRow(row);
            //id = DF_out.addObject();
            DF_out.setValue(MDL_IMAGE,fnImg,id);
            //DF_out.setValue(MDL_ANGLE_COMPARISON,output, id);
        }

        i++;
    }
    angular_distance /= i;
    shift_distance /=i;

    if (fillOutput)
    {
        DF_out.write(fn_out + ".xmd");
        Histogram1D hist;
        compute_hist(vec_diff, hist, 0, 180, 180);
        hist.write(fn_out + "_vec_diff_hist.txt",MDL_ANGLE_DIFF,MDL_COUNT);
        compute_hist(shift_diff, hist, 20);
        hist.write(fn_out + "_shift_diff_hist.txt",MDL_SHIFT_DIFF,MDL_COUNT);
        if (verbose==2)
        {
            compute_hist(rot_diff, hist, 100);
            hist.write(fn_out + "_rot_diff_hist.txt");
            compute_hist(tilt_diff, hist, 100);
            hist.write(fn_out + "_tilt_diff_hist.txt");
            compute_hist(psi_diff, hist, 100);
            hist.write(fn_out + "_psi_diff_hist.txt");
            compute_hist(X_diff, hist, 20);
            hist.write(fn_out + "_X_diff_hist.txt");
            compute_hist(Y_diff, hist, 20);
            hist.write(fn_out + "_Y_diff_hist.txt");
        }
    }

    std::cout << "Global angular distance = " << angular_distance << std::endl;
    std::cout << "Global shift   distance = " << shift_distance   << std::endl;
}

void ProgAngularDistance::computeWeights()
{
	MetaData DF1sorted, DF2sorted, DFweights;
	MDLabel label=MDL::str2Label(idLabel);
	DF1sorted.sort(DF1,label);
	DF2sorted.sort(DF2,label);
	std::vector<MDLabel> labels;
	labels.push_back(label);
	labels.push_back(MDL_ANGLE_ROT);
	labels.push_back(MDL_ANGLE_TILT);
	labels.push_back(MDL_ANGLE_PSI);
	labels.push_back(MDL_SHIFT_X);
	labels.push_back(MDL_SHIFT_Y);
	labels.push_back(MDL_FLIP);
    DF1sorted.keepLabels(labels);
    DF2sorted.keepLabels(labels);
    MDLabel angleDiffLabel, shiftDiffLabel, weightLabel;
	switch (set)
	{
	case 1:
		angleDiffLabel=MDL_ANGLE_DIFF;
		shiftDiffLabel=MDL_SHIFT_DIFF;
		weightLabel=MDL_WEIGHT_JUMPER;
		break;
	case 2:
		angleDiffLabel=MDL_ANGLE_DIFF2;
		shiftDiffLabel=MDL_SHIFT_DIFF2;
		weightLabel=MDL_WEIGHT_JUMPER2;
		break;
	case 0:
		angleDiffLabel=MDL_ANGLE_DIFF0;
		shiftDiffLabel=MDL_SHIFT_DIFF0;
		weightLabel=MDL_WEIGHT_JUMPER0;
		break;
	}

    // for(MDIterator __iter(__md); __iter.hasNext(); __iter.moveNext())
    MDIterator iter1(DF1sorted), iter2(DF2sorted);
    std::vector< Matrix1D<double> > ang1, ang2;
    Matrix1D<double> rotTiltPsi(5);
	size_t currentId;
	bool anotherImageIn2=iter2.hasNext();
	size_t id1, id2;
	bool mirror;
    while (anotherImageIn2)
    {
    	ang1.clear();
    	ang2.clear();

    	// Take current id
    	DF2sorted.getValue(label,currentId,iter2.objId);

    	// Grab all the angles in DF2 associated to this id
    	bool anotherIteration=false;
    	do
    	{
    		DF2sorted.getValue(label,id2,iter2.objId);
			anotherIteration=false;
    		if (id2==currentId)
    		{
				DF2sorted.getValue(MDL_ANGLE_ROT,XX(rotTiltPsi),iter2.objId);
				DF2sorted.getValue(MDL_ANGLE_TILT,YY(rotTiltPsi),iter2.objId);
				DF2sorted.getValue(MDL_ANGLE_PSI,ZZ(rotTiltPsi),iter2.objId);
				DF2sorted.getValue(MDL_SHIFT_X,VEC_ELEM(rotTiltPsi,3),iter2.objId);
				DF2sorted.getValue(MDL_SHIFT_Y,VEC_ELEM(rotTiltPsi,4),iter2.objId);
				DF2sorted.getValue(MDL_FLIP,mirror,iter2.objId);
				if (mirror)
				{
					double rotp, tiltp, psip;
					Euler_mirrorY(XX(rotTiltPsi),YY(rotTiltPsi),ZZ(rotTiltPsi), rotp, tiltp, psip);
					XX(rotTiltPsi)=rotp;
					YY(rotTiltPsi)=tiltp;
					ZZ(rotTiltPsi)=psip;
				}
				ang2.push_back(rotTiltPsi);
				iter2.moveNext();
				if (iter2.hasNext())
					anotherIteration=true;
    		}
    	} while (anotherIteration);

    	// Advance Iter 1 to catch Iter 2
    	double N=0, cumulatedDistance=0, cumulatedDistanceShift=0;
    	size_t newObjId=0;
    	if (iter1.objId>0)
    	{
			DF1sorted.getValue(label,id1,iter1.objId);
			while (id1<currentId && iter1.hasNext())
			{
				iter1.moveNext();
				if (iter1.hasNext())
					DF1sorted.getValue(label,id1,iter1.objId);
			}

			// If we are at the end of DF1, then we did not find id1 such that id1==currentId
			if (!iter1.hasNext())
				break;

			// Grab all the angles in DF1 associated to this id
			anotherIteration=false;
			do
			{
				DF1sorted.getValue(label,id1,iter1.objId);
				anotherIteration=false;
				if (id1==currentId)
				{
					DF1sorted.getValue(MDL_ANGLE_ROT,XX(rotTiltPsi),iter1.objId);
					DF1sorted.getValue(MDL_ANGLE_TILT,YY(rotTiltPsi),iter1.objId);
					DF1sorted.getValue(MDL_ANGLE_PSI,ZZ(rotTiltPsi),iter1.objId);
					DF1sorted.getValue(MDL_SHIFT_X,VEC_ELEM(rotTiltPsi,3),iter1.objId);
					DF1sorted.getValue(MDL_SHIFT_Y,VEC_ELEM(rotTiltPsi,4),iter1.objId);
					DF1sorted.getValue(MDL_FLIP,mirror,iter1.objId);
					if (mirror)
					{
						double rotp, tiltp, psip;
						Euler_mirrorY(XX(rotTiltPsi),YY(rotTiltPsi),ZZ(rotTiltPsi), rotp, tiltp, psip);
						XX(rotTiltPsi)=rotp;
						YY(rotTiltPsi)=tiltp;
						ZZ(rotTiltPsi)=psip;
					}
					ang1.push_back(rotTiltPsi);
					iter1.moveNext();
					if (iter1.hasNext())
						anotherIteration=true;
				}
			} while (anotherIteration);

			// Process both sets of angles
			cumulatedDistance=0;
			N=0;
			for (size_t i=0; i<ang2.size(); ++i)
			{
				const Matrix1D<double> &anglesi=ang2[i];
				double rot2=XX(anglesi);
				double tilt2=YY(anglesi);
				double psi2=ZZ(anglesi);
				double bestDistance=1e38, bestShiftDistance=1e38;
				for (size_t j=0; j<ang1.size(); ++j)
				{
					const Matrix1D<double> &anglesj=ang1[j];
					double rot1=XX(anglesj);
					double tilt1=YY(anglesj);
					double psi1=ZZ(anglesj);
					double dist = SL.computeDistance(rot1, tilt1, psi1, rot2, tilt2, psi2, false,
													 check_mirrors, object_rotation);
					if (dist<bestDistance)
					{
						bestDistance=dist;
						bestShiftDistance=fabs(VEC_ELEM(anglesi,3)-VEC_ELEM(anglesj,3))+fabs(VEC_ELEM(anglesi,4)-VEC_ELEM(anglesj,4));
					}
				}
				bestShiftDistance*=0.5;
				if (bestDistance<360)
				{
					cumulatedDistance+=bestDistance;
					cumulatedDistanceShift+=bestShiftDistance;
					N++;
				}
			}
			newObjId=DFweights.addObject();
			DFweights.setValue(label,currentId,newObjId);
    	}
    	else
    		N=0;

    	if (N>0)
    	{
			double meanDistance=cumulatedDistance/ang2.size();
			DFweights.setValue(angleDiffLabel,meanDistance,newObjId);
			double meanDistanceShift=cumulatedDistanceShift/ang2.size();
			DFweights.setValue(shiftDiffLabel,meanDistanceShift,newObjId);
			ang1.clear();
			ang2.clear();
    	}
    	else
    		if (newObjId>0)
    		{
				DFweights.setValue(angleDiffLabel,-1.0,newObjId);
				DFweights.setValue(shiftDiffLabel,-1.0,newObjId);
				ang1.clear();
				ang2.clear();
    		}
        anotherImageIn2=iter2.hasNext();
    }
    if (ang2.size()>0)
    {
		size_t newObjId=DFweights.addObject();
		DFweights.setValue(label,currentId,newObjId);
		DFweights.setValue(angleDiffLabel,-1.0,newObjId);
		DFweights.setValue(shiftDiffLabel,-1.0,newObjId);
    }

    // If there are more images in MD1 than in MD2, set the last images to 0
    while (iter2.hasNext())
    {
		size_t newObjId=DFweights.addObject();
		DFweights.setValue(label,currentId,newObjId);
		DFweights.setValue(angleDiffLabel,-1.0,newObjId);
		DFweights.setValue(shiftDiffLabel,-1.0,newObjId);
		iter2.moveNext();
    }

    // Calculate the deviation with respect to angleDiff=0 of the angular distances
    std::vector<double> angleDistances, shiftDistances;
	DFweights.getColumnValues(angleDiffLabel,angleDistances);
	DFweights.getColumnValues(shiftDiffLabel,shiftDistances);

    double sigma2=0, sigma2D=0, N=0;
    for (size_t i=0; i<angleDistances.size(); ++i)
    {
    	double d=angleDistances[i];
    	if (d>0)
		{
    		sigma2+=d*d;
        	d=shiftDistances[i];
        	sigma2D+=d*d;
    		++N;
		}
    }
    sigma2/=N;
    sigma2=std::max(minSigma*minSigma,sigma2);
    std::cout << "Sigma of angular distances=" << sqrt(sigma2) << std::endl;
    sigma2D/=N;
    sigma2D=std::max(minSigmaD*minSigmaD,sigma2D);
    std::cout << "Sigma of shift distances=" << sqrt(sigma2D) << std::endl;

    // Adjust the jumper weights according to a Gaussian
    double isigma2=-0.5/sigma2;
    double isigma2D=-0.5/sigma2D;
    FOR_ALL_OBJECTS_IN_METADATA(DFweights)
    {
    	double d;
    	double weight=1.0;
		DFweights.getValue(angleDiffLabel,d,__iter.objId);
		if (d>0)
			weight=exp(d*d*isigma2);
		else
			weight=0.5;

		if (minSigmaD>0)
		{
			DFweights.getValue(shiftDiffLabel,d,__iter.objId);
			if (d>0)
				weight*=exp(d*d*isigma2D);
			else
				weight*=0.5;
		}
		DFweights.setValue(weightLabel,weight,__iter.objId);
    }

    // Transfer these weights to the DF2 metadata
    MetaData DF2weighted;
	if (DF2.containsLabel(angleDiffLabel))
		DF2.removeLabel(angleDiffLabel);
	if (DF2.containsLabel(angleDiffLabel))
		DF2.removeLabel(angleDiffLabel);
	if (DF2.containsLabel(weightLabel))
		DF2.removeLabel(weightLabel);
    DF2weighted.join1(DF2,DFweights,label,INNER);
    FOR_ALL_OBJECTS_IN_METADATA(DF2weighted)
    {
    	double d;
    	DF2weighted.getValue(angleDiffLabel,d,__iter.objId);
    	if (d<0)
    	{
    		// DF2weighted.setValue(MDL_ENABLED,-1,__iter.objId);
    		DF2weighted.setValue(angleDiffLabel,0.0,__iter.objId);
    	}
    	DF2weighted.getValue(shiftDiffLabel,d,__iter.objId);
    	if (d<0)
    	{
    		// DF2weighted.setValue(MDL_ENABLED,-1,__iter.objId);
    		DF2weighted.setValue(shiftDiffLabel,0.0,__iter.objId);
    	}
    }
    DF2weighted.removeDisabled();
    DF2weighted.write(fn_out+"_weights.xmd");
}
