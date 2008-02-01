/***************************************************************************
 *
 * Authors:    Slavica Jonic                slavica.jonic@a3.epfl.ch (2004)
 *             Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
 *
 * Biomedical Imaging Group, EPFL.
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

#include "nma_alignment.h"

#include <iostream>
#include <data/args.h>
#include <data/selfile.h>
#include <data/docfile.h>

// Empty constructor =======================================================
Prog_nma_alignment_prm::Prog_nma_alignment_prm()
{
    currentImg=NULL;
    each_image_produces_an_output = true;
}

// Read arguments ==========================================================
void Prog_nma_alignment_prm::read(int argc, char **argv)
{
    Prog_parameters::read(argc, argv);
    fnPDB=getParameter(argc,argv,"-pdb");
    fnModeList=getParameter(argc,argv,"-modes");
    sampling_rate=textToFloat(getParameter(argc,argv,"-sampling_rate"));
    fnOut=getParameter(argc,argv,"-o");
    produce_side_info();
}

// Show ====================================================================
void Prog_nma_alignment_prm::show()
{
    Prog_parameters::show();
    std::cout << "PDB:           " << fnPDB         << std::endl
              << "Mode list:     " << fnModeList    << std::endl
	      << "Sampling rate: " << sampling_rate << std::endl
	      << "Output:        " << fnOut         << std::endl
    ;
}

// usage ===================================================================
void Prog_nma_alignment_prm::usage()
{
    Prog_parameters::usage();
    std::cerr << "   -pdb <PDB filename>     : PDB Model to compute NMA\n"
              << "   -modes <filename>       : File with a list of mode filenames\n"
	      << "   -sampling_rate <Ts>     : in Angstroms/pixel\n"
	      << "   -o <output filename>    : File for the assignment\n"
    ;
}

// Produce side information ================================================
const Prog_nma_alignment_prm *global_NMA_prog;

void Prog_nma_alignment_prm::produce_side_info()
{
    // If dont_modify_header
    if (dont_modify_header)
        each_image_produces_an_output = false;

    // Read list of modes
    std::ifstream fhModeList;
    fhModeList.open(fnModeList.c_str());
    if (!fhModeList)
       REPORT_ERROR(1,(std::string)"Cannot open "+fnModeList+" for reading");
    while (!fhModeList.eof())
    {
        std::string modeName;
	fhModeList >> modeName;
	if (modeName!="") modeList.push_back(modeName);
    }
    fhModeList.close();
    
    // Get the size of the images in the selfile
    SelFile SFin;
    SFin.read(fn_in);
    int expYdim, expXdim;
    SFin.ImgSize(expYdim, expXdim);
    imgSize=expYdim;
    
    // Set the pointer of the program to this object
    global_NMA_prog=this;
}

// Create deformed PDB =====================================================
FileName Prog_nma_alignment_prm::createDeformedPDB(
    const Matrix1D<double> &trial, int pyramidLevel) const
{
    std::string command;
    FileName fnRandom;
    fnRandom.init_random(20);
    command=(std::string)"./movemode.pl "+fnPDB+" "+modeList[0]+" "+
        floatToString((float)trial(0)*1000)+" > inter; mv -f inter deformedPDB_"+
	fnRandom+".pdb";
    system(command.c_str());
    for (int i=1; i<XSIZE(trial)-5; ++i)
    {
        command=(std::string)"./movemode.pl "+fnRandom+" "+modeList[i]+" "+
            floatToString((float)trial(i)*1000)+" > inter; mv -f inter deformedPDB_"+
	    fnRandom+".pdb";
        system(command.c_str());
    }
    command=(std::string)"xmipp_convert_pdb2vol"+
       " -i deformedPDB_"+fnRandom+".pdb"+
       " -sampling_rate "+floatToString((float)sampling_rate)+
       " -size "+integerToString(ROUND(imgSize))+
       " -quiet";
    system(command.c_str());

    if (pyramidLevel!=0)
    {
        command=(std::string)"xmipp_scale_pyramid"+
           " -i deformedPDB_"+fnRandom+".vol"+
           " -reduce -levels "+integerToString(pyramidLevel)+
           " -quiet";
        system(command.c_str());
    }

    return fnRandom;
}

// Perform complete search =================================================
void Prog_nma_alignment_prm::performCompleteSearch(
    const FileName &fnRandom, int pyramidLevel, int pyramidLevelCont,
    Matrix1D<double> &trial) const
{
    std::string command;
    
    // Prepare image and volume
    if (pyramidLevel!=0)
    {
    	// Reduce the image
    	command=(std::string)"xmipp_scale_pyramid -i "+
	   currentImg->name()+" -o downimg_"+fnRandom+".xmp "+
	   "-reduce -levels "+integerToString(pyramidLevel)+
	   " -quiet";
	system(command.c_str());
    }
    else
    {
        // Make links
    	command=(std::string)"ln -s "+
	   currentImg->name()+" downimg_"+fnRandom+".xmp";
	system(command.c_str());
    }
    command=(std::string)"xmipp_selfile_create "+
       "downimg_"+fnRandom+".xmp > selfile_"+fnRandom+".sel";
    system(command.c_str());
    
    // Perform complete search
    command=(std::string)"xmipp_angular_discrete_assign"+
        " -i selfile_"+fnRandom+".sel"+
	" -ref deformedPDB_"+fnRandom+".vol"+
	" -oang angledisc_"+fnRandom+".txt"+
	" -proj_step 10 -psi_step 10 "+
	" -max_shift_change "+integerToString(ROUND((double)imgSize/
	   (20.0*pow(2.0,(double)pyramidLevel))))+
	" -quiet";
    system(command.c_str());
    
    // Pickup results
    DocFile DF;
    DF.read("angledisc_"+fnRandom+".txt");
    DF.adjust_to_data_line();
    trial(XSIZE(trial)-5)=DF(0);
    trial(XSIZE(trial)-4)=DF(1);
    trial(XSIZE(trial)-3)=DF(2);
    trial(XSIZE(trial)-2)=DF(3)*pow(2.0,(double)pyramidLevel-pyramidLevelCont);
    trial(XSIZE(trial)-1)=DF(4)*pow(2.0,(double)pyramidLevel-pyramidLevelCont);
    
    // Scale shifts if necessary
    if (pyramidLevel!=pyramidLevelCont)
    {
        // Change the selfile
        command=(std::string)"xmipp_selfile_create "+
            currentImg->name()+" > selfile_"+fnRandom+".sel";
        system(command.c_str());

        // Apply header
        DF.set(3,trial(XSIZE(trial)-2));
        DF.set(4,trial(XSIZE(trial)-1));
        DF.write();
	command=(std::string)"xmipp_header_assign"+
	    " -i angledisc_"+fnRandom+".txt"+
	    " -o selfile_"+fnRandom+".sel"+
	    " -force -quiet";
	system(command.c_str());
    }
}

// Continuous assignment ===================================================
double Prog_nma_alignment_prm::performContinuousAssignment(
    const FileName &fnRandom, int pyramidLevel,
    int pyramidLevelDiscrete, Matrix1D<double> &trial) const
{
    std::string command;
    if (currentStage==1)
        command=(std::string)"xmipp_selfile_create "+
           "downimg_"+fnRandom+".xmp > selfile_"+fnRandom+".sel";
    else
        command=(std::string)"xmipp_selfile_create "+
           currentImg->name()+" > selfile_"+fnRandom+".sel";
    system(command.c_str());

    // Perform alignment
    command=(std::string)"xmipp_angular_continuous_assign"+
        " -i selfile_"+fnRandom+".sel"+
	" -ref deformedPDB_"+fnRandom+".vol"+
	" -oang anglecont_"+fnRandom+".txt"+
	" -quiet";
    system(command.c_str());
    
    // Pick up results
    DocFile DF;
    DF.read("anglecont_"+fnRandom+".txt");
    DF.adjust_to_data_line();
    trial(XSIZE(trial)-5)=DF(0);
    trial(XSIZE(trial)-4)=DF(1);
    trial(XSIZE(trial)-3)=DF(2);
    trial(XSIZE(trial)-2)=DF(3)*pow(2.0,(double)pyramidLevel);
    trial(XSIZE(trial)-1)=DF(4)*pow(2.0,(double)pyramidLevel);
    
    return DF(5);
}

// Compute fitness =========================================================
double Prog_nma_alignment_prm::computeFitness(Matrix1D<double> &trial) const
{
    int pyramidLevelDisc=2;
    int pyramidLevelCont=(currentStage==1)?2:0;

    FileName fnRandom=createDeformedPDB(trial,pyramidLevelCont);
    if (currentStage==1)
    {
        performCompleteSearch(fnRandom,pyramidLevelDisc,pyramidLevelCont,trial);
    }
    else
    {
        // Set the header of the image to the best parameters of Stage 1
	ImageXmipp img(currentImg->name());
        img.clear_fFlag_flag();
	double rot =bestStage1(XSIZE(bestStage1)-5);
	double tilt=bestStage1(XSIZE(bestStage1)-4);
	double psi =bestStage1(XSIZE(bestStage1)-3);
	int xshift =bestStage1(XSIZE(bestStage1)-2);
	int yshift =bestStage1(XSIZE(bestStage1)-1);
        img.set_eulerAngles(rot, tilt, psi);
        img.set_originOffsets(xshift, yshift);
	img.write();
    }
    double fitness=performContinuousAssignment(fnRandom,pyramidLevelCont,
        pyramidLevelDisc,trial);
    std::string command=(std::string)"rm -f *"+fnRandom+"*";
    system(command.c_str());
    return fitness;
}

// Assign NMA and alignment parameters =====================================
double wrapperFitness(double *prm)
{
    Matrix1D<double> trial;
    trial.initZeros(global_NMA_prog->modeList.size()+5);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(trial)
        trial(i)=prm[i+1];
    double fitness=global_NMA_prog->computeFitness(trial);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(trial)
        prm[i+1]=trial(i);
    std::cout << "Trial=" << trial.transpose() << " ---> " << fitness << std::endl;
    return fitness;
}

Matrix1D<double> Prog_nma_alignment_prm::assignParameters(const ImageXmipp &img)
{
    Matrix1D<double> parameters;
    parameters.initZeros(modeList.size()+5);
    currentImg=&img;

    Matrix1D<double> steps(XSIZE(parameters));
    double fitness;
    int iter;
    steps.init_constant(1);
    for (int i=1; i<=5; ++i)
        steps(XSIZE(steps)-i)=0;
    currentStage=1;
    powellOptimizer(parameters, 1, XSIZE(steps), &wrapperFitness,
        0.01, fitness, iter, steps, true);
    bestStage1=parameters;
    currentStage=2;
    powellOptimizer(parameters, 1, XSIZE(steps), &wrapperFitness,
        0.1, fitness, iter, steps, true);
    parameters.resize(XSIZE(parameters)+1);
    parameters(XSIZE(parameters)-1)=fitness;
    listAssignments.push_back(parameters);
    img_names.push_back(img.name());
    return parameters;
 }

// Finish computations======================================================
void Prog_nma_alignment_prm::finish_processing()
{
    int p = listAssignments.size();
    DocFile DF;
    DF.reserve(p + 1);
    std::string comment="Headerinfo columns: ";
    for (int i=0; i<modeList.size(); ++i)
       comment+=modeList[i]+" ";
    comment+="Predicted_Rot Predicted_Tilt Predicted_Psi Predicted_ShiftX Predicted_ShiftY Cost";
    DF.append_comment(comment);
    for (int i = 0; i < p; i++)
    {
        DF.append_comment(img_names[i]);
        DF.append_data_line(listAssignments[i]);
    }
    DF.write(fnOut);
}
