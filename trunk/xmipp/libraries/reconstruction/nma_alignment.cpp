/***************************************************************************
 *
 * Authors:    Slavica Jonic                slavica.jonic@impmc.jussieu.fr
 *             Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
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
#include <data/metadata.h>
#include "../../external/condor/ObjectiveFunction.h"
#include "../../external/condor/Solver.h"
#include "../../external/condor/tools.h"
#include "../../external/condor/Vector.h"

// Empty constructor =======================================================
Prog_nma_alignment_prm::Prog_nma_alignment_prm()
{
    MPIversion = false;
    currentImg=NULL;
    each_image_produces_an_output = true;
}

// Read arguments ==========================================================
void Prog_nma_alignment_prm::read(int argc, char **argv)
{
    Prog_parameters::read(argc, argv);
    fnPDB=getParameter(argc,argv,"-pdb");
    fnModeList=getParameter(argc,argv,"-modes");
    scale_defamp=textToFloat(getParameter(argc,argv,"-deformation_scale"));
    sampling_rate=textToFloat(getParameter(argc,argv,"-sampling_rate"));
    symmetry=getParameter(argc,argv,"-sym","c1");
    fnmask=getParameter(argc,argv,"-mask","");
    fnOut=getParameter(argc,argv,"-oang");
    useFixedGaussian = checkParameter(argc, argv, "-fixed_Gaussian");
    if (useFixedGaussian)
        sigmaGaussian=textToFloat(getParameter(argc,argv,"-fixed_Gaussian","-1"));
    do_centerPDB = checkParameter(argc, argv, "-centerPDB");
    
    if (!MPIversion) produce_side_info();
}

// Show ====================================================================
void Prog_nma_alignment_prm::show()
{
    Prog_parameters::show();
    std::cout << "PDB:                " << fnPDB            << std::endl
              << "Mode list:          " << fnModeList       << std::endl
              << "Amplitude scale:    " << scale_defamp     << std::endl
	      << "Sampling rate:      " << sampling_rate    << std::endl
	      << "Symmetry:           " << symmetry         << std::endl
              << "Mask:               " << fnmask           << std::endl
              << "Use fixed Gaussian: " << useFixedGaussian << std::endl
	      << "Sigma of Gaussian:  " << sigmaGaussian    << std::endl
              << "Center PDB:         " << do_centerPDB     << std::endl
	      << "Output:             " << fnOut            << std::endl
    ;
}

// usage ===================================================================
void Prog_nma_alignment_prm::usage()
{
    Prog_parameters::usage();
    std::cerr << "   -pdb <PDB filename>                : PDB Model to compute NMA\n"
              << "   -modes <filename>                  : File with a list of mode filenames\n"
              << "   -deformation_scale                 : Scaling factor to scale deformation amplitude\n"
	      << "   -sampling_rate <Ts>                : in Angstroms/pixel\n"
	      << "  [-sym]                              : Symmetry file or point group\n"
              << "  [-mask]                             : Mask\n"
              << "  [-fixed_Gaussian <std>]             : For pseudo atoms fixed_Gaussian must be used. The the standard deviation <std> may or may not be used\n"
              << "  [-centerPDB]                        : Center the PDB structure\n"
	      << "   -oang <output filename>            : File for the assignment\n"
    ;
}

// Produce side information ================================================
const Prog_nma_alignment_prm *global_NMA_prog;

void Prog_nma_alignment_prm::produce_side_info(int rank)
{
    // Read list of modes
    MetaData SFmodelist;
    SFmodelist.read(fnModeList);
    FileName tempname;
    int i=0;
    FOR_ALL_OBJECTS_IN_METADATA(SFmodelist){      
       SFmodelist.getValue(MDL_IMAGE,tempname);
       modeList.push_back(tempname);
       //std::cout << "MODE NAME=" << tempname  << " " << std::endl;
       i++;
    }
    
    // Get the size of the images in the selfile
    MetaData SFin;
    SFin.read(fn_in);
    //FileName tempname;
    SFin.getValue(MDL_IMAGE,tempname);
    ImageXmipp tempimage;
    tempimage.read(tempname);
    
    imgSize=YSIZE(tempimage());
    int imgnumber=SFin.size();

//std::cout << "IMAGESIZE=" << imgSize  << "" << "NAME" << tempname << std::endl;
    
    // Set the pointer of the program to this object
    global_NMA_prog=this;
        
    rangen=rank;
}

// Create deformed PDB =====================================================
FileName Prog_nma_alignment_prm::createDeformedPDB(int pyramidLevel) const
{
    std::string command;
    FileName fnRandom;

    fnRandom.init_random(19);  
    
    fnRandom=fnRandom+integerToString(rangen);  //the last character is the rank
    
    command=(std::string)"xmipp_move_along_NMAmode "+fnPDB+" "+modeList[0]+" "+
        floatToString((float)trial(0)*scale_defamp)+" > inter"+fnRandom+"; mv -f inter"+fnRandom+" deformedPDB_"+
	fnRandom+".pdb";
    system(command.c_str());
  
    for (int i=1; i<XSIZE(trial)-5; ++i)
    {
	command=(std::string)"xmipp_move_along_NMAmode deformedPDB_"+fnRandom+".pdb "+modeList[i]+" "+
            floatToString((float)trial(i)*scale_defamp)+" > inter"+fnRandom+"; mv -f inter"+fnRandom+" deformedPDB_"+
	    fnRandom+".pdb";

        system(command.c_str());	
    }

    command=(std::string)"xmipp_convert_pdb2vol"+
       " -i deformedPDB_"+fnRandom+".pdb"+
       " -size "+integerToString(ROUND(imgSize))+
       " -sampling_rate " + floatToString((float)sampling_rate) +
       " >& /dev/null";
    
    if (do_centerPDB) command+=(std::string)" -centerPDB "; 

    if (useFixedGaussian)
    {
       if (sigmaGaussian<0) 
       {
          command+=(std::string)" -fixed_Gaussian ";
       }
       else 
       {
          command+=(std::string)" -fixed_Gaussian " + floatToString((float)sigmaGaussian)+ " -intensityColumn Bfactor ";
       }
    }


    system(command.c_str());

//std::cout << "COMMAND=" << command  << std::endl;

    command=(std::string)"xmipp_fourier_filter"+
       " -i deformedPDB_"+fnRandom+".vol"+
       " -sampling "+floatToString((float)sampling_rate)+
       " -low_pass 15 -fourier_mask raised_cosine 0.1 >& /dev/null";

    system(command.c_str());

    if (pyramidLevel!=0)
    {
        command=(std::string)"xmipp_scale_pyramid"+
           " -i deformedPDB_"+fnRandom+".vol"+
           " -reduce -levels "+integerToString(pyramidLevel)+
           " -quiet";

        system(command.c_str());
    }

//std::cout << "COMMAND=" << command  << std::endl;

    return fnRandom;
}

// Perform complete search =================================================
void Prog_nma_alignment_prm::performCompleteSearch(
    const FileName &fnRandom, int pyramidLevel) const
{
    std::string command;
    
    // Reduce the image
    command=(std::string)"xmipp_scale_pyramid -i "+
       currentImg->name()+" -o downimg_"+fnRandom+".xmp "+
       "-reduce -levels "+integerToString(pyramidLevel)+
       " -quiet";

    system(command.c_str());

//std::cout << "COMMAND=" << command  << std::endl;

    command=(std::string)"xmipp_selfile_create "+
       "downimg_"+fnRandom+".xmp > selfile_"+fnRandom+".sel";

    system(command.c_str());
    
    command=(std::string)"mkdir ref" + fnRandom;
    system(command.c_str());	
        
    command = (std::string)"-i deformedPDB_"+fnRandom+".vol"+ " -o ref" + fnRandom + "/ref" + fnRandom +
        "_ -sampling_rate 5 -sym " + symmetry;
    command=(std::string)
        " xmipp_angular_project_library "+command + " -quiet";

    system(command.c_str());
    
//std::cout << "COMMAND=" << command  << std::endl;

    command=(std::string)"xmipp_selfile_create \"ref" + fnRandom + "/ref*.xmp\" > ref"+fnRandom+"_.sel";
    system(command.c_str());
    
    command=(std::string)"mv ref" + fnRandom + "/*.doc .";
    system(command.c_str());

    if (fnmask != "") 
    {
        command=(std::string) " xmipp_mask -i ref"+fnRandom+"_.sel -mask "+fnmask;

        system(command.c_str());
    }
		
    command=(std::string)" xmipp_header_extract -i selfile_"+fnRandom+".sel -o docexp"+fnRandom+".txt";

    system(command.c_str());

    command=(std::string)" xmipp_angular_discrete_assign"+
        " -i docexp"+fnRandom+".txt" +
	" -ref ref"+fnRandom+"_.sel"+
	" -oang angledisc_"+fnRandom+".txt"+
	" -psi_step 5 "+
	" -max_shift_change "+integerToString(ROUND((double)imgSize/
	   (10.0*pow(2.0,(double)pyramidLevel))))+
	" -5D -sym " + symmetry + " -quiet";

    system(command.c_str()); 
//std::cout << "COMMAND=" << command  << std::endl; 

}

// Continuous assignment ===================================================
double Prog_nma_alignment_prm::performContinuousAssignment(
    const FileName &fnRandom, int pyramidLevel) const
{
    std::string command;
    if (pyramidLevel==0)
    {
        // Make links
        command=(std::string)"ln -sf "+
           currentImg->name()+" downimg_"+fnRandom+".xmp";
        system(command.c_str());
    }

    // Perform alignment
    command=(std::string)"xmipp_angular_continuous_assign"+
        " -ang angledisc_"+fnRandom+".txt"+
	" -ref deformedPDB_"+fnRandom+".vol"+
	" -oang anglecont_"+fnRandom+".txt"+
	" -quiet";
    system(command.c_str());

//   std::cout << "COMMAND=" << command  << std::endl;
    
    // Pick up results
    MetaData DF;
    DF.read("anglecont_"+fnRandom+".txt");
    DF.getValue(MDL_ANGLEROT,trial(XSIZE(trial)-5));
    DF.getValue(MDL_ANGLETILT,trial(XSIZE(trial)-4));
    DF.getValue(MDL_ANGLEPSI,trial(XSIZE(trial)-3));
    DF.getValue(MDL_SHIFTX,trial(XSIZE(trial)-2));trial(XSIZE(trial)-2)*=pow(2.0,(double)pyramidLevel);
    DF.getValue(MDL_SHIFTY,trial(XSIZE(trial)-1));trial(XSIZE(trial)-1)*=pow(2.0,(double)pyramidLevel);
    double tempvar;
    DF.getValue(MDL_COST,tempvar);
    return tempvar;
}

void Prog_nma_alignment_prm::update_bestfit(double fitness,int dim) const
{  
    if (fitness < fitness_min(0))
    { 
	fitness_min(0) = fitness;
	for (int i=0; i<dim+5; i++)
        {
            trial_best(i)=trial(i);
        }
    }
}

// Compute fitness =========================================================
double ObjFunc_nma_alignment::eval(Vector X, int *nerror)
{
    int dim=global_NMA_prog->modeList.size();
    
    for (int i=0; i<dim; i++)
    {
        global_NMA_prog->trial(i)=X[i];
    }

    int pyramidLevelDisc=1;
    int pyramidLevelCont=(global_NMA_prog->currentStage==1)?1:0;
    
    FileName fnRandom=global_NMA_prog->createDeformedPDB(pyramidLevelCont);
    if (global_NMA_prog->currentStage==1)
    {
         global_NMA_prog->performCompleteSearch(fnRandom,pyramidLevelDisc);
    }
    else
    {
        double rot, tilt, psi, xshift, yshift;
        MetaData DF;
        
	rot = global_NMA_prog->bestStage1(XSIZE(global_NMA_prog->bestStage1)-5);
	tilt = global_NMA_prog->bestStage1(XSIZE(global_NMA_prog->bestStage1)-4);
	psi = global_NMA_prog->bestStage1(XSIZE(global_NMA_prog->bestStage1)-3);
	xshift = global_NMA_prog->bestStage1(XSIZE(global_NMA_prog->bestStage1)-2);
	yshift = global_NMA_prog->bestStage1(XSIZE(global_NMA_prog->bestStage1)-1);
       
	    DF.addObject();
	    DF.setValue(MDL_IMAGE,"downimg_"+fnRandom+".xmp");
	    DF.setValue(MDL_ANGLEROT,rot);
	    DF.setValue(MDL_ANGLETILT,tilt);
	    DF.setValue(MDL_ANGLEPSI,psi);
	    DF.setValue(MDL_SHIFTX,xshift);
	    DF.setValue(MDL_SHIFTY,yshift);
	    
        DF.write("angledisc_"+fnRandom+".txt");

    }
    double fitness=global_NMA_prog->performContinuousAssignment(fnRandom,pyramidLevelCont);

    
   
    std::string command=(std::string)"rm -rf *"+fnRandom+"* &";
    system(command.c_str());
    
    //std::cout << "Trial=" << global_NMA_prog->trial.transpose() << " ---> " << fitness << std::endl;
    
    global_NMA_prog->update_bestfit(fitness,dim);
    
    return fitness;
}

ObjFunc_nma_alignment::ObjFunc_nma_alignment(int _t, int _n)
{

}

void Prog_nma_alignment_prm::assignParameters(ImageXmipp &img)
{
    
    FileName imgname=img.name();

//std::cout << "NAME=" << imgname  << std::endl;
    
    //double rhoStart=1e-0, rhoEnd=1e-4;
    double rhoStart=1e-0, rhoEnd=1e-3;
    
    int niter=1000;

    ObjectiveFunction *of;

    int dim=modeList.size();

    parameters.initZeros(dim+5);
    currentImg=&img;

    trial.initZeros(dim+5);
    trial_best.initZeros(dim+5);

    fitness_min.initZeros(1);
    fitness_min(0)=1000000.0;

    currentStage=1;

    of=new ObjFunc_nma_alignment(1,dim);
    
    of->xStart.setSize(dim);
    for (int i=0; i<dim; i++) of->xStart[i]=0.;
    
    strcpy(of->name,("OF1_"+integerToString(rangen)).c_str());
    of->setSaveFile();

    CONDOR(rhoStart, rhoEnd, niter, of);
    of->printStats();
        FILE *ff=fopen(("res1_"+integerToString(rangen)+".txt").c_str(),"w");
    fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
    fclose(ff);

    double fitness=of->valueBest;
    //std::cout << "Best fitness = " << fitness << std::endl;
    double *dd=of->xBest;
    /*for (int i=0; i<dim; i++)
    {
        std::cout << "Best deformations = " << dd[i] << std::endl;
    }*/

    for (int i=0; i<dim+5; i++)
    {
        parameters(i)=trial_best(i);
        trial(i)=trial_best(i);
    }
    
    bestStage1=trial;
    
    delete of;
    
    currentStage=2;
    fitness_min(0)=1000000.0;
    
    of=new ObjFunc_nma_alignment(1,dim);

    of->xStart.setSize(dim);
    for (int i=0; i<dim; i++) of->xStart[i]=parameters(i);
    strcpy(of->name,("OF2_"+integerToString(rangen)).c_str());
    of->setSaveFile();
    
    rhoStart=1e-3, rhoEnd=1e-4;
    //rhoStart=1e-3, rhoEnd=1e-3;
    
    CONDOR(rhoStart, rhoEnd, niter, of);
    of->printStats();
    ff=fopen(("res2_"+integerToString(rangen)+".txt").c_str(),"w");
    fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
    fclose(ff);

    fitness=of->valueBest;
    //std::cout << "Best fitness = " << fitness << std::endl;
    dd=of->xBest;
    /*for (int i=0; i<dim; i++)
    {
        std::cout << "Best deformations = " << dd[i] << std::endl;
    }*/

    for (int i=0; i<dim+5; i++)
    {
        trial(i)=trial_best(i);
    }

    for (int i=dim; i<dim+5; i++)
    {
        parameters(i-dim)=trial_best(i);
    }

    for (int i=0; i<dim; i++)
    {
        parameters(5+i)=trial_best(i)*scale_defamp;
    }

    img.read(imgname);
   
    parameters.resize(XSIZE(parameters)+1);
    parameters(XSIZE(parameters)-1)=fitness_min(0);
    listAssignments.push_back(parameters);
    img_names.push_back(imgname);
           
    	DF_out.addObject();
    	DF_out.setValue(MDL_IMAGE,imgname);
    	DF_out.setValue(MDL_ANGLEROT,parameters(0));
    	DF_out.setValue(MDL_ANGLETILT,parameters(1));
    	DF_out.setValue(MDL_ANGLEPSI,parameters(2));
    	DF_out.setValue(MDL_SHIFTX,parameters(3));
    	DF_out.setValue(MDL_SHIFTY,parameters(4));
    	     
        std::vector<double> vectortemp;
        for (int j = 5; j < 5+dim; j++)
        {
        	vectortemp.push_back(parameters(j));
        }
        DF_out.setValue(MDL_NMA,vectortemp);
        DF_out.setValue(MDL_COST,parameters(5+dim));
        
        DF_out.write(fnOut+integerToString(rangen));
 
    delete of;
 }

// Finish computations======================================================
void Prog_nma_alignment_prm::finish_processing()
{  
        int p = listAssignments.size();
        MetaData DF;
               
        for (int i = 0; i < p; i++)
        {
        	DF.addObject();
        	DF.setValue(MDL_IMAGE,img_names[i]);
        	DF.setValue(MDL_ANGLEROT,listAssignments[i](0));
        	DF.setValue(MDL_ANGLETILT,listAssignments[i](1));
        	DF.setValue(MDL_ANGLEPSI,listAssignments[i](2));
        	DF.setValue(MDL_SHIFTX,listAssignments[i](3));
        	DF.setValue(MDL_SHIFTY,listAssignments[i](4));
        	
        	int xsz=XSIZE(listAssignments[i]);
            std::vector<double> vectortemp;
            for (int j = 5; j < xsz-1; j++)
            {
            	vectortemp.push_back(listAssignments[i](j));
            }
            DF.setValue(MDL_NMA,vectortemp);
            DF.setValue(MDL_COST,listAssignments[i](xsz-1));
        }
        DF.write(fnOut);
}
