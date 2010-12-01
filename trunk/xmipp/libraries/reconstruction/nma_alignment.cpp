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
#include <data/metadata_extension.h>

#include <iostream>
#include <data/args.h>
#include <unistd.h>
#include <sys/stat.h>
#include "../../external/condor/Solver.h"
#include "../../external/condor/tools.h"

// Empty constructor =======================================================
ProgNmaAlignment::ProgNmaAlignment()
{
    rangen = 0;
    currentImgName="";
    each_image_produces_an_output = true;
}

// Params definition ============================================================
void ProgNmaAlignment::defineParams()
{
    XmippMetadataProgram::defineParams();
    addParamsLine("   -pdb <PDB_filename>                : PDB Model to compute NMA");
    addParamsLine("   -oang <output_filename>            : File for the assignment");
    addParamsLine("   -modes <filename>                  : File with a list of mode filenames");
    addParamsLine("  [-deformation_scale <s=1>]          : Scaling factor to scale deformation amplitude");
    addParamsLine("  [-sampling_rate <Ts=1>]             : in Angstroms/pixel");
    addParamsLine("  [-mask <m=\"\">]                    : Mask");
    addParamsLine("  [-gaussian_Fourier <s=0.5>]         : Weighting sigma in Fourier space");
    addParamsLine("  [-gaussian_Real    <s=0.5>]         : Weighting sigma in Real space");
    addParamsLine("  [-zerofreq_weight  <s=0.>]          : Zero-frequency weight");
    addParamsLine("  [-centerPDB]                        : Center the PDB structure");
    addParamsLine("  [-filterVol <cutoff=15.>]           : Flter the volume from the PDB structure. Default cut-off is 15 A.");
    addParamsLine("  [-fixed_Gaussian <std=-1>]          : For pseudo atoms fixed_Gaussian must be used.");
    addParamsLine("                                      : Default standard deviation <std> is read from PDB file.");
}

// Read arguments ==========================================================
void ProgNmaAlignment::readParams()
{
    XmippMetadataProgram::readParams();
    fnPDB = getParam("-pdb");
    fnOut = getParam("-oang");
    fnModeList = getParam("-modes");
    scale_defamp = getDoubleParam("-deformation_scale");
    sampling_rate = getDoubleParam("-sampling_rate");
    fnmask = getParam("-mask");
    gaussian_DFT_sigma = getDoubleParam( "-gaussian_Fourier");
    gaussian_Real_sigma = getDoubleParam( "-gaussian_Real");
    weight_zero_freq = getDoubleParam( "-zerofreq_weight");
    do_centerPDB = checkParam("-centerPDB");
    do_FilterPDBVol = checkParam("-filterVol");
    if (do_FilterPDBVol)
        cutoff_LPfilter = getDoubleParam("-filterVol");
    useFixedGaussian = checkParam("-fixed_Gaussian");
    if (useFixedGaussian)
        sigmaGaussian = getDoubleParam("-fixed_Gaussian");

}

// Show ====================================================================
void ProgNmaAlignment::show()
{
    XmippMetadataProgram::show();
    std::cout << "PDB:                 " << fnPDB               << std::endl
    << "Output:              " << fnOut               << std::endl
    << "Mode list:           " << fnModeList          << std::endl
    << "Amplitude scale:     " << scale_defamp        << std::endl
    << "Sampling rate:       " << sampling_rate       << std::endl
    << "Mask:                " << fnmask              << std::endl
    << "Gaussian Fourier:    " << gaussian_DFT_sigma  << std::endl
    << "Gaussian Real:       " << gaussian_Real_sigma << std::endl
    << "Zero-frequency weight:"<< weight_zero_freq    << std::endl
    << "Center PDB:          " << do_centerPDB        << std::endl
    << "Filter PDB volume    " << do_FilterPDBVol     << std::endl
    << "Use fixed Gaussian:  " << useFixedGaussian    << std::endl
    << "Sigma of Gaussian:   " << sigmaGaussian       << std::endl
    ;
}



// Produce side information ================================================
ProgNmaAlignment *global_NMA_prog;

void ProgNmaAlignment::createWorkFiles()
{
  MetaData mdTodo, mdDone;
  mdTodo = mdIn;
  if (exists("nmaDone.xmd"))
  {
    mdDone.read("nmaDone.xmd");
    mdTodo.subtraction(mdDone, MDL_IMAGE);
  }
  else//if not exists create metadata only with headers
  {
    mdDone.addLabel(MDL_IMAGE);
    mdDone.addLabel(MDL_ENABLED);
    mdDone.addLabel(MDL_ANGLEROT);
    mdDone.addLabel(MDL_ANGLETILT);
    mdDone.addLabel(MDL_ANGLEPSI);
    mdDone.addLabel(MDL_SHIFTX);
    mdDone.addLabel(MDL_SHIFTY);
    mdDone.addLabel(MDL_NMA);
    mdDone.addLabel(MDL_COST);
    mdDone.write("nmaDone.xmd");
  }
  mdIn = mdTodo;
}

void ProgNmaAlignment::preProcess()
{
    // Read list of modes
    MetaData SFmodelist(fnModeList);
    FileName tempname;
    FOR_ALL_OBJECTS_IN_METADATA(SFmodelist)
    {
        SFmodelist.getValue(MDL_IMAGE,tempname);
        modeList.push_back(tempname);
    }
    // Get the size of the images in the selfile
    ImgSize(mdIn, imgSize);
    // Set the pointer of the program to this object
    global_NMA_prog = this;
    //create some neededs files
    createWorkFiles();
}

// Create deformed PDB =====================================================
FileName ProgNmaAlignment::createDeformedPDB(int pyramidLevel) const
{
    std::string command;
    FileName fnRandom;

    fnRandom.initRandom(19);

    fnRandom=fnRandom+integerToString(rangen);  //the last character is the rank

    command=(std::string)"xmipp_move_along_NMAmode "+fnPDB+" "+modeList[0]+" "+
            floatToString((float)trial(0)*scale_defamp)+" > inter"+fnRandom+"; mv -f inter"+fnRandom+" deformedPDB_"+
            fnRandom+".pdb";
    system(command.c_str());

    for (int i=1; i<VEC_XSIZE(trial)-5; ++i)
    {
        command=(std::string)"xmipp_move_along_NMAmode deformedPDB_"+fnRandom+".pdb "+modeList[i]+" "+
                floatToString((float)trial(i)*scale_defamp)+" > inter"+fnRandom+"; mv -f inter"+fnRandom+" deformedPDB_"+
                fnRandom+".pdb";

        system(command.c_str());
    }

    command=(std::string)"xmipp_convert_pdb2vol"+
            " -i deformedPDB_"+fnRandom+".pdb"+
            " -size "+integerToString(ROUND(imgSize))+
            " -sampling_rate " + floatToString((float)sampling_rate);

    if (do_centerPDB)
        command+=(std::string)" -centerPDB ";

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
    command+=" >& /dev/null";

    system(command.c_str());

    if (do_FilterPDBVol)
    {
        command=(std::string)"xmipp_fourier_filter"+
                " -i deformedPDB_"+fnRandom+".vol"+
                " -sampling "+floatToString((float)sampling_rate)+
                " -low_pass " + floatToString((float)cutoff_LPfilter) + " -fourier_mask raised_cosine 0.1 >& /dev/null";
        system(command.c_str());
    }

    if (pyramidLevel!=0)
    {
        command=(std::string)"xmipp_scale_pyramid"+
                " -i deformedPDB_"+fnRandom+".vol"+
                " -reduce -levels "+integerToString(pyramidLevel)+
                " -v 0";

        system(command.c_str());
    }

    return fnRandom;
}

// Perform complete search =================================================
void ProgNmaAlignment::performCompleteSearch(
    const FileName &fnRandom, int pyramidLevel) const
{
    std::string command;

    // Reduce the image
    FileName fnDown=(std::string)"downimg_"+fnRandom+".xmp";
    if (pyramidLevel!=0)
    {
    	Image<double> I, Ireduced;
    	I.read(currentImgName);
    	reduceBSpline(BSPLINE3,Ireduced(),I());
    	Ireduced.write(fnDown);
    }
    else
        link(currentImgName.c_str(),fnDown.c_str());

    mkdir(((std::string)"ref" + fnRandom).c_str(),S_IRWXU);

    command = (std::string)"xmipp_angular_project_library -i deformedPDB_"+fnRandom+".vol"+
              " -o ref" + fnRandom + "/ref" + fnRandom+".stk"
              " --sampling_rate 25 -v 0";
    FileName fnRefSel=(std::string)"ref" + fnRandom + "/ref" + fnRandom+".doc";
    system(command.c_str());

    if (fnmask != "")
    {
        command=(std::string) " xmipp_mask -i "+fnRefSel+" -mask "+fnmask;
        system(command.c_str());
    }

    // Perform alignment
    command=(std::string)" xmipp_angular_discrete_assign"+
            " -i downimg_"+fnRandom+".xmp"+
            " -ref "+fnRefSel+
            " -o angledisc_"+fnRandom+".txt"+
            " -psi_step 5 "+
            " -max_shift_change "+integerToString(ROUND((double)imgSize/
                                                  (10.0*pow(2.0,(double)pyramidLevel))))+
            " -search5D -v 0";
    system(command.c_str());
}

// Continuous assignment ===================================================
double ProgNmaAlignment::performContinuousAssignment(
    const FileName &fnRandom, int pyramidLevel) const
{
    // Perform alignment
    std::string command=(std::string)"xmipp_angular_continuous_assign"+
            " -i angledisc_"+fnRandom+".txt"+
            " -ref deformedPDB_"+fnRandom+".vol"+
            " -o anglecont_"+fnRandom+".txt"+
            " -gaussian_Fourier " + floatToString((float)gaussian_DFT_sigma) +
            " -gaussian_Real " + floatToString((float)gaussian_Real_sigma) +
            " -zerofreq_weight " + floatToString((float)weight_zero_freq) +
            " -v 0";
    system(command.c_str());

    // Pick up results
    MetaData DF;
    DF.read("anglecont_"+fnRandom+".txt");
    DF.getValue(MDL_ANGLEROT,trial(VEC_XSIZE(trial)-5));
    DF.getValue(MDL_ANGLETILT,trial(VEC_XSIZE(trial)-4));
    DF.getValue(MDL_ANGLEPSI,trial(VEC_XSIZE(trial)-3));
    DF.getValue(MDL_SHIFTX,trial(VEC_XSIZE(trial)-2));
    trial(VEC_XSIZE(trial)-2)*=pow(2.0,(double)pyramidLevel);
    DF.getValue(MDL_SHIFTY,trial(VEC_XSIZE(trial)-1));
    trial(VEC_XSIZE(trial)-1)*=pow(2.0,(double)pyramidLevel);
    double tempvar;
    DF.getValue(MDL_COST,tempvar);
    return tempvar;
}

void ProgNmaAlignment::updateBestFit(double fitness,int dim)
{
    if (fitness < fitness_min(0))
    {
        fitness_min(0) = fitness;
        trial_best = trial;
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

        rot = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-5);
        tilt = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-4);
        psi = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-3);
        xshift = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-2);
        yshift = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-1);

        DF.addObject();
        FileName fnDown = (std::string)"downimg_"+fnRandom+".xmp";
        DF.setValue(MDL_IMAGE,fnDown);
        DF.setValue(MDL_ENABLED,1);
        DF.setValue(MDL_ANGLEROT,rot);
        DF.setValue(MDL_ANGLETILT,tilt);
        DF.setValue(MDL_ANGLEPSI,psi);
        DF.setValue(MDL_SHIFTX,xshift);
        DF.setValue(MDL_SHIFTY,yshift);

        DF.write((std::string)"angledisc_"+fnRandom+".txt");
        link(global_NMA_prog->currentImgName.c_str(),fnDown.c_str());
    }
    double fitness=global_NMA_prog->performContinuousAssignment(fnRandom,pyramidLevelCont);

    std::string command = (std::string)"rm -rf *"+fnRandom+"* &";
    //std::cerr << "cleanning: " << command << std::endl;
    system(command.c_str());

    //std::cout << "Trial=" << global_NMA_prog->trial.transpose() << " ---> " << fitness << std::endl;

    global_NMA_prog->updateBestFit(fitness,dim);

    return fitness;
}

ObjFunc_nma_alignment::ObjFunc_nma_alignment(int _t, int _n)
{
}

void ProgNmaAlignment::processImage(const FileName &fnImg, const FileName &fnImgOut, long int objId)
{
    double rhoStart=1e-0, rhoEnd=1e-3;

    int niter=1000;

    ObjectiveFunction *of;

    int dim=modeList.size();

    parameters.initZeros(dim+5);
    currentImgName=fnImg;

    trial.initZeros(dim+5);
    trial_best.initZeros(dim+5);

    fitness_min.initZeros(1);
    fitness_min(0)=1000000.0;

    currentStage=1;
    //std::cerr << std::endl << "DEBUG: ===== Processing image " << fnImg << " at stage: " << currentStage << std::endl;
    of=new ObjFunc_nma_alignment(1,dim);

    of->xStart.setSize(dim);
    for (int i=0; i<dim; i++)
        of->xStart[i]=0.;

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

    bestStage1 = trial = parameters = trial_best;

    delete of;

    currentStage = 2;
    //std::cerr << std::endl << "DEBUG: ===== Processing image " << fnImg << " at stage: " << currentStage << std::endl;
    fitness_min(0) = 1000000.0;

    of=new ObjFunc_nma_alignment(1,dim);

    of->xStart.setSize(dim);
    for (int i=0; i<dim; i++)
        of->xStart[i]=parameters(i);
    strcpy(of->name,("OF2_"+integerToString(rangen)).c_str());
    of->setSaveFile();

    rhoStart=1e-3, rhoEnd=1e-4;

    CONDOR(rhoStart, rhoEnd, niter, of);
    of->printStats();
    ff=fopen(("res2_"+integerToString(rangen)+".txt").c_str(),"w");
    fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
    fclose(ff);

    fitness=of->valueBest;
    std::cout << "Best fitness = " << fitness << std::endl;
    dd=of->xBest;
    for (int i=0; i<dim; i++)
    {
        std::cout << "Best deformations = " << dd[i] << std::endl;
    }

    trial = trial_best;

    for (int i=dim; i<dim+5; i++)
    {
        parameters(i-dim)=trial_best(i);
    }

    for (int i=0; i<dim; i++)
    {
        parameters(5+i)=trial_best(i)*scale_defamp;
    }

    parameters.resize(VEC_XSIZE(parameters)+1);
    parameters(VEC_XSIZE(parameters)-1)=fitness_min(0);

    writeImageParameters(fnImg);
    delete of;
}

void ProgNmaAlignment::writeImageParameters(const FileName &fnImg)
{
  MetaData md;
  md.addObject();
  md.setValue(MDL_IMAGE,fnImg);
  md.setValue(MDL_ENABLED,1);
  md.setValue(MDL_ANGLEROT,parameters(0));
  md.setValue(MDL_ANGLETILT,parameters(1));
  md.setValue(MDL_ANGLEPSI,parameters(2));
  md.setValue(MDL_SHIFTX,parameters(3));
  md.setValue(MDL_SHIFTY,parameters(4));

  int dim=modeList.size();
  std::vector<double> vectortemp;
  for (int j = 5; j < 5+dim; j++)
  {
      vectortemp.push_back(parameters(j));
  }

  md.setValue(MDL_NMA,vectortemp);
  md.setValue(MDL_COST,parameters(5+dim));

  md.append("nmaDone.xmd");
}

// Finish computations======================================================
void ProgNmaAlignment::postProcess()
{
//    int p = assignments.size();
//    MetaData DF;
//
//    for (int i = 0; i < p; i++)
//    {
//        DF.addObject();
//        DF.setValue(MDL_IMAGE,img_names[i]);
//        DF.setValue(MDL_ENABLED,1);
//        DF.setValue(MDL_ANGLEROT,assignments[i](0));
//        DF.setValue(MDL_ANGLETILT,assignments[i](1));
//        DF.setValue(MDL_ANGLEPSI,assignments[i](2));
//        DF.setValue(MDL_SHIFTX,assignments[i](3));
//        DF.setValue(MDL_SHIFTY,assignments[i](4));
//
//        int xsz=VEC_XSIZE(assignments[i]);
//        std::vector<double> vectortemp;
//        for (int j = 5; j < xsz-1; j++)
//        {
//            vectortemp.push_back(assignments[i](j));
//        }
//        DF.setValue(MDL_NMA,vectortemp);
//        DF.setValue(MDL_COST,assignments[i](xsz-1));
//    }
//    DF.write(fnOut);
}
