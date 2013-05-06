/***************************************************************************
 *
 * Authors:    Slavica Jonic                slavica.jonic@impmc.upmc.fr
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
#include "data/metadata_extension.h"
#include "program_extension.h"

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include "../../external/condor/Solver.h"
#include "../../external/condor/tools.h"

// Empty constructor =======================================================
ProgNmaAlignment::ProgNmaAlignment() {
	rangen = 0;
	resume = false;
	currentImgName = "";
	each_image_produces_an_output = false;
	produces_an_output = true;
	progVolumeFromPDB = new ProgPdbConverter();
	projMatch = false;
}

ProgNmaAlignment::~ProgNmaAlignment() {
	delete progVolumeFromPDB;
}

// Params definition ============================================================
void ProgNmaAlignment::defineParams() {
	addUsageLine("Align images with an atomic or pseudo-atomic (from EM volume) structure by computing deformation amplitudes along normal modes, three Euler angles, and two in-plane shifts");
        defaultComments["-i"].clear();
	defaultComments["-i"].addComment("Metadata with image filenames");
	defaultComments["-o"].clear();
	defaultComments["-o"].addComment("Metadata with output Euler angles, shifts, and deformation amplitudes");
	XmippMetadataProgram::defineParams();
	addParamsLine("   --pdb <PDB_filename>                : Reference atomic or pseudo-atomic structure in PDB format");
	addParamsLine("  [--odir <outputDir=\".\">]           : Output directory");
	addParamsLine("  [--resume]                           : Resume processing");
	addParamsLine("==Generation of the deformed reference volumes==");
	addParamsLine("   --modes <filename>                  : File with a list of mode filenames");
	addParamsLine("  [--sampling_rate <Ts=1>]             : Pixel size in Angstroms");
	addParamsLine("  [--filterVol <cutoff=15.>]           : Filter the volume after deforming. If this option is used, the default cut-off is 15 A.");
	addParamsLine("  [--centerPDB]                        : Center the PDB structure");
	addParamsLine("  [--fixed_Gaussian <std=-1>]          : For pseudo-atoms fixed_Gaussian must be used.");
	addParamsLine("                                       : Default standard deviation <std> is read from the PDB file.");
	addParamsLine("==Combined elastic and rigid-body alignment==");
	addParamsLine("  [--trustradius_scale <s=1>]          : Positive scaling factor to scale the initial trust region radius");
	addParamsLine("  [--mask <m=\"\">]                    : 2D masking  of the projections of the deformed volume");
	addParamsLine("  [--projMatch]                        : Real-space projection matching instead of the default combined global (wavelet-space) and local (Fourier central slice) projection matching");
	addParamsLine("                                       :+Note that wavelet-based method needs the image size to be power of 2");
	addParamsLine("  [--discrAngStep <ang=10>]            : Angular sampling step for computing the reference projections (this is the minimum step for real-space projection matching while the minimum step does not exist in the default method)");
	addParamsLine("  [--gaussian_Fourier <s=0.5>]         : Sigma of Gaussian weigthing in Fourier space (parameter of central-slice method)");
	addParamsLine("  [--gaussian_Real    <s=0.5>]         : Sigma of Gaussian weigthing in real space for spline interpolation in Fourier space (parameter of central-slice method)");
	addParamsLine("  [--zerofreq_weight  <s=0.>]          : Zero-frequency weight (parameter of central-slice method)");
	addExampleLine("xmipp_nma_alignment -i images.sel --pdb 2tbv.pdb --modes modelist.xmd --trustradius_scale 1.2 --sampling_rate 3.2 -o output.xmd --resume");
}

// Read arguments ==========================================================
void ProgNmaAlignment::readParams() {
	XmippMetadataProgram::readParams();
	fnPDB = getParam("--pdb");
	fnOutDir = getParam("--odir");
	fnModeList = getParam("--modes");
	resume = checkParam("--resume");
	trustradius_scale = abs(getDoubleParam("--trustradius_scale"));
	sampling_rate = getDoubleParam("--sampling_rate");
	fnmask = getParam("--mask");
	gaussian_DFT_sigma = getDoubleParam("--gaussian_Fourier");
	gaussian_Real_sigma = getDoubleParam("--gaussian_Real");
	weight_zero_freq = getDoubleParam("--zerofreq_weight");
	do_centerPDB = checkParam("--centerPDB");
	do_FilterPDBVol = checkParam("--filterVol");
	if (do_FilterPDBVol)
		cutoff_LPfilter = getDoubleParam("--filterVol");
	useFixedGaussian = checkParam("--fixed_Gaussian");
	//if (useFixedGaussian)
	sigmaGaussian = getDoubleParam("--fixed_Gaussian");
	projMatch = checkParam("--projMatch");
	discrAngStep = getDoubleParam("--discrAngStep");
}

// Show ====================================================================
void ProgNmaAlignment::show() {
	XmippMetadataProgram::show();
	std::cout
            << "Output directory:     " << fnOutDir << std::endl
	        << "PDB:                  " << fnPDB << std::endl
			<< "Resume:               " << resume << std::endl
			<< "Mode list:            " << fnModeList << std::endl
			<< "Trust-region scale:   " << trustradius_scale << std::endl
			<< "Pixel size:           " << sampling_rate << std::endl
			<< "Mask:                 " << fnmask << std::endl
			<< "Center PDB:           " << do_centerPDB << std::endl
			<< "Filter PDB volume:    " << do_FilterPDBVol << std::endl
			<< "Use pseudo-atoms:     " << useFixedGaussian << std::endl
			<< "Pseudo-atom sigma:    " << sigmaGaussian << std::endl
			<< "Projection matching:  " << projMatch << std::endl
			<< "AngularSamplingStep:  " << discrAngStep << std::endl
	                << "Gaussian Fourier:     " << gaussian_DFT_sigma << std::endl
			<< "Gaussian Real:        " << gaussian_Real_sigma << std::endl
			<< "Zero-frequency weight:" << weight_zero_freq << std::endl;
}

// Produce side information ================================================
ProgNmaAlignment *global_nma_prog;

void ProgNmaAlignment::createWorkFiles() {
	MetaData *pmdIn = getInputMd();
	MetaData mdTodo, mdDone;
	mdTodo = *pmdIn;
	FileName fn(fnOutDir+"/nmaDone.xmd");
	if (fn.exists() && resume) {
		mdDone.read(fn);
		mdTodo.subtraction(mdDone, MDL_IMAGE);
	} else //if not exists create metadata only with headers
	{
		mdDone.addLabel(MDL_IMAGE);
		mdDone.addLabel(MDL_ENABLED);
		mdDone.addLabel(MDL_ANGLE_ROT);
		mdDone.addLabel(MDL_ANGLE_TILT);
		mdDone.addLabel(MDL_ANGLE_PSI);
		mdDone.addLabel(MDL_SHIFT_X);
		mdDone.addLabel(MDL_SHIFT_Y);
		mdDone.addLabel(MDL_NMA);
		mdDone.addLabel(MDL_COST);
		mdDone.write(fn);
	}
	*pmdIn = mdTodo;
}

void ProgNmaAlignment::preProcess() {
	MetaData SF(fnModeList);
	numberOfModes = SF.size();
	// Get the size of the images in the selfile
	imgSize = xdimOut;
	// Set the pointer of the program to this object
	global_nma_prog = this;
	//create some neededs files
	createWorkFiles();
}

void ProgNmaAlignment::finishProcessing() {
	XmippMetadataProgram::finishProcessing();
	rename((fnOutDir+"/nmaDone.xmd").c_str(), fn_out.c_str());
}

// Create deformed PDB =====================================================
FileName ProgNmaAlignment::createDeformedPDB(int pyramidLevel) const {
	String program;
	String arguments;
	FileName fnRandom;
	fnRandom.initUniqueName(nameTemplate,fnOutDir);
	const char * randStr = fnRandom.c_str();

	program = "xmipp_pdb_nma_deform";
	arguments = formatString(
			"--pdb %s -o %s_deformedPDB.pdb --nma %s --deformations ",
			fnPDB.c_str(), randStr, fnModeList.c_str());
	for (size_t i = 0; i < VEC_XSIZE(trial) - 5; ++i)
		arguments += floatToString(trial(i)) + " ";
	runSystem(program, arguments, false);

	program = "xmipp_volume_from_pdb";
	arguments = formatString(
			"-i %s_deformedPDB.pdb --size %i --sampling %f -v 0", randStr,
			imgSize, sampling_rate);

	if (do_centerPDB)
		arguments.append(" --centerPDB ");

	if (useFixedGaussian) {
		arguments.append(" --fixed_Gaussian ");
		if (sigmaGaussian >= 0)
			arguments += formatString("%f --intensityColumn Bfactor",
					sigmaGaussian);
	}
	progVolumeFromPDB->read(arguments);
	progVolumeFromPDB->tryRun();

	if (do_FilterPDBVol) {
		program = "xmipp_transform_filter";
		arguments = formatString(
						"-i %s_deformedPDB.vol --sampling %f --fourier low_pass %f  -v 0",
						randStr, sampling_rate, cutoff_LPfilter);
		runSystem(program, arguments, false);
	}

	if (pyramidLevel != 0) {
		Image<double> I;
		FileName fnDeformed = formatString("%s_deformedPDB.vol",randStr);
		I.read(fnDeformed);
		selfPyramidReduce(BSPLINE3, I(), pyramidLevel);
		I.write(fnDeformed);
	}

	return fnRandom;
}

// Perform complete search =================================================
void ProgNmaAlignment::performCompleteSearch(const FileName &fnRandom,
		int pyramidLevel) const {
	String program;
	String arguments;
	const char * randStr = fnRandom.c_str();

	// Reduce the image
	FileName fnDown = formatString("%s_downimg.xmp", fnRandom.c_str());
	if (pyramidLevel != 0) {
		Image<double> I;
		I.read(currentImgName);
		selfPyramidReduce(BSPLINE3, I(), pyramidLevel);
		I.write(fnDown);
	} else
		link(currentImgName.c_str(), fnDown.c_str());

	mkdir((fnRandom+"_ref").c_str(), S_IRWXU);

	double angSampling=2*RAD2DEG(atan(1.0/((double) imgSize / pow(2.0, (double) pyramidLevel+1))));
	angSampling=std::max(angSampling,discrAngStep);
	program = "xmipp_angular_project_library";
	arguments = formatString(
			"-i %s_deformedPDB.vol -o %s_ref/ref.stk --sampling_rate %f -v 0",
			randStr, randStr, angSampling);
	if (projMatch)
		arguments +=formatString(
						" --compute_neighbors --angular_distance -1 --experimental_images %s_downimg.xmp", randStr);

	runSystem(program, arguments, false);

	const char * refSelStr = formatString("%s_ref/ref.doc", randStr).c_str();
	const char * refStkStr = formatString("%s_ref/ref.stk", randStr).c_str();

	if (fnmask != "") {
		program = "xmipp_transform_mask";
		arguments = formatString("-i %s --mask binary_file %s", refSelStr,
				fnmask.c_str());
		runSystem(program, arguments, false);
	}

	// Perform alignment
	if (!projMatch) {
		program = "xmipp_angular_discrete_assign";
		arguments = formatString(
						"-i %s_downimg.xmp --ref %s -o %s_angledisc.xmd --psi_step 5 --max_shift_change %d --search5D -v 0",
						randStr, refSelStr, randStr, (int)round((double) imgSize / (10.0 * pow(2.0, (double) pyramidLevel))));
	} else {
		program = "xmipp_angular_projection_matching";
		arguments =	formatString(
				        "-i %s_downimg.xmp --ref %s -o %s_angledisc.xmd --search5d_step 1 --max_shift %d -v 0",
				        randStr, refStkStr, randStr, (int)round((double) imgSize / (10.0 * pow(2.0, (double) pyramidLevel))));
	}
	runSystem(program, arguments, false);
}

// Continuous assignment ===================================================
double ProgNmaAlignment::performContinuousAssignment(const FileName &fnRandom,
		int pyramidLevel) const {
	// Perform alignment
	const char * randStr = fnRandom.c_str();
	String fnResults=formatString("%s_anglecont.xmd", randStr);
	bool costSource=true;
	if (!projMatch) {
		String program = "xmipp_angular_continuous_assign";
		String arguments =
				formatString(
						"-i %s_angledisc.xmd --ref %s_deformedPDB.vol -o %s --gaussian_Fourier %f --gaussian_Real %f --zerofreq_weight %f -v 0",
						randStr, randStr, fnResults.c_str(), gaussian_DFT_sigma,
						gaussian_Real_sigma, weight_zero_freq);
		runSystem(program, arguments, false);
	}
	else
	{
		costSource=false;
		if (currentStage==1)
			fnResults=formatString("%s_angledisc.xmd", randStr);
		else
		{
			mkdir((fnRandom+"_ref").c_str(), S_IRWXU);

			String program = "xmipp_angular_project_library";
			double angSampling=RAD2DEG(atan(1.0/((double) imgSize / pow(2.0, (double) pyramidLevel+1))));
			angSampling=std::max(angSampling,discrAngStep);
			String arguments = formatString(
					"-i %s_deformedPDB.vol -o %s_ref/ref.stk --sampling_rate %f -v 0",
					randStr, randStr,angSampling);
			arguments +=formatString(
							" --compute_neighbors --angular_distance %f --experimental_images %s_angledisc.xmd",
							2.5*angSampling,randStr);
			runSystem(program, arguments, false);

			const char * refStkStr = formatString("%s_ref/ref.stk", randStr).c_str();
			program = "xmipp_angular_projection_matching";
			arguments =	formatString(
							"-i %s_downimg.xmp --ref %s -o %s --search5d_step 1 --max_shift %d -v 0",
							randStr, refStkStr, fnResults.c_str(), (int)round((double) imgSize / (10.0 * pow(2.0, (double) pyramidLevel))));
			runSystem(program, arguments, false);
		}
	}

	// Pick up results
	MetaData DF(fnResults);
	MDRow row;
	DF.getRow(row, DF.firstObject());
	row.getValue(MDL_ANGLE_ROT, trial(VEC_XSIZE(trial) - 5));
	row.getValue(MDL_ANGLE_TILT, trial(VEC_XSIZE(trial) - 4));
	row.getValue(MDL_ANGLE_PSI, trial(VEC_XSIZE(trial) - 3));
	row.getValue(MDL_SHIFT_X, trial(VEC_XSIZE(trial) - 2));
	trial(VEC_XSIZE(trial) - 2) *= pow(2.0, (double) pyramidLevel);
	row.getValue(MDL_SHIFT_Y, trial(VEC_XSIZE(trial) - 1));
	trial(VEC_XSIZE(trial) - 1) *= pow(2.0, (double) pyramidLevel);
	double tempvar;
	if (!costSource) {
		row.getValue(MDL_MAXCC, tempvar);
		tempvar = -tempvar;
	} else
		row.getValue(MDL_COST, tempvar);
	return tempvar;
}

void ProgNmaAlignment::updateBestFit(double fitness, int dim) {
	if (fitness < fitness_min(0)) {
		fitness_min(0) = fitness;
		trial_best = trial;
	}
}

// Compute fitness =========================================================
double ObjFunc_nma_alignment::eval(Vector X, int *nerror) {
	int dim = global_nma_prog->numberOfModes;

	for (int i = 0; i < dim; i++) {
		global_nma_prog->trial(i) = X[i];
	}

	int pyramidLevelDisc = 1;
	int pyramidLevelCont = (global_nma_prog->currentStage == 1) ? 1 : 0;

	FileName fnRandom = global_nma_prog->createDeformedPDB(pyramidLevelCont);
	const char * randStr = fnRandom.c_str();

	if (global_nma_prog->currentStage == 1) {
		global_nma_prog->performCompleteSearch(fnRandom, pyramidLevelDisc);
	} else {
		double rot, tilt, psi, xshift, yshift;
		MetaData DF;

		rot = global_nma_prog->bestStage1(
				VEC_XSIZE(global_nma_prog->bestStage1) - 5);
		tilt = global_nma_prog->bestStage1(
				VEC_XSIZE(global_nma_prog->bestStage1) - 4);
		psi = global_nma_prog->bestStage1(
				VEC_XSIZE(global_nma_prog->bestStage1) - 3);
		xshift = global_nma_prog->bestStage1(
				VEC_XSIZE(global_nma_prog->bestStage1) - 2);
		yshift = global_nma_prog->bestStage1(
				VEC_XSIZE(global_nma_prog->bestStage1) - 1);

		size_t objId = DF.addObject();
		FileName fnDown = formatString("%s_downimg.xmp", randStr);
		DF.setValue(MDL_IMAGE, fnDown, objId);
		DF.setValue(MDL_ENABLED, 1, objId);
		DF.setValue(MDL_ANGLE_ROT, rot, objId);
		DF.setValue(MDL_ANGLE_TILT, tilt, objId);
		DF.setValue(MDL_ANGLE_PSI, psi, objId);
		DF.setValue(MDL_SHIFT_X, xshift, objId);
		DF.setValue(MDL_SHIFT_Y, yshift, objId);

		DF.write(formatString("%s_angledisc.xmd", randStr));
		link(global_nma_prog->currentImgName.c_str(), fnDown.c_str());
	}
	double fitness = global_nma_prog->performContinuousAssignment(fnRandom,
			pyramidLevelCont);

	runSystem("rm", formatString("-rf %s* &", randStr));

	global_nma_prog->updateBestFit(fitness, dim);

	return fitness;
}

ObjFunc_nma_alignment::ObjFunc_nma_alignment(int _t, int _n) {
}

void ProgNmaAlignment::processImage(const FileName &fnImg,
		const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut) {
	static size_t imageCounter = 0;
	++imageCounter;

	double rhoStart=trustradius_scale*250.;
        double rhoEnd=trustradius_scale*50.;

	int niter=10000;

	ObjectiveFunction *of;

	int dim = numberOfModes;

	parameters.initZeros(dim + 5);
	currentImgName = fnImg;
	sprintf(nameTemplate, "_node%d_img%lu_XXXXXX", rangen, (long unsigned int)imageCounter);

	trial.initZeros(dim + 5);
	trial_best.initZeros(dim + 5);

	fitness_min.initZeros(1);
	fitness_min(0) = 1000000.0;

	currentStage = 1;
#ifdef DEBUG
	std::cerr << std::endl << "DEBUG: ===== Node: " << rangen
	<<" processing image " << fnImg <<"(" << objId << ")"
	<< " at stage: " << currentStage << std::endl;
#endif
	of = new ObjFunc_nma_alignment(1, dim);

	of->xStart.setSize(dim);
	for (int i = 0; i < dim; i++)
		of->xStart[i] = 0.;

#ifdef DEBUG
	strcpy(of->name,("OF1_"+integerToString(rangen)).c_str());
	of->setSaveFile();
#endif

	CONDOR(rhoStart, rhoEnd, niter, of);
#ifdef DEBUG
	of->printStats();
	FILE *ff = fopen(("res1_"+integerToString(rangen)+".xmd").c_str(),"w");
	fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
	fclose(ff);
#endif

	bestStage1 = trial = parameters = trial_best;

	delete of;

	currentStage = 2;
#ifdef DEBUG
	std::cerr << std::endl << "DEBUG: ===== Node: " << rangen
	<<" processing image " << fnImg <<"(" << objId << ")"
	<< " at stage: " << currentStage << std::endl;
#endif

	fitness_min(0) = 1000000.0;

	of = new ObjFunc_nma_alignment(1, dim);

	of->xStart.setSize(dim);
	for (int i = 0; i < dim; i++)
		of->xStart[i] = parameters(i);
#ifdef DEBUG
	strcpy(of->name,("OF2_"+integerToString(rangen)).c_str());
	of->setSaveFile();
#endif

	rhoStart=trustradius_scale*50.-1., rhoEnd=0.5;
	CONDOR(rhoStart, rhoEnd, niter, of);
#ifdef DEBUG
	of->printStats();
	ff=fopen(("res2_"+integerToString(rangen)+".xmd").c_str(),"w");
	fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
	fclose(ff);
#endif

#ifdef DEBUG
	std::cout << "Best fitness = " << fitness << std::endl;
	for (int i=0; i<dim; i++)
	{
		std::cout << "Best deformations = " << dd[i] << std::endl;
	}
#endif

	trial = trial_best;

	for (int i = dim; i < dim + 5; i++) {
		parameters(i - dim) = trial_best(i);
	}

	for (int i = 0; i < dim; i++) {
		parameters(5 + i) = trial_best(i);
	}

	parameters.resize(VEC_XSIZE(parameters) + 1);
	parameters(VEC_XSIZE(parameters) - 1) = fitness_min(0);

	writeImageParameters(fnImg);
	delete of;
}

void ProgNmaAlignment::writeImageParameters(const FileName &fnImg) {
	MetaData md;
	size_t objId = md.addObject();
	md.setValue(MDL_IMAGE, fnImg, objId);
	md.setValue(MDL_ENABLED, 1, objId);
	md.setValue(MDL_ANGLE_ROT, parameters(0), objId);
	md.setValue(MDL_ANGLE_TILT, parameters(1), objId);
	md.setValue(MDL_ANGLE_PSI, parameters(2), objId);
	md.setValue(MDL_SHIFT_X, parameters(3), objId);
	md.setValue(MDL_SHIFT_Y, parameters(4), objId);

	int dim = numberOfModes;
	std::vector<double> vectortemp;

	for (int j = 5; j < 5 + dim; j++) {
		vectortemp.push_back(parameters(j));
	}

	md.setValue(MDL_NMA, vectortemp, objId);
	md.setValue(MDL_COST, parameters(5 + dim), objId);

	md.append(fnOutDir+"/nmaDone.xmd");
}
