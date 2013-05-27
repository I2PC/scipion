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

#include "flexible_alignment.h"
#include "data/metadata_extension.h"
#include "program_extension.h"
#include "libraries/reconstruction/pdb_nma_deform.h"

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include "../../external/bilib/headers/messagedisplay.h"
#include "../../external/bilib/headers/error.h"
#include "../../external/bilib/configs.h"

// Empty constructor =======================================================
ProgFlexibleAlignment::ProgFlexibleAlignment()
{
    rangen = 0;
    resume = false;
    currentImgName = "";
    each_image_produces_an_output = false;
    produces_an_output = true;
}

ProgFlexibleAlignment::~ProgFlexibleAlignment()
{
    //delete progVolumeFromPDB;
}

// Params definition ============================================================
void ProgFlexibleAlignment::defineParams()
{
    addUsageLine("Compute deformation parameters according to a set of NMA modes");
    defaultComments["-o"].clear();
    defaultComments["-o"].addComment("Metadata with output alignment and deformations");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --pdb <PDB_filename>                : PDB Model to compute NMA");
    addParamsLine("  [--odir <outputDir=\".\">]           : Output directory");
    addParamsLine("  [--resume]                           : Resume processing");
    addParamsLine("==Generation of the deformed volumes==");
    addParamsLine("   --modes <filename>                  : File with a list of mode filenames");
    addParamsLine("  [--defampsampling <s=200>]           : Deformation sampling");
    addParamsLine("  [--maxdefamp <s=2000>]               : Maximum deformation amplitude");
    addParamsLine("  [--translsampling <s=2>]             : Translational sampling");
    addParamsLine("  [--maxtransl <max=10>]               : Maximum translation");
    addParamsLine("  [--sampling_rate <Ts=1>]             : in Angstroms/pixel");
    addParamsLine("  [--filterVol <cutoff=15.>]           : Filter the volume after deforming. Default cut-off is 15 A.");
    addParamsLine("  [--centerPDB]                        : Center the PDB structure");
    addParamsLine("  [--fixed_Gaussian <std=-1>]          : For pseudo atoms fixed_Gaussian must be used.");
    addParamsLine("                                       : Default standard deviation <std> is read from PDB file.");
    addParamsLine("==Angular assignment and mode detection==");
    addParamsLine("  [--mask <m=\"\">]                    : 2D Mask applied to the reference images of the deformed volume");
    addParamsLine("                                       :+Note that wavelet assignment needs the input images to be of a size power of 2");
    addParamsLine("  [--minAngularSampling <ang=3>]       : Minimum angular sampling rate");
    addParamsLine("  [--gaussian_Real    <s=0.5>]         : Weighting sigma in Real space");
    addParamsLine("  [--zerofreq_weight  <s=0.>]          : Zero-frequency weight");
    addParamsLine("  [--sigma    <s=10>]                  : Sigma");
    addParamsLine("  [--max_iter  <N=60>]                 : Maximum number of iterations");
    addExampleLine("xmipp_nma_alignment -i images.sel --pdb 2tbv.pdb --modes modelist.xmd --sampling_rate 6.4 -o output.xmd --resume");
}

// Read arguments ==========================================================
void ProgFlexibleAlignment::readParams()
{
    XmippMetadataProgram::readParams();
    fnPDB = getParam("--pdb");
    fnOutDir = getParam("--odir");
    fnModeList = getParam("--modes");
    resume = checkParam("--resume");
    maxdefamp = getDoubleParam("--maxdefamp");
    defampsampling = getDoubleParam("--defampsampling");
    translsampling = getDoubleParam("--translsampling");
    maxtransl = getDoubleParam("--maxtransl");
    sampling_rate = getDoubleParam("--sampling_rate");
    fnmask = getParam("--mask");
    gaussian_Real_sigma = getDoubleParam("--gaussian_Real");
    weight_zero_freq = getDoubleParam("--zerofreq_weight");
    do_centerPDB = checkParam("--centerPDB");
    do_FilterPDBVol = checkParam("--filterVol");
    if (do_FilterPDBVol)
        cutoff_LPfilter = getDoubleParam("--filterVol");
    useFixedGaussian = checkParam("--fixed_Gaussian");
    //if (useFixedGaussian)
    sigmaGaussian = getDoubleParam("--fixed_Gaussian");
    minAngularSampling = getDoubleParam("--minAngularSampling");
    sigma = getDoubleParam("--sigma");
    max_no_iter = getIntParam("--max_iter");
}

// Show ====================================================================
void ProgFlexibleAlignment::show()
{
    XmippMetadataProgram::show();
    std::cout
    << "Output directory:     " << fnOutDir << std::endl
    << "PDB:                  " << fnPDB << std::endl
    << "Resume:               " << resume << std::endl
    << "Mode list:            " << fnModeList << std::endl
    << "Deformation sampling: " << defampsampling << std::endl
    << "Maximum amplitude:    " << maxdefamp << std::endl
    << "Transl. sampling:     " << translsampling << std::endl
    << "Max. Translation:     " << maxtransl << std::endl
    << "Sampling rate:        " << sampling_rate << std::endl
    << "Mask:                 " << fnmask << std::endl
    << "Center PDB:           " << do_centerPDB << std::endl
    << "Filter PDB volume     " << do_FilterPDBVol << std::endl
    << "Use fixed Gaussian:   " << useFixedGaussian << std::endl
    << "Sigma of Gaussian:    " << sigmaGaussian << std::endl
    << "minAngularSampling:   " << minAngularSampling << std::endl
    << "Gaussian Real:        " << gaussian_Real_sigma << std::endl
    << "Zero-frequency weight:" << weight_zero_freq << std::endl
    << "Sigma:                " << sigma << std::endl
    << "Max. Iter:            " << max_no_iter << std::endl
    ;
}

// Produce side information ================================================
ProgFlexibleAlignment *global_flexible_prog;

void ProgFlexibleAlignment::createWorkFiles()
{
    MetaData *pmdIn = getInputMd();
    MetaData mdTodo, mdDone;
    mdTodo = *pmdIn;
    FileName fn(fnOutDir+"/nmaDone.xmd");
    if (fn.exists() && resume)
    {
        mdDone.read(fn);
        mdTodo.subtraction(mdDone, MDL_IMAGE);
    }
    else //if not exists create metadata only with headers
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

void ProgFlexibleAlignment::preProcess()
{
    MetaData SF(fnModeList);
    numberOfModes = SF.size();
    SF.getColumnValues(MDL_NMA_MODEFILE,modeList);

    // Get the size of the images in the selfile
    imgSize = xdimOut;
    // Set the pointer of the program to this object
    global_flexible_prog = this;
    //create some neededs files
    createWorkFiles();
}

void ProgFlexibleAlignment::finishProcessing()
{
    XmippMetadataProgram::finishProcessing();
    rename((fnOutDir+"/nmaDone.xmd").c_str(), fn_out.c_str());
}

// Create deformed PDB =====================================================
FileName ProgFlexibleAlignment::createDeformedPDB()
{
    String program;
    String arguments;
    FileName fnRandom;
    fnRandom.initUniqueName(nameTemplate,fnOutDir);
    const char * randStr = fnRandom.c_str();
    std::cout << "210" << randStr << std::endl;
    program = "xmipp_pdb_nma_deform";
    arguments = formatString(
                    "--pdb %s -o %s_deformedPDB.pdb --nma %s --deformations ",
                    fnPDB.c_str(), randStr, fnModeList.c_str());
    for (size_t i = 5; i < VEC_XSIZE( trial ); i++)
        arguments += floatToString(trial(i)) + " ";
    runSystem(program, arguments, false);


    return fnRandom;
}

void  ProjectionRefencePoint(Matrix1D<double>  &Parameters,
                             int      dim,
                             double   *R,
                             double   *Tr,
                             MultidimArray<double>  &proj_help_test,
                             MultidimArray<double>  &P_esp_image,
                             //double   *S_mu,
                             int      Xwidth,
                             int      Ywidth,
                             double   sigma)
{
    double *coord_gaussian;
    double *ksi_v,*coord_img;
    double *ro_ksi_v, *ro_coord_img, *ro_coord_gaussian;
    //double S_mu;
    int    psi_max = (int)(sqrt(3)*Ywidth/2);
    int    kx,ky;
    double    a0,a1,a2;
    double centre_Xwidth, centre_Ywidth;
    double sum2,hlp;

    proj_help_test.initZeros(Xwidth,Ywidth);
    ksi_v    = (double*)malloc( (size_t) 4L * sizeof(double));//01
    for (int i=0; i<dim+5; i++)
    {
        global_flexible_prog->trial(i) = Parameters(i);
        //std::cout << "344Parameters("<<i<<") = " <<Parameters(i) << std::endl;
    }
    std::string command;
    String program;
    String arguments;
    FileName fnRandom;
    std::cout << 310<<Parameters<< std::endl;
    fnRandom.initUniqueName(global_flexible_prog->nameTemplate,global_flexible_prog->fnOutDir);
    const char * randStr = fnRandom.c_str();
    //std::cout << 312 << randStr <<std::endl;

    program = "xmipp_pdb_nma_deform";
    arguments = formatString(
                    "--pdb %s -o %s_deformedPDB.pdb --nma %s --deformations ",
                    global_flexible_prog->fnPDB.c_str(), randStr, global_flexible_prog->fnModeList.c_str());
    for (size_t i = 5; i < VEC_XSIZE(global_flexible_prog->trial) ; i++){
        float aaa=global_flexible_prog->scdefamp*Parameters(i);
        arguments += floatToString(aaa) + " ";
    }
    runSystem(program, arguments, false);
    std::cout << "304 arguments " << arguments << std::endl;


    String  deformed_pdb = formatString("%s_deformedPDB.pdb",randStr);
    centre_Xwidth = double(Xwidth-1)/2.0;
    centre_Ywidth = double(Ywidth-1)/2.0;
    //std::cout << "377  centre_Xwidth = " <<centre_Xwidth  << std::endl;
    Matrix1D<double> limit0(3), limitF(3), centerOfMass(3);
    const char *intensityColumn = " ";
    computePDBgeometry(global_flexible_prog->fnPDB, centerOfMass, limit0, limitF, intensityColumn);
    centerOfMass = (limit0 + limitF)/2;
    //std::cout << "c383enterOfMass = " <<centerOfMass  << std::endl;
    std::ifstream fh_deformedPDB;
    fh_deformedPDB.open(deformed_pdb.c_str());

    if (!fh_deformedPDB)
        REPORT_ERROR(ERR_UNCLASSIFIED, (std::string)"Prog_PDBPhantom_Parameters::protein_geometry:"
                     "Cannot open " + deformed_pdb + " for reading");

    // Process all lines of the filem+".pdb").c_str());
    if (!fh_deformedPDB)
        REPORT_ERROR(ERR_UNCLASSIFIED, (std::string)"Prog_PDBPhantom_Parameters::protein_geometry:"
                     "Cannot open " + deformed_pdb + " for reading");

    //ksi_v    = (double*)malloc( (size_t) 4L * sizeof(double));//02
    ro_ksi_v    = (double*)malloc( (size_t) 4L * sizeof(double));//03

    ksi_v[0] = 0.0;
    ksi_v[1] = 0.0;
    ksi_v[3] = 1.0;

    coord_gaussian = (double*)malloc( (size_t) 4L * sizeof(double));//04
    ro_coord_gaussian = (double*)malloc( (size_t) 4L * sizeof(double));//05

    coord_gaussian[3] = 1.0;

    coord_img = (double*)malloc( (size_t) 4L * sizeof(double));//06
    ro_coord_img = (double*)malloc( (size_t) 4L * sizeof(double));//07

    coord_img[2] = 0.0;
    coord_img[3] = 1.0;
    double x,y,z;
    std::string kind;
    std::string line;
    std::string atom_type;
    int ttt =0;
    std::cout << "351 Parameters " << Parameters <<std::endl;

    while (!fh_deformedPDB.eof())
    {
        ttt++;
        // Read an ATOM line
        getline(fh_deformedPDB, line);
        if (line == "")
            continue;
        kind = line.substr(0,4);
        if (kind != "ATOM" && kind !="HETA")
            continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        atom_type = line.substr(13,2);
        x = textToFloat(line.substr(30,8));
        y = textToFloat(line.substr(38,8));
        z = textToFloat(line.substr(46,8));
        coord_gaussian[0] = ( x - centerOfMass(0) ) / (global_flexible_prog->sampling_rate);
        coord_gaussian[1] = ( y - centerOfMass(1) ) / (global_flexible_prog->sampling_rate);
        coord_gaussian[2] = ( z - centerOfMass(2) ) / (global_flexible_prog->sampling_rate);

        for(int ksi=-psi_max;ksi<=psi_max;ksi++)
        {
            ksi_v[2]=(double)ksi;
            MatrixMultiply( R, ksi_v, ro_ksi_v,4L, 4L, 1L);

            for(ky=0;ky<Ywidth;ky++)
            {
                coord_img[1] = (double)ky - centre_Ywidth;
                for(kx=0;kx<Xwidth;kx++)
                {
                    coord_img[0] = (double)kx - centre_Xwidth;
                    MatrixMultiply( Tr, coord_img, ro_coord_img,4L, 4L, 1L);

                    a0 = ro_coord_img[0]+ro_ksi_v[0]+coord_gaussian[0];
                    a1 = ro_coord_img[1]+ro_ksi_v[1]+coord_gaussian[1];
                    a2 = ro_coord_img[2]+ro_ksi_v[2]+coord_gaussian[2];
                    proj_help_test(kx-Xwidth/2,ky-Ywidth/2) += (double)exp( - (a0*a0+a1*a1+a2*a2)/(sigma*sigma) );

                }

            }

        }

    }//std::cout << "465 TTT=  "<<ttt << std::endl;
    std::cout << "400 Parameters " << Parameters <<std::endl;
    Image<double> test1;
    test1()=proj_help_test;

    // Close file
    fh_deformedPDB.close();
    /*command=(std::string)"rm -f deformedPDB_"+ fnRandom+".pdb";
    system(command.c_str());*/
    std::cout << "409 Parameters " << Parameters <<std::endl;
    // To calculate the value of cost function
    for(sum2 = 0.0, kx=0;kx<Xwidth;kx++)
    {
        for(ky=0;ky<Ywidth;ky++)
        {
            hlp   = (double)(proj_help_test(kx-Xwidth/2,ky-Ywidth/2)-P_esp_image(kx-Xwidth/2,ky-Ywidth/2));
            hlp   = (double)hlp*hlp;
            sum2 += (double)hlp;
        }
    }
    std::cout << "420 Parameters " << Parameters <<std::endl;
    global_flexible_prog->costfunctionvalue = (double)sum2;
    /* Calcule the real projection*/
    std::cout << "503 *costfunctionvalue=" << global_flexible_prog->costfunctionvalue<< std::endl;
    std::cout << "424 Parameters " << Parameters <<std::endl;

    free(ksi_v);               std::cout << "486" << std::endl;
    free(ro_ksi_v);            std::cout << "486" << std::endl;
    free(coord_gaussian);      std::cout << "487" << std::endl;
    free(ro_coord_gaussian);   std::cout << "487" << std::endl;
    free(coord_img);           std::cout << "488" << std::endl;
    free(ro_coord_img);        std::cout << "488" << std::endl;


    //return (!ERROR);
}

/* ------------------------------------------------------------------------- */
/* Calculate the values of partial valur of P_mu_image                       */
/* ------------------------------------------------------------------------- */
int partialpfunction(Matrix1D<double>  &Parameters,
                     Matrix1D<double>  &centerOfMass,
                     double            *R,
                     double            *Tr,
                     double            *DR0,
                     double            *DR1,
                     double            *DR2,
                     MultidimArray<double>  &DP_Rz1,
                     MultidimArray<double>  &DP_Ry,
                     MultidimArray<double>  &DP_Rz2,
                     MultidimArray<double>  &DP_x,
                     MultidimArray<double>  &DP_y,
                     MultidimArray<double>  &DP_q,
                     //double            *cost,
                     MultidimArray<double>  &P_mu_image,
                     MultidimArray<double>  &P_esp_image,
                     int               Xwidth,
                     int               Ywidth)
{
    int     psi_max = (int)(sqrt(3)*128/(global_flexible_prog->sampling_rate));
    double  help, a0,a1,a2;
    double  *help_v,*coord_gaussian,*coord_img;
    double  *ModeValues;
    int     Line_number = 0;
    int     kx,ky;
    int    dim = global_flexible_prog-> numberOfModes;
    double  centre_Xwidth, centre_Ywidth;
    centre_Xwidth = (double)(Xwidth - 1)/2.0;
    centre_Ywidth = (double)(Ywidth - 1)/2.0;

    DP_Rz1.initZeros(Xwidth,Ywidth);
    DP_Ry. initZeros(Xwidth,Ywidth);
    DP_Rz2.initZeros(Xwidth,Ywidth);
    DP_x.  initZeros(Xwidth,Ywidth);
    DP_y.  initZeros(Xwidth,Ywidth);
    DP_q.  initZeros(dim,Xwidth,Ywidth);

    std::ifstream ModeFile;
    //std::string modefilename = modeList[0];
    std::string modefilename = global_flexible_prog->modeList[0];
    ModeFile.open(modefilename.c_str());
    if(ModeFile.fail())
        REPORT_ERROR(ERR_UNCLASSIFIED, (std::string) modefilename + " for reading");
    std::string line;
    while(getline(ModeFile,line))
    {
        Line_number++;
    }
    ModeFile.close();
    ModeValues = (double*)malloc( (size_t) 3*Line_number*dim  * sizeof(double));//08
    //std::string line;
    std::string x,y,z;



    for (int i=0; i<dim;i++)
    {
        modefilename = global_flexible_prog->modeList[i];
        ModeFile.open(modefilename.c_str());
        int n =0;
        while(getline(ModeFile,line))
        {
            x = line.substr( 3,10);
            y = line.substr(14,10);
            z = line.substr(27,10);
            ModeValues[i*3*Line_number + n*3 + 0] = atof( x.c_str());
            ModeValues[i*3*Line_number + n*3 + 1] = atof( y.c_str());
            ModeValues[i*3*Line_number + n*3 + 2] = atof( z.c_str());
            n++;
        }
        ModeFile.close();
    }
    //clear line;
    
    for (int i=0; i<dim+5; i++)
    {
        global_flexible_prog->trial(i) = Parameters(i);
        std::cout << "partialpfunctionParameters(" << i << ")" << Parameters(i) << std::endl;
    }

    FileName fnRandom = global_flexible_prog->createDeformedPDB();
    std::ifstream fh_deformedPDB;
    fh_deformedPDB.open((fnRandom+ "_deformedPDB.pdb").c_str());
    if (!fh_deformedPDB)
        REPORT_ERROR(ERR_UNCLASSIFIED, (std::string)"Prog_PDBPhantom_Parameters::protein_geometry:" "Cannot open " + fnRandom+"_deformedPDB.pdb" + " for reading");

    // Process all lines of the file
    help_v    = ( double* )malloc( (size_t) 4L * sizeof(double));
    /*if (help_v == (double *)NULL)
{
        WRITE_ERROR(partialpfunction, "ERROR - Not enough memory for help_v");
        return(ERROR);
}*/
    help_v[0] = 0.0;
    help_v[1] = 0.0;
    help_v[3] = 1.0;
    std::cout << " 586help_v[3] = " << help_v[3]  << std::endl;
    coord_gaussian     = (double*)malloc( (size_t) 4L * sizeof(double));
    /*if (coord_gaussian == (double *)NULL)
{
        WRITE_ERROR(partialpfunction, "ERROR - Not enough memory for help_v");
        free(help_v);
        return(ERROR);
}*/
    coord_gaussian[3] = 1.0;

    coord_img = (double*)malloc( (size_t) 4L * sizeof(double));
    /*if (coord_img == (double *)NULL)
{
        WRITE_ERROR(partialpfunction, "ERROR - Not enough memory for help_v");
        free(help_v);
        free(coord_gaussian);
        return(ERROR);
}*/
    coord_img[2] = 0.0;
    coord_img[3] = 1.0;
    std::cout << " 606coord_img[3] = " << coord_img[3] << std::endl;
    int k=0;
    //std::string line;
    std::string kind;
    while (!fh_deformedPDB.eof())
    {
        // Read an ATOM line
        //getline(fh_deformedPDB,line);
        getline(fh_deformedPDB,line);
        if (line == "")
            continue;
        kind = line.substr(0,4);
        if (kind != "ATOM" )
            continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        std::string atom_type = line.substr(13,2);
        double x = textToFloat(line.substr(30,8));
        double y = textToFloat(line.substr(38,8));
        double z = textToFloat(line.substr(46,8));

        // Correct position
        coord_gaussian[0] = ( x - centerOfMass(0) ) / global_flexible_prog->sampling_rate;
        coord_gaussian[1] = ( y - centerOfMass(1) ) / global_flexible_prog->sampling_rate;
        coord_gaussian[2] = ( z - centerOfMass(2) ) / global_flexible_prog->sampling_rate;

        //MatrixMultiply( Tr, coord_gaussian, coord_gaussian,4L, 4L, 1L);
        if (MatrixMultiply( Tr, coord_gaussian, coord_gaussian,4L, 4L, 1L) == ERROR)
        {
            WRITE_ERROR(partialpfunction, "Error returned by MatrixMultiply");
            free(help_v);
            return(ERROR);
        }


        for(int ksi=-psi_max;ksi<=psi_max;ksi++)
        {
            coord_img[3] = (double)ksi;
            for(kx=0;kx<Xwidth;kx++)
            {
                for(ky=0;ky<Ywidth;ky++)
                {

                    coord_img[0] = (double)kx - centre_Xwidth;
                    coord_img[1] = (double)ky - centre_Ywidth;

                    //MatrixTimesVector( Tr, coord_img, help_v,4L, 4L);
                    if (MatrixTimesVector( Tr, coord_img, help_v,4L, 4L) == ERROR)
                    {
                        WRITE_ERROR(partialpfunction, "Error returned by MatrixMultiply");
                        free(help_v);
                        free(coord_gaussian);
                        free(coord_img);
                        return(ERROR);
                    }
                    a0 = help_v[0]+coord_gaussian[0];//std::cout << " 691 a0 = " << a0 ;
                    a1 = help_v[1]+coord_gaussian[1];//std::cout << " 692 a1 = " << a1 ;
                    a2 = help_v[2]+coord_gaussian[2];//std::cout << " 693 a2 = " << a2 ;

                    help =(double)exp( -(a0*a0+a1*a1+a2*a2)/(global_flexible_prog->sigma * global_flexible_prog->sigma) );                    //std::cout << " 693 help = " << help ;
                    //std::cout << " 691 sigma = " << (global_flexible_prog->sigma * global_flexible_prog->sigma) <<std::endl ;
                    //help              += P_mu_image(kx-Xwidth,ky-Ywidth);
                    //P_mu_image(kx-Xwidth,ky-Ywidth) =  help;

                    //MatrixTimesVector( DR0, coord_img, help_v,4L, 4L);
                    if (MatrixTimesVector( DR0, coord_img, help_v,4L, 4L) == ERROR)
                    {
                        WRITE_ERROR(partialpfunction, "Error returned by MatrixMultiply");
                        free(help_v);
                        free(coord_gaussian);
                        free(coord_img);
                        return(ERROR);
                    }

                    DP_Rz1(kx,ky) += help * (a0*help_v[0]+a1*help_v[1]+a2*help_v[2]);
                    if (MatrixTimesVector( DR1, coord_img, help_v,4L, 4L) == ERROR)
                    {
                        WRITE_ERROR(partialpfunction, "Error returned by MatrixMultiply");
                        free(help_v);
                        free(coord_gaussian);
                        free(coord_img);
                        return(ERROR);
                    }
                    DP_Ry(kx,ky) += help * (a0*help_v[0]+a1*help_v[1]+a2*help_v[2]);
                    if (MatrixTimesVector( DR2, coord_img, help_v,4L, 4L) == ERROR)
                    {
                        WRITE_ERROR(partialpfunction, "Error returned by MatrixMultiply");
                        free(help_v);
                        free(coord_gaussian);
                        free(coord_img);
                        return(ERROR);
                    }
                    DP_Rz2(kx,ky) += help * (a0*help_v[0]+a1*help_v[1]+a2*help_v[2]);
                    DP_x(kx,ky)   += help * (a0*R[0]+a1*R[1]+a2*R[2]); //global_flexible_prog->sctrans * DP_x(kx,ky)
                    DP_y(kx,ky)   += help * (a0*R[4]+a1*R[5]+a2*R[6]); //global_flexible_prog->sctrans * DP_y(kx,ky)
                    for(int i = 0; i < dim; i++)
                    {
                        DP_q(i,kx,ky) += global_flexible_prog->scdefamp * help * (a0*ModeValues[i*3*Line_number +  k*3] + a1*ModeValues[i*3*Line_number + k*3 + 1] + a2*ModeValues[i*3*Line_number + k*3 + 2]);
                    }
                }
            }
        }
        k++;
    }
    std::cout << " 714 help = " << help << std::endl;
    fh_deformedPDB.close();
    /*std::string command=(std::string)"rm -f deformedPDB_"+ fnRandom+".pdb";
    system(command.c_str());
    command=(std::string)"rm -f downimg_"+fnRandom+".xmp";
    system(command.c_str());*/
    return(!ERROR);
}/*end of partialpfunction*/


/* ------------------------------------------------------------------------- */
/* Gradient and Hessian at pixel                                             */
/* ------------------------------------------------------------------------- */
void gradhesscost_atpixel(
    double *Gradient,
    double *Hessian,
    double *helpgr,
    double difference)
{
    int      trialSize = VEC_XSIZE(global_flexible_prog->trial);

    for (int i=0;i<trialSize;i++)
    {
        Gradient[i] += difference * helpgr[i];
        for(int j=0;j<=i;j++)
        {
            Hessian[i*trialSize + j] += helpgr[j]*helpgr[i];
        }
    }
    //return(!ERROR);
}/* End of gradhesscost_atpixel */

/* ------------------------------------------------------------------------- */
/* Calculate the values of Gradient and Hessian                              */
/* ------------------------------------------------------------------------- */
// int Prog_flexali_gauss_prm::return_gradhesscost(
int return_gradhesscost(
              Matrix1D<double>  &centerOfMass,
              double *Gradient,
              double *Hessian,
              Matrix1D<double>  &Parameters,
              int     dim,
              MultidimArray<double>  &rg_projimage,
              MultidimArray<double>  &P_esp_image,
              int    Xwidth,
              int    Ywidth)
          {
              double  phi,theta,psi,x0,y0;
              double  *Rz1,*Ry,*Rz2,*R,*DRz1,*DRy,*DRz2,*DR0,*DR1,*DR2,*Tr;
              double  *hlp,*helpgr;
              int     i,j;
              double  difference;
              double  SinPhi,CosPhi,SinPsi,CosPsi,SinTheta,CosTheta;
              int     trialSize = VEC_XSIZE(global_flexible_prog->trial);
              double  lambda = 1000.0;
              double  sigmalocal = global_flexible_prog->sigma;
              int     half_Xwidth = Xwidth/2;
              int     half_Ywidth = Ywidth/2;

              //global_flexible_prog->costfunctionvalue = 0.0;


              helpgr = (double *)malloc((size_t)(dim +5) * sizeof(double));
              for (int i = 0; i < dim + 5; i++)
              {
                  global_flexible_prog->trial(i) = Parameters(i);
              }
              phi   = Parameters(0);
              theta = Parameters(1);
              psi   = Parameters(2);
              x0    = Parameters(3);
              y0    = Parameters(4);
              std::cout << "813 phi   = "   << phi   << std::endl;
              std::cout << "814 theta = "   << theta << std::endl;
              std::cout << "815 psi   = "   << psi   << std::endl;
              std::cout << "816 x0    = "   << x0    << std::endl;
              std::cout << "817 y0    = "   << y0    << std::endl;
              //defamp = (double *)malloc((size_t) dim * sizeof(double));
              //if (defamp == (double *)NULL)
              //{
              //    WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for defamp");
              //    return(ERROR);
              //}
              SinPhi   = sin(phi);
              CosPhi   = cos(phi);
              SinPsi   = sin(psi);
              CosPsi   = cos(psi);
              SinTheta = sin(theta);
              CosTheta = cos(theta);
              Rz1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Rz1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Rz1");
                  return(ERROR);
              }

              Ry = (double *)malloc((size_t) 16L * sizeof(double));
              if (Ry == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Ry");
                  free(Rz1);
                  return(ERROR);
              }
              std::cout << " 819 phi = " << phi << std::endl;
              Rz2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Rz2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Rz2");
                  free(Rz1);
                  free(Ry);
                  return(ERROR);
              }
              std::cout << " 829 phi = " << phi << std::endl;
              if (GetIdentitySquareMatrix(Rz2, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }
              std::cout << " 839 phi = " << phi << std::endl;
              hlp    =  Rz2;
              *hlp++ =  CosPsi;
              *hlp   =  SinPsi;
              hlp    += (std::ptrdiff_t)3L;
              *hlp++ =  - SinPsi;
              *hlp   =  CosPsi;

              if (GetIdentitySquareMatrix(Rz1, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }

              hlp    =  Rz1;
              *hlp++ =  CosPhi;
              *hlp   =  SinPhi;
              hlp    += (std::ptrdiff_t)3L;
              *hlp++ =  - SinPhi;
              *hlp   =  CosPhi;
              std::cout << " 863 phi = " << phi << std::endl;
              if (GetIdentitySquareMatrix(Ry, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");
                  //free(defamp);
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }

              hlp  =  Ry;
              *hlp =  CosTheta;
              hlp  += (std::ptrdiff_t)2L;
              *hlp =  - SinTheta;
              hlp  += (std::ptrdiff_t)6L;
              *hlp =  SinTheta;
              hlp  += (std::ptrdiff_t)2L;
              *hlp =  CosTheta;

              R = (double *)malloc((size_t) 16L * sizeof(double));
              if (R == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for R");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }



              if (multiply_3Matrices(Rz2, Ry, Rz1, R, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  return(ERROR);
              }

              DRz1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DRz1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DRz1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  return(ERROR);
              }

              DRy = (double *)malloc((size_t) 16L * sizeof(double));
              if (DRy == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DRy");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  return(ERROR);
              }

              DRz2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DRz2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DRz2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  return(ERROR);
              }

              for (i = 0L; i < 16L; i++)
              {
                  DRz2[i] = 0.0;
                  DRz1[i] = 0.0;
                  DRy[i] = 0.0;
              }

              hlp    =  DRz2;
              *hlp++ =  - SinPsi;
              *hlp   =  CosPsi;
              hlp    += (std::ptrdiff_t)3L;
              *hlp++ =  - CosPsi;
              *hlp   =  - SinPsi;

              hlp    =  DRz1;
              *hlp++ =  - SinPhi;
              *hlp   =  CosPhi;
              hlp    += (std::ptrdiff_t)3L;
              *hlp++ =  - CosPhi;
              *hlp   =  - SinPhi;

              hlp  =  DRy;
              *hlp =  - SinTheta;
              hlp  += (std::ptrdiff_t)2L;
              *hlp =  - CosTheta;
              hlp  += (std::ptrdiff_t)6L;
              *hlp =  CosTheta;
              hlp  += (std::ptrdiff_t)2L;
              *hlp =  - SinTheta;
              for (i = 0L; i < 16L; i++)
              {
                  std::cout << "DRz2["<<i<<"] = " << DRz2[i];
                  std::cout << "DRz1["<<i<<"] = " << DRz1[i];
                  std::cout << "DRy ["<<i<<"] = " << DRy[i]  << std::endl;
              }
              std::cout << " 975 phi = " << phi << std::endl;
              DR0 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DR0 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DR0");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  return(ERROR);
              }
              std::cout << " 990 phi = " << phi << std::endl;
              DR1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DR1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DR1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  return(ERROR);
              }
              std::cout << " 1006 phi = " << phi << std::endl;
              DR2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DR2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DR2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  return(ERROR);
              }

              std::cout << " 1024 phi = " << phi << std::endl;
              if (multiply_3Matrices(Rz2, Ry, DRz1, DR0, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  //free(defamp);
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }

              std::cout << " 1043 phi = " << phi << std::endl;
              if (multiply_3Matrices(Rz2, DRy, Rz1, DR1, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }


              if (multiply_3Matrices(DRz2, Ry, Rz1, DR2, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }

              Tr = (double *)malloc((size_t) 16L * sizeof(double));
              if (Tr == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Tr");
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }

              if (GetIdentitySquareMatrix(Tr, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");

                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(Tr);
                  return(ERROR);
              }

              std::cout << " 1256 phi = " << phi << std::endl;
              hlp  =  Tr;
              hlp  += (std::ptrdiff_t)3L;
              *hlp =  - x0; //global_flexible_prog->sctrans * x0
              hlp  += (std::ptrdiff_t)5L;
              *hlp =  - y0; //global_flexible_prog->sctrans * y0
std::cout << "1060 Parameters =" << Parameters << "*cost = "<< global_flexible_prog->costfunctionvalue << std::endl;
              ProjectionRefencePoint(Parameters,dim,R,Tr,rg_projimage,P_esp_image,Xwidth,Ywidth,sigmalocal);
              std::cout << "1062 Parameters =" << Parameters << std::endl;
              /*if (ProjectionRefencePoint(dim,R,Tr,P_mu_image,cost,Xwidth,Ywidth,sigma) == ERROR)
          {
               WRITE_ERROR(return_gradhesscost, "ProjectionRefencePoint");
               //free(defamp);
               free(DP_q);
               free(DP_y);
               free(DP_x);
               free(DP_Rz2);
               free(DP_Ry);
               free(DP_Rz1);
               free(Rz1);
               free(Ry);
               free(Rz2);
               free(R);
               free(DRz1);
               free(DRy);
               free(DRz2);
               free(DR0);
               free(DR1);
               free(DR2);
               free(Tr);
               return(ERROR);global_flexible_prog->
               }*/

              MultidimArray<double> DP_Rx(Xwidth,Ywidth),DP_Ry(Xwidth,Ywidth),DP_Rz2(Xwidth,Ywidth);
              MultidimArray<double> DP_x(Xwidth,Ywidth),DP_y(Xwidth,Ywidth);
              MultidimArray<double> DP_q(dim,Xwidth,Ywidth);
              if (partialpfunction(Parameters,centerOfMass,R,Tr,DR0,DR1,DR2,DP_Rx,DP_Ry,DP_Rz2,DP_x,DP_y,DP_q,rg_projimage,P_esp_image,Xwidth,Ywidth) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by partialpfunction");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(Tr);
                  return(ERROR);
              }

              helpgr = (double *)malloc((size_t)((long)trialSize) * sizeof(double));
              if (helpgr == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for helpgr");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(Tr);
                  return(ERROR);
              }


              for( i=0;i<trialSize;i++)
              {
                  Gradient[i] = 0.0;
                  for( j=0;j<trialSize;j++)
                  {
                      Hessian[i*trialSize+j] = 0.0;
                  }
              }
              std::cout << " 1352 phi = " << phi << std::endl;
              for(int kx=0;kx<Xwidth;kx++)
                  for(int ky=0;ky<Ywidth;ky++)
                  {
                      difference = rg_projimage(kx-half_Xwidth,ky-half_Ywidth) - P_esp_image(kx-half_Xwidth,ky-half_Ywidth);
                      //std::cout << " 1399 difference = " << difference << std::endl;
                      helpgr[0] = DP_Rx(kx,ky); //std::cout << "  helpgr= " << helpgr[0];
                      helpgr[1] = DP_Ry(kx,ky); //std::cout << "||" << helpgr[1];
                      helpgr[2] = DP_Rz2(kx,ky);//std::cout << "||" << helpgr[2];
                      helpgr[3] = DP_x(kx,ky);  //std::cout << "||" << helpgr[3];
                      helpgr[4] = DP_y(kx,ky);  //std::cout << "||" << helpgr[4];
                      for(j=0;j<dim;j++)
                      {
                          helpgr[5+j] =DP_q(j,kx,ky);//std::cout << "||" << helpgr[5+j];
                      }
                      gradhesscost_atpixel(Gradient,Hessian,helpgr,difference);
                      /*if (gradhesscost_atpixel(Gradient,Hessian,helpgr,difference) == ERROR)
                  {
                       WRITE_ERROR(return_gradhesscost, "Error returned by partialpfunction");
                       //free(defamp);
                       free(DP_q);
                       free(DP_y);
                       free(DP_x);
                       free(DP_Rz2);
                       free(DP_Ry);
                       free(DP_Rz1);
                       free(Rz1);
                       free(Ry);
                       free(Rz2);
                       free(R);
                       free(DRz1);
                       free(DRy);
                       free(DRz2);
                       free(DR0);
                       free(DR1);
                       free(DR2);
                       free(Tr);
                       return(ERROR);
                       }*/
                  }
              for (i=0;i<trialSize;i++)
              {
                  for(j=0;j<trialSize;j++)
                      std::cout << Hessian[i*trialSize+j] << "||";
                  std::cout << std::endl;
              }

              std::cout << " 1392 phi = " << phi << std::endl;

              for (int i=0;i<trialSize;i++)
              {
                  Hessian[i*(trialSize) + i] += lambda*Hessian[i*trialSize + i] ;
                  if(i+1<trialSize)
                      for(int j=i+1;j<trialSize;j++)
                      {
                          Hessian[i*trialSize + j]= Hessian[j*trialSize + i];
                      }
              }
              for (i=0;i<trialSize;i++)
              {
                  for(j=0;j<trialSize;j++)
                      std::cout << Hessian[i*trialSize+j] << "||";
                  std::cout << std::endl;
              }

              for(i=0;i<trialSize;i++)
              {
                  std::cout << "Gradient = " << Gradient[i] << std::endl;
              }


              std::cout << " 1446 phi = " << phi << std::endl;
              free(helpgr);
              free(Rz1);
              std::cout << " 1454 phi = " << phi << std::endl;
              free(Ry);
              std::cout << " 1455 phi = " << phi << std::endl;
              free(Rz2);
              std::cout << " 1456 phi = " << phi << std::endl;
              free(R);
              std::cout << " 1457 phi = " << phi << std::endl;
              free(DRz1);
              std::cout << " 1458 phi = " << phi << std::endl;
              free(DRy);
              std::cout << " 1459 phi = " << phi << std::endl;
              free(DRz2);
              std::cout << " 1460 phi = " << phi << std::endl;
              free(DR0);
              std::cout << " 1461 phi = " << phi << std::endl;
              free(DR1);
              std::cout << " 1462 phi = " << phi << std::endl;
              free(DR2);
              std::cout << " 1463 phi = " << phi << std::endl;
              free(Tr);
              std::cout << " 1464 phi = " << phi << std::endl;
              return(!ERROR);
          }/* End of return_gradhesscost */

          /* ------------------------------------------------------------------------- */
          /* Optimizer                                                                 */
          /* ------------------------------------------------------------------------- */
          int levenberg_cst2(
              MultidimArray<double>  &lc_P_mu_image,
              MultidimArray<double>  &P_esp_image,
              Matrix1D<double>  &centerOfMass,
              double *beta,
              double *alpha,
              //double *cost,
              Matrix1D<double>  &Parameters,
              double OldCost,
              double *lambda,
              double LambdaScale,
              long   *iter,
              double tol_angle,
              double tol_shift,
              double tol_defamp,
              int    *IteratingStop,
              size_t Xwidth,
              size_t Ywidth)
          {

#ifndef DBL_EPSILON
#define DBL_EPSILON 2e-16
#endif
              double   *da, epsilon = DBL_EPSILON;
              double   *u,   *v, *w;
              double   *t;
              double   wmax, thresh;
              int      i, j;
              long     ma = (long) VEC_XSIZE(global_flexible_prog->trial);
              double   hlp;
              double   *costchanged=NULL,max_defamp;
              int      dim = global_flexible_prog->numberOfModes;
              double   *a;
              int     Status = !ERROR;
              std::cout << "1563 *lambda =" << *lambda << "  lambda ="<<*lambda<<"  LambdaScale ="<<LambdaScale<< std::endl;
              std::cout << "1470 Parameters" << Parameters << std::endl;
              std::cout << tol_angle<<tol_shift<<tol_defamp<<Xwidth<<Ywidth <<std::endl;

              costchanged  = (double *)malloc((size_t)(1L) * sizeof(double));
              a  = (double *)malloc((size_t)ma * sizeof(double));
              /*if (a == (double *)NULL)
              {
                  WRITE_ERROR(levenberg_cst2, "ERROR - Not enough memory for a in levenberg_cst2");
                  return(ERROR);
              }*/
              for(i=0;i<ma;i++)
              {
                  a[i]=Parameters(i);
                  std::cout << "a[i]" << a[i]<<std::endl;
              }


              da = (double *)malloc((size_t)ma * sizeof(double));
              /*if (da == (double *)NULL)
              {
                  WRITE_ERROR(levenberg_cst2, "ERROR - Not enough memory for da in levenberg_cst2");
                  free(a);
                  return(ERROR);
              }*/
              u = (double *)malloc((size_t)(ma * ma) * sizeof(double));
              /*if (u == (double *)NULL)
              {
                  free(a);
                  free(da);
                  WRITE_ERROR(levenberg_cst2, "ERROR - Not enough memory for u in levenberg_cst2");
                  return(ERROR);
              }*/
              v = (double *)malloc((size_t)(ma * ma) * sizeof(double));
              /*if (v == (double *)NULL)
              {
                  free(a);
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst2, "ERROR - Not enough memory for v in levenberg_cst2");
                  return(ERROR);
              }*/
              w = (double *)malloc((size_t)ma * sizeof(double));
              /*if (w == (double *)NULL)
              {
                  free(a);
                  free(v);
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst2, "ERROR - Not enough memory for w in levenberg_cst2");
                  return(ERROR);
              }*/

              std::cout << "1547 cost" << global_flexible_prog->costfunctionvalue_cst << std::endl;
              t = u;
              for (i = 0L; (i < ma); t += (std::ptrdiff_t)(ma + 1L), i++)
              {
                  for (j = 0L; (j < ma); alpha++, j++)
                  {
                      *u++ = -*alpha;/*std::cout << " u= " << *u;*/
                  }

                  *t *= 1.0 + *lambda;
              }
              u -= (std::ptrdiff_t)(ma * ma);
              alpha -= (std::ptrdiff_t)(ma * ma);
              for(i=0L;(i<ma*ma);i++)
              {
                  std::cout <<",  u"<<i<< " = " << u[i] ;
              }
              

              SingularValueDecomposition(u, ma, ma, w, v, SVDMAXITER, &Status);

              /*if (SingularValueDecomposition(u, ma, ma, w, v, SVDMAXITER, &Status) == ERROR)
              {
                  free(a);
                  std::cout << "1571" << std::endl;
                  free(w);
                  std::cout << "1572" << std::endl;
                  free(v);
                  std::cout << "1573" << std::endl;
                  free(u);
                  std::cout << "1574" << std::endl;
                  free(da);
                  std::cout << "1575" << std::endl;
                  WRITE_ERROR(levenberg_cst2, "ERROR - Unable to perform svdcmp in levenberg_cst2");
                  return(ERROR);
              }*/
              wmax = 0.0;
              std::cout << "1571" << std::endl;
              t = w + (std::ptrdiff_t)ma;
              while (--t >= w)
                  if (*t > wmax)
                      wmax = *t;
              thresh = epsilon * wmax;
              w += (std::ptrdiff_t)ma;
              j = ma;
              while (--j >= 0L)
              {
                  if (*--w < thresh)
                  {
                      *w = 0.0;
                      for (i = 0; (i < ma); i++)
                      {
                          u[i * ma + j] = 0.0;
                          v[i * ma + j] = 0.0;
                      }
                  }
              }
              

              SingularValueBackSubstitution(u, w, v, ma, ma, beta, da, &Status);

              /*if (SingularValueBackSubstitution(u, w, v, ma, ma, beta, da, &Status) == ERROR)
              {
                  free(a);
                  free(w);
                  free(v);
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst2, "ERROR - Unable to perform svbksb in levenberg_cst2");
                  return(ERROR);
              }*/
              std::cout << "1600" << std::endl;

              v  =  (double *)memcpy(v, a, (size_t)ma * sizeof(double));
              t  =  v + (std::ptrdiff_t)ma;
              a  += (std::ptrdiff_t)ma;
              da += (std::ptrdiff_t)ma;
              while (--t >= v)
              {
                  da--;
                  *t = *--a;
                  *t += *da;
              }
              std::cout << "1612" << std::endl;
              //Matrix1D<double>  Parameters((int)ma);
              for(i=0;i<ma;i++)
                  Parameters(i) = v[i];

              return_gradhesscost(centerOfMass,w, u,  Parameters,dim,lc_P_mu_image,P_esp_image,Xwidth,Ywidth);
              costchanged[0]=global_flexible_prog->costfunctionvalue;

              /*if (return_gradhesscost(centerOfMass,w, u,  Parameters,dim,lc_P_mu_image,P_esp_image,Xwidth,Ywidth) == ERROR)
              {
                  WRITE_ERROR(levenberg_cst2, "Error returned by total_gradhesscost");
                  free(a);
                  free(w);
                  free(v);
                  free(u);
                  free(da);
                  return(ERROR);
              }*/

              std::cout << "1633 Parameters" << Parameters<<std::endl;
              std::cout << "w" << w[0] << "||" << w[1] << "||" << w[2] << "||" << w[3] << "||" << w[4] << "||" << w[5] << "||" << w[6] << "||" << std::endl;
              std::cout << "u = " ;
              for (i=0;i<ma;i++)
              {
                  for(j=0;j<ma;j++)
                      std::cout << u[i*ma+j] << "||";
                  std::cout << std::endl;
              }

              for(max_defamp=0.0, i=0L; i< dim;i++)
              {
                  hlp = fabs(a[i+5] - v[i+5]);
                  if (hlp > max_defamp)
                      max_defamp = fabs(a[i+5] - v[i+5]);
              }

              (*iter)++;
              if (costchanged[0] < OldCost)
              {
                  if ((fabs(a[0] - v[0]) < tol_angle) && (fabs(a[1] - v[1]) < tol_angle) && (fabs(a[2] - v[2]) < tol_angle))
                  {
                      if ((fabs(a[3] - v[3]) < tol_shift) && (fabs(a[4] - v[4]) < tol_shift))
                      {
                          if(max_defamp < tol_defamp)
                          {
                              *IteratingStop = 1;
                          }
                      }
                  }
              }

              if (global_flexible_prog->costfunctionvalue_cst == -1.0)
              {
                  if (costchanged[0] < OldCost)
                  {
                      for (i = 0L; (i < ma); i++)
                      {
                          for (j = 0L; (j < ma); j++)
                              alpha[i * ma + j] = u[i * ma + j];
                          beta[i] = w[i];
                          a[i] = v[i];
                      }
                      global_flexible_prog->costfunctionvalue_cst = costchanged[0];
                  }

                 /* free(w);
                  free(u);
                  free(da);
                  free(v);*/
                  return(!ERROR);
              }

              for (i = 0L; (i < ma); i++)
              {
                  for (j = 0L; (j < ma); j++)
                      alpha[i * ma + j] = u[i * ma + j];
                  beta[i] = w[i];
                  a[i] = v[i];
              }
              global_flexible_prog->costfunctionvalue_cst = costchanged[0];


#ifndef DBL_MIN
#define DBL_MIN 1e-26
#endif
#ifndef DBL_MAX
#define DBL_MAX 1e+26
#endif

              std::cout << "1774 cost" << global_flexible_prog->costfunctionvalue_cst<<std::endl;

              if (costchanged[0] < OldCost)
              {
                  if (*lambda > DBL_MIN)
                  {
                      *lambda /= LambdaScale;
                  }
                  else
                  {
                      *IteratingStop = 1;
                  }
              }
              else
              {
                  if (*lambda < DBL_MAX)
                  {
                      *lambda *= LambdaScale;
                  }
                  else
                  {
                      *IteratingStop = 1;
                  }
              }
              std::cout << "1797 *cost" << global_flexible_prog->costfunctionvalue_cst<<std::endl;
              /*free(w);
              std::cout << "1798 *cost" << global_flexible_prog->costfunctionvalue_cst<<std::endl;
              free(v);
              std::cout << "1799 cost" << global_flexible_prog->costfunctionvalue_cst<<std::endl;
              free(u);
              std::cout << "1800 cost" << global_flexible_prog->costfunctionvalue_cst<<std::endl;
              free(da);
              std::cout << "1801 cost" << global_flexible_prog->costfunctionvalue_cst<<std::endl;*/
              return(!ERROR);
          } /* End of levenberg_cst2 */



          /* Registration ------------------------------------------------------------ */
          //cstregistrationcontinuous(centerOfMass,cost,Parameters,P_mu_image,P_esp_image,Xwidth,Ywidth);
          int  cstregistrationcontinuous(
              Matrix1D<double>  &centerOfMass,
              //double            *cost,
              Matrix1D<double>  &Parameters,
              MultidimArray<double>  &cst_P_mu_image,
              MultidimArray<double>  &P_esp_image,
              size_t            Xwidth,
              size_t            Ywidth)
          {
              const double    epsilon = DBL_EPSILON;
              int             Status = !ERROR, DoDesProj, IteratingStop, FlagMaxIter;
              long            MaxIter, MaxIter1, MaxIter2, iter;
              long            SizeIm, MaxNumberOfFailures, SatisfNumberOfSuccesses, nSuccess, nFailure;
              double          *pntr_time, *pntr_FailureIter, *pntr_RedftCompProj, *pntr_ImdftCompProj;
              double          LambdaScale=2., OldCost, tol_angle, tol_shift,tol_defamp;
              double          OneIterInSeconds,  *Gradient, *Hessian;
              time_t          time1, time2, *tp1 = NULL, *tp2 = NULL;
              long            dim = (long) global_flexible_prog->numberOfModes;
              long            MaxNoIter,MaxNoFailure,SatisfNoSuccess;

              double lambda=1000.;

              DoDesProj        = 0;
              MaxNoIter        = global_flexible_prog->max_no_iter;
              MaxNoFailure     = (long)(0.3 * MaxNoIter);
              SatisfNoSuccess  = (long)(0.7 * MaxNoIter);
              MaxIter  = MaxNoIter;
              MaxIter1 = MaxIter - 1L;

              MaxNumberOfFailures     = MaxNoFailure;
              SatisfNumberOfSuccesses = SatisfNoSuccess;
              tol_angle  = 0.0;
              tol_shift  = 0.0;
              tol_defamp = 0.0;

              std::cout << "cost = " << global_flexible_prog->costfunctionvalue_cst << std::endl;

              Gradient = (double *)malloc((size_t) (dim + 5L) * sizeof(double));
              if (Gradient == (double *)NULL)
              {
                  WRITE_ERROR(cstregistrationcontinuous, "ERROR - Not enough memory for Gradient");
                  return(ERROR);
              }

              Hessian = (double *)malloc((size_t) (dim+5L)*(dim+5L) * sizeof(double));
              if (Hessian == (double *)NULL)
              {
                  WRITE_ERROR(cstregistrationcontinuous, "ERROR - Not enough memory for Hessian");
                  free(Gradient);
                  return(ERROR);
              }

              if (DoDesProj && (global_flexible_prog->currentStage == 2))
              {
		  Parameters(0)=10;
                  Parameters(1)=0;
                  Parameters(2)=0;
                  Parameters(3)=0;
                  Parameters(4)=0;
                  Parameters(5)=0.5;
                  //Parameters(6)=0;
              }

              std::cout << "1819Parameters = " << Parameters << std::endl;
              if (return_gradhesscost(centerOfMass,Gradient, Hessian, Parameters,
                                      dim,cst_P_mu_image,P_esp_image,Xwidth,Ywidth) == ERROR)
              {
                  WRITE_ERROR(cstregistrationcontinuous, "Error returned by return_gradhesscost");
                  //free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  return(ERROR);
              }
              std::cout << "1935 cost = " << global_flexible_prog->costfunctionvalue_cst << std::endl;
              time2 = time(tp2);
              OneIterInSeconds = difftime(time2, time1);
              if (DoDesProj && (global_flexible_prog->currentStage == 2))
              {
                  Image<double> Itemp;
                  //Itemp().setXmippOrigin();
                  Itemp()=cst_P_mu_image;
          	  Itemp.write("test_refimage.spi");

                  free(Gradient);
                  free(Hessian);
                  return(!ERROR);
              }
              std::cout << "1864time2 = " << time2 << std::endl;
              if ((MaxIter == 0L) && (!DoDesProj))
              {
                  free(Gradient);
                  free(Hessian);
                  return(!ERROR);
              }
              nSuccess = 0L;
              nFailure = 0L;
              OldCost   = global_flexible_prog->costfunctionvalue_cst;
              iter = - 1L;
              IteratingStop = 0;
              FlagMaxIter = (MaxIter != 1L);
              if (!FlagMaxIter)
            	  global_flexible_prog->costfunctionvalue_cst = - 1.0;
              std::cout << "1984 *cost = " << global_flexible_prog->costfunctionvalue_cst << std::endl;

              do
              {
                  time1 = time(tp1);

                  std::cout << "1992 *cost = " << global_flexible_prog->costfunctionvalue_cst << " cost = " << global_flexible_prog->costfunctionvalue_cst << std::endl;
                  std::cout << "cost = " << global_flexible_prog->costfunctionvalue_cst << std::endl;
                  if (levenberg_cst2(cst_P_mu_image,P_esp_image,centerOfMass,Gradient, Hessian,  Parameters,OldCost, &lambda, LambdaScale, &iter, tol_angle, tol_shift, tol_defamp,&IteratingStop,Xwidth,Ywidth) == ERROR)
                  {
                      WRITE_ERROR(cstregistrationcontinuous, "Error returned by levenberg_cst2");
                      //free(Parameters);
                      free(Gradient);
                      free(Hessian);
                      return(ERROR);
                  }
                  std::cout << "2002 *cost = " << global_flexible_prog->costfunctionvalue_cst  << " cost = " << global_flexible_prog->costfunctionvalue_cst  << std::endl;
                  time2 = time(tp2);

                  if (global_flexible_prog->costfunctionvalue_cst < OldCost)
                  {
                      OldCost = global_flexible_prog->costfunctionvalue_cst;
                      nSuccess++;
                      pntr_FailureIter++;
                      if (nSuccess >= SatisfNumberOfSuccesses)
                      {
                          break;
                      }
                      if (IteratingStop)
                      {
                          break;
                      }
                  }
                  else
                  {
                      nFailure++;
                      *pntr_FailureIter++ = (iter + 1L);
                  }

              }
              while ((nFailure <= MaxNumberOfFailures) && (iter < MaxIter1) && FlagMaxIter);
              std::cout << "1950 OldCost = " << OldCost << std::endl;

              free(Gradient);
              free(Hessian);

              return(!ERROR);
          }

          void ProgFlexibleAlignment::performCompleteSearch(int pyramidLevel)
          {

              int       dim = numberOfModes;
              int       ModMaxdeDefamp, ModpowDim,help;
              double    SinPhi,CosPhi,SinPsi,CosPsi;
              double    phi,theta,psi;
              double    *Rz1,*Ry,*Rz2,*hlp;
              double    S_muMin = 1e30;

              double    *R,*Tr;
              double    x0,y0;
              //double    *cost;
              costfunctionvalue = 0.0;
              Matrix1D<double> Parameters(dim+5);
              Matrix1D<double> limit0(3), limitF(3), centerOfMass(3);
              const char *intensityColumn = "Bfactor";
              computePDBgeometry(fnPDB, centerOfMass, limit0, limitF, intensityColumn);
              centerOfMass = (limit0 + limitF)/2;

              FileName fnRandom;
              fnRandom.initUniqueName(nameTemplate,fnOutDir);
              const char * randStr = fnRandom.c_str();
              std::string command;

          	// Reduce the image
          	fnDown = formatString("%s_downimg.xmp", fnRandom.c_str());
          	if (pyramidLevel != 0) {
          		Image<double> I;
          		I.read(currentImgName);
          		selfPyramidReduce(BSPLINE3, I(), pyramidLevel);
          		I.write(fnDown);
          	} 
                /*else
          		link(currentImgName.c_str(), fnDown.c_str());*/

std::cout << "2074" << fnDown << std::endl;

              Image<double> imgtemp;
              imgtemp.read(fnDown);


              imgtemp().setXmippOrigin();
              //P_esp_image(imgtemp);
              MultidimArray<double>  P_mu_image;
              P_mu_image.resize(imgtemp());
              int Xwidth = XSIZE(imgtemp());
              int Ywidth = YSIZE(imgtemp());
              double    reduce_rate = 1;
              std::cout << "2116 reduce_rate = "<< reduce_rate << std::endl;

              //command=(std::string)"rm -f "+ fnDown;
              //system(command.c_str());

              ModMaxdeDefamp = (int)floor(maxdefamp/defampsampling) + 1;
              ModpowDim      = (int)pow(ModMaxdeDefamp,dim);

              Rz1 = (double *)malloc((size_t) 16L * sizeof(double));
              /*if (Rz1 == (double *)NULL)
          {
                  WRITE_ERROR(performCompleteSearch, "ERROR - Not enough memory for Rz1");
                  //free(Q);
                  return(ERROR);
          }*/

              Ry = (double *)malloc((size_t) 16L * sizeof(double));
              /*if (Ry == (double *)NULL)
          {
                 WRITE_ERROR(performCompleteSearch, "ERROR - Not enough memory for Ry");
                  //free(Q);
                  free(Rz1);
                  return(ERROR);
          }*/

              Rz2 = (double *)malloc((size_t) 16L * sizeof(double));
              /*if (Rz2 == (double *)NULL)
          {
                  WRITE_ERROR(performCompleteSearch, "ERROR - Not enough memory for Rz2");
                  //free(Q);
                  free(Rz1);
                  free(Ry);
                  return(ERROR);
          }*/

              R = (double *)malloc((size_t) 16L * sizeof(double));
              /*if (R == (double *)NULL)
          {
                  WRITE_ERROR(performCompleteSearch, "ERROR - Not enough memory for R");
                  //free(Q);
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
          }*/

              Tr = (double *)malloc((size_t) 16L * sizeof(double));
              /*if (R == (double *)NULL)
          {
              WRITE_ERROR(performCompleteSearch, "ERROR - Not enough memory for R");
                  //free(Q);
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(Tr);
                  return(ERROR);
          }*/


              for(phi = 0; phi < 360; phi += minAngularSampling)
              {
                  Parameters(0) = phi;
                  trial(0) = phi;
                  SinPhi = sin(phi/180*PI);
                  CosPhi = cos(phi/180*PI);


                  GetIdentitySquareMatrix(Ry, 4L);
                  /*if (GetIdentitySquareMatrix(Ry, 4L) == ERROR)
              {
                      WRITE_ERROR(performCompleteSearch, "Error returned by GetIdentitySquareMatrix");
                      //free(Q);
                      free(Rz1);
                      free(Ry);
                      free(Rz2);
                      free(Tr);
                      return(ERROR);
              }*/

                  hlp = Ry;
                  *hlp = CosPhi;
                  hlp += (std::ptrdiff_t)2L;
                  *hlp = - SinPhi;
                  hlp += (std::ptrdiff_t)6L;
                  *hlp = SinPhi;
                  hlp += (std::ptrdiff_t)2L;
                  *hlp = CosPhi;

                  for(theta = 0; theta < 180; theta += minAngularSampling )
                  {
                      Parameters(1) = theta;
                      trial(1)  =  theta;

                      GetIdentitySquareMatrix(Rz1, 4L);
                      /*if (GetIdentitySquareMatrix(Rz1, 4L) == ERROR)
                  {
                           WRITE_ERROR(performCompleteSearch, "Error returned by GetIdentitySquareMatrix");
                           //free(Q);
                           free(Rz1);
                           free(Ry);
                           free(Rz2);
                           free(Tr);
                           return(ERROR);
                  }*/

                      hlp = Rz1;
                      *hlp++ = CosPhi;
                      *hlp = SinPhi;
                      hlp += (std::ptrdiff_t)3L;
                      *hlp++ = - SinPhi;
                      *hlp = CosPhi;


                      for(psi = 0; psi < 360; psi += minAngularSampling)
                      {
                          Parameters(2) = psi;
                          trial(2) =  psi;
                          SinPsi   =  sin(psi/180*PI);
                          CosPsi   =  cos(psi/180*PI);


                          GetIdentitySquareMatrix(Rz2, 4L);
                          /*if (GetIdentitySquareMatrix(Rz2, 4L) == ERROR)
                      {
                              WRITE_ERROR(performCompleteSearch, "Error returned by GetIdentitySquareMatrix");
                              //free(Q);
                              free(Rz1);
                              free(Ry);
                              free(Rz2);
                              free(Tr);
                              return(ERROR);
                      }*/

                          hlp = Rz2;
                          *hlp++ = CosPsi;
                          *hlp = SinPsi;
                          hlp += (std::ptrdiff_t)3L;
                          *hlp++ = - SinPsi;
                          *hlp = CosPsi;


                          multiply_3Matrices(Rz2, Ry, Rz1, R, 4L, 4L, 4L, 4L);
                          /*if (multiply_3Matrices(Rz2, Ry, Rz1, R, 4L, 4L, 4L, 4L) == ERROR)
                      {
                               WRITE_ERROR(performCompleteSearch, "Error returned by multiply_3Matrices");
                               free(Rz1);
                               free(Ry);
                               free(Rz2);
                               free(R);
                               free(Tr);
                               return(ERROR);
                      }*/

                          for(x0 = 0; x0 <= maxtransl; x0 += translsampling)
                          {
                              Parameters(3) = x0/reduce_rate;
                              trial(3)  =  x0;
                              for(y0 = 0; y0 <= maxtransl; y0 += translsampling)
                              {
                                  Parameters(4) = y0/reduce_rate;
                                  trial(4)  =  y0;


                                  GetIdentitySquareMatrix(Tr, 4L);
                                  /*if (GetIdentitySquareMatrix(Tr, 4L) == ERROR)
                              {
                                      WRITE_ERROR(performCompleteSearch, "Error returned by GetIdentitySquareMatrix");
                                      //free(Q);
                                      free(Rz1);
                                      free(Ry);
                                      free(Rz2);
                                      free(Tr);
                                      return(ERROR);
                              }*/
                                  hlp = Tr;
                                  hlp += (std::ptrdiff_t)3L;
                                  *hlp = - x0;
                                  hlp += (std::ptrdiff_t)5L;
                                  *hlp = - y0;


                                  MatrixMultiply( R, Tr,Tr,4L, 4L, 4L);
                                  /*if (MatrixMultiply( R, Tr,Tr,4L, 4L, 4L) == ERROR)
                              {
                                       WRITE_ERROR(performCompleteSearch, "Error returned by multiply_3Matrices");
                                       free(Rz1);
                                       free(Ry);
                                       free(Rz2);
                                       free(R);
                                       free(Tr);
                                       return(ERROR);
                              }*/


                                  for(int i=0;i< ModpowDim; i++)
                                  {
                                      help =  i;
                                      for(int j=dim-1; j>=0;j--)
                                      {
                                          Parameters(j+5) = (double) (help % ModMaxdeDefamp) * defampsampling;
                                          trial(j+5) = (double) (help % ModMaxdeDefamp) * defampsampling;
                                          help = (int)floor(help/ModMaxdeDefamp);

                                      }

                                      ProjectionRefencePoint(Parameters,dim,R,Tr,P_mu_image,imgtemp(),Xwidth,Ywidth,sigma);
                                      /*if (ProjectionRefencePoint(Parameters,dim,R,Tr,P_mu_image,imgtemp(),cost,Xwidth,Ywidth,sigma) == ERROR)
                                  {
                                           WRITE_ERROR(return_gradhesscost, "Error returned by ProjectionRefencePoint");
                                           free(Rz1);
                                           free(Ry);
                                           free(Rz2);
                                           free(R);
                                           free(Tr);
                                           return(ERROR);
                                  }*/

                                      if (costfunctionvalue < S_muMin)
                                      {//2
                                          S_muMin  =  costfunctionvalue;
                                          for(int k = 0; k< dim + 5 ; k++ )
                                          {
                                              trial_best(k) = trial(k);
                                          }

                                      }//2--2213
                                      std::cout << "trial" << trial <<"  cost=" << costfunctionvalue <<std::endl;
                                      std::cout << "trial_best" << trial_best <<"  cost_best=" <<S_muMin <<std::endl;
                                  }//3--
                              }
                          }
                      }
                  }

              }

              std::cout << "trial" << trial <<costfunctionvalue <<std::endl;
              std::cout << "trial_best" << trial_best <<S_muMin <<std::endl;
              for (int i=0; i<dim+5; i++)
              {
                  trial(i)=trial_best(i);
                  parameters(i)=trial_best(i);
              }

              //if (pyramidLevel!=0)
              //{
              //    std::string command=(std::string)"rm -f " + fnDown;
              //   system(command.c_str());
              //}

              MetaData DF_out_discrete;
              size_t id=DF_out_discrete.addObject();
              DF_out_discrete.setValue(MDL_IMAGE,currentImgName,id);
              std::vector<double> aux;
              aux.resize(VEC_XSIZE(parameters));
              FOR_ALL_ELEMENTS_IN_MATRIX1D(parameters)
              aux[i]=parameters(i);
              DF_out_discrete.setValue(MDL_NMA,aux,id);
              DF_out_discrete.write("result_angular_discrete.xmd");
              //command=(std::string)"rm -f "+ fnDown;
              //system(command.c_str());
          }

          // Continuous assignment ===================================================
          double ProgFlexibleAlignment::performContinuousAssignment(int pyramidLevel) 
          {
              int    dim = numberOfModes;
              //double *cost;
              std::cout << "Bfactor001" << std::endl;
              costfunctionvalue = 0.0;
              costfunctionvalue_cst = 0.0;
              //std::cout << "cost = " << costfunctionvalue << std::endl;

              Matrix1D<double> Parameters;

              /*for (int i=0; i<dim+5; i++)
              {
                  Parameters(i) = trial(i);
                  std::cout << "performContinuousAssignment_Parameters(" << i << ")" << Parameters(i) << std::endl;
              }*/
              Parameters = trial_best;

              std::cout << "pyramidLevel  " << pyramidLevel << std::endl;
              

              std::string command;
              if (pyramidLevel==0)
              {
                  // Make links
                  /*command=formatString(
                		  "ln -sf %s %s",
                		  currentImgName.c_str(), fnDown.c_str());*/

                  FileName fnRandom;
                  fnRandom.initUniqueName(nameTemplate,fnOutDir);
                  const char * randStr = fnRandom.c_str();
                  fnDown = formatString("%s_downimg.xmp", fnRandom.c_str());

                  Image<double> I;
          	  I.read(currentImgName);
                  selfPyramidReduce(BSPLINE3, I(), pyramidLevel);
                  I.write(fnDown);
          	

                  //std::cout << "command.c_str()  " << command.c_str() << std::endl;
                  //system(command.c_str());
                  std::cout << "1986  " <<pyramidLevel<< std::endl;
              }
              std::cout << "19860  " <<pyramidLevel<< std::endl;

              //std::string arguments;
              //arguments = formatString("%s_downimg.xmp",fnRandom.c_str());
              Image<double> imgtemp;
              imgtemp.read(fnDown);
              std::cout << "1991" << std::endl;

              imgtemp().setXmippOrigin();
              MultidimArray<double>  P_mu_image,P_esp_image;
              P_esp_image = imgtemp();
              P_mu_image.resize(imgtemp());
              int Xwidth = XSIZE(imgtemp());
              int Ywidth = YSIZE(imgtemp());
              //double P_mu_image[Xwidth*Xwidth];
              Matrix1D<double> limit0(3), limitF(3), centerOfMass(3);
              const char *intensityColumn = "Bfactor";
              //Matrix1D<double> test=centerOfMass;
              std::cout << "Bfactor001" << std::endl;
              computePDBgeometry(fnPDB, centerOfMass, limit0, limitF, intensityColumn);

              std::cout << "Bfactor" << std::endl;
              centerOfMass = (limit0 + limitF)/2;

              cstregistrationcontinuous(centerOfMass,Parameters,P_mu_image,P_esp_image,Xwidth,Ywidth);
              std::cout << "centerOfMass_Bfactor003" << std::endl;
              /* Insert code for continious alignment */

              std::cout << "trial " << trial << std::endl;

              trial(3) *= pow(2.0, (double) pyramidLevel);
	      trial(4) *= pow(2.0, (double) pyramidLevel);

              trial_best=trial;

              std::cout << "trial_best " << trial_best << std::endl;


              /*command=formatString("rm -f %s ", fnDown.c_str());
              system(command.c_str());*/

              std::cout << "rm fndown " << std::endl;

              double outcost = costfunctionvalue_cst;
              
              std::cout << "outcost" << outcost << std::endl;

              return outcost;
          }


// Compute fitness =========================================================
          double ProgFlexibleAlignment::eval()
          {
              int dim = numberOfModes;

              int pyramidLevelDisc = 1;
              int pyramidLevelCont = (currentStage == 1) ? 1 : 0;

              if (currentStage == 1)
              {
                  performCompleteSearch(pyramidLevelDisc);
              }
              else
              {

                  //link(currentImgName.c_str(), fnDown.c_str());
              }
              double fitness = performContinuousAssignment(pyramidLevelCont);

              //std::cout << "Trial=" << trial_best.transpose() << " ---> " << fitness << std::endl;

              std::cout << "Fitness" << std::endl;
              return fitness;
          }


          void ProgFlexibleAlignment::processImage(const FileName &fnImg,
                  const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
          {
              static size_t imageCounter = 0;
              ++imageCounter;

              double scdefamp = 1000;

              int dim = numberOfModes;

              parameters.initZeros(dim + 5);
              currentImgName = fnImg;
              sprintf(nameTemplate, "_node%d_img%ld_XXXXXX", rangen, (long int)imageCounter);

              trial.initZeros(dim + 5);
              trial_best.initZeros(dim + 5);

              currentStage = 1;
#ifdef DEBUG

              std::cerr << std::endl << "DEBUG: ===== Node: " << rangen
              <<" processing image " << fnImg <<"(" << objId << ")"
              << " at stage: " << currentStage << std::endl;
#endif

              double fitness=eval();
              bestStage1 = trial = parameters = trial_best;

              currentStage = 2;
#ifdef DEBUG

              std::cerr << std::endl << "DEBUG: ===== Node: " << rangen
              <<" processing image " << fnImg <<"(" << objId << ")"
              << " at stage: " << currentStage << std::endl;
#endif

              fitness=eval();

std::cout << "step1" << std::endl;

              trial = trial_best;
              parameters = trial_best;

std::cout << "step2" << std::endl;

              parameters.resize(VEC_XSIZE(parameters) + 1);
              parameters(VEC_XSIZE(parameters) - 1) = fitness;

std::cout << "step3" << std::endl;

              writeImageParameters(fnImg);

          }

          void ProgFlexibleAlignment::writeImageParameters(const FileName &fnImg)
          {
              
std::cout << "parameters0" << parameters(0) << std::endl;
std::cout << "parameters1" << parameters(1) << std::endl;
std::cout << "parameters2" << parameters(2) << std::endl;
std::cout << "parameters3" << parameters(3) << std::endl;
std::cout << "parameters4" << parameters(4) << std::endl;

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

std::cout << "step4" << std::endl;

              for (int j = 5; j < 5 + dim; j++)
              {
                  vectortemp.push_back(parameters(j));
              }

std::cout << "step5" << std::endl;

              md.setValue(MDL_NMA, vectortemp, objId);

std::cout << "step6" << std::endl;

              md.setValue(MDL_COST, parameters(5 + dim), objId);

std::cout << "step7" << std::endl;

              md.append(fnOutDir+"/nmaDone.xmd");

std::cout << "step8" << std::endl;
          }
