/***************************************************************************
 *
 * Authors:
 *
 * Roberto Marabini
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

#include "angular_project_library.h"

/* Empty constructor ------------------------------------------------------- */
ProgAngularProjectLibrary::ProgAngularProjectLibrary()
{
    /** sampling object 1 by default*/
    mysampling.setSampling(1);
    Vshears=NULL;
    Vfourier=NULL;

}

ProgAngularProjectLibrary::~ProgAngularProjectLibrary()
{
    delete Vshears;
    delete Vfourier;
}

/* Read parameters --------------------------------------------------------- */
void ProgAngularProjectLibrary::readParams()
{
    input_volume = getParam("-i");
    output_file = getParam("-o");
    output_file_root = output_file.withoutExtension();
    fn_sym = getParam("--sym");
    fn_sym_neigh=checkParam("--sym_neigh")?getParam("--sym_neigh"):fn_sym;
    sampling = getDoubleParam("--sampling_rate");
    psi_sampling = getDoubleParam("--psi_sampling");
    max_tilt_angle = getDoubleParam("--max_tilt_angle");
    min_tilt_angle = getDoubleParam("--min_tilt_angle");
    angular_distance_bool = checkParam("--angular_distance");
    angular_distance=0.;
    if(angular_distance_bool)
    {
        FnexperimentalImages = getParam("--experimental_images");
        angular_distance = getDoubleParam("--angular_distance");
    }
    compute_closer_sampling_point_bool= checkParam("--closer_sampling_points");
    if(compute_closer_sampling_point_bool)
        FnexperimentalImages = getParam("--experimental_images");
    if (STR_EQUAL(getParam("--method"), "real_space"))
        projType = REALSPACE;
    if (STR_EQUAL(getParam("--method"), "shears"))
        projType = SHEARS;
    if (STR_EQUAL(getParam("--method"), "fourier"))
    {
        projType = FOURIER;
        paddFactor = getDoubleParam("--method", 1);
        maxFrequency = getDoubleParam("--method", 2);
        String degree = getParam("--method", 3);
        if (degree == "nearest")
            BSplineDeg = NEAREST;
        else if (degree == "linear")
            BSplineDeg = LINEAR;
        else if (degree == "bspline")
            BSplineDeg = BSPLINE3;
        else
            REPORT_ERROR(ERR_ARG_BADCMDLINE, "The interpolation kernel can be : nearest, linear, bspline");
    }

    //NOTE perturb in computed after the even sampling is computes
    //     and max tilt min tilt applied
    perturb_projection_vector=getDoubleParam("--perturb");
    compute_neighbors_bool=checkParam("--compute_neighbors");
    remove_points_far_away_from_experimental_data_bool=checkParam("--near_exp_data");
    if(remove_points_far_away_from_experimental_data_bool)
        FnexperimentalImages = getParam("--experimental_images");
    fn_groups = getParam("--groups");
    only_winner = checkParam("--only_winner");
}

/* Usage ------------------------------------------------------------------- */
void ProgAngularProjectLibrary::defineParams()
{
    addUsageLine("Create a gallery of projections from a volume");
    addUsageLine("+see [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Angular_project_library_v3][here]] for further information about this program.");

    addParamsLine("   -i <input_volume_file>       : Input Volume");
    addParamsLine("   -o <output_file_name>        : stack with output files");
    addParamsLine("  [--sym <symmetry=c1>]         : Symmetry to define sampling ");
    addParamsLine("                                : One of the 17 possible symmetries in");
    addParamsLine("                                : single particle electron microscopy");
    addParamsLine("  [--sampling_rate <Ts=5>]      : Distance in degrees between sampling points");
    addParamsLine("==+Extra parameters==");
    addParamsLine("  [--sym_neigh <symmetry>]      : symmetry used to define neighbors, by default");
    addParamsLine("                                : same as sym");
    addParamsLine("  [--psi_sampling <psi=360>]    : sampling in psi, 360 -> no sampling in psi");
    addParamsLine("  [--max_tilt_angle <tmax=180>] : maximum tilt angle in degrees");
    addParamsLine("  [--min_tilt_angle <tmin=0>]   : minimum tilt angle in degrees");
    addParamsLine("  [--experimental_images <docfile=\"\">] : doc file with experimental data");
    addParamsLine("  [--angular_distance <ang=20>]     : Do not search a distance larger than...");
    addParamsLine("  requires --experimental_images;");
    addParamsLine("  [--closer_sampling_points]    : create doc file with closest sampling points");
    addParamsLine("  requires --experimental_images;");
    addParamsLine("  [--near_exp_data]             : remove points far away from experimental data");
    addParamsLine("  requires --experimental_images;");
    addParamsLine("  [--compute_neighbors]         : create doc file with sampling point neighbors");
    addParamsLine("  requires --angular_distance;");
    addParamsLine("  [--method <method=fourier>]              : Projection method");
    addParamsLine("        where <method>");
    addParamsLine("                real_space                    : Makes projections by ray tracing in real space");
    addParamsLine("                fourier <pad=1> <maxfreq=0.25> <interp=bspline> : Takes a central slice in Fourier space");
    addParamsLine("                                              : pad controls the padding factor, by default, the padded volume is");
    addParamsLine("                                              : the same than the original volume. ");
    addParamsLine("                                              : maxfreq is the maximum frequency for the pixels and by default ");
    addParamsLine("                                              : pixels with frequency more than 0.25 are not considered.");
    addParamsLine("                                              : interp is the method for interpolation and the values can be: ");
    addParamsLine("                                              : nearest:          Nearest Neighborhood  ");
    addParamsLine("                                              : linear:           Linear  ");
    addParamsLine("                                              : bspline:          Cubic BSpline  ");
    addParamsLine("  [--perturb <sigma=0.0>]       : gaussian noise projection unit vectors ");
    addParamsLine("                                : a value=sin(sampling_rate)/4  ");
    addParamsLine("                                : may be a good starting point ");
    addParamsLine("  [--groups <selfile=\"\">]     : selfile with groups");
    addParamsLine("  [--only_winner]               : if set each experimental");
    addParamsLine("                                : point will have a unique neighbor");

    addExampleLine("Sample at 2 degrees and use c6 symmetry:", false);
    addExampleLine("xmipp_angular_project_library -i in.vol -o out.stk --sym c6 --sampling_rate 2");

}

/* Show -------------------------------------------------------------------- */
void ProgAngularProjectLibrary::show()
{
    if (!verbose)
        return;
    std::cout << "output input_volume root:  " << input_volume << std::endl
    << "output files root:         " << output_file_root << std::endl
    << "Sampling rate:             " << sampling    << std::endl
    << "symmetry sampling group:   " << fn_sym << std::endl
    << "symmetry neighbohound group: " << fn_sym_neigh << std::endl
    << "max_tilt_angle:            " << max_tilt_angle << std::endl
    << "min_tilt_angle:            " << min_tilt_angle << std::endl
    << "psi_sampling:              " << psi_sampling << std::endl
    << "compute_neighbors:         " << compute_neighbors_bool << std::endl
    << "only_winner:               " << only_winner << std::endl
    << "verbose:                   " << verbose << std::endl
    << "projection method:         ";

    if (projType == FOURIER)
    {
        std::cout << " fourier " <<std::endl;
        std::cout << "     pad factor: "   << paddFactor <<std::endl;
        std::cout << "     maxFrequency: " << maxFrequency <<std::endl;
        std::cout << "     interpolator: ";
        if (BSplineDeg == NEAREST)
            std::cout << " nearest" <<std::endl;
        else if (BSplineDeg == LINEAR)
            std::cout << " linear" <<std::endl;
        else if (BSplineDeg == BSPLINE3)
            std::cout << " bspline" <<std::endl;
    }
    else if (projType == REALSPACE)
        std::cout << " realspace " <<std::endl;

    if (angular_distance_bool)
        std::cout << "angular_distance:          " << angular_distance << std::endl;
    if (FnexperimentalImages.size() > 0)
        std::cout << "experimental_images:       " << FnexperimentalImages << std::endl;
    std::cout << "compute_closer_sampling_point_bool:" << compute_closer_sampling_point_bool
    << std::endl;
    if (perturb_projection_vector!=0)
        std::cout << "perturb_projection_vector: " << perturb_projection_vector << std::endl;
}

void ProgAngularProjectLibrary::project_angle_vector (int my_init, int my_end, bool verbose)
{
    Projection P;
    FileName fn_proj;
    double rot,tilt,psi;
    int mySize;
    int numberStepsPsi = 1;

    mySize=my_end-my_init+1;
    if (psi_sampling < 360)
    {
        numberStepsPsi = (int) (359.99999/psi_sampling);
        mySize *= numberStepsPsi;
    }

    if (verbose)
        init_progress_bar(mySize);
    int myCounter=0;


    for (double mypsi=0;mypsi<360;mypsi += psi_sampling)
        for (int i=0;i<my_init;i++)
            myCounter++;

//    if (shears && XSIZE(inputVol())!=0 && VShears==NULL)
//        VShears=new RealShearsInfo(inputVol());
    if (projType == SHEARS && XSIZE(inputVol())!=0 && Vshears==NULL)
        Vshears=new RealShearsInfo(inputVol());
    if (projType == FOURIER && XSIZE(inputVol())!=0 && Vfourier==NULL)
        Vfourier=new FourierProjector(inputVol(),
        		                      paddFactor,
        		                      maxFrequency,
        		                      BSplineDeg);

    for (double mypsi=0;mypsi<360;mypsi += psi_sampling)
    {
        for (int i=my_init;i<=my_end;i++)
        {
            if (verbose)
                progress_bar(i-my_init);
            psi= mypsi+ZZ(mysampling.no_redundant_sampling_points_angles[i]);
            tilt=      YY(mysampling.no_redundant_sampling_points_angles[i]);
            rot=       XX(mysampling.no_redundant_sampling_points_angles[i]);

//            if (shears)
//                projectVolume(*VShears, P, Ydim, Xdim, rot,tilt,psi);
//            else
//                projectVolume(inputVol(), P, Ydim, Xdim, rot,tilt,psi);
            if (projType == SHEARS)
                projectVolume(*Vshears, P, Ydim, Xdim,   rot, tilt, psi);
            else if (projType == FOURIER)
                projectVolume(*Vfourier, P, Ydim, Xdim,  rot, tilt, psi);
            else if (projType == REALSPACE)
                projectVolume(inputVol(), P, Ydim, Xdim, rot, tilt, psi);


            P.setEulerAngles(rot,tilt,psi);
            P.setDataMode(_DATA_ALL);
            P.write(output_file,(size_t) (numberStepsPsi * i + mypsi +1),true,WRITE_REPLACE);
        }
    }
    if (verbose)
        progress_bar(mySize);
}

/* Run --------------------------------------------------------------------- */
void ProgAngularProjectLibrary::run()
{
    /////////////////////////////
    // PreRun for all nodes but not for all works
    /////////////////////////////
    //only rank 0
	mysampling.verbose=verbose;
    show();
    //all ranks
    mysampling.setSampling(sampling);
    srand ( time(NULL) );
    //process the symmetry file
    //only checks symmetry and set pg_order and pg_group, no memory allocation
    if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(ERR_VALUE_INCORRECT,
                     (std::string)"Invalid symmetry" +  fn_sym);
    if(perturb_projection_vector!=0)
    {
        int my_seed;
        my_seed=rand();
        // set noise deviation and seed
        mysampling.setNoise(perturb_projection_vector,my_seed);
    }
    if(angular_distance_bool!=0)
        mysampling.setNeighborhoodRadius(angular_distance);//irrelevant
    //true -> half_sphere
    mysampling.computeSamplingPoints(false,max_tilt_angle,min_tilt_angle);
    //only rank 0
    //mysampling.createSymFile(fn_sym,symmetry, sym_order);
    //all nodes
    mysampling.SL.readSymmetryFile(fn_sym);
    //store symmetry matrices, this is faster than computing them each time
    mysampling.fillLRRepository();
    //mpi_barrier here
    //all working nodes must read symmetry file
    //and experimental docfile if apropiate
    //symmetry_file = symmetry + ".sym";
    //SL.readSymmetryFile(symmetry_file)
    // We first sample The  whole sphere
    // Then we remove point redundant due to sampling symmetry
    // use old symmetry, this is geometric does not use L_R
    mysampling.removeRedundantPoints(symmetry, sym_order);

    //=========================
    //======================
    //recompute symmetry with neigh symmetry
    // If uncomment neighbour are not OK. BE CAREFUL
#define BREAKSIMMETRY
#ifdef BREAKSIMMETRY
    if (!mysampling.SL.isSymmetryGroup(fn_sym_neigh, symmetry, sym_order))
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         (std::string)"Invalid neig symmetry" +  fn_sym_neigh);
        mysampling.SL.readSymmetryFile(fn_sym_neigh);
        mysampling.fillLRRepository();
#endif
#undef BREAKSIMMETRY
        //precompute product between symmetry matrices and experimental data
    if (FnexperimentalImages.size() > 0)
        mysampling.fillExpDataProjectionDirectionByLR(FnexperimentalImages);

    //remove points not close to experimental points, only for no symmetric cases
    if (FnexperimentalImages.size() > 0 &&
        remove_points_far_away_from_experimental_data_bool)
    {
        // here we remove points no close to experimental data, neight symmetry must be use
        mysampling.removePointsFarAwayFromExperimentalData();
    }
    if(compute_closer_sampling_point_bool)
    {
        //find sampling point closer to experimental point (only 0) and bool
        //and save docfile with this information
        // use neight symmetry
        mysampling.findClosestSamplingPoint(FnexperimentalImages,output_file_root);
    }
    //only rank 0
    //write docfile with vectors and angles
    mysampling.createAsymUnitFile(output_file_root);
    //all nodes
    //If there is no reference available exit
    try
    {
        inputVol.read(input_volume);
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        exit(0);
    }
    inputVol().setXmippOrigin();
    Xdim = XSIZE(inputVol());
    Ydim = YSIZE(inputVol());

    if (compute_neighbors_bool)
    {
        // new symmetry
        mysampling.computeNeighbors(only_winner);
        mysampling.saveSamplingFile(output_file_root,false);
    }
    //release some memory
    mysampling.exp_data_projection_direction_by_L_R.clear();

    unlink(output_file.c_str());

    //mpi master should divide doc in chuncks
    //in this serial program there is a unique chunck
    //angle information is in
    //mysampling.no_redundant_sampling_points_vector[i]
    //Run for all works
    project_angle_vector(0,
                         mysampling.no_redundant_sampling_points_angles.size()-1,verbose);

    //only rank 0 create sel file
    MetaData  mySFin, mySFout;
    FileName fn_temp;
    mySFin.read(output_file_root+"_angles.doc");
    size_t myCounter=0;
    size_t id;
    int ref;
    for (double mypsi=0;mypsi<360;mypsi += psi_sampling)
    {
        FOR_ALL_OBJECTS_IN_METADATA(mySFin)
        {
            double x,y,z, rot, tilt, psi;
            mySFin.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
            mySFin.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
            mySFin.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
            mySFin.getValue(MDL_X,x,__iter.objId);
            mySFin.getValue(MDL_Y,y,__iter.objId);
            mySFin.getValue(MDL_Z,z,__iter.objId);
            mySFin.getValue(MDL_REF,ref,__iter.objId);
            fn_temp.compose( ++myCounter,output_file);
            id = mySFout.addObject();
            mySFout.setValue(MDL_IMAGE,fn_temp,id);
            mySFout.setValue(MDL_ENABLED,1,id);
            mySFout.setValue(MDL_ANGLE_ROT,rot,id);
            mySFout.setValue(MDL_ANGLE_TILT,tilt,id);
            mySFout.setValue(MDL_ANGLE_PSI,psi+mypsi,id);
            mySFout.setValue(MDL_X,x,id);
            mySFout.setValue(MDL_Y,y,id);
            mySFout.setValue(MDL_Z,z,id);
            mySFout.setValue(MDL_SCALE,1.0,id);
            mySFout.setValue(MDL_REF,ref,id);
        }
    }
    mySFout.setComment("x,y,z refer to the coordinates of the unitary vector at direction given by the euler angles");
    mySFout.write(output_file_root+".doc");
    unlink((output_file_root+"_angles.doc").c_str());

    if (fn_groups!="")
        createGroupSamplingFiles();
}

void ProgAngularProjectLibrary::createGroupSamplingFiles(void)
{

    //#define DEBUGTIME
#ifdef  DEBUGTIME
    time_t start,end;
    double time_dif;
    time (&start);
#endif

    //load txt file
    mysampling.readSamplingFile(output_file_root,false);
#ifdef  DEBUGTIME

    time (&end);
    time_dif = difftime (end,start);
    start=end;
    printf ("re-read entire sampling file after %.2lf seconds\n", time_dif );
#endif

    StringVector blockList;
    getBlocksInMetaDataFile(fn_groups,blockList);
    FileName fn_temp, fn_exp;
    FileName my_output_file_root;
    MetaData SFBlock;

    fn_exp = FnexperimentalImages.removeBlockName();
    int igrp=1;
    for (StringVector::iterator it= blockList.begin();
         it!=blockList.end(); it++,igrp++)
    {
        my_output_file_root.compose(output_file_root + "_group",igrp,"");
        std::cerr<<"Writing group sampling file "<< my_output_file_root<<std::endl;

        fn_temp.compose(*it,fn_exp);
        SFBlock.read(fn_temp);
        if (SFBlock.size() > 0)//Do we really need this check?
            //I guess so since user may have supplied a particular
            //defocus classification. ROB
        {
            mysampling.fillExpDataProjectionDirectionByLR(fn_temp);//SFBlock@fn_groups
            if(compute_closer_sampling_point_bool)
            {
                //find sampling point closer to experimental point (only 0) and bool
                //and save docfile with this information
                mysampling.findClosestSamplingPoint(fn_temp,my_output_file_root);
            }

            //save saveSamplingFile
            if (compute_neighbors_bool)
            {
                mysampling.computeNeighbors(only_winner);
                mysampling.saveSamplingFile(my_output_file_root,false);
            }
        }
    }
#ifdef  DEBUGTIME
    time (&end);
    time_dif = difftime (end,start);
    start=end;
    printf ("Written all group sampling files after %.2lf seconds\n", time_dif );
#endif


}
