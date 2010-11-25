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
    mysampling.SetSampling(1);
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
    shears = checkParam("--shears");
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
    addUsageLine("Example of use: Sample at 2 degrees and use c6 symmetry");
    addUsageLine("   xmipp_angular_project_library -i in.vol -o out -sym c6 -sampling_rate 2");
    addParamsLine("   -i <input_volume_file>           : Input Volume");
    addParamsLine("   -o <root_file_name>         : Root for output files");
    addParamsLine("  [--sym <symmetry=c1>]         : Symmetry to define sampling ");
    addParamsLine("                               : One of the 17 possible symmetries in");
    addParamsLine("                               : single particle electron microscopy");
    addParamsLine("  [--sampling_rate <Ts=5>]          : Distance in degrees between sampling points");
    addParamsLine("==+Extra parameters==");
    addParamsLine("  [--sym_neigh <symmetry>]      : symmetry used to define neighbors, by default");
    addParamsLine("                               : same as sym");
    addParamsLine("  [--psi_sampling <psi=360>]        : sampling in psi, 360 -> no sampling in psi");
    addParamsLine("  [--max_tilt_angle <tmax=91>]      : maximum tilt angle in degrees");
    addParamsLine("  [--min_tilt_angle <tmin=-91>]     : minimum tilt angle in degrees");
    addParamsLine("  [--experimental_images <docfile=\"\">] : doc file with experimental data");
    addParamsLine("  [--angular_distance <ang=20>]     : Do not search a distance larger than...");
    addParamsLine("  requires --experimental_images;");
    addParamsLine("  [--closer_sampling_points]    : create doc file with closest sampling points");
    addParamsLine("  requires --experimental_images;");
    addParamsLine("  [--near_exp_data]             : remove points far away from experimental data");
    addParamsLine("  requires --experimental_images;");
    addParamsLine("  [--compute_neighbors]         : create doc file with sampling point neighbors");
    addParamsLine("  requires --angular_distance;");
    addParamsLine("  [--shears]                    : use projection shears to generate projections");
    addParamsLine("  [--perturb <sigma=0.0>]       : gaussian noise projection unit vectors ");
    addParamsLine("                         : a value=sin(sampling_rate)/4  ");
    addParamsLine("                         : may be a good starting point ");
    addParamsLine("  [--groups <selfile=\"\">]     : selfile with groups");
    addParamsLine("  [--only_winner]               : if set each experimental");
    addParamsLine("                               : point will have a unique neighbor");
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
    << "shears:                    " << shears << std::endl
    ;
    if (angular_distance_bool)
        std::cout << "angular_distance:          " << angular_distance << std::endl;
    if (FnexperimentalImages.size() > 0)
        std::cout << "experimental_images:       " << FnexperimentalImages << std::endl;
    std::cout << "compute_closer_sampling_point_bool:" << compute_closer_sampling_point_bool
    << std::endl;
    if (perturb_projection_vector!=0)
        std::cout << "perturb_projection_vector: " << perturb_projection_vector << std::endl;
}

void
ProgAngularProjectLibrary::project_angle_vector(
    int my_init, int my_end, bool verbose)
{
    Projection P;
    FileName fn_proj;
    double rot,tilt,psi;
    int mySize;
    mySize=my_end-my_init+1;
    if (psi_sampling < 360)
        mySize *= (int) (359.99999/psi_sampling);
    if (verbose)
        init_progress_bar(mySize);
    int myCounter=0;


    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
        for (int i=0;i<my_init;i++)
            myCounter++;

    if (shears && XSIZE(inputVol())!=0)
    {
        prepareStructVolume(inputVol(),VShears);
        inputVol.clear();
    }
    //fn_proj.compose(output_file_root, ++myCounter,"xmp");

    unlink(output_file.c_str());

    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
    {
        for (int i=my_init;i<=my_end;i++)
        {
            if (verbose)
                progress_bar(i-my_init);
            psi= mypsi+ZZ(mysampling.no_redundant_sampling_points_angles[i]);
            tilt=      YY(mysampling.no_redundant_sampling_points_angles[i]);
            rot=       XX(mysampling.no_redundant_sampling_points_angles[i]);

            if (shears)
                project_Volume(VShears, P, Ydim, Xdim,rot,tilt,psi);
            else
                project_Volume(inputVol(), P, Ydim, Xdim,rot,tilt,psi);

            //fn_proj.compose(output_file_root, ++myCounter,"xmp");
            P.write(output_file,-1,true,WRITE_APPEND);
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
    //all ranks
    //only rank 0

    show();
    //all ranks
    mysampling.SetSampling(sampling);
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
        mysampling.SetNoise(perturb_projection_vector,my_seed);
    }
    if(angular_distance_bool!=0)
        mysampling.SetNeighborhoodRadius(angular_distance);//irelevant
    //true -> half_sphere
    mysampling.Compute_sampling_points(false,max_tilt_angle,min_tilt_angle);
    //only rank 0
    //mysampling.create_sym_file(fn_sym,symmetry, sym_order);
    //all nodes
    mysampling.SL.read_sym_file(fn_sym);
    //store symmetry matrices, this is faster than computing them each time
    mysampling.fill_L_R_repository();
    //mpi_barrier here
    //all working nodes must read symmetry file
    //and experimental docfile if apropiate
    //symmetry_file = symmetry + ".sym";
    //SL.read_sym_file(symmetry_file)
    // We first sample The  whole sphere
    // Then we remove point redundant due to sampling symmetry
    // use old symmetry, this is geometric does not use L_R
    mysampling.remove_redundant_points(symmetry, sym_order);

    //=========================
    //======================
    //recompute symmetry with neigh symmetry
    if (!mysampling.SL.isSymmetryGroup(fn_sym_neigh, symmetry, sym_order))
        REPORT_ERROR(ERR_VALUE_INCORRECT,
                     (std::string)"Invalid neig symmetry" +  fn_sym_neigh);
    mysampling.SL.read_sym_file(fn_sym_neigh);
    mysampling.fill_L_R_repository();
    //precompute product between symmetry matrices and experimental data
    if (FnexperimentalImages.size() > 0)
        mysampling.fill_exp_data_projection_direction_by_L_R(FnexperimentalImages);

    //remove points not close to experimental points, only for no symmetric cases
    if (FnexperimentalImages.size() > 0 &&
        remove_points_far_away_from_experimental_data_bool)
    {
        // here we remove points no close to experimental data, neight symmetry must be use
        mysampling.remove_points_far_away_from_experimental_data();
    }
    if(compute_closer_sampling_point_bool)
    {
        //find sampling point closer to experimental point (only 0) and bool
        //and save docfile with this information
        // use neight symmetry
        mysampling.find_closest_sampling_point(FnexperimentalImages,output_file_root);
    }
    //only rank 0
    //write docfil with vectors and angles
    mysampling.create_asym_unit_file(output_file_root);
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
        mysampling.compute_neighbors(only_winner);
        mysampling.save_sampling_file(output_file_root,false);
    }
    //release some memory
    mysampling.exp_data_projection_direction_by_L_R.clear();
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
    int myCounter=-1;
    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
    	FOR_ALL_OBJECTS_IN_METADATA(mySFin)
    	{
            double x,y,z, rot, tilt, psi;
            mySFin.getValue(MDL_ANGLEROT,rot);
            mySFin.getValue(MDL_ANGLETILT,tilt);
            mySFin.getValue(MDL_ANGLEPSI,psi);
            mySFin.getValue(MDL_X,x);
            mySFin.getValue(MDL_Y,y);
            mySFin.getValue(MDL_Z,z);
            fn_temp.compose( ++myCounter,output_file);
            mySFout.addObject();
            mySFout.setValue(MDL_IMAGE,fn_temp);
            mySFout.setValue(MDL_ENABLED,1);
            mySFout.setValue(MDL_ANGLEROT,rot);
            mySFout.setValue(MDL_ANGLETILT,tilt);
            mySFout.setValue(MDL_ANGLEPSI,psi+mypsi);
            mySFout.setValue(MDL_X,x);
            mySFout.setValue(MDL_Y,y);
            mySFout.setValue(MDL_Z,z);
        }
    mySFout.write(output_file_root+".doc");
    unlink((output_file_root+"_angles.doc").c_str());
}

void ProgAngularProjectLibrary::createGroupSamplingFiles(void)
{

    //#define DEBUGTIME
#ifdef  DEBUGTIME
    #include <ctime>

    time_t start,end;
    double time_dif;
    time (&start);
#endif

    //load txt file
    mysampling.read_sampling_file(output_file_root,false);
#ifdef  DEBUGTIME

    time (&end);
    time_dif = difftime (end,start);
    start=end;
    printf ("re-read entire sampling file after %.2lf seconds\n", time_dif );
#endif

    MetaData mySF(fn_groups);
    FileName fn_temp;
    FileName my_output_file_root;
    int igrp = 0;
    FOR_ALL_OBJECTS_IN_METADATA(mySF)
    {
        igrp++;
        mySF.getValue(MDL_IMAGE,fn_temp);
        if (fn_temp=="")
            break;
        my_output_file_root.compose(output_file_root + "_group",igrp,"");
        std::cerr<<"Writing group sampling file "<< my_output_file_root<<std::endl;

        if (fn_temp.size() > 0)
        {
            mysampling.fill_exp_data_projection_direction_by_L_R(fn_temp);
            if(compute_closer_sampling_point_bool)
            {
                //find sampling point closer to experimental point (only 0) and bool
                //and save docfile with this information
                mysampling.find_closest_sampling_point(fn_temp,my_output_file_root);
            }

            //save save_sampling_file
            if (compute_neighbors_bool)
            {
                mysampling.compute_neighbors(only_winner);
                mysampling.save_sampling_file(my_output_file_root,false);
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
