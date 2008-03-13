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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "create_projection_library.h"

/* Empty constructor ------------------------------------------------------- */
Prog_create_projection_library_Parameters::Prog_create_projection_library_Parameters()
{
    /** sampling object 1 by default*/
    mysampling.SetSampling(1);
}


/* Read parameters --------------------------------------------------------- */
void Prog_create_projection_library_Parameters::read(int argc, char **argv)
{
    input_volume = getParameter(argc, argv, "-i");
    output_file_root = getParameter(argc, argv, "-o");
    fn_sym = getParameter(argc, argv, "-sym","c1");
//    symmetry = getParameter(argc, argv, "-sym","c1");
//    sym_order = textToInteger(getParameter(argc, argv, "-sym_order", "1"));
    sampling = textToFloat(getParameter(argc, argv, "-sampling_rate", "5"));
    psi_sampling = textToFloat(getParameter(argc, argv, "-psi_sampling", "360"));
    max_tilt_angle = textToFloat(getParameter(argc, argv, "-max_tilt_angle","91"));
    min_tilt_angle = textToFloat(getParameter(argc, argv, "-min_tilt_angle","-91"));
    angular_distance_bool = checkParameter(argc, argv,"-angular_distance");
    angular_distance=0.;
    if(angular_distance_bool)
    {
         if(!checkParameter(argc, argv, "-experimental_images"))
         {
             std::cerr << "Docfile with experimental images euler angles is missing" << std::endl;
             exit(0);
         } 
            FnexperimentalImages = getParameter(argc, argv, "-experimental_images","");
            angular_distance = textToFloat(getParameter(argc, argv,"-angular_distance"));
    }
    compute_closer_sampling_point_bool= checkParameter(argc, argv,"-closer_sampling_points");
    if(compute_closer_sampling_point_bool)
    {
         if(!checkParameter(argc, argv, "-experimental_images"))
         {
             std::cerr << "Docfile with experimental images euler angles is missing" << std::endl;
             exit(0);
         } 
            FnexperimentalImages = getParameter(argc, argv, "-experimental_images","");
    }
    quiet = checkParameter(argc, argv,"-quiet");
    //NOTE perturb in computed after the even sampling is computes
    //     and max tilt min tilt applied
    perturb_projection_vector=textToFloat(getParameter(argc,argv,"-perturb","0"));       
    compute_neighbors_bool=checkParameter(argc, argv,"-compute_neighbors");
    remove_points_far_away_from_experimental_data_bool=
                   checkParameter(argc, argv,"-near_exp_data");
    if(remove_points_far_away_from_experimental_data_bool)
    {
         if(!checkParameter(argc, argv, "-experimental_images"))
         {
             std::cerr << "Docfile with experimental images euler angles is missing" << std::endl;
             exit(0);
         } 
            FnexperimentalImages = getParameter(argc, argv, "-experimental_images","");
    }
    if (angular_distance_bool==false && compute_neighbors_bool==true)
         {
             std::cerr << "If -compute_neighbors requires -angular_distance" << std::endl;
             exit(0);
         } 
        
}

/* Usage ------------------------------------------------------------------- */
void Prog_create_projection_library_Parameters::usage()
{
    std::cerr << "create_projection_library\n"
    << "   -i input_volume             : Input Volume\n"
    << "   -o root_file_name           : Root for output files\n"
    << "  [-sym cn]   :One of the 17 possible symmetries in\n"
    << "                                single particle electronmicroscopy\n"
    << "                                i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, dnh, t, td, th, o, oh, i, ih\n"
    << "                               : where n may change from 1 to 99\n"
    << "  [-sampling_rate 5]           : Distance in degrees between sampling points\n"
    << "  [-psi_sampling 360]          : sampling in psi, 360 -> no sampling in psi\n"
    << "  [-max_tilt_angle  91]        : maximum tilt angle in degrees\n"
    << "  [-min_tilt_angle -91]        : minimum tilt angle in degrees\n"
    << "  [-experimental_images \"\"]  : doc file with experimental data\n"
    << "  [-angular_distance 20]       : do not search a distance larger than...\n"
    << "  [-closer_sampling_points]    : create doc file with closest sampling points\n"
    << "  [-compute_neighbors]         : create doc file with sampling point neighbors\n"
    << "  [-quiet]                     : do not show messages\n"
    << "  [-perturb default=0.0]       : gaussian noise projection unit vectors \n"
    << "			         a value=sin(sampling_rate)/4  \n"
    << "			         may be a good starting point \n"
    << "\n"
    << "Example of use: Sample at 2 degrees and use c6 symmetry\n"
    << "   xmipp_create_projection_library -i in.vol -o out "
    << "    -symmetry cn -sym_order 6 -sampling_rate 2\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_create_projection_library_Parameters::show()
{
    if (quiet) return;    
    std::cout << "output input_volume root:  " << input_volume << std::endl
              << "output files root:         " << output_file_root << std::endl
              << "Sampling rate:             " << sampling    << std::endl
              << "symmetry group:            " << fn_sym << std::endl
              << "max_tilt_angle:            " << max_tilt_angle << std::endl
              << "min_tilt_angle:            " << min_tilt_angle << std::endl
              << "psi_sampling:              " << psi_sampling << std::endl
              << "compute_neighbors:         " << compute_neighbors_bool << std::endl
              << "quiet:                     " << quiet << std::endl
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
Prog_create_projection_library_Parameters::project_angle_vector(
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
        
    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
    {
       for (int i=my_init;i<=my_end;i++)
       {    
           if (verbose)
               progress_bar(i-my_init);
           psi= mypsi+ZZ(mysampling.no_redundant_sampling_points_angles[i]);
           tilt=      YY(mysampling.no_redundant_sampling_points_angles[i]);
           rot=       XX(mysampling.no_redundant_sampling_points_angles[i]);

           project_Volume(inputVol(), P, Ydim, Xdim,rot,tilt,psi);

           fn_proj.compose(output_file_root, myCounter++,"xmp");
           P.write(fn_proj);
       }
    }
    if (verbose)
        progress_bar(mySize);

}



/* Run --------------------------------------------------------------------- */
void Prog_create_projection_library_Parameters::run()
{ 
    #define DEBUGTIME
    #ifdef  DEBUGTIME
    #include <ctime>
    
    time_t start,end;
    double time_dif;
    time (&start);

    #endif
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
    if (!mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
         REPORT_ERROR(3005, (std::string)"create_projection_library::run Invalid symmetry" +  fn_sym);
    if(perturb_projection_vector!=0)
        {
        int my_seed;
        my_seed=rand();
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
    //mpi_barrier here
    //all working nodes must read symmetry file
    //and experimental docfile if apropiate
    //symmetry_file = symmetry + ".sym";
    //SL.read_sym_file(symmetry_file)
    mysampling.remove_redundant_points(symmetry, sym_order);
    //remove points not close to experimental points, only for no symmetric cases
    #ifdef  DEBUGTIME
    time (&end);
    time_dif = difftime (end,start); start=end;
    printf ("remove_redundant_points after %.2lf seconds\n", time_dif );
    #endif
    if (FnexperimentalImages.size() > 0 && 
        remove_points_far_away_from_experimental_data_bool)
        {	
        mysampling.remove_points_far_away_from_experimental_data(FnexperimentalImages);
        #ifdef  DEBUGTIME
        time (&end);
        time_dif = difftime (end,start); start=end;
        printf ("remove_points_far_away_from_experimental_data after %.2lf seconds\n", time_dif );
        #endif
        }
    if(compute_closer_sampling_point_bool)
	    {
	    //find sampling point closer to experimental point (only 0) and bool
	    //and save docfile with this information
	    mysampling.find_closest_sampling_point(FnexperimentalImages,output_file_root);
        #ifdef  DEBUGTIME
        time (&end);
        time_dif = difftime (end,start); start=end;
        printf ("find_closest_sampling_point after %.2lf seconds\n", time_dif );
        #endif
        }
    //only rank 0
    mysampling.create_asym_unit_file(output_file_root);
    //all nodes
    inputVol.read(input_volume);
    inputVol().setXmippOrigin();
    Xdim = XSIZE(inputVol());
    Ydim = YSIZE(inputVol());
    #ifdef  DEBUGTIME
    time (&end);
    time_dif = difftime (end,start); start=end;
    printf ("compute_neighbors before %.2lf seconds\n", time_dif );
    #endif
    if (compute_neighbors_bool)
        {
	    mysampling.compute_neighbors();
	    #ifdef  DEBUGTIME
	    time (&end);
	    time_dif = difftime (end,start); start=end;
	    printf ("compute_neighbors after %.2lf seconds\n", time_dif );
	    #endif
	    mysampling.save_sampling_file(output_file_root);
	    }
    //mpi master should divide doc in chuncks
    //in this serial program there is a unique chunck
    //angle information is in
    //mysampling.no_redundant_sampling_points_vector[i]
exit(0);
    //Run for all works
    project_angle_vector(0,
                 mysampling.no_redundant_sampling_points_angles.size()-1,!quiet);
	#ifdef  DEBUGTIME
	time (&end);
	time_dif = difftime (end,start); start=end;
	printf ("project_angle_vector after %.2lf seconds\n", time_dif );
	#endif
                 
    //only rank 0 create sel file
    SelFile  mySF;
    FileName fn_temp;
    int myCounter=0;
    
    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
       for (int i=0;i<=mysampling.no_redundant_sampling_points_angles.size()-1;i++)
       {    
        fn_temp.compose(output_file_root, myCounter++,"xmp");
        mySF.insert(fn_temp);
       }
    fn_temp=output_file_root+".sel";   
    mySF.write(fn_temp);         
}

