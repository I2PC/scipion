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
    symmetry = getParameter(argc, argv, "-symmetry","cn");
    sym_order = textToInteger(getParameter(argc, argv, "-sym_order", "1"));
    sampling = textToFloat(getParameter(argc, argv, "-sampling_rate", "5"));
    psi_sampling = textToFloat(getParameter(argc, argv, "-psi_sampling", "360"));
    max_tilt_angle = textToFloat(getParameter(argc, argv, "-max_tilt_angle","91"));
    min_tilt_angle = textToFloat(getParameter(argc, argv, "-min_tilt_angle","-91"));
    angular_distance_bool = checkParameter(argc, argv,"-angular_distance");
    if(angular_distance_bool)
    {
         if(!checkParameter(argc, argv, "-experimental_images"))
         {
             cerr << "Docfile with experimental images euler angles is missing" << endl;
             exit(0);
         } 
            FnexperimentalImages = getParameter(argc, argv, "-experimental_images","");
            angular_distance = textToFloat(getParameter(argc, argv,"-angular_distance"));
    }
}

/* Usage ------------------------------------------------------------------- */
void Prog_create_projection_library_Parameters::usage()
{
    cerr << "create_projection_library\n"
    << "   -i input_volume             : Input Volume\n"
    << "   -o root_file_name           : Root for output files\n"
    << "  [-symmetry cn]   :One of the 17 possible symmetries in\n"
    << "                                single particle electronmicroscopy\n"
    << "                                i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, dnh, t, td, th, o, oh, i, ih\n"
    << "  [-sym_order 1]               : For infinite groups symmetry order\n"
    << "  [-sampling_rate 5]           : Distance in degrees between sampling points\n"
    << "  [-psi_sampling 360]          : sampling in psi, 360 -> no sampling in psi\n"
    << "  [-max_tilt_angle  91]        : maximum tilt angle in degrees\n"
    << "  [-min_tilt_angle -91]        : minimum tilt angle in degrees\n"
    << "  [-experimental_images \"\"]  : doc file with experimental data\n"
    << "  [-angular_distance 20]       : do not search a distance larger than...\n"
    << "\n"
    << "Example of use: Sample at 2 degrees and use c6 symmetry\n"
    << "   xmipp_create_projection_library -i in.vol -o out "
    << "    -symmetry cn -sym_order 6 -sampling_rate 2\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_create_projection_library_Parameters::show()
{
    
    cout
    << "output input_volume root:  " << input_volume << endl
    << "output files root:         " << output_file_root << endl
    << "Sampling rate:             " << sampling    << endl
    << "symmetry group:            " << symmetry << endl
    << "symmetry order:            " << sym_order << endl
    << "max_tilt_angle:            " << max_tilt_angle << endl
    << "min_tilt_angle:            " << min_tilt_angle << endl
    << "psi_sampling:              " << psi_sampling << endl
    ;
    if(angular_distance_bool)
        cout << "angular_distance:          " << angular_distance << endl
             << "experimental_images:       " << FnexperimentalImages << endl
        ;
}



void
Prog_create_projection_library_Parameters::project_angle_vector(int
my_init, int my_end, bool verbose)
{
    Projection P;
    FileName fn_proj;
    double rot,tilt,psi;
    int mySize;
    mySize=my_end-my_init+1;
    if(psi_sampling < 360)
       mySize *= (int) (359.99999/psi_sampling);
    if(verbose)
       init_progress_bar(mySize);
    int myCounter=0;
    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
       for (int i=0;i<my_init;i++)
         myCounter++;
         
    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
    {
       for (int i=my_init;i<=my_end;i++)
       {    
           if(verbose)
               progress_bar(i-my_init);
           psi= mypsi+ZZ(close_points_angles[i]);
           tilt=      YY(close_points_angles[i]);
           rot=       XX(close_points_angles[i]);

    //progress bar     
           project_Volume(inputVol(), P, Ydim, Xdim,rot,tilt,psi);

           fn_proj.compose(output_file_root, myCounter++,"xmp");
           P.write(fn_proj);
       }
    }
    if(verbose)
        progress_bar(mySize);

}

#ifdef NEVEDEFINED
//UNTESTED CODE
void
Prog_create_projection_library_Parameters::createDocFileWithExperimentalImages(void)
{   
    SelFile         mySF
    DocFile         myDF;
    Matrix1D<double> docline;
    FileName        fn_img,fn_doc;
    headerXmipp     nyhead;

    if( (symmetry!="cn" || sym_order!=1) && () )
    {
        cerr << "experimental_images option only available for symmetry "
             << "cn and order 1" 
             <<endl;
        exit(0);     
    }  
    
    if(SFexperimentalImages!="")
    {
        docline.initZeros(3);

        mySF.read(DFexperimentalImages);
        mySF.go_beginning(); 

        myDF.reserve(mySF.ImgNo());
        myDF.append_comment("Headerinfo columns: rot (1) , tilt (2), psi (3)");

        while (!mySF.eof())
        {
            fn_img = mySF.NextImg();
            nyhead.read(fn_img);
            docline(0) = head.Phi();
            docline(1) = head.Theta();
            docline(2) = head.Psi();
            myDF.append_comment(fn_img);
            myDF.append_data_line(docline);
        }
        fn_doc = DFexperimentalImages
        fn_doc.remove_extension("sel");   
        fn_doc.add_extension("exp_doc");
        myDF.write(.substitute_extension("sel", "doc"));
    }
}
#endif

void
Prog_create_projection_library_Parameters::remove_points_not_close_to_experimental_points(void)
{   
    if((symmetry!="cn" || sym_order!=1) && angular_distance_bool==true)
    {
        cerr << "experimental_images option only available for symmetry "
             << "cn and order 1" 
             <<endl;
        exit(0);     
    }  
    //return new docfile if not restrictions
    if(!angular_distance_bool)
    {
        for (int i = 0; i < mysampling.no_redundant_sampling_points_angles.size(); i++)
        {
         close_points_angles.push_back(mysampling.no_redundant_sampling_points_angles[i]);
        }    
    }
    else
    {

        DocFile          DFi;
        //init both df
        DFi.read(FnexperimentalImages);//experiemntal points
        double ref_tilt,ref_rot;
        double img_tilt,img_rot;
        for (int i = 0; i < mysampling.no_redundant_sampling_points_angles.size(); i++)
        {
            ref_rot  = XX(mysampling.no_redundant_sampling_points_angles[i]);
            ref_tilt = YY(mysampling.no_redundant_sampling_points_angles[i]);
            // act_rot_range is tilt-angle dependent!
            DFi.go_first_data_line();//experimental data
            while (!DFi.eof())
            {
                img_rot=DFi(0);
                img_tilt=DFi(1);
                // I knw this comparison is not accurate
                //but is fast and we do not need much precission
                if ( ABS(realWRAP(img_rot-ref_rot,-180.,180.)) <= angular_distance &&
                     ABS(realWRAP(img_tilt-ref_tilt,-180.,180.)) <= angular_distance )
                {
                     close_points_angles.push_back(mysampling.sampling_points_angles[i]);
                     break;
                }
            DFi.next();
            }
         }      
    }                     
}
/* Run --------------------------------------------------------------------- */
void Prog_create_projection_library_Parameters::run()
{
    /////////////////////////////
    // PreRun for all nodes but not for all works
    /////////////////////////////
    //only rank 0
    show();
    //all ranks
    mysampling.SetSampling(sampling);
    //mysampling.SetNeighborhoodRadius(0.);//irelevant
    //true -> half_sphere
    mysampling.Compute_sampling_points(false,max_tilt_angle,min_tilt_angle);
    //only rank 0
    mysampling.create_sym_file(symmetry, sym_order);
    //mpi_barrier here
    //all working nodes must read symmetry file
    //and experimental docfile if apropiate
    //symmetry_file = symmetry + ".sym";
    //SL.read_sym_file(symmetry_file)
    mysampling.remove_redundant_points(symmetry, sym_order);
    //remove points not close to experimental points, only for no symmetric cases
    //if(DFexperimentalImages!="")
    remove_points_not_close_to_experimental_points();
    //only rank 0
    mysampling.create_asym_unit_file(output_file_root);
    //all nodes
    inputVol.read(input_volume);
    inputVol().setXmippOrigin();
    Xdim = XSIZE(inputVol());
    Ydim = YSIZE(inputVol());
    //mysampling.compute_neighbors();
    
    //mpi maaster should divide doc in chuncks
    //in this serial program there is a unique chunck
    //angle information is in
    //mysampling.no_redundant_sampling_points_vector[i]

    //Run for all works
    project_angle_vector(0,
                 close_points_angles.size()-1);
                 
    //only rank 0 create sel file
    SelFile  mySF;
    FileName fn_temp;
    int myCounter=0;
    
    for (int mypsi=0;mypsi<360;mypsi += psi_sampling)
       for (int i=0;i<=close_points_angles.size()-1;i++)
       {    
        fn_temp.compose(output_file_root, myCounter++,"xmp");
        mySF.insert(fn_temp);
       }
    fn_temp=output_file_root+".sel";   
    mySF.write(fn_temp);         
}

