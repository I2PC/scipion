/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
#ifndef _XMIPPPROG_HH
#   define _XMIPPPROG_HH

#include "xmippImages.hh"
#include "xmippVolumes.hh"

/**@name Programs
    Here you will find some routines which might help you to write programs
    with a common structure where only a small part changes. These routines
    make the heavy work for you and you only have to supply the changing
    part from one to another.
*/
//@{
/** Virtual Program parameters class.
    Program parameters must be a class inheriting from this class. 
    Here you are an example \Ref{SF_main}*/
class Prog_parameters {
public:
   /// Output extension
   FileName oext;
   /// Output root
   FileName oroot;
   /// Output file
   FileName fn_out;
   /// Input file
   FileName fn_in;
   /// This flag for application of the transformation as stored in the header
   bool     apply_geo;
   /// Use this flag for not writing at every image
   bool     each_image_produces_an_output;
   /// Use this flag for not producing a time bar
   bool     allow_time_bar;
public:
   /// Empty constructor
   Prog_parameters() {oroot=oext=fn_out=fn_in="";
      each_image_produces_an_output=TRUE; allow_time_bar=TRUE; apply_geo=FALSE;}
   /// Read the basic parameters defined for this class
   virtual void read(int argc, char **argv);
   /// Show these parameters
   virtual void show();
   /// Usage
   virtual void usage();
   /// Get input size
   void get_input_size(int &Zdim, int &Ydim, int &Xdim);
   /// Get the number of images to process
   int get_images_to_process();
};

#define IMAGE2IMAGE     1
#define FOURIER2FOURIER 2
#define IMAGE2FOURIER   3
#define FOURIER2IMAGE   4
#define IMAGE2FILE      5
#define FILE2FILE       6

/** Main for SelFiles and individual images without parameters.
    This routine implements the main of a program which reads either an
    image/volume or a selfile, processes each image/volume separately
    (this is the part
    you have to supply, together with the extra parameters needed by the
    program), for instance, add some noise,
    and then save the results. The output images can be the same input ones,
    or new ones (in case of an image, the user must supply the output name,
    and for selection files, he must supply an output extension). If an error
    occurs, the routine shows the error message and exits 1, ie, the program
    is aborted.
    
    Here you have an example of programming and calling use.
    \begin{verbatim}
add_noise.cc:    
    
    #include <XmippData/xmippTypes.hh> 

    class Add_noise_parameters: public Prog_parameters {
    public:
       double noise_avg;
       double noise_stddev;
       void read(int argc, char **argv)
          {noise_avg=AtoF(get_param(argc,argv,"-noise_avg","0"));
           noise_stddev=AtoF(get_param(argc,argv,"-noise_stddev"));}
      void show()
          {cout << "Noise avg=" << noise_avg << endl
                << "Noise stddev=" << noise_stddev << endl;}
      void usage()
         {cerr << "   -noise_avg <avg>         : Gaussian noise average\n"
               << "   -noise_stddev <stddev>   : Gaussian standard deviation\n";}
    };
    
    bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
      Add_noise_parameters *eprm=(Add_noise_parameters *) prm;
      img().add_noise(eprm->noise_avg, eprm->noise_stddev);
      return true; // Return true if the image is sucessfully treated
                // and false if it is not. This information is used
                // to create the output selfile
    }

    bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
      Add_noise_parameters *eprm=(Add_noise_parameters *) prm;
      vol().add_noise(eprm->noise_avg, eprm->noise_stddev);
      return true;
    }

    void main (int argc, char **argv) {
       Add_noise_parameters prm;
       SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
    }
    \end{verbatim}
    
    You would have to compile it
    \begin{verbatim}
    XmippCompile add_noise.cc
    \end{verbatim}
    
    And now use it!!
    \begin{verbatim}
    1: add_noise -i g0ta00001.xmp                   -noise_stddev 0.1
    2: add_noise -i g0ta00001.xmp -o g0ta00001.pmx  -noise_stddev 0.1
    3: add_noise -i g0t.sel                         -noise_stddev 0.1
    4: add_noise -i g0t.sel       -oext pmx         -noise_stddev 0.1
    \end{verbatim}    

    Example 1 changes the angles of a single image and store the result
    in the same image, while example 2 stores it in a different one.
    Example 3 changes the angles of a set of images modifying them.
    Example 4 creates a new set of images with extension .pmx. Additionally
    a selection file is created called g0tpmx.sel

    The operation mode determines how the program works. There are five
    valid modes: IMAGE2IMAGE, FOURIER2FOURIER, IMAGE2FOURIER, FOURIER2IMAGE,
    IMAGE2FILE, FILE2FILE.
    
    If the input and output types are the same the image is rewritten in
    the input argument of the function
    \\ Ex: process_img(ImageXmipp &i, Prog_parameters *prm)
    \\Otherwise, there must be different input (1st) and output arguments (2nd)
    \\ Ex: process_img(FourierImageXmipp &input, ImageXmipp &output,
          Prog_parameters *prm)
    \\ Ex: process_img(ImageXmipp &input, const FileName &fn_out,
          Prog_parameters *prm)

    If the program mode is FILE2FILE then the process function must be of
    the form: bool process(const FileName &fn_in, const FileName &fn_out, 
    const Prog_parameters *prm); and it must return true if the file
    was successfully processed. This function must be provided as the
    process_img function to SF_main, the process_vol can be set to NULL.
*/
void SF_main(int argc, char **argv, Prog_parameters *prm,
   void *process_img,
   void *process_vol, int operation_mode=IMAGE2IMAGE);
//@}
#endif
