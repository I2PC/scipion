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
#ifndef _PROG_PROJECTION_HH
#  define _PROG_PROJECTION_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippProjection.hh>
#include "../phantom.hh"
#include "../projection.hh"
#include "Prog_project_crystal.hh"

/**@name Projection program */
//@{
/* Projection Program Parameters ------------------------------------------- */
/** Parameter class for the project program */
class Prog_Project_Parameters {
public:
   /// Filename with the \Ref{Projection_Parameters}.
   FileName fn_proj_param;
   /// Selection file with all projections
   FileName fn_sel_file;
   /** Filename with the special crystal parameters
      (\Ref{Crystal_Projection_Parameters}) */
   FileName fn_crystal;
   /// Only create angles, do not project
   bool only_create_angles;

   #define TELL_SHOW_ANGLES 0x1
   /** Debugging variable.
       This is a bitwise flag with the following valid labels:
       \\TELL_SHOW_ANGLES: the program shows the angles for each image.*/
   int tell;
public:
   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message. */
   void read(int argc, char **argv);

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage();
};

/* Angular range ----------------------------------------------------------- */
/** Angular range definition.
    This structure is used to store an angular range, ie, from a given
    ang0 to angF with a number of samples in between. The range can
    be DETERMINISTIC (all samples are equally distributed in the range),
    RANDOM (all samples are picked uniformly random in the range)
    or RANDOM BY GROUPS (for instance, given a set of 3 Euler angles,
    we may have two groups with different psi angles random by group,
    ie, the psi angle is random but inside the group the angle remains
    constant, see the following table for an example).
    \begin{verbatim}
         ROT      TILT      PSI
       -------  --------  -------
          0        0       34.56
          0       90       34.56
          0      180       34.56
          0      270       34.56

          0        0      132.87
          0       90      132.87
          0      180      132.87
          0      270      132.87
    \end{verbatim}
*/
struct Angle_range {
   /// initial angular value
   double    ang0;
   /// final angular value
   double    angF;
   /// No. of samples
   int      samples;

   #define ANGLE_RANGE_DETERMINISTIC 0
   #define ANGLE_RANGE_RANDOM_GROUPS 1
   #define ANGLE_RANGE_RANDOM        2
   #define ANGLE_EVENLY              3
   /** Kind of range.
       The kind of range can be any of these three:
       ANGLE_RANGE_DETERMINISTIC, ANGLE_RANGE_RANDOM_GROUPS,
       ANGLE_RANGE_RANDOM, ANGLE_EVENLY */
   int      randomness;
   /// Mean of the noise that must be added to the definition of the angle
   double    Navg;
   /// Stddev of the noise that must be added to the definition of the angle
   double    Ndev;
};

/* Projection parameters --------------------------------------------------- */
/** Projecting parameters.
    This class reads a set of projection parameters in a file (see
    \URL[Project]{../../../Applications/Src/Project/Help/project.html}
    for more information about the file structure) and extract the
    useful information from it. This class allows also to write in
    the same file format, for doing so you must first introduce the
    right values in the class fields, and the call to the procedure
    write */
class Projection_Parameters {
public:
   /**@name Filenames and size */
   //@{
   /** Phantom filename.
       It can be a Xmipp volume or a mathematically defined phantom. */
   FileName fn_phantom;
   /// Starting name for all projections
   string   fn_projection_seed;
   /// First projection number. By default, 1.
   int      starting;
   /// Extension for projection filenames. This is optional
   string   fn_projection_extension;
   
   /// Projection Xdim
   int      proj_Xdim;
   /// Projection Ydim
   int      proj_Ydim;
   
   /// Debugging level. See \Ref{Prog_Project_Parameters::tell}
   int tell;
   //@}

   /**@name Angle Mode: Giving ranges */
   //@{
   /// Enable angle range mode (0 or 1)
   int enable_angle_range;
   /// Rotational angle range
   Angle_range rot_range;
   /// Tilting angle range
   Angle_range tilt_range;
   /// Psi angle range
   Angle_range psi_range;
   //@}

   /**@name Angle mode: Using a Document File */
   //@{
   /// Document filename
   FileName fn_angle;
   /// First number in the document file is "rot","tilt" or "psi"
   string   ang1;
   /// Second number in the document file is "rot","tilt" or "psi"
   string   ang2;
   /// Third number in the document file is "rot","tilt" or "psi"
   string   ang3;
   //@}

   /**@name Noise.
      Noise to the pixels is introduced by adding Gaussian noise, while
      to the center is done by translating the image by a random
      quantity */
   //@{
   /// Bias to be applied to each pixel grey value */
   double    Npixel_avg;
   /// Standard deviation of the noise to be added to each pixel grey value
   double    Npixel_dev;

   /// Bias to apply to the image center
   double    Ncenter_avg;
   /// Standard deviation of the image center
   double    Ncenter_dev;
   //@}
public:
   /** From Program Parameters.
       This function loads the Projection Parameters from the parameters
       given to the program (PROJECT). */
   void from_prog_params(const Prog_Project_Parameters &prog_prm);

   /** Read projection parameters from a file.
       An exception is thrown if the file is not found or any of the
       parameters is not found in the right place.*/
   void read(FileName fn_proj_param) _THROW;

   /** Write projection parameters to a file.
       The projection parameters are written into a file wth the same
       structure as always. If the file cannot be openned for output
       an exception is thrown. */
   void write(FileName fn_proj_param) _THROW;
};

/** Project program Side information.
    This class contains side necessary information for the Project program.
    This information can be obtained from the parameters and is basically
    the Xmipp volume or phantom description plus a flag saying which
    of the two is valid. */
class PROJECT_Side_Info {
public:
   /// Document File for the projecting angles. Order: rot, tilt, psi
   DocFile        DF;     
   /// Projecting from a voxel volume?
   int            voxel_mode;
   /// Phantom Xmipp volume
   VolumeXmipp    phantom_vol;
   /// Phantom mathematical description
   Phantom        phantom_descr;
public:
   /** Produce Project Side information.
       This function produce the side information from the project
       program parameters. Basically it loads the phantom, sets
       the phantom mode to voxel or mathemattical description and
       generates or read the projection angles.
       
       Projection from voxel volumes is not implemented yet and
       an exception is thrown. */
   void produce_Side_Info(const Projection_Parameters &prm) _THROW;
};

/* Assigning angles -------------------------------------------------------- */
/** Assign angles from the projection parameters to a Document file.
    Given a set of projection parameters this function returns a 
    document file with a set of angles according to the specifications.
    The order in the output document file is rotational, tilting and
    psi angle.

    The assignment can be done from another document file (with any
    angle order) or internally generated according to the ranges defined
    in the parameters.
    
    The total number of angles is returned. The Document File is cleared,
    the first key in the document file is the starting key of the
    projection set. The current line of the document file is set to
    the beginning of the file.
*/
int PROJECT_Assign_angles(DocFile &DF, const Projection_Parameters &prm);

/* Effectively project ----------------------------------------------------- */
/** Effectively project.
    This is the routine which effectively projects, it needs the projection
    parameters and the side information, ie, the loaded phantom and list
    of angles from which project. The Projection field will keep
    at the end the last projection, this is useful in case you want
    to project only one image, although it is also written to disk.
    The returned number is the total number of projections generated.
    A selection file with all images is also returned.*/
int PROJECT_Effectively_project(const Projection_Parameters &prm,
    PROJECT_Side_Info &side, const Crystal_Projection_Parameters &prm_crystal,
    Projection &proj, SelFile &SF);

/* Main routine ------------------------------------------------------------ */
/** Main Project routine.
    Generate a set of projections given the projection parameters.
    This is the main projecting routine. This function generates a set
    of projections according to the projecting parameters defined (see
    \Ref{Prog_Project_Parameters} to know more about how to specify
    everything). The projections are written to disk.
    
    The Projection field will keep
    at the end the last projection, this is useful in case you want
    to project only one image, although it is also written to disk.
    The returned number is the total number of projections generated.
    A selection file with all images is also returned (and saved if any
    name has been given in the parameters).*/
int ROUT_project(Prog_Project_Parameters &prm, Projection &proj, SelFile &SF);
//@}
#endif
