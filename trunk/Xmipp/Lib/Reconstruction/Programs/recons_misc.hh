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
#ifndef _RECONS_MISC_HH
#  define _RECONS_MISC_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippSelFiles.hh>
#include "../symmetries.hh"
#include "../projection.hh"
#include "../grids.hh"
#include "Basic_art.hh"

/**@name Reconstruction Miscellanea */
//@{
/**@name Projection sorting */
//@{
/** Reconstruction information.
   This structure contains information for all projections which are
   going to participate in the reconstruction.
   This structure has also got
   information for the symmetry implementation. If there is any symmetry
   then an entry in this table is created using the same projection name
   but different symmetry matrices (only the matrix index is annotated
   in this structure). The Euler angles stored for the symmetrized image
   are the final ones, ie, the original Euler angles symmetrized according
   to the symmetry matrix. Then the symmetry identificator kept in this
   structure is only used to keep some track of what matrix was used to
   symmetrize.
*/
struct Recons_info {
   /// Projection filename
   FileName fn_proj;
   /// CTF filename
   FileName fn_ctf;
   /// Rotational angle
   float  rot;
   /// Tilting angle
   float  tilt;
   /// Psi angle
   float  psi;
   /** Symmetry number.
       This number express to which symmetry matrix this projection
       is related to (-1: without symmetry, 0: using symmetry matrix 0,
       1: using symmetry matrix 1 ...) */
   int    sym;
}; 

/** Build from a Selection File and a Symmetry List. 
    The result is stored in the Recons_info array which should point
    to NULL when it is not initialized. */
void build_recons_info(SelFile &selfile, SelFile &selctf, const FileName &fn_ctf,
   const SymList &SL, Recons_info * &IMG_Inf, bool do_not_use_symproj);

//@}

/* ------------------------------------------------------------------------- */
/**@name Variability analysis */
//@{
/** Variability structure */
class VariabilityClass {
public:
   typedef enum {VAR_none, VAR_measuring, VAR_analyzing} t_VAR_status;
   t_VAR_status VAR_state;
   int Zoutput_volume_size;
   int Youtput_volume_size;
   int Xoutput_volume_size;
   Basic_ART_Parameters *prm;

   /// Vector of training vectors
   vector < matrix3D<double> > VA;

   /// Number of updates so far
   int N;   

   /// Constructor
   VariabilityClass(Basic_ART_Parameters *_prm,
      int _Zoutput_volume_size, int _Youtput_volume_size,
      int _Xoutput_volume_size);

   /** Start a new ART iteration. */
   void newIteration();
   
   /** Update data with a new volume.
       The update volume is set to zeros after this function */
   void newUpdateVolume(GridVolume *ptr_vol_out, Projection &read_proj);

   /** Finish analysis. */
   void finishAnalysis();
};
//@}

/* ------------------------------------------------------------------------- */
/**@name POCS */
//@{
/** POCS structure */
class POCSClass {
public:
   typedef enum {POCS_measuring, POCS_use, POCS_lowering, POCS_N_measure,
       POCS_N_use} t_POCS_status;
   t_POCS_status POCS_state;
   double POCS_avg;
   double POCS_stddev;
   double POCS_min;
   double POCS_max;
   double POCS_mean_error;
   double POCS_max_error;
   double POCS_global_mean_error;
   int POCS_freq;
   int POCS_i;
   int POCS_vec_i;
   int POCS_used;
   int POCS_N;
   int Zoutput_volume_size;
   int Youtput_volume_size;
   int Xoutput_volume_size;
   bool apply_POCS;
   matrix1D<double> POCS_errors;
   Basic_ART_Parameters *prm;

   /// Constructor
   POCSClass(Basic_ART_Parameters *_prm,
      int _Zoutput_volume_size, int _Youtput_volume_size,
      int _Xoutput_volume_size);
   
   /// Start New ART iteration
   void newIteration();

   /// Start new Projection
   void newProjection();

   /// Apply
   void apply(GridVolume &vol_blobs, int it, int images);
};
//@}

//@}

#endif
