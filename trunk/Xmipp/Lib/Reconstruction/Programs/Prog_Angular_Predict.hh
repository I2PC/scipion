/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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
#ifndef _PROG_ANGULAR_PREDICT
   #define _PROG_ANGULAR_PREDICT

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippProgs.hh>
#include <Classification/xmippPC.hh>
#include "Prog_Angular_Project.hh"
#include "Prog_Angular_Distance.hh"
#include "Prog_Angular_Denoise.hh"
#include "../symmetries.hh"

#include <map>
#include <algorithm>

/**@name Angular Predict */
//@{
/** Angular Predict parameters. */
class Prog_angular_predict_prm: public Prog_parameters {
public:
   /** Root filename from Angular Reference */
   FileName fn_ref;
   /** Selfile with the reference projections */
   FileName fn_refsel;
   /** Root filename from Angular Project */
   FileName fn_pcaproj;
   /** List of PCAs to use for projection */
   matrix1D<double> PCA_list;
   /** Filename with the angles. As provided by the movement file
       provided by xmipp_project */
   FileName fn_ang;
   /** Filename for the output angles.*/
   FileName fn_out_ang;
   /** Filename for the symmetry file */
   FileName fn_sym;
   /** Maximum projection change.
       -1 means all are allowed. */
   double max_proj_change;
   /** Maximum psi change.
       -1 means all are allowed */
   double max_psi_change;
   /** Psi step */
   double psi_step;
   /** Maximum shift change */
   double max_shift_change;
   /** Shift step */
   double shift_step;
   /** Threshold for discarding images.
       If it is 20%, only the 20% of the images will be kept each round. */
   double th_discard;
   /** Scale for the denoising */
   int s_denoise;
   /** Check mirrors.
      If 1 then mirror versions of the experimental images are also explored.*/
   int check_mirrors;
   /** Way to pick views.
       0 maximum correlation of the first group.
       1 average of the most populated group.
       2 maximum correlation of the most populated group. */
   int pick;
   #define TELL_ROT_TILT  1
   #define TELL_PSI_SHIFT 2
   #define TELL_OPTIONS   4
   /** Show level.*/
   int tell;

public:
   // Vector with the rotational angles to learn
   vector<double> rot;
   // Vector with the tilting angles to learn
   vector<double> tilt;
   // Visible space
   Prog_angular_denoise_prm Angular_denoise;
   // PCA projector
   Prog_angular_project_prm PCA_projector;
   /* Vector with all projections of the reference images
      in the selected PCAs. */
   vector<xmippCTVectors *> PCA_ref;
   // Index of the current processed image
   int current_img;
   // Vector of predicted rotational angles
   vector<double> predicted_rot;
   // Vector of predicted tilting angles
   vector<double> predicted_tilt;
   // Vector of predicted psi angles
   vector<double> predicted_psi;
   // Vector of predicted shiftX
   vector<double> predicted_shiftX;
   // Vector of predicted shiftY
   vector<double> predicted_shiftY;
   // Vector of predicted corr
   vector<double> predicted_corr;
   // Parameters for computing distances
   Prog_angular_distance_prm distance_prm;
public:
   /// Read argument from command line
   void read(int argc, char **argv) _THROW;

   /// Show
   void show();

   /// Usage
   void usage();

   /** Produce side info.
       Read the PCASet file and the Mask_no.
       An exception is thrown if any of the files is not found*/
   void produce_side_info() _THROW;

   /** Build candidate list.
       Build a candidate list with all possible reference projections
       which are not further than the maximum allowed change from
       the given image. 
       
       The list size is the total number of reference images. For each
       image the list is true if it is still a candidate.*/
   void build_ref_candidate_list(const ImageXmipp &I,
      vector<bool> &candidate_list, vector<double> &cumulative_corr,
      vector<double> &sumxy, vector<double> &sumxx, vector<double> &sumyy);

   /** Refine candidate list via correlation.
       Given the projection of the image at hand on a certain PCA
       space, and the list of alive candidates, this function
       correlates the input image with all alive candidates and
       leave to pass on ly th% of the images.
       
       m is the PCA being studied.*/
   void refine_candidate_list_with_correlation(int m, xmippVector &PCA,
      xmippCTVectors &TS, vector<bool> &candidate_list,
      vector<double> &cumulative_corr,
      vector<double> &sumxy, vector<double> &sumxx, vector<double> &sumyy,
      double &dim, double th=50);

   /** Evaluate candidates by correlation. The evaluation is returned in
       candidate_rate. Furthermore, this function returns the threshold for
       passing in the "score exam", a 7*/
   double evaluate_candidates(const vector<double> &vscore,
      const vector<int> &candidate_idx, vector<double> &candidate_rate,
      double weight);
   
   /** Group views.
       The input images are supposed to be ordered by rate.
       The groups are also sorted by rate. */
   void group_views(const vector<double> &vrot,
      const vector<double> &vtilt, const vector<double> &vpsi,
      const vector<int> &best_idx, const vector<int> &candidate_idx,
      vector< vector<int> > &groups);

   /** Pick the best image from the groups.
       If method == 0 it takes the maximum of the first group (the
       one with best rate). If method==1, it takes the maximum
       of the most populated group. */
   int pick_view(int method,
      vector< vector<int> > &groups,
      vector<double> &vscore,
      vector<double> &vrot,
      vector<double> &vtilt,
      vector<double> &vpsi,
      const vector<int> &best_idx,
      const vector<int> &candidate_idx, const vector<double> &candidate_rates);

   /** Predict rotational and tilting angles.
       The function returns the two assigned angles and the corresponding
       correlation. The index of the best matching reference image is also
       returned. The function predict shift and psi angle calls this
       one for evaluating each possible combination.*/
   double predict_rot_tilt_angles(ImageXmipp &I,
      double &assigned_rot, double &assigned_tilt, int &best_ref_idx) _THROW;

   /** Predict angles and shift.
       This function searches in the shift-psi space and for each combination
       it correlates with the whole reference set. */
   double predict_angles(ImageXmipp &I,
      double &assigned_shiftX, double &assigned_shiftY,
      double &assigned_rot, double &assigned_tilt, double &assigned_psi);

   /** Finish processing.
       Close all output files. */
   void finish_processing();
};
//@}
#endif
