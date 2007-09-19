/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>
#include <reconstruction/symmetrize.h>
#include <vector>

typedef struct Wedgelist {int num; double th0; double thF;} wedgelist;

#define FOR_ALL_MODELS() for (int refno=0;refno<nr_ref; refno++)
#define FOR_ALL_ROT() for (int irot=0; irot<nr_rot; irot++ ) 
#define FOR_ALL_TILT() for (int itilt=0; itilt<nr_tilt; itilt++ ) 
#define FOR_ALL_PSI() for (int ipsi=0; ipsi<nr_psi; ipsi++ ) 
#define FOR_ALL_LIMITED_TRANSLATIONS() for (int itrans=0; itrans<nr_trans; itrans++)
#define FOR_ALL_DEFOCUS_GROUPS() for (int ifocus=0; ifocus<nr_focus; ifocus++) 
#define SIGNIFICANT_WEIGHT_LOW 1e-4

/**@name MLalign3D */
//@{
/** MLalign3D parameters. */
class Prog_MLalign3D_prm {
public:
  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_ref,fn_root,fn_frac,fn_sig,fn_doc,fn_misalign;
  FileName fn_wlist,fn_sym,fn_solv,fn_solv2, fn_mask;
  /** sigma-value for origin offsets */
  double sigma_offset;
  /** Vector containing estimated fraction for each model */
  vector<double> alpha_k;
  /** Flag whether to fix estimates for model fractions */
  bool fix_fractions;
  /** Flag whether to fix estimate for sigma of noise */
  bool fix_sigma_noise;
  /** Flag whether to fix estimate for sigma of origin offsets */
  bool fix_sigma_offset;
  /** Starting iteration */
  int istart;
  /** Number of iterations to be performed */
  int Niter;
  /** dimension of the images, and sqrt(3)*dimension (for rotated wedges)  */
  int dim, bigdim;
  /** Number of steps to sample rot, tilt, psi and limited translations */
  int nr_rot, nr_rot_tilt, nr_tilt, nr_psi, nr_trans;
  /** Range and step size for angular searches */
  double rot0,  rotF,  rot_step;
  double tilt0, tiltF, tilt_step;
  double psi0,  psiF,  psi_step;
  /** Number of reference images */
  int nr_ref;
  /** Verbose level:
      1: gives progress bar (=default)
      0: gives no output to screen at all */
  int verb;
  // SelFile images (working, test and reference set)
  SelFile SF, SFr;
  // Vector for images to hold current references
  vector <Matrix3D<double> > Iref, Iold;
  /** Maximum-shift for conventional LSQ refinement */
  double max_shift;
  // For all tomograms: angles, offsets and wedge parameters
  vector<double> img_rot,img_tilt,img_psi,img_xoff,img_yoff,img_zoff,img_th0,img_thF,img_wednr;
  // Matrices for calculating PDF of (in-plane) translations
  Matrix3D<double> pdf_trans, Mr2;
  /** Number for which limited translation is zero */
  int zero_trans;
  /** Matrix of sigma2-values  */
  Matrix3D<double> Msigma2;
  /** Offsets for limited translations */
  vector<Matrix1D<double> > Vtrans;
  /** number of different wedges */
  int nr_wedge;
  Matrix3D<double> smooth_edge_mask, corr_mask;
  Matrix3D<int> outside_mask, mask;
  /* wedgelist */
  vector<wedgelist> wedges;
  /* Symmetry information */
  SymList SL;
  /* Use CCF mode instead of ML */
  bool ccf_mode;
  /* sigma noise */
  double sigma_noise2;
  /* store as double: dim x dim x dim */
  double dim3;
  /* for smoothing */
  double theta, theta_step, theta0;
  /* also for smoothing */
  int nr_img;
  /* for CCF-mode only: subtract current volume from average */
  bool do_subtract_current;


public:
  /// Read arguments from command line
  void read(int argc, char **argv);
  
  /// Show
  void show();

  /// Usage
  void usage();

  /// Extended Usage
  void extended_usage();

  /// Read input images all in memory
  void produce_Side_info();

  /// Calculate prior probability density funtion of the translations
  void calculate_pdf_trans();

  /// Calculate constrained CCF as Beck et al., Science 2005
  void CCF_integrate(Matrix3D<double> &Mimg, Matrix2D<double> &A_img, 
		     vector<Matrix3D<double> > &wsum_Mimgs,
		     vector<Matrix3D<double> > &wsum_Mwedge,
		     double &th0, double &thF, vector<double> &sumw, double &maxccf, 
		     int &opt_refno, double &opt_rot, double &opt_tilt, double &opt_psi,
		     double &opt_xoff, double &opt_yoff, double &opt_zoff);

  /// Calculate weighted ML averages for new model and new model
  /// parameters using real-space likelihood functions
  void ML_integrate(Matrix3D<double> &Mimg, Matrix2D<double> &A_img, 
		    vector<Matrix3D<double> > &wsum_Mimgs,
		    vector<Matrix3D<double> > &wsum_Mwedge,
		    double &wsum_sigma_noise2,  double &wsum_sigma_offset,  
		    double &th0, double &thF, 
		    vector<double> &sumw, double &LL, double &fracweight, 
		    int &opt_refno, double &opt_rot, double &opt_tilt, double &opt_psi,
		    double &opt_xoff, double &opt_yoff, double &opt_zoff);

  /// Integrate over all experimental images
  void ML_sum_over_all_images(SelFile &SF, vector<Matrix3D<double> > &Iref, 
			      double &LL, double &sumcorr, DocFile &DFo, 
			      vector<Matrix3D<double> > &wsum_Mref,
			      vector<Matrix3D<double> > &wsum_Mwedge,
			      double &wsum_sigma_noise2,
                              double &wsum_sigma_offset, vector<double> &sumw);

  /// Update all model parameters
  void update_parameters(vector<Matrix3D<double> > &wsum_Mref,
			 vector<Matrix3D<double> > &wsum_Mwedge, 
			 double &wsum_sigma_noise2, double &wsum_sigma_offset, 
			 vector<double> &sumw, double &sumcorr, double &sumw_allrefs, int iter);

  // Solvent flattening
  void solvent_flattening(FileName &fn_solvent);

  // Planned misalignment
  void misalign(Matrix2D<double> &A_img, DocFile DFmis, FileName fn_img);

  // Calculate symmetrized tomogram and its corresponding missing wedge
  void symmetrize_tomogram(Matrix3D<double> &Min, Matrix3D<double> &Mout, 
			   Matrix3D<double> &Mwedge, 
			   SymList &SL, Matrix2D<double> A, 
			   double th0, double thF, bool do_inverse=true);

  /// Write out reference images, selfile and logfile
  void write_output_files(const int iter, SelFile &SF, DocFile &DF, 
			  double &sumw_allrefs, vector<double> &sumw, 
			  double &LL, double &avecorr);


};				    
//@}
