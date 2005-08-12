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
#include <XmippData/xmippFFT.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippFilters.hh>
#include <XmippData/xmippMasks.hh>
#include <vector>

#define FOR_ALL_MODELS() for (int refno=0;refno<n_ref; refno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<nr_psi; ipsi++ ) 
#define FOR_ALL_FLIPS() for (int iflip=0; iflip<nr_flip; iflip++)
#define FOR_ALL_LIMITED_TRANSLATIONS() for (int itrans=0; itrans<nr_trans; itrans++)
#define FOR_ALL_DEFOCUS_GROUPS() for (int ifocus=0; ifocus<nr_focus; ifocus++) 
#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLANGLE 1.75

/**@name MLalign2D */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm {
public:
  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_ref,fn_root,fn_frac,fn_sig,fn_cv;
  /** If true: use least-squares instead of maximum likelihood target */
  bool LSQ_rather_than_ML;
  /** Command line */
  string cline;
  /** Sigma value for expected pixel noise */
  double sigma_noise;
  /** sigma-value for origin offsets */
  double sigma_offset;
  /** Vector containing estimated fraction for each model */
  vector<double> alpha_k;
  /** Vector containing estimated fraction for mirror of each model */
  vector<double> mirror_fraction;
  /** Flag for checking mirror images of all references */
  bool do_mirror;
  /** Flag whether to fix estimates for model fractions */
  bool fix_fractions;
  /** Flag whether to fix estimate for sigma of origin offset */
  bool fix_sigma_offset;
  /** Flag whether to fix estimate for sigma of noise */
  bool fix_sigma_noise;
  /** Starting iteration */
  int istart;
  /** Number of iterations to be performed */
  int Niter;
  /** dimension of the images */
  int dim;
  /** Number of steps to sample in-plane rotation in 90 degrees */
  int nr_psi;
  /** Number of operations in "flip-array" (depending on do_mirror) */
  int nr_flip;
  /** Sampling rate for in-plane rotation */
  float psi_step;
  /** Number of reference images */
  int n_ref;
  /** Sum of squared amplitudes of the references */
  vector<double> A2;
  /** Verbose level:
      1: gives progress bar (=default)
      0: gives no output to screen at all */
  int verb;
  /** Stopping criterium */
  double eps;
  // Write out document file with orientations & models that give max. probability for each image
  bool write_docfile;
  // Write out selection files according to model assignments
  bool write_selfiles;
  // Write out images and selfile after each iteration
  bool write_intermediate;
  // SelFile images (working, test and reference set)
  SelFile SF, SFcv, SFr;
  // vector for flipping (i.e. 90/180-degree rotations) matrices
  vector<matrix2D<double> > F;
  // Vector for images to hold current references
  vector <ImageXmipp> Iref, Iold;
  /** Maximum-shift for conventional LSQ refinement */
  double max_shift;
  // Vector to store optimal origin offsets (for CC-precalculation)
  vector<vector<matrix1D<double> > > imgs_offsets;
  // For fourier mode: optimal origin offsets
  vector<double > offset_x,offset_y;
  // Matrices for calculating PDF of (in-plane) translations
  matrix2D<double> P_phi, Mr2;
  // Fast mode
  bool fast_mode;
  // Fast mode
  double C_fast;
  // Fast mode
  double Paccept_fast;
  /** If true: use fourier instead of real space maximum likelihood target */
  bool fourier_mode;
  /** If true: use cross-validation in fourier_mode */
  bool do_cv;
  /** Very crude origin pixel atrifact correction */
  bool do_esthetics;
  /** Number of limited translations */
  int nr_trans;
  /** Number for which limited translation is zero */
  int zero_trans;
  /** For fourier-mode: matrix of sigma2-values (one for each defocuss-group  */
  vector<matrix2D<double> > Msigma2;
  /** Offsets for limited translations */
  vector<matrix1D<double> > Vtrans;
  /** number of defocus groups */
  int nr_focus;
  matrix2D<int> resol_mask;

public:
  /// Read arguments from command line
  void read(int argc, char **argv) _THROW;
  
  /// Show
  void show();

  /// Usage
  void usage(bool ML3D=false);

  /// Extended Usage
  void extended_usage(bool ML3D=false);

  /// Read input images all in memory
  void produce_Side_info() _THROW;

  /// Read document file with defocus groups
  void read_defocus_groups() _THROW;

  /// Calculate initial sigma2 from average power spectrum of the
  /// experimental images
  void estimate_initial_sigma2() _THROW;

  /// Generate initial references from random subset averages
  void generate_initial_references() _THROW;

  /// Calculate probability density distribution for in-plane transformations
  void calculate_pdf_phi() _THROW;

  /// Fill vector of matrices with all rotations of reference
  void rotate_reference(vector<ImageXmipp> &Iref, bool &also_real_space, vector <vector< matrix2D<double> > > &Mref,
			vector <vector< matrix2D<complex<double> > > > &Fref) _THROW;

  /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
  void reverse_rotate_reference(vector <vector< matrix2D<complex<double> > > > &Fnew, 
				vector <vector< matrix2D<double> > > &Mnew, bool &real_space, 
				vector<matrix2D<double> > &Mref) _THROW;
 
  /// Pre-calculate which model and phi have significant probabilities without taking translations into account!
  void preselect_significant_model_phi(matrix2D<double> &Mimg, vector<matrix1D<double> > &offsets, 
				       vector <vector< matrix2D<double > > > &Mref, 
				       matrix2D<int> &Msign) _THROW;

  /// Calculate weighted ML averages for new model and new model
  /// parameters using fourier-space likelihood functions
  void ML_integrate_FS_model_phi(matrix2D<double> &Mimg, 
				 vector <vector< matrix2D<complex<double> > > > &Fref,
				 vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs,
				 matrix2D<double > &sigma2, matrix2D<double > &Mwsum_sigma2, 
				 vector<double> &sumw, vector<double> &sumw_mirror, 
				 double &LL, double &fracweight, int &opt_refno, double &opt_psi, 
				 double &opt_xoff, double &opt_yoff, 
				 double &sumw_cv, bool &cv_flag) _THROW;

  /// Calculate weighted ML averages for new model and new model parameters 
  void ML_integrate_model_phi_trans(matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
				    matrix2D<int> &Msignificant,
				    vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs, 
				    double &wsum_sigma_noise, double &wsum_sigma_offset, 
				    vector<double> &sumw, vector<double> &sumw_mirror, 
				    double &LL, double &fracweight, int &opt_refno, double &opt_psi, 
				    matrix1D<double> &opt_offsets, vector<matrix1D<double> > &opt_offsets_ref) _THROW;

  /// Calculate LSQ averages for new model and new model parameters 
  void LSQ_search_model_phi_trans(matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
				  double &max_shift, 
				  vector <vector< matrix2D<double> > > &Msum_imgs, 
				  vector<double> &sumw, vector<double> &sumw_mirror, 
				  double &minSQ, int &opt_refno, double &opt_psi, 
				  matrix1D<double> &opt_offsets) _THROW;

  /// Integrate over all experimental images
  void ML_sum_over_all_images(SelFile &SF, vector<ImageXmipp> &Iref, 
			      double &LL, double &sumcorr, DocFile &DFo, 
			      vector<matrix2D<double> > &wsum_Mref,
			      double &wsum_sigma_noise, vector<matrix2D<double> > &Mwsum_sigma2, 
                              double &sumw_cv, double &wsum_sigma_offset, 
			      vector<double> &sumw, vector<double> &sumw_mirror,
			      vector<int> &count_defocus) _THROW;

  /// Update all model parameters
  void update_parameters(vector<matrix2D<double> > &wsum_Mref,
			 double &wsum_sigma_noise, vector<matrix2D<double> > &Mwsum_sigma2, 
			 double &sumw_cv, double &wsum_sigma_offset,
			 vector<double> &sumw, vector<double> &sumw_mirror, 
			 double &sumcorr, double &sumw_allrefs,
			 vector<int> &count_defocus) ;

  /// check convergence
  bool check_convergence(vector<double> &conv);

  /// Output some parameters to screen
  void output_to_screen(int &iter, double &sumcorr, double &LL);

  /// Write out reference images, selfile and logfile
  void write_output_files(const int iter, SelFile &SF, DocFile &DF, 
			  double &sumw_allrefs, double &LL, double &avecorr, vector<double> &conv) _THROW;


};				    
//@}
