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
#include <XmippData/xmippGeometry.hh>
#include <XmippData/xmippFilters.hh>
#include <XmippData/xmippMasks.hh>
#include <Reconstruction/CTF.hh>
#include <vector>

#define FOR_ALL_MODELS() for (int refno=0;refno<n_ref; refno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<nr_psi; ipsi++ ) 
#define FOR_ALL_FLIPS() for (int iflip=0; iflip<nr_flip; iflip++)
#define FOR_ALL_LIMITED_TRANSLATIONS() for (int itrans=0; itrans<nr_trans; itrans++)
#define FOR_ALL_DEFOCUS_GROUPS() for (int ifocus=0; ifocus<nr_focus; ifocus++) 
#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLVALUE 1e-4
#define SMALLANGLE 1.75

/**@name MLalign2D */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm {
public:
  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_ref,fn_root,fn_frac,fn_sig,fn_doc,fn_ctf;
  /** If true: use maximum cross-correlation instead of maximum likelihood target */
  bool maxCC_rather_than_ML;
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
  int dim, dim2;
  /** Number of steps to sample in-plane rotation in 90 degrees */
  int nr_psi;
  /** Number of operations in "flip-array" (depending on do_mirror) */
  int nr_flip;
  /** Sampling rate for in-plane rotation */
  float psi_step;
  /** Total degrees in FOR_ALL_ROTATIONS */
  double psi_max;
  /** Total number of no-mirror rotations in FOR_ALL_FLIPS */
  int nr_nomirror_flips;
  /** Number of reference images */
  int n_ref;
  /** Total number of experimental images */
  int nr_exp_images;
  /** Sum of squared amplitudes of the references */
  vector<double> A2;
  /** Verbose level:
      1: gives progress bar (=default)
      0: gives no output to screen at all */
  int verb;
  /** Stopping criterium */
  double eps;
  /** Write out document file with orientations & models that give
      max. probability for each image */
  bool write_docfile;
  /** Write out selection files according to model assignments */
  bool write_selfiles;
  /** Write out images and selfile after each iteration */
  bool write_intermediate;
  /** SelFile images (working, test and reference set) */
  SelFile SF, SFr;
  /** vector for flipping (i.e. 90/180-degree rotations) matrices */
  vector<matrix2D<double> > F;
  /** Vector for images to hold references (new & old) */
  vector <ImageXmipp> Iref, Iold;
  /** Matrices for calculating PDF of (in-plane) translations */
  matrix2D<double> P_phi, Mr2;
  /** Fast mode */
  bool fast_mode;
  /** Fast mode */
  double C_fast;
  /** Very crude origin pixel atrifact correction */
  bool do_esthetics;
  /** Maximum shift to be trusted */
  double max_shift;
  /** Limit translational searches */
  bool limit_trans;
  /** Number of limited translations */
  int nr_trans;
  /** Number for which limited translation is zero */
  int zero_trans;
  /** Offsets for limited translations */
  vector<matrix1D<double> > Vtrans;
  /** Start all optimal offsets from zero values */
  bool zero_offsets;
  /** Limited search range for origin offsets */
  double search_shift;
  /** Limit orientational searches */
  bool limit_rot;
  /** Limited search range for projection directions */
  double search_rot;
  /** Save memory options */
  bool save_mem1, save_mem2, save_mem3;
  /** Vectors to store old phi, theta, xoff and yoff for all images */
  vector<float> imgs_oldphi, imgs_oldtheta, imgs_oldxoff, imgs_oldyoff;


public:
  /// Read arguments from command line
  void read(int argc, char **argv);
  
  /// Show
  void show(bool ML3D=false);

  /// Usage
  void usage();

  /// Extended Usage
  void extended_usage(bool ML3D=false);

  /// Setup lots of stuff
  void produce_Side_info();

  /// Read reference images in memory & set offset vectors
  /// (This produce_side_info is Selfile-dependent!)
  void produce_Side_info2();

  /// Read and write optimal translations to disc 
  /// (not to store them all in memory)
  void write_offsets(FileName fn, vector<double> &data);
  void read_offsets(FileName fn, vector<double> &data);

  /// Generate initial references from random subset averages
  void generate_initial_references();

  /// Calculate probability density distribution for in-plane transformations
  void calculate_pdf_phi();

  /// Fill vector of matrices with all rotations of reference
  void rotate_reference(vector<ImageXmipp> &Iref, 
			bool fill_real_space, 
			bool fill_fourier_space, 
			vector <vector< matrix2D<double> > > &Mref,
			vector <vector< matrix2D<complex<double> > > > &Fref);

  /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
  void reverse_rotate_reference(vector <vector< matrix2D<complex<double> > > > &Fnew, 
				vector <vector< matrix2D<double> > > &Mnew, bool real_space, 
				vector<matrix2D<double> > &Mref);

  /// Calculate which references have projection directions close to
  /// phi and theta
  void preselect_directions(float &phi, float &theta,
			    vector<double> &pdf_directions);

  /// Pre-calculate which model and phi have significant probabilities 
  /// without taking translations into account!
  void preselect_significant_model_phi(matrix2D<double> &Mimg, vector<double> &offsets, 
				       vector <vector< matrix2D<double > > > &Mref, 
				       matrix2D<int> &Msignificant,
				       vector<double > &pdf_directions);

  // Calculate the FT of a translated matrix using a phase shift in
  // Fourier space
  void Fourier_translate2D(const matrix2D<complex<double> > &Fimg, 
			   int focus, matrix1D<double> &trans,
			   matrix2D<complex<double> > &Fimg_shift);

  // If not determined yet: search optimal offsets using maxCC
  // Then for all optimal translations, calculate all translated FTs 
  // for each of the flipped variants
  void calculate_fourier_offsets(matrix2D<double> &Mimg, int focus,
				 vector <vector< matrix2D<complex<double> > > > &Fref,
				 matrix2D<double> &ctf, vector<double> &offsets,
				 vector<vector<matrix2D<complex<double> > > > &Fimg_trans,
				 matrix2D<int> &Moffsets, matrix2D<int> &Moffsets_mirror);

  // Calculate translated matrices for all limited translations
  // for each of the flipped variants
  void calculate_realspace_offsets(matrix2D<double> &Mimg, vector<double > &offsets,
				   vector<double > &pdf_directions,
				   vector<vector<matrix2D<double> > > &Mimg_trans,
				   matrix2D<int> &Moffsets, matrix2D<int> &Moffsets_mirror);

  // ML-integration over limited translations,
  // and with -fast way of selection significant rotations
  void ML_integrate_locally(matrix2D<double> &Mimg, 
			    vector <vector< matrix2D<double> > > &Mref, 
			    vector <vector< matrix2D<double> > > &Mwsum_imgs, 
			    double &wsum_sigma_noise, double &wsum_sigma_offset, 
			    vector<double> &sumw, vector<double> &sumw_mirror, 
			    double &LL, double &fracweight, 
			    int &opt_refno, double &opt_psi, 
			    matrix1D<double> &opt_offsets, 
			    vector<double> &opt_offsets_ref,
			    vector<double> &pdf_directions);

  /// ML-integration over all (or -fast) translations
  void ML_integrate_complete(matrix2D<double> &Mimg, 
			     vector <vector< matrix2D<complex<double> > > > &Fref, 
			     matrix2D<int> &Msignificant,
			     vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs, 
			     double &wsum_sigma_noise, double &wsum_sigma_offset, 
			     vector<double> &sumw, vector<double> &sumw_mirror, 
			     double &LL, double &fracweight, int &opt_refno, double &opt_psi, 
			     matrix1D<double> &opt_offsets, vector<double> &opt_offsets_ref,
			     vector<double > &pdf_directions);

  /// Calculate maxCC averages for new model and new model parameters 
  void maxCC_search_complete(matrix2D<double> &Mimg, 
			     vector <vector< matrix2D<complex<double> > > > &Fref, 
			     vector <vector< matrix2D<double> > > &Mref, 
			     double &max_shift, 
			     vector <vector< matrix2D<double> > > &Msum_imgs, 
			     vector<double> &sumw, vector<double> &sumw_mirror, 
			     double &minSQ, int &opt_refno, double &opt_psi, 
			     matrix1D<double> &opt_offsets,
			     vector<double> &pdf_directions);

  /// Integrate over all experimental images
  void ML_sum_over_all_images(SelFile &SF, vector<ImageXmipp> &Iref, 
			      double &LL, double &sumcorr, DocFile &DFo, 
			      vector<matrix2D<double> > &wsum_Mref,
			      double &wsum_sigma_noise, double &wsum_sigma_offset, 
			      vector<double> &sumw, vector<double> &sumw_mirror);

  /// Update all model parameters
  void update_parameters(vector<matrix2D<double> > &wsum_Mref,
			 double &wsum_sigma_noise, double &wsum_sigma_offset, 
			 vector<double> &sumw, vector<double> &sumw_mirror, 
			 double &sumcorr, double &sumw_allrefs);

  /// check convergence
  bool check_convergence(vector<double> &conv);

  /// Output some parameters to screen
  void output_to_screen(int &iter, double &sumcorr, double &LL);

  /// Write out reference images, selfile and logfile
  void write_output_files(const int iter, SelFile &SF, DocFile &DF, DocFile &DFo,
			  double &sumw_allrefs, double &LL, double &avecorr, 
			  vector<double> &conv);

};
//@}
