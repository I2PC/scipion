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
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<(nr_psi); ipsi++ ) 
#define FOR_ALL_FLIPS() for (int iflip=0; iflip<nr_flip; iflip++)
#define SIGNIFICANT_WEIGHT_LOW 1e-8

/**@name MLalign2D */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm {
public:
  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_ref,fn_root,fn_frac;
  /** If true: use least-squares instead of maximum likelihood target */
  bool LSQ_rather_than_ML;
  /** Command line */
  string cline;
  /** Sigma value for expected pixel noise */
  double sigma_noise;
  /** sigma-value for origin offsets */
  double sigma_offset;
  /** vector for flag of active for all models */
  vector<bool> active_model;
  /** Vector containing estimated fraction for each model */
  vector<double> alpha_k;
  /** Vector containing estimated fraction for mirror of each model */
  vector<double> mirror_fraction;
  /** Flag for checking mirror images of all references */
  bool do_mirror;
  /** Flag whether to apply shifts in header of 2D-images */
  bool apply_shifts;
  /** Flag whether to fix estimates for model fractions */
  bool fix_fractions;
  /** Flag whether to fix estimate for sigma of origin offset */
  bool fix_sigma_offset;
  /** Flag whether to fix estimate for sigma of noise */
  bool fix_sigma_noise;
  /** Maximum-shift for conventional LSQ refinement */
  double max_shift;
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
  /** Number of pixels in one image */
  double npix;
  /** Sum of squared amplitudes of the references (all scaled to first one) */
  double A2;
  /** Verbose level:
      (1=default) gives progress bar & info 
      (0) gives no output to screen at all */
  int verb;
  /** Stopping criterium */
  double eps;
  // Sometimes (with very noisy data) artifacts at the origin pixel (0,0) are introduced
  // Probably this is caused by the applied scaling in Fourier space
  // The following hidden parameter serves to correct the origin pixel by averaging over its neighbouring 4 pixels
  bool do_esthetics;
  // Write out document file with orientations & models that give max. probability for each image
  bool write_docfile;
  // Write out images and selfile after each iteration
  bool write_intermediate;

public:
  // SelFile images (working and reference set)
  SelFile SF, SFr;
  // vector for flipping (i.e. 90/180-degree rotations) matrices
  vector<matrix2D<double> > F;

  // Vector for images to hold current references
  vector <ImageXmipp> Iref, Iold;

public:
  /// Read argument from command line
  void read(int argc, char **argv) _THROW;
  
  /// Read input images all in memory
  void read_all_input_in_memory() _THROW;

  /// Generate initial references from random subset averages
  void generate_initial_references() _THROW;

  /// Show
  void show();

  /// Usage
  void usage();

  /// Calculate probability density distribution for in-plane transformations
  void calculate_pdf_phi(double &sigma_offset, matrix2D<double> &P_phi, matrix2D<double> &Mr2) _THROW;

  /// Fill vector of matrices with all rotations of reference
  void rotate_reference(vector<ImageXmipp> &Iref, vector<bool> &active_model, 
			vector <vector< matrix2D<complex<double> > > > &Fref) _THROW;

  /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
  void reverse_rotate_reference(vector <vector< matrix2D<complex<double> > > > &Fnew, 
				vector<bool> &active_model, vector<matrix2D<double> > &Mref) _THROW;

  /// Calculate weighted averages for new model and new model parameters 
  void ML_integrate_phi_one_image(matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
				  vector<double> &P_model, matrix2D<double> &P_phi, matrix2D<double> &Mr2,
				  vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs, 
				  double &wsum_sigma_noise, double &wsum_sigma_offset, vector<double> &sumw, vector<double> &sumw_mirror, 
				  double &LL, double &maxcorr, 
				  int opt_refno, double opt_psi, double opt_xoff, double opt_yoff) _THROW;

  /// Integrate over all experimental images
  void ML_sum_over_all_images(SelFile &SF, vector<bool> &active_model, vector<ImageXmipp> &Iref, 
			  matrix2D<double> &P_phi, matrix2D<double> &Mr2, vector<double> &P_model, vector<double> &mirror_fraction, 
			  double &LL, double &avecorr, DocFile &DFo, 
			  vector<matrix2D<double> > &wsum_Mref,
			  double &wsum_sigma_noise, double &wsum_sigma_offset, vector<double> &sumw, vector<double> &sumw_mirror) _THROW;

  /// Write out reference images, selfile and logfile
  void write_output_files(const int iter, SelFile &SF, DocFile &DF, 
			  double &sumw_allrefs, double &LL, double &avecorr, vector<double> &conv) _THROW;

  /// Main routine
  void MLalign2D() _THROW;

};				    
//@}
