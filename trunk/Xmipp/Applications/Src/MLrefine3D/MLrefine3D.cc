/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)   
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

/* INCLUDES ---------------------------------------------------------------- */
#include <Reconstruction/Programs/Prog_Refine3d.hh> 

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {

  int                         iter,c,nn,imgno,volno,opt_refno,converged=0;;
  double                      LL,sumw_allrefs,convv,sumcorr;
  vector<double>              conv;
  double                      aux,wsum_sigma_noise, wsum_sigma_offset;
  vector<matrix2D<double> >   wsum_Mref;
  vector<ImageXmipp>          Ireg;
  vector<double>              sumw,sumw_cv,sumw_mirror;
  matrix2D<double>            P_phi,Mr2;
  FileName                    fn_doc,fn_sel,fn_iter,fn_tmp;
  matrix1D<double>            oneline(0);
  DocFile                     DFo,DFf;
  SelFile                     SFa;
 
  Prog_Refine3d_prm           prm;
  Prog_MLalign2D_prm          ML2D_prm;

  // Get input parameters
  try {

    // Read command line
    prm.read(argc,argv);
    fn_iter=prm.fn_root+"_it";
    prm.show();
    // Write starting volume(s) to disc with correct name for iteration loop
    prm.remake_SFvol(prm.istart-1,true);

    // Read MLalign2D-stuff
    ML2D_prm.read(argc,argv);
    if (!check_param(argc,argv,"-psi_step")) ML2D_prm.psi_step=prm.angular;
    ML2D_prm.fn_root=prm.fn_root;
    ML2D_prm.fast_mode=true;
    ML2D_prm.do_mirror=true;
    ML2D_prm.write_docfile=true;
    ML2D_prm.write_selfiles=false;
    ML2D_prm.fn_ref=prm.fn_root+"_lib.sel";
    // Project volume and read lots of stuff into memory
    prm.project_reference_volume(ML2D_prm.SFr);
    ML2D_prm.SF=prm.SF;

    ML2D_prm.produce_Side_info();
    ML2D_prm.produce_Side_info2();    
    ML2D_prm.show(true);

    // Initialize some stuff
    DFf.reserve(2*ML2D_prm.SFr.ImgNo()+4);
    for (int refno=0; refno<ML2D_prm.n_ref; refno++) conv.push_back(-1.);
    ML2D_prm.Iold.clear(); // To save memory

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}
    
  try {

    // Loop over all iterations
    iter=prm.istart;
    while (!converged && iter<=prm.Niter ) {

      if (prm.verb>0) {
	cerr        << "--> 3D-EM volume refinement:  iteration " << iter <<" of "<< prm.Niter<<endl;
	prm.fh_hist << "--> 3D-EM volume refinement:  iteration " << iter <<" of "<< prm.Niter<<endl;
      }

      DFo.reserve(2*prm.SF.ImgNo()+1);
      DFo.clear();
      if (ML2D_prm.maxCC_rather_than_ML) 
	DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");
      else 
	DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");

      // Pre-calculate pdfs
      if (!ML2D_prm.maxCC_rather_than_ML) ML2D_prm.calculate_pdf_phi();

      // Integrate over all images
      ML2D_prm.ML_sum_over_all_images(prm.SF,ML2D_prm.Iref,LL,sumcorr,DFo,wsum_Mref,
				      wsum_sigma_noise,wsum_sigma_offset,sumw,sumw_mirror); 
	
      // Update model parameters
      ML2D_prm.update_parameters(wsum_Mref,wsum_sigma_noise,wsum_sigma_offset,
				 sumw,sumw_mirror,sumcorr,sumw_allrefs); 

      if (ML2D_prm.write_intermediate) 
	ML2D_prm.write_output_files(iter,SFa,DFf,DFo,sumw_allrefs,LL,sumcorr,conv);
      else ML2D_prm.output_to_screen(iter,sumcorr,LL);
      if (ML2D_prm.maxCC_rather_than_ML) 
	prm.fh_hist << " Average maxCC = "<<sumcorr<<endl;
      else
	prm.fh_hist << " LL = "<<LL<<" sigma_noise= "<<ML2D_prm.sigma_noise<<endl;

      // Reconstruct new volumes from the reference images
      for (volno=0; volno<prm.Nvols; volno++)
	prm.reconstruction(argc, argv, iter, volno);

      // Update the reference volume selection file
      prm.remake_SFvol(iter,false);

      // Mask (and filter) the volumes
      prm.post_process_volumes(argc, argv);

      // Check convergence 
      if (prm.check_convergence(iter)) {
	converged=1;
	if (prm.verb>0) cerr <<"--> Optimization converged!"<<endl;
      } 

      // Re-project volumes
      if (!converged) {
	prm.project_reference_volume(ML2D_prm.SFr);
	// Read new references from disc (I could just as well keep them in memory, maybe...)
	ML2D_prm.SFr.go_beginning();
	c=0;
	while(!ML2D_prm.SFr.eof()) {
	  ML2D_prm.Iref[c].read(ML2D_prm.SFr.NextImg(),false,false,false,false);
	  ML2D_prm.Iref[c]().set_Xmipp_origin();
	  c++;
	}
      }

      iter++;
    } // end loop iterations

    if (!converged && prm.verb>0) cerr <<"--> Optimization was stopped before convergence was reached!"<<endl;
  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}
}




