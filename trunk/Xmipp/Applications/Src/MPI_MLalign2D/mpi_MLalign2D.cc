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
#include <Reconstruction/Programs/Prog_MLalign2D.hh> 
#include <mpi.h>


/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {

  // For parallelization
  int rank, size, num_img_tot;

  int c,nn,imgno,opt_refno,iaux;
  double LL,sumw_allrefs,convv,sumcorr,sumw_cv;
  bool converged;
  vector<double> conv;
  double aux,wsum_sigma_noise, wsum_sigma_offset;
  vector<matrix2D<double > > wsum_Mref;
  vector<ImageXmipp> Ireg;
  vector<double> sumw,sumw_mirror;
  matrix2D<double> P_phi,Mr2,Maux;
  vector<matrix2D<double> > Mwsum_sigma2;
  vector<int> count_defocus;
  FileName fn_img,fn_tmp;
  matrix1D<double> oneline(0);
  DocFile DFo,DFf;
  SelFile SFo,SFa;
    
  Prog_MLalign2D_prm prm;

  // Init Parallel interface		
  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
   
  // Get input parameters
  try {
    
    // Read subsequently to avoid problems in restart procedure
    for (int proc=0; proc<size; proc++) {    
      if (proc==rank) prm.read(argc,argv);
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // Create references from random subset averages, or read them from selfile
    if (prm.fn_ref=="") {
      if (prm.n_ref!=0) {
	if (rank==0) {
	  prm.generate_initial_references();
	} else {
	  prm.fn_ref=prm.fn_root+"_it";
	  prm.fn_ref.compose(prm.fn_ref,0,"sel");
	}
	MPI_Barrier(MPI_COMM_WORLD);
      } else {
	REPORT_ERROR(1,"Please provide -ref or -nref");
      }
    }

    if (prm.fourier_mode && rank==0) prm.estimate_initial_sigma2();
    MPI_Barrier(MPI_COMM_WORLD);

    // Select only relevant part of selfile for this rank
    prm.SF.mpi_select_part(rank,size,num_img_tot);

    prm.produce_Side_info();

    if (rank==0) prm.show();
    else  prm.verb=0;

  } catch (Xmipp_error XE) {if (rank==0) {cout << XE; prm.usage();} MPI_Finalize(); exit(1);} 

    
  try {

    Maux.resize(prm.dim,prm.dim);
    Maux.set_Xmipp_origin();
    DFo.reserve(2*prm.SF.ImgNo()+1);
    DFf.reserve(2*prm.SFr.ImgNo()+4);
    SFa.reserve(prm.Niter*prm.n_ref);
    SFa.clear();

  // Loop over all iterations
    for (int iter=prm.istart; iter<=prm.Niter; iter++) {

      if (prm.verb>0) cerr << "  multi-reference refinement:  iteration " << iter <<" of "<< prm.Niter<<endl;

      for (int refno=0;refno<prm.n_ref; refno++) prm.Iold[refno]()=prm.Iref[refno]();

      DFo.clear();
      if (rank==0) 
	if (prm.LSQ_rather_than_ML)
	  DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");
	else
	  DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");

      // Pre-calculate pdfs
      if (!prm.LSQ_rather_than_ML) prm.calculate_pdf_phi();

      // Integrate over all images
      prm.ML_sum_over_all_images(prm.SF,prm.Iref,LL,sumcorr,DFo, 
				 wsum_Mref,wsum_sigma_noise,Mwsum_sigma2,sumw_cv,
				 wsum_sigma_offset,sumw,sumw_mirror,count_defocus); 

      // Here MPI_allreduce of all wsums,LL and sumcorr !!!
      MPI_Allreduce(&LL,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      LL=aux;
      MPI_Allreduce(&sumcorr,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      sumcorr=aux;
      MPI_Allreduce(&wsum_sigma_noise,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      wsum_sigma_noise=aux;
      MPI_Allreduce(&wsum_sigma_offset,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      wsum_sigma_offset=aux;
      MPI_Allreduce(&sumw_cv,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      sumw_cv=aux;
      for (int ifocus=0;ifocus<prm.nr_focus;ifocus++) {
	MPI_Allreduce(MULTIDIM_ARRAY(Mwsum_sigma2[ifocus]),MULTIDIM_ARRAY(Maux),
		      MULTIDIM_SIZE(Mwsum_sigma2[ifocus]),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	Mwsum_sigma2[ifocus]=Maux;
	MPI_Allreduce(&count_defocus[ifocus],&iaux,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	count_defocus[ifocus]=iaux;
      }
      for (int refno=0;refno<prm.n_ref; refno++) { 
	MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]),MULTIDIM_ARRAY(Maux),
		      MULTIDIM_SIZE(wsum_Mref[refno]),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	wsum_Mref[refno]=Maux;
	MPI_Allreduce(&sumw[refno],&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	sumw[refno]=aux;
	MPI_Allreduce(&sumw_mirror[refno],&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	sumw_mirror[refno]=aux;
      }

      // Update model parameters
      prm.update_parameters(wsum_Mref,wsum_sigma_noise, Mwsum_sigma2,
			    sumw_cv,wsum_sigma_offset,
			    sumw,sumw_mirror,sumcorr,sumw_allrefs,
			    count_defocus);    

      // Check convergence 
      converged=prm.check_convergence(conv);

      if (rank==0) {
	if (prm.write_intermediate)
	  prm.write_output_files(iter,SFa,DFf,sumw_allrefs,LL,sumcorr,conv);
	else prm.output_to_screen(iter,sumcorr,LL);
      }
      
      if (converged) {
	if (prm.verb>0) cerr <<" Optimization converged!"<<endl;
	break;
      }

    } // end loop iterations


    // All nodes write out temporary DFo and SFo
    if (prm.write_docfile) {
      fn_img.compose(prm.fn_root,rank,"tmpdoc");
      DFo.write(fn_img);
    }
    if (prm.write_selfiles) {
      for (int refno=0;refno<prm.n_ref; refno++) { 
	DFo.go_beginning();
	if (rank==0) DFo.next();
	SFo.clear();
	for (int n=0; n<DFo.dataLineNo(); n++ ) {
	  fn_img=((DFo.get_current_line()).get_text()).erase(0,3);
	  DFo.adjust_to_data_line();
	  if ((refno+1)==(int)DFo(5)) SFo.insert(fn_img,SelLine::ACTIVE);
	  DFo.next();
	}
	fn_tmp=prm.fn_root+"_ref";
	fn_tmp.compose(fn_tmp,refno+1,"");
	fn_tmp.compose(fn_tmp,rank,"tmpsel");
	SFo.write(fn_tmp);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0) {
      // Write out converged structures
      prm.write_output_files(-1,SFa,DFf,sumw_allrefs,LL,sumcorr,conv);
      
      if (prm.write_docfile) {
	// Write out docfile with optimal transformation & references
	DFo.clear();
	for (int rank2=0; rank2<size; rank2++) {
	  fn_img.compose(prm.fn_root,rank2,"tmpdoc");
	  int ln=DFo.LineNo();
	  DFo.append(fn_img);
	  DFo.locate(DFo.get_last_key());
	  DFo.next();
	  DFo.remove_current();
	  system(((string)"rm -f "+fn_img).c_str());
	}
	fn_tmp=prm.fn_root+".doc";
	DFo.write(fn_tmp);
      }
      if (prm.write_selfiles) {
	// Also write out selfiles of all experimental images, classified according to optimal reference image
	for (int refno=0;refno<prm.n_ref; refno++) { 
	  SFo.clear();
	  for (int rank2=0; rank2<size; rank2++) {
	    fn_img=prm.fn_root+"_ref";
	    fn_img.compose(fn_img,refno+1,"");
	    fn_img.compose(fn_img,rank2,"tmpsel");
	    SFo.append(fn_img);
	    system(((string)"rm -f "+fn_img).c_str());
	  }
	  fn_tmp=prm.fn_root+"_ref";
	  fn_tmp.compose(fn_tmp,refno+1,"sel");
	  SFo.write(fn_tmp);
	}
      }
    }

  } catch (Xmipp_error XE) {if (rank==0) {cout << XE; prm.usage();} MPI_Finalize(); exit(1);}


  MPI_Finalize();	
  return 0;

}




