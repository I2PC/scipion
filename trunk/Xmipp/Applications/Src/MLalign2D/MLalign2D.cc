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

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {

  int c,nn,imgno,opt_refno;
  double LL,sumw_allrefs,convv,sumcorr;
  bool converged;
  vector<double> conv;
  double aux,wsum_sigma_noise, wsum_sigma_offset;
  vector<matrix2D<double > > wsum_Mref;
  vector<ImageXmipp> Ireg;
  vector<double> sumw,sumw_mirror;
  matrix2D<double> P_phi,Mr2,Maux;
  FileName fn_img,fn_tmp;
  matrix1D<double> oneline(0);
  DocFile DFo,DFf;
  SelFile SFo,SFa;

  Prog_MLalign2D_prm prm;

  // Get input parameters
  try {
    prm.read(argc,argv);
    // Create references from random subset averages, or read them from selfile
    if (prm.fn_ref=="") {
      if (prm.n_ref!=0) {
	prm.generate_initial_references();
      } else {
	REPORT_ERROR(1,"Please provide -ref or -nref");
      }
    }
    prm.produce_Side_info();
    prm.show();

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}
    
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

      conv.clear();
      DFo.clear();
      if (prm.LSQ_rather_than_ML) 
	DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");
      else 
	DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Pmax/sumP (8)");

      // Pre-calculate pdfs
      if (!prm.LSQ_rather_than_ML) prm.calculate_pdf_phi();

      // Integrate over all images
      prm.ML_sum_over_all_images(prm.SF,prm.Iref, 
				 LL,sumcorr,DFo, 
				 wsum_Mref,wsum_sigma_noise,wsum_sigma_offset,sumw,sumw_mirror); 

      // Update model parameters
      sumw_allrefs=0.;
      for (int refno=0;refno<prm.n_ref; refno++) {
	if (sumw[refno]>0.) {
	  prm.Iref[refno]()=wsum_Mref[refno];
	  prm.Iref[refno]()/=sumw[refno];
	  prm.Iref[refno].weight()=sumw[refno];
	  sumw_allrefs+=sumw[refno];
	  if (prm.do_esthetics) MAT_ELEM(prm.Iref[refno](),0,0)=
			      (MAT_ELEM(prm.Iref[refno](),1,0)+MAT_ELEM(prm.Iref[refno](),0,1)+
			       MAT_ELEM(prm.Iref[refno](),-1,0)+MAT_ELEM(prm.Iref[refno](),0,-1))/4;
	  if (!prm.fix_fractions) prm.alpha_k[refno]=sumw[refno]/prm.SF.ImgNo();
	  if (!prm.fix_fractions) prm.mirror_fraction[refno]=sumw_mirror[refno]/sumw[refno];
	} else {
	  prm.Iref[refno].weight()=0.;
	  prm.Iref[refno]().init_zeros();
	  prm.alpha_k[refno]=0.;
	  prm.mirror_fraction[refno]=0.;
	}
      }
      if (!prm.fix_sigma_offset) prm.sigma_offset=sqrt(wsum_sigma_offset/(2*sumw_allrefs));
      if (!prm.fix_sigma_noise)  prm.sigma_noise=sqrt(wsum_sigma_noise/(sumw_allrefs*prm.dim*prm.dim));
                                                                                     
      sumcorr/=sumw_allrefs;

      // Check convergence 
      converged=true;
      for (int refno=0;refno<prm.n_ref; refno++) { 
	if (prm.Iref[refno].weight()>0.) {
	  Maux=mul_elements(prm.Iold[refno](),prm.Iold[refno]());
	  convv=1/(Maux.compute_avg());
	  Maux=prm.Iold[refno]()-prm.Iref[refno]();
	  Maux=mul_elements(Maux,Maux);
	  convv*=Maux.compute_avg();
	  conv.push_back(convv);
	  if (convv>prm.eps) converged=false;
	} else {
	  conv.push_back(-1.);
	}
      }

      if (prm.write_intermediate) {
	prm.write_output_files(iter,SFa,DFf,sumw_allrefs,LL,sumcorr,conv);
      } else {
	// Output current parameter values to screen 
	if (prm.verb>0) { 
	  if (prm.LSQ_rather_than_ML) cout <<"  iter "<<iter<<" <CC>= "+FtoA(sumcorr,10,5);
	  else {
	    cout <<"  iter "<<iter<<" noise= "<<FtoA(prm.sigma_noise,10,7)<<" offset= "<<FtoA(prm.sigma_offset,10,7);
	    cout <<"  LL= "<<LL<<" <Pmax/sumP>= "<<sumcorr<<endl;
	    cout <<"  Model  fraction  mirror-fraction "<<endl;
	    for (int refno=0;refno<prm.n_ref; refno++)  
	      cout <<"  "<<ItoA(refno+1,5)<<" "<<FtoA(prm.alpha_k[refno],10,7)<<" "<<FtoA(prm.mirror_fraction[refno],10,7)<<endl;
	  }
	}
      }
      
      if (converged) {
	if (prm.verb>0) cerr <<" Optimization converged!"<<endl;
	break;
      }

    } // end loop iterations


    // Write out converged structures
    prm.write_output_files(-1,SFa,DFf,sumw_allrefs,LL,sumcorr,conv);
      
    if (prm.write_docfile) {
      // Write out docfile with optimal transformation & references
      fn_img=prm.fn_root+".doc";
      DFo.write(fn_img);
    }
    if (prm.write_selfiles) {
      // Also write out selfiles of all experimental images, classified according to optimal reference image
      for (int refno=0;refno<prm.n_ref; refno++) { 
	DFo.go_beginning();
	SFo.clear();
	for (int n=0; n<DFo.dataLineNo(); n++ ) {
	  DFo.next();
	  fn_img=((DFo.get_current_line()).get_text()).erase(0,3);
	  DFo.adjust_to_data_line();
	  if ((refno+1)==(int)DFo(5)) SFo.insert(fn_img,SelLine::ACTIVE);
	}
	fn_tmp=prm.fn_root+"_ref";
	fn_tmp.compose(fn_tmp,refno+1,"sel");
	SFo.write(fn_tmp);
      }
    }

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}
}




