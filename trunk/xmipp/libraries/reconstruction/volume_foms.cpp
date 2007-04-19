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

#include "volume_foms.h"
#include "phantom.h"

#include <data/histogram.h>
#include <data/docfile.h>

// Some nice definitions ...................................................
#define L    VOLMATRIX(*vol_label)
#define P    VOLMATRIX(*vol_phantom)
#define R    VOLMATRIX(*vol_recons)
#define M    VOLMATRIX(*vol_mask)
#define pi   VOLVOXEL (*vol_phantom  ,k,i,j)
#define ri   VOLVOXEL (*vol_recons   ,k,i,j)
#define li   VOLVOXEL (*vol_label    ,k,i,j)
#define mi   VOLVOXEL (*vol_mask     ,k,i,j)
#define l_r  VOLVOXEL (*vol_label,   (int)ZZ(r), (int)YY(r), (int)XX(r))
#define p_r  VOLVOXEL (*vol_phantom, (int)ZZ(r), (int)YY(r), (int)XX(r))
#define r_r  VOLVOXEL (*vol_recons,  (int)ZZ(r), (int)YY(r), (int)XX(r))
#define m_r  VOLVOXEL (*vol_mask,    (int)ZZ(r), (int)YY(r), (int)XX(r))

typedef struct coord {
   double x;
   double y;
   double z;
};


/* ------------------------------------------------------------------------- */
/* Compute Number of voxels in each feature                                  */
/* ------------------------------------------------------------------------- */
void compute_voxels_in_feat(Volume *vol_label,
   matrix1D<double> &feat_voxels) {
   int sel_feat;

   // Choose dimensions for output vector
   int label_no=(int)(*vol_label)().compute_max();
   feat_voxels.resize(label_no+1); // 0, 1, ..., FeatNo()

   FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*vol_label)) {
      sel_feat = (int) VOLVOXEL(*vol_label,k,i,j);
      if (sel_feat>=0)
         VEC_ELEM(feat_voxels,sel_feat)++;
   }
}

/* ------------------------------------------------------------------------- */
/* Show voxels in feature                                                    */
/* ------------------------------------------------------------------------- */
void show_voxels_in_feat(const Volume *vol_phantom,
   const Volume *vol_recons, const Volume *vol_label,
   int selected_feat, ostream &out) {
   out << "#Voxels for feature " << selected_feat;

   FOR_ALL_ELEMENTS_IN_MATRIX3D(L) {
      if (li==selected_feat) {
         out << "(z=" << k << ",y=" << i << ",x=" << j << ")"
             << "Phantom=" << pi
             << " Recons=" << ri << endl;
      }
   }
}

void show_voxels_in_feat(const Volume *vol_phantom,
   const Volume *vol_recons, const Volume *vol_label,
   const Phantom &phantom_descr,
   int selected_feat, ostream &out) {
   matrix1D<double> r(3);
   int true_feat=ABS(selected_feat);
   if        (selected_feat>0) {
      out << "#Voxels for feature " << selected_feat << endl;
   } else if (selected_feat<0) {
      out << "#Voxel for border of feature  " << selected_feat << endl;
   } else
      out << "#Voxels for background\n";

   out << "#  distance   Phantom      Recons      (x,y,z)       relative\n";

   FOR_ALL_ELEMENTS_IN_MATRIX3D(L) {
      if (li==selected_feat) {
         if (selected_feat!=0) {
            XX(r)=(double)j-XX(phantom_descr(true_feat)->Center);
            YY(r)=(double)i-YY(phantom_descr(true_feat)->Center);
            ZZ(r)=(double)k-ZZ(phantom_descr(true_feat)->Center);
         } else {
            XX(r)=(double)j;
            YY(r)=(double)i;
            ZZ(r)=(double)k;
         }
         out.width(10); out << r.module() << " ";
         out.width(10); out << pi << " ";
         out.width(10); out << ri << " ";
         out.width(3);  out << " (" << j << ",";
         out.width(3);  out << i << ",";
         out.width(3);  out << k << ")";
                        out << " (" << r.transpose() << ")" << endl;
      }
   }
}

/* ------------------------------------------------------------------------- */
/* Compute structural consistency FOMs for a feature                         */
/* ------------------------------------------------------------------------- */
void compute_sc_FOMs(
   const Volume *vol_phantom, const Volume *vol_recons,
   const Volume *vol_label, const Volume *vol_mask, int sel_feat,
   double &scL2_FOM, double &scL1_FOM, double &scmu_FOM, double &scdev_FOM,
   double &scrange_FOM, double &sccorr_FOM, double &scinf_FOM, bool tell) {

   double L2_error    =0, L1_error   =0;
   double phantom_sum =0, recons_sum =0, phantom_recons_sum=0;
   double phantom_sum2=0, recons_sum2=0;
   int    no_samples  =0;
   double Mp=-1e30, mp=1e30;
   double Mr=-1e30, mr=1e30;

   FOR_ALL_ELEMENTS_IN_MATRIX3D(L) {
      if (!mi) continue; // If it is not in the mask continue
      int f=(int)li;
      if (sel_feat==-1 || f==sel_feat) {
         double diff    = pi-ri;
         L2_error     += diff*diff*0.25;
         L1_error     += ABS(diff);
         phantom_sum  += pi;
         recons_sum   += ri;
         phantom_sum2 += pi*pi;
         recons_sum2  += ri*ri;
	 phantom_recons_sum += pi*ri;
         no_samples++;
         Mp=MAX(Mp,pi); Mr=MAX(Mr,ri);
         mp=MIN(mp,pi); mr=MIN(mr,ri);
      }
   }
   double phantom_avg    = (no_samples>0)? phantom_sum/no_samples:0;
   double recons_avg     = (no_samples>0)? recons_sum /no_samples:0;
   double phantom_stddev = (no_samples>0)?
      sqrt(phantom_sum2/no_samples - phantom_avg*phantom_avg):0;
   double recons_stddev  = (no_samples>0)?
      sqrt(recons_sum2 /no_samples - recons_avg *recons_avg):0;
   double covariance     = (no_samples>0)?
      (phantom_recons_sum-phantom_avg*recons_sum-recons_avg*phantom_sum+
         phantom_avg*recons_avg*no_samples)/no_samples:0;

   if (no_samples>0) {
      scL2_FOM    = (double)(1-1.0/no_samples*L2_error);
      scL1_FOM    = (double)(1-0.5/no_samples*L1_error);
      scmu_FOM    = (double)(1-0.5*ABS(phantom_avg    - recons_avg));
      scdev_FOM   = (double)(1-    ABS(phantom_stddev - recons_stddev));
      scrange_FOM = (double)(1-0.5*(ABS(Mp-Mr) + ABS(mp-mr)));
      if (phantom_stddev!=0 && recons_stddev!=0)
         sccorr_FOM  = covariance/(phantom_stddev*recons_stddev);
      else sccorr_FOM=-1;
   } else {
      scL2_FOM = scL1_FOM = scmu_FOM = scdev_FOM = scrange_FOM =
         sccorr_FOM = -1.0;
   }

   // Compute mutual information FOM
   histogram1D hist_phantom, hist_recons;
   histogram2D hist_phantom_recons;
   hist_phantom.init(mp,Mp,200);
   hist_recons.init (mr,Mr,200);
   hist_phantom_recons.init(mp,Mp,200,mr,Mr,200);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(L) {
      if (!mi) continue; // If it is not in the mask continue
      int f=(int)li;
      if (sel_feat==-1 || f==sel_feat) {
         hist_phantom.insert_value(pi);
         hist_recons. insert_value(ri);
	 hist_phantom_recons.insert_value(pi,ri);
      }
   }
   // Make them probabilities
   double phantom_information=0;
   double recons_information=0;
   double phantom_recons_information=0;
   if (no_samples!=0) {
      hist_phantom /= no_samples;
      hist_recons  /= no_samples;
      hist_phantom_recons /= no_samples;
      FOR_ALL_ELEMENTS_IN_MATRIX1D(hist_phantom)
	 if (hist_phantom(i)!=0)
            phantom_information+=(-hist_phantom(i)*log(hist_phantom(i)));
      FOR_ALL_ELEMENTS_IN_MATRIX1D(hist_recons)
	 if (hist_recons(i)!=0)
            recons_information+=(-hist_recons(i)*log(hist_recons(i)));
      FOR_ALL_ELEMENTS_IN_MATRIX2D(hist_phantom_recons)
	 if (hist_phantom_recons(i,j)!=0)
            phantom_recons_information+=
	       (-hist_phantom_recons(i,j)*log(hist_phantom_recons(i,j)));
      if (phantom_information!=0 && recons_information!=0)
         scinf_FOM=phantom_recons_information/(phantom_information*recons_information);
   } else scinf_FOM=-1;

   // Show process
   if (tell&0x4) { // This is the flag SHOW_PROCESS of Prog_evaluate.h
      cout << "   Selected feature: ";
      switch (sel_feat) {
         case -1: cout << "whole volume\n"; break;
         case  0: cout << "background\n"; break;
         default: cout << sel_feat << endl; break;
      }
      cout << "   L2_error:         " << L2_error/no_samples << endl;
      cout << "   L1_error:         " << L1_error/no_samples << endl;
      cout << "   phantom average:  " << phantom_avg << endl;
      cout << "   recons average:   " << recons_avg << endl;
      cout << "   phantom stddev:   " << phantom_stddev << endl;
      cout << "   recons stddev:    " << recons_stddev << endl;
      cout << "   phantom range:    [" << mp << "," << Mp << "]\n";
      cout << "   recons range:     [" << mr << "," << Mr << "]\n";
      cout << "   covariance:       " << covariance << endl
           << "   phantom inf:      " << phantom_information << endl
	   << "   recons inf:       " << recons_information << endl
	   << "   mutual inf:       " << phantom_recons_information << endl
           << endl;
   }
}

/* ------------------------------------------------------------------------- */
/* Compute Histogram based FOMs                                              */
/* ------------------------------------------------------------------------- */
//#define DEBUG
void compute_hs_FOMs(Volume *vol_phantom,
   Volume *vol_recons, const Volume *vol_label, const Volume *vol_mask,
   int sel_feat, const Phantom &phantom_descr,
   int back_mode, double back_param,
   double &hsmu_FOM, double &hsbr_FOM, double &hsdt_FOM, double &hsvr_FOM,
   int tell, const FileName &fn_histog) {

   // Mean separability variables ..........................................
   double mpf=0, mpb=0, mpo=0; // mean and variance of the phantom inner
   double vpf=0, vpb=0, vpo=0; // foreground, background and outer foreground
   double mrf=0, mrb=0, mro=0; // mean and variance of the reconstruction inner
   double vrf=0, vrb=0, vro=0; // foreground, background and outer foreground
   int   Nf =0, Nb =0, No =0; // Number of voxels
   double zpi,    zri;         // Normalized inference variable for the phantom
                              // and the reconstruction, between
   double zpo,    zro;         // Normalized inference variable for the phantom
                              // and the reconstruction, between
   double phantom_confidence_i, recons_confidence_i;
   double phantom_confidence_o, recons_confidence_o;

   // Histogram detection error variables .................................
   histogram1D Hpf, Hpb;    // Phantom fore and background histograms
   histogram1D Hrf, Hrb;    // Reconstruction for and background histograms
   double min, max;
   (*vol_phantom)().compute_double_minmax(min,max);
      Hpf.init(min,max,100);
      Hpb.init(min,max,100);
   (*vol_recons)().compute_double_minmax(min,max);
      Hrf.init(min,max,100);
      Hrb.init(min,max,100);
   double phantom_dt_err, recons_dt_err;

   // Vertical resolution variables ........................................
   int compute_vrFOM=phantom_descr(sel_feat)->Type=="dcy";
   double mp1=0, mp2=0, mp3=0, vp1=0, vp2=0, vp3=0;
   double mr1=0, mr2=0, mr3=0, vr1=0, vr2=0, vr3=0;
   double Nvr=0;             // Number of points in the planes
   double vrp, vrr;          // Vertical resolution for phantom and recons

   DCylinder *DC;
   double z1,z2,z3;
   if (compute_vrFOM) {
      DC=(DCylinder *) phantom_descr(sel_feat);
      z1=ROUND(ZZ(DC->Center)+DC->separation/2+DC->height/2);
      z2=ROUND(ZZ(DC->Center)-DC->separation/2-DC->height/2);
      z3=ROUND(ZZ(DC->Center));
      #ifdef DEBUG
         cout << "z1=" << z1 << endl;
         cout << "z2=" << z2 << endl;
         cout << "z3=" << z3 << endl;
      #endif
   }

   // Other variables ......................................................
   matrix1D<double> r(3), corner1(3), corner2(3), aux1(3), aux2(3);
   Feature *Back       = phantom_descr(sel_feat)->background(back_mode,back_param);
   Feature *Inner_fore = phantom_descr(sel_feat)->scale(1.0/3.0);
   Feature *Mid_fore   = phantom_descr(sel_feat)->scale(2.0/3.0);

   // Classify every voxel .................................................
   // foreground, background or nowhere (border or another feature)
   Back->corners(vol_recons,corner1,corner2);
   #ifdef DEBUG
      cout << "Foreground " << phantom_descr(sel_feat);
      cout << "Background " << Back;
      cout << "Inner  foreground " << Inner_fore;
      cout << "Middle foreground " << Mid_fore;
      cout << "Corner1 " << corner1.transpose() << endl;
      cout << "Corner2 " << corner2.transpose() << endl;
   #endif
   FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) {
      #ifdef DEBUG
         cout << "Studying " << r.transpose();
      #endif
      if (!m_r) {
         #ifdef DEBUG
            cout << "    Not in the mask\n";
         #endif
         continue; // If it is not in the mask continue
      }
      int totally_inside_foreground =
         phantom_descr(sel_feat)->voxel_inside(r)==8;
         #ifdef DEBUG
            cout << "   It's totally inside the foreground\n";
         #endif

      // Check if inside background
      if (Back->voxel_inside(r)==8) {
         // Is it purely in background?
         if (l_r==0) {
            #ifdef DEBUG
               cout << "   It's purely in the background\n";
            #endif
            mpb += p_r; vpb += p_r*p_r;
            mrb += r_r; vrb += r_r*r_r;
            Nb++;
            Hpb.insert_value(p_r);
            Hrb.insert_value(r_r);
         // touching the inner foreground?
         } else if (Inner_fore->voxel_inside(r)>0) {
            #ifdef DEBUG
               cout << "   It's in the inner foreground\n";
            #endif
            mpf += p_r; vpf += p_r*p_r;
            mrf += r_r; vrf += r_r*r_r;
            Nf++;
         // Touching the outer foreground
         } else if (Mid_fore->voxel_inside(r)!=8 && totally_inside_foreground) {
            #ifdef DEBUG
               cout << "   It's in the outer foreground\n";
            #endif
            mpo += p_r; vpo += p_r*p_r;
            mro += r_r; vro += r_r*r_r;
            No++;
         }

         // Check if it purely in the foreground for the histograms
         if (totally_inside_foreground) {
            Hpf.insert_value(p_r);
            Hrf.insert_value(r_r);
         }
      }

      // Check if for vrFOM
      if (compute_vrFOM && ABS(ZZ(r)-z2)<XMIPP_EQUAL_ACCURACY &&
          totally_inside_foreground) {
         #ifdef DEBUG
            cout << "   Adding it to hsvr\n";
         #endif
                   mp2 += p_r; vp2 += p_r*p_r; mr2 += r_r; vr2 += r_r*r_r;
         ZZ(r)=z1; mp1 += p_r; vp1 += p_r*p_r; mr1 += r_r; vr1 += r_r*r_r;
         ZZ(r)=z3; mp3 += p_r; vp3 += p_r*p_r; mr3 += r_r; vr3 += r_r*r_r;
         ZZ(r)=z2;
         Nvr++;
      }
   }

   // Finish computing the means and variances .............................
   if (Nf!=0) {mpf/=Nf; mrf/=Nf; vpf=vpf/Nf-mpf*mpf; vrf=vrf/Nf-mrf*mrf;}
   else {cout << "   Watch out, there is no voxel in the inner foreground\n";}
   if (Nb!=0) {mpb/=Nb; mrb/=Nb; vpb=vpb/Nb-mpb*mpb; vrb=vrb/Nb-mrb*mrb;}
   else {cout << "   Watch out, there is no voxel in the inside background\n";}
   if (No!=0) {mpo/=No; mro/=No; vpo=vpo/No-mpo*mpo; vro=vro/No-mro*mro;}
   else {cout << "   Watch out, there is no voxel in the outer foreground\n";}

   // Compute mean separability ............................................
   // Compute inference confidence for phantom
   // We are interested in integral {from -zr to zr} N(0,1) dz
   // Because of the erfc function definition we have to change the z value
   if (Nb!=0 && Nf!=0) {
      zri= ABS(mrf-mrb)/sqrt(vrf/Nf+vrb/Nb);
      recons_confidence_i=erf(zri/sqrt(2.0));

      if (vpf+vpb==0) {zpi=-1; hsmu_FOM=zri; phantom_confidence_i=1;}
      else {
         zpi=ABS(mpf-mpb)/sqrt(vpf/Nf+vpb/Nb);
         phantom_confidence_i=erf(sqrt(2.0)*zpi);
         // hsmu_FOM=1-ABS(phantom_confidence-recons_confidence);
         hsmu_FOM=zri-zpi;
      }
   } else hsmu_FOM=-1;

   // Compute border separability ..........................................
   if (No!=0 && Nb!=0) {
      zro= ABS(mro-mrb)/sqrt(vro/No+vrb/Nb);
      recons_confidence_i=erf(zro/sqrt(2.0));

      if (vpo+vpb==0) {zpo=-1; hsbr_FOM=zro; phantom_confidence_o=1;}
      else {
         zpo= ABS(mpo-mpb)/sqrt(vpo/No+vpb/Nb);
         phantom_confidence_o=erf(sqrt(2.0)*zpo);
         hsbr_FOM=zro-zpo;
      }
   } else hsbr_FOM=-1;

   // Compute detectability error ..........................................
   phantom_dt_err=detectability_error(Hpb,Hpf);
   recons_dt_err =detectability_error(Hrb,Hrf);
   hsdt_FOM=1-ABS(phantom_dt_err-recons_dt_err);

   // Compute vertical resolution FOM ......................................
   if (compute_vrFOM && Nvr!=0) {
      mp1 /= Nvr; mp2 /= Nvr; mp3 /= Nvr;
      mr1 /= Nvr; mr2 /= Nvr; mr3 /= Nvr;
      vp1 = vp1/Nvr-mp1*mp1; vp2 = vp2/Nvr-mp2*mp2; vp3 = vp3/Nvr-mp3*mp3;
      vr1 = vr1/Nvr-mr1*mr1; vr2 = vr2/Nvr-mr2*mr2; vr3 = vr3/Nvr-mr3*mr3;

      #ifdef DEBUG
         cout << "mr1=" << mr1 << "\nmr2=" << mr2 << "\nmr3=" << mr3 << endl
              << "vr1=" << vr1 << "\nvr2=" << vr2 << "\nvr3=" << vr3 << endl
              << "mp1=" << mp1 << "\nmp2=" << mp2 << "\nmp3=" << mp3 << endl
              << "vp1=" << vp1 << "\nvp2=" << vp2 << "\nvp3=" << vp3 << endl
              << "Nvr=" << Nvr;
      #endif
      if (vr1+vr2+4*vr3!=0) {
         vrr=(mr1+mr2-2*mr3)/sqrt(vr1+vr2+4*vr3);
         if (vp1+vp2+4*vp3!=0) {
            vrp=(mp1+mp2-2*mp3)/sqrt(vp1+vp2+4*vp3);
            hsvr_FOM=vrr/vrp;
         } else
            hsvr_FOM=vrr;
      } else {
         cout << "   Watch out, there is no deviation in reconstruction \n"
              << "   for vertical resolution FOM\n";
         hsvr_FOM=-1;
      }
   } else {
      cout << "   Watch out, there is no voxel for vertical resolution FOM\n";
      hsvr_FOM=-1;
   }

   // Show results .........................................................
   if (tell&0x4) { // This is the flag SHOW_PROCESS of Prog_evaluate.h
      cout << "Foreground " << phantom_descr(sel_feat);
      cout << "Background " << Back;
      cout << "Inner  foreground " << Inner_fore;
      cout << "Middle foreground " << Mid_fore;
      printf("   Selected feature: %d\n",sel_feat);
      printf("      Inner Foreground\n");
      printf("         phantom avg(): % 8.5f   var(): % 8.5f\n",mpf,vpf);
      printf("         recons  avg(): % 8.5f   var(): % 8.5f\n",mrf,vrf);
      printf("         Number of voxels: %d\n",Nf);
      printf("      Outer Foreground\n");
      printf("         phantom avg(): % 8.5f   var(): % 8.5f\n",mpo,vpo);
      printf("         recons  avg(): % 8.5f   var(): % 8.5f\n",mro,vro);
      printf("         Number of voxels: %d\n",No);
      printf("      Background\n");
      printf("         phantom avg(): % 8.5f   var(): % 8.5f\n",mpb,vpb);
      printf("         recons  avg(): % 8.5f   var(): % 8.5f\n",mrb,vrb);
      printf("         Number of voxels: %d\n",Nb);
      printf("      Inner-Background Separability\n");
      printf("         z(phantom): % 12.5f confidence: %9.7f\n",zpi,
         phantom_confidence_i);
      printf("         z(recons) : % 12.5f confidence: %9.7f\n",zri,
         recons_confidence_i);
      printf("      Border Separability\n");
      printf("         z(phantom): % 12.5f confidence: %9.7f\n",zpo,
         phantom_confidence_o);
      printf("         z(recons) : % 12.5f confidence: %9.7f\n",zro,
         recons_confidence_o);
      printf("      Old Hot spot detectability\n");
      printf("         z(phantom): % 9.5f\n",
         sqrt((double)Nb*Nf)*ABS(mpf-mpb)/sqrt(Nb*vpb+Nf*vpf));
      printf("         z(recons): % 9.5f\n",
         sqrt((double)Nb*Nf)*ABS(mrf-mrb)/sqrt(Nb*vrb+Nf*vrf));
      if (compute_vrFOM) {
         printf("      Vertical resoulution\n");
         printf("         p1.avg(): % 8.5f   r1.avg(): % 8.5f\n",mp1,mr1);
         printf("         p2.avg(): % 8.5f   r2.avg(): % 8.5f\n",mp2,mr2);
         printf("         p3.avg(): % 8.5f   r3.avg(): % 8.5f\n",mp3,mr3);
         printf("         p1.var(): % 8.5f   r1.var(): % 8.5f\n",vp1,vr1);
         printf("         p2.var(): % 8.5f   r2.var(): % 8.5f\n",vp2,vr2);
         printf("         p3.var(): % 8.5f   r3.var(): % 8.5f\n",vp3,vr3);
      }
      printf("      Detectability error\n");
      printf("         phantom dt error: %10.8f\n",phantom_dt_err);
      printf("         recons  dt error: %10.8f\n",recons_dt_err);
      printf("\n");
   }

   if (tell&0x8) { // This is the flag SAVE_HISTOGRAMS of Prog_evaluate.h
      DocFile DF;
      matrix1D<double> aux(6);
      DF.reserve(XSIZE(Hpf)+2);
      DF.append_comment("Histograms");
      DF.append_comment("Format: <phantom value> <Back> <Fore> <recons value> <Back> <Fore>");
      for (int i=STARTINGX(Hpf); i<FINISHINGX(Hpf); i++) {
          double pd; Hpf.index2val(i,pd);
          double rd; Hrf.index2val(i,rd);
          VEC_ELEM(aux,0)=pd;
          VEC_ELEM(aux,1)=Hpb(i);
          VEC_ELEM(aux,2)=Hpf(i);
          VEC_ELEM(aux,3)=rd;
          VEC_ELEM(aux,4)=Hrb(i);
          VEC_ELEM(aux,5)=Hrf(i);
          DF.append_data_line(aux);
      }
      if (fn_histog=="stdout") cout << DF;
      else                     DF.write(fn_histog);
   }

   // Free extra features
   delete Back;
   delete Inner_fore;
   delete Mid_fore;

}
#undef DEBUG

/* ------------------------------------------------------------------------- */
/* Compute Directional FOMs                                                  */
/* ------------------------------------------------------------------------- */
void compute_dr_FOMs(const Volume *vol_phantom, const Volume *vol_recons,
   const Volume *vol_mask,
   double rot, double tilt,
   Image *img_histog, int no_steps, double max_background,
   double &drrt_FOM, int tell, const FileName &fn_radon) {

   // Auxiliar variables ...................................................
   matrix3D<double> rot_phantom, rot_recons, rot_mask;
   matrix1D<double> RT_phantom, RT_recons, RTno;
   histogram1D      hist_recons;

   // Rotate volumes .......................................................
   // to align the given direction with Z axis
   if (rot!=0 || tilt!=0) {
      Euler_rotate(P,rot,tilt, 0.0F, rot_phantom);
      Euler_rotate(R,rot,tilt, 0.0F, rot_recons );
      Euler_rotate(M,rot,tilt, 0.0F, rot_mask );
   } else {
      rot_phantom=P;
      rot_recons=R;
      rot_mask=M;
   }

   // Initialise output image and Radon Transform ..........................
   (*img_histog)().resize(rot_recons.SliNo(),no_steps);
   (*img_histog)().startingY()=rot_recons.startingZ();
   compute_hist(rot_recons,hist_recons,max_background,
      rot_recons.compute_max(), no_steps);

   RT_recons.resize(rot_recons.SliNo());
   RT_recons.startingX()=rot_recons.startingZ();
   RT_phantom=RT_recons;
   RTno.init_zeros(RT_recons);

   // Run over the volume ..................................................
   FOR_ALL_ELEMENTS_IN_MATRIX3D(rot_phantom) {
      if (!mi) continue; // If it is not in the mask continue

      // Radon Transforms
      VEC_ELEM(RT_phantom,k) += VOL_ELEM(rot_phantom,k,i,j);
      VEC_ELEM(RT_recons ,k) += VOL_ELEM(rot_recons, k,i,j);
      VEC_ELEM(RTno,k)++;

      // Histogram Image
      int hh; hist_recons.val2index(VOL_ELEM(rot_recons,k,i,j),hh);
      if (hh!=-1) IMGPIXEL(*img_histog,k,hh)++;
   }

   // Compute Radon FOM ....................................................
   int N=0; double sum=0;
   for (int i=STARTINGX(RT_recons); i<=FINISHINGX(RT_recons); i++)
      if (RTno(i)!=0) {
         // Finish computing the Radon Transform
         RT_phantom(i) /= RTno(i);
         RT_recons(i)  /= RTno(i);

         // Sum for the drrt_FOM
         sum += ABS(RT_phantom(i)-RT_recons(i));
         N++;
      }
   if (N==0) drrt_FOM=-1;
   else      drrt_FOM=1-0.5*sum/N;

   // Show results .........................................................
   if (tell&0x4) { // This is the flag SHOW_PROCESS of Prog_evaluate.h
      DocFile DF;
      matrix1D<double> aux(4);
      DF.reserve(XSIZE(RT_recons)+2);
      DF.append_comment("Radon Transforms\n");
      DF.append_comment("Format: phantom recons diff no.voxels");
      for (int i=STARTINGX(RT_recons); i<=FINISHINGX(RT_recons); i++) {
           VEC_ELEM(aux,0)=VEC_ELEM(RT_phantom,i);
           VEC_ELEM(aux,1)=VEC_ELEM(RT_recons,i);
           VEC_ELEM(aux,2)=VEC_ELEM(RT_phantom,i)-VEC_ELEM(RT_recons,i);
           VEC_ELEM(aux,3)=VEC_ELEM(RTno,i);
           DF.append_data_line(aux);
      }
      if (fn_radon=="stdout") cout << DF;
      else                    DF.write(fn_radon);
   }
}

/* ------------------------------------------------------------------------- */
/* Compute Distance map                                                      */
/* ------------------------------------------------------------------------- */
#define di        VOLVOXEL(*vol_distance,k,i,j)
void compute_distance_map(const Volume *vol_label, const Phantom &label,
   const Volume *vol_mask, Volume *vol_distance) {
   // Variables ............................................................
   // Some vectors to perform the distance operations
   struct coord    auxcoord;
   vector<coord>   border;
   matrix1D<double> r(3), v_dist(3);

   // Search for all voxels in the border ..................................
   FOR_ALL_ELEMENTS_IN_MATRIX3D(L) {
      if (li==0) continue;
      if (li<0 ||
            VOLVOXEL(*vol_label,k  ,i  ,j-1)==0 ||
            VOLVOXEL(*vol_label,k  ,i  ,j+1)==0 ||
            VOLVOXEL(*vol_label,k  ,i-1,j  )==0 ||
            VOLVOXEL(*vol_label,k  ,i-1,j-1)==0 ||
            VOLVOXEL(*vol_label,k  ,i-1,j+1)==0 ||
            VOLVOXEL(*vol_label,k  ,i+1,j  )==0 ||
            VOLVOXEL(*vol_label,k  ,i+1,j-1)==0 ||
            VOLVOXEL(*vol_label,k  ,i+1,j-1)==0 ||
            VOLVOXEL(*vol_label,k-1,i  ,j  )==0 ||
            VOLVOXEL(*vol_label,k-1,i  ,j-1)==0 ||
            VOLVOXEL(*vol_label,k-1,i  ,j+1)==0 ||
            VOLVOXEL(*vol_label,k-1,i-1,j  )==0 ||
            VOLVOXEL(*vol_label,k-1,i-1,j-1)==0 ||
            VOLVOXEL(*vol_label,k-1,i-1,j+1)==0 ||
            VOLVOXEL(*vol_label,k-1,i+1,j  )==0 ||
            VOLVOXEL(*vol_label,k-1,i+1,j-1)==0 ||
            VOLVOXEL(*vol_label,k-1,i+1,j+1)==0 ||
            VOLVOXEL(*vol_label,k+1,i  ,j  )==0 ||
            VOLVOXEL(*vol_label,k+1,i  ,j-1)==0 ||
            VOLVOXEL(*vol_label,k+1,i  ,j+1)==0 ||
            VOLVOXEL(*vol_label,k+1,i-1,j  )==0 ||
            VOLVOXEL(*vol_label,k+1,i-1,j-1)==0 ||
            VOLVOXEL(*vol_label,k+1,i-1,j+1)==0 ||
            VOLVOXEL(*vol_label,k+1,i+1,j  )==0 ||
            VOLVOXEL(*vol_label,k+1,i+1,j-1)==0 ||
            VOLVOXEL(*vol_label,k+1,i+1,j+1)==0 ) {
         auxcoord.x=j; auxcoord.y=i; auxcoord.z=k;
         border.push_back(auxcoord);
      }
   }
   int bmax=border.size();

   // Compute distance map .................................................
   (*vol_distance)().init_zeros(L);

   init_progress_bar(ZSIZE(L));
   for (int k=STARTINGZ(L); k<=FINISHINGZ(L); k++) {
      progress_bar(k-STARTINGZ(L));
      for (int i=STARTINGY(L); i<=FINISHINGY(L); i++)
         for (int j=STARTINGX(L); j<=FINISHINGX(L); j++) {
            if (li<0) continue; // If it is in the border, continue
            if (!mi) {di=-1; continue;} // If it is not in the mask, continue

            // Search in all voxels in the border ..........................
            di=MAXFLOAT;
            for (int b=0; b<bmax; b++) {
                double dist=sqrt((k-border[b].z)*(k-border[b].z)+
                                (i-border[b].y)*(i-border[b].y)+
                                (j-border[b].x)*(j-border[b].x));
                if (dist<di) di=dist;
            }
         }
   }
   progress_bar(ZSIZE(L));
}

/* ------------------------------------------------------------------------- */
/* Compute Distance based FOMs                                               */
/* ------------------------------------------------------------------------- */
void compute_ds_FOMs(const Volume *vol_phantom, const Volume *vol_recons,
   const Volume *vol_label,
   const Volume *vol_distance, double &dsbl_FOM, double &dsad_FOM) {
// Some extra variables ....................................................
   int Nb=0;
   double blurr_d_1=0;
   double appear_d=0;
   double error;

// Compute errors ..........................................................
   FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*vol_label)) {
      if (li==0 && di!=-1) {
         error=(pi-ri)*(pi-ri)/4;
         blurr_d_1 +=1/di*error;
         appear_d  +=di*error;
         Nb++;
      }
   }

// Compute FOMS ............................................................
   /* CO: dsbl_FOM=1/(1.0/Nb*1/2*blurr_d_1); */
   dsbl_FOM=(double)(1-1.0/Nb*blurr_d_1);
   dsad_FOM=(double)(1-appear_d/Nb);
}

/* ------------------------------------------------------------------------- */
/* Show shape                                                                */
/* ------------------------------------------------------------------------- */
void show_shape(const Volume *vol_phantom, const Volume *vol_recons,
   const Volume *vol_label,
   const Volume *vol_distance, const FileName &fn_shape) {
   DocFile DF;
   matrix1D<double> aux(3);

   double r;

   int N=0;
   FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*vol_label))
      if (di!=-1) N++;
   DF.reserve(N+1);

   DF.append_comment("format: r phantom recons");

   FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*vol_label)) {
      if (di==-1) continue;
      if      (li< 0) r=0;
      else if (li==0) r= di;
      else            r=-di;
      VECTOR_R3(aux,r,pi,ri);
      DF.append_data_line(aux);
   }

   // Save shape
   if (fn_shape=="stdout") cout << DF;
   else DF.write(fn_shape);

   DF.clear();
}

/* ------------------------------------------------------------------------- */
/* Compute resolution                                                        */
/* ------------------------------------------------------------------------- */// Based on Spider
#ifdef NEVER_DEFINED
void compute_resolution(VolumeXmipp &vol_phantom,
   VolumeXmipp &vol_recons, double &resolution) {
   int tell=1;
   // If the phantom is not saved, save it
   FileName ext=vol_recons.name().get_extension();
   FileName original_phantom_name=vol_phantom.name();
   vol_phantom.rename((string)"superfeo."+ext);
   vol_phantom.write();

   // Generate and run spider batch
   ofstream spider_batch;
   spider_batch.open(((string)"b01."+ext).c_str());
   if (!spider_batch)
      REPORT_ERROR(1,"Compute resol:: Cannot open file for Spider batch");
   spider_batch
      << "rf 3\n"
      << vol_phantom.name().without_extension() << endl
      << vol_recons .name().without_extension() << endl
      << "1" << endl
      << "0.001 1000" << endl
      << "C" << endl
      << "90" << endl
      << "5" << endl
      << "resolution" << endl
      << "en\n"
   ;
   spider_batch.close();
   char *spider_prog=getenv("SPIDER");
   if (spider_prog==NULL)
       REPORT_ERROR(1,"Compute resol:: The environment variable SPIDER is not set");
   system(((string)spider_prog+" "+ext+" b01").c_str());
   vol_phantom.rename(original_phantom_name);

   // Process output file
   DocFile DF;
   try {
      DF.read((string)"resolution."+ext);
   } catch (Xmipp_error XE) {resolution=-1; return;}
   bool look_for_crossing=FALSE, found=FALSE;
   int i=1, imax=DF.get_last_key();
   while (!found && i<=imax) {
      if (!look_for_crossing) look_for_crossing=DF(i,2)>DF(i,3);
      else                    found=DF(i,2)<DF(i,3);
      i++;
   }
   i--;

   // Compute exact resolution
   if (found) {
      double x0=DF(i-1,0), y0=DF(i-1,2)-DF(i-1,3);
      double xF=DF(i  ,0), yF=DF(i  ,2)-DF(i  ,3);
      resolution=x0-y0*(xF-x0)/(yF-y0);

      //system(((string)"rm resolution."+ext).c_str());
      system("rm LOG* results* b01*");
      system(((string)"rm superfeo."+ext).c_str());
   } else resolution=-1;

   // Show results .........................................................
   if (tell&0x4) // This is the flag SHOW_PROCESS of Prog_evaluate.h
      system(((string)"cat resolution."+ext).c_str());

//   system(((string)"rm resolution."+ext).c_str());
}

// Based on Bsoft
void compute_resolution(VolumeXmipp &vol_phantom,
  VolumeXmipp &vol_recons, double &resolution) {
   // Prepare volumes to be evaluated
   if (vol_phantom.name()!="")
      system(((string)"ln -s "+vol_phantom.name()+" superfeo.spi").c_str());
   else vol_phantom.write("superfeo.spi");
   if (vol_recons.name()!="")
      system(((string)"ln -s "+vol_recons.name()+" superfeo2.spi").c_str());
   else vol_recons.write("superfeo2.spi");

   // Run bresolve of bsoft to compute resolution
   system("bresolve -s 1 -m superfeo.spi superfeo2.spi > superfeo3");

   // Process output file
   ifstream fh_resol;
   fh_resol.open("superfeo3");
   if (!fh_resol)
      REPORT_ERROR(1,
         "compute_resolution: Cannot open results file from Bresolve");
   try {
      fh_resol >> resolution;
      fh_resol >> resolution;
   } catch (...) {
      cerr << "compute_resolution: Error reading resolution, ignoring it\n";
      resolution=-1;
   }
   fh_resol.close();

   // Remove intermidiate files
   system("rm superfeo.spi superfeo2.spi superfeo3");
}
#endif

// Based on Xmipp
void compute_resolution(VolumeXmipp &vol_phantom,
  VolumeXmipp &vol_recons, double &resolution) {
  matrix1D<double> frequency, FSC;
  resolution=compute_FSC(vol_phantom,vol_recons,1,frequency,FSC);
}

#ifdef NEVER_DEFINED
// Based on Bsoft
double compute_FSC(VolumeXmipp &vol_phantom,
   VolumeXmipp &vol_recons, double sampling_rate,
   matrix1D<double> &frequency, matrix1D<double> &FSC) {
   double resolution=-1;
   frequency.clear();
   FSC.clear();

   // Prepare volumes to be evaluated
   if (vol_phantom.name()!="")
      system(((string)"ln -s "+vol_phantom.name()+" superfeo.spi").c_str());
   else vol_phantom.write("superfeo.spi");
   if (vol_recons.name()!="")
      system(((string)"ln -s "+vol_recons.name()+" superfeo2.spi").c_str());
   else vol_recons.write("superfeo2.spi");

   // Run bresolve of bsoft to compute resolution
   string command=(string)"bresolve -s "+FtoA(sampling_rate)+
      " -v4 -m superfeo.spi superfeo2.spi > superfeo3";
   system(command.c_str());

   // Process output file
   ifstream fh_resol;
   fh_resol.open("superfeo3");
   if (!fh_resol)
      REPORT_ERROR(1,
         "compute_FSC: Cannot open results file from Bresolve");
   try {
      string line;

      // Find the line that says radius
      bool radius_found=false;
      do {
         if (!getline(fh_resol,line))
            REPORT_ERROR(1,
               "compute_FSC: unexpected end of file from Bresolve");
         if (line.find("Radius")==0) radius_found=true;
      } while (!radius_found);

      // Copy until the line that says Fourier
      bool fourier_found=false;
      vector<double> auxfreq, auxFSC;
      vector<string> tokens;
      do {
         if (!getline(fh_resol,line))
            REPORT_ERROR(1,
               "compute_FSC: unexpected end of file from Bresolve");
         if (line.find("Fourier")==0) {fourier_found=true; continue;}

         // It is a valid line
         tokenize(line,tokens);
         double f  =AtoF(tokens[1]);
         double fsc=AtoF(tokens[4]);
         auxfreq.push_back(f);
         auxFSC. push_back(fsc);
         if (resolution==-1 && fsc<0.5) resolution=f;
      } while (!fourier_found);

      // Now copy
      frequency.resize(auxfreq.size());
      FSC.resize(auxfreq.size());
      FOR_ALL_ELEMENTS_IN_MATRIX1D(frequency) {
         frequency(i)=auxfreq[i];
         FSC      (i)=auxFSC [i];
      }
   } catch (Xmipp_error XE) {
      cerr << XE << endl
           << "compute_resolution: Error reading resolution, ignoring it\n";
   }
   fh_resol.close();

   // Fix a bug of Bresolve
   if (FSC(0)>1) FSC(0)=1;

   // Remove intermidiate files
   system("rm superfeo.spi superfeo2.spi superfeo3");
   return resolution;
}
#endif

// Based on Xmipp
double compute_FSC(VolumeXmipp &vol_phantom,
   VolumeXmipp &vol_recons, double sampling_rate,
   matrix1D<double> &frequency, matrix1D<double> &FSC) {
   double resolution=-1;
   frequency.clear();
   FSC.clear();

   // Prepare volumes to be evaluated
   if (vol_phantom.name()!="")
      system(((string)"ln -s "+vol_phantom.name()+" superfeo.vol").c_str());
   else vol_phantom.write("superfeo.vol");
   if (vol_recons.name()!="")
      system(((string)"ln -s "+vol_recons.name()+" superfeo2.vol").c_str());
   else vol_recons.write("superfeo2.vol");

   // Run bresolve of bsoft to compute resolution
   string command=(string)"xmipp_resolution -sam "+FtoA(sampling_rate)+
      " -ref superfeo.vol -i superfeo2.vol";
   system(command.c_str());

   // Process output file
   ifstream fh_resol;
   fh_resol.open("superfeo2.vol.frc");
   if (!fh_resol)
      REPORT_ERROR(1,
         "compute_FSC: Cannot open results file from xmipp_resolution");
   try {
      string line;
      // Skip first line
      getline(fh_resol,line);

      // Read the rest
      vector<double> auxfreq, auxFSC;
      vector<string> tokens;
      while (getline(fh_resol,line)) {
         tokenize(line,tokens);
         double f  =AtoF(tokens[0]);
         double fsc=AtoF(tokens[1]);
         auxfreq.push_back(f);
         auxFSC. push_back(fsc);
         if (resolution==-1 && fsc<0.5) resolution=f;
      }

      // Now copy
      frequency.resize(auxfreq.size());
      FSC.resize(auxfreq.size());
      FOR_ALL_ELEMENTS_IN_MATRIX1D(frequency) {
         frequency(i)=auxfreq[i];
         FSC      (i)=auxFSC [i];
      }
   } catch (Xmipp_error XE) {
      cerr << XE << endl
           << "compute_resolution: Error reading resolution, ignoring it\n";
   }
   fh_resol.close();

   // Remove intermidiate files
   system("rm superfeo.vol superfeo2.vol superfeo2.vol.frc superfeo2.vol.dpr");
   return resolution;
}
