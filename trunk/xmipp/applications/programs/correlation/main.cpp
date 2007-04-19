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

#include <data/progs.h>
#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>

class Similarity_parameters: public Prog_parameters {
public:
  FileName    fn_ref,fn_msk;
  ImageXmipp  refI,MI;
  VolumeXmipp refV,MV;
  matrix3D<int> mask3D;
  matrix2D<int> mask2D;
  bool usemask,docc,doeu,domi,doco;

public:
  Similarity_parameters() {}
  void read(int argc, char **argv) {
    Prog_parameters::read(argc,argv);
    try {
      usemask=false;
      docc=true; doeu=true; domi=true;
      fn_ref=get_param(argc,argv,"-ref");
      if (Is_ImageXmipp(fn_ref)) {
	refI.read(fn_ref,false,false,apply_geo);
	refI().set_Xmipp_origin();
	fn_msk=get_param(argc,argv,"-mask","");
	if (fn_msk!="") {
	  usemask=true;
	  MI.read(fn_msk,false,false,apply_geo);
	  MI().set_Xmipp_origin();
	  mask2D.resize(MI());
	  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(MI()) {
	    dMij(mask2D,i,j)=(int)dMij(MI(),i,j);
	  }
	  MI().core_deallocate();
	}
      } else if (Is_VolumeXmipp(fn_ref)) {
	refV.read(fn_ref);
	refV().set_Xmipp_origin();
	fn_msk=get_param(argc,argv,"-mask","");
	if (fn_msk!="") {
	  usemask=true;
	  MV.read(fn_msk);
	  MV().set_Xmipp_origin();
	  mask3D.resize(MV());
	  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(MV()) {
	    dVkij(mask3D,k,i,j)=(int)dVkij(MV(),k,i,j);
	  }
	  MV().core_deallocate();
	}
      } else REPORT_ERROR(1,"Reference is not an image or a volume");
      if (check_param(argc,argv,"-co")) {domi=false; doeu=false; docc=false;}
      if (check_param(argc,argv,"-cc")) {domi=false; doeu=false; doco=false;}
      if (check_param(argc,argv,"-mi")) {docc=false; doeu=false; doco=false;}
      if (check_param(argc,argv,"-eu")) {domi=false; docc=false; doco=false;}

   } catch (Xmipp_error XE) {cout << XE; usage(); exit(1);}
  }

  void show() {
    cout << "Reference file: " << fn_ref << endl;
    if (usemask) cout << "mask file: " << fn_msk << endl;
    Prog_parameters::show();
    cout <<endl;
  }
  void usage() {
    cerr << "   -ref <input file>        : Filename for reference image/volume \n";
    Prog_parameters::usage();
    cerr << "  [-mask <input mask>]      : Restrict similarity calculation to region within the mask\n";
    cerr << "  [-co ]                    : Only calculate correlation (i.e. signal product) \n";
    cerr << "  [-cc ]                    : Only calculate cross-correlation coefficient \n";
    cerr << "  [-eu ]                    : Only calculate euclidian distance \n";
    cerr << "  [-mi ]                    : Only calculate mutual information\n";
  }
};


void process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Similarity_parameters *eprm=(Similarity_parameters *) prm;

   double co,cc,eu,mi;
   if (!eprm->usemask) {
     if (eprm->doco) co=correlation(eprm->refI(),img());
     if (eprm->docc) cc=correlation_index(eprm->refI(),img());
     if (eprm->doeu) eu=euclidian_distance(eprm->refI(),img());
     if (eprm->domi) mi=mutual_information(eprm->refI(),img());
   } else {
     if (eprm->doco) co=correlation(eprm->refI(),img(),&eprm->mask2D);
     if (eprm->docc) cc=correlation_index(eprm->refI(),img(),&eprm->mask2D);
     if (eprm->doeu) eu=euclidian_distance(eprm->refI(),img(),&eprm->mask2D);
     if (eprm->domi) mi=mutual_information(eprm->refI(),img(),0,0,&eprm->mask2D);
   }

   cout << img.name()<<": ";
   if (eprm->doco) cout << " co= "<<co;
   if (eprm->docc) cout << " cc= "<<cc;
   if (eprm->doeu) cout << " eu= "<<eu;
   if (eprm->domi) cout << " mi= "<<mi;
   cout <<endl;
}


void process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Similarity_parameters *eprm=(Similarity_parameters *) prm;

   double co,cc,eu,mi;
   if (!eprm->usemask) {
     if (eprm->doco) co=correlation(eprm->refV(),vol());
     if (eprm->docc) cc=correlation_index(eprm->refV(),vol());
     if (eprm->doeu) eu=euclidian_distance(eprm->refV(),vol());
     if (eprm->domi) mi=mutual_information(eprm->refV(),vol());
   } else {
     if (eprm->doco) co=correlation(eprm->refV(),vol(),&eprm->mask3D);
     if (eprm->docc) cc=correlation_index(eprm->refV(),vol(),&eprm->mask3D);
     if (eprm->doeu) eu=euclidian_distance(eprm->refV(),vol(),&eprm->mask3D);
     if (eprm->domi) mi=mutual_information(eprm->refV(),vol(),0,0,&eprm->mask3D);
   }
   cout << vol.name()<<": ";
   if (eprm->docc) cout << " co= "<<co;
   if (eprm->docc) cout << " cc= "<<cc;
   if (eprm->doeu) cout << " eu= "<<eu;
   if (eprm->domi) cout << " mi= "<<mi;
   cout <<endl;

}


int main (int argc, char **argv) {
   Similarity_parameters prm;
   prm.allow_time_bar=false;
   prm.each_image_produces_an_output=false;
   prm.apply_geo=true;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}


