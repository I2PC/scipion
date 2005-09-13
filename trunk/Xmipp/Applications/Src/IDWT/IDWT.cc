/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippWavelets.hh>

class IDWT_parameters: public Prog_parameters {
   string IDWT_type;
public:
   void read(int argc, char **argv) {
      Prog_parameters::read(argc,argv);
      IDWT_type=get_param(argc,argv,"-type","DAUB12");
      if      (IDWT_type=="DAUB4")  set_DWT_type(DAUB4);
      else if (IDWT_type=="DAUB12") set_DWT_type(DAUB12);
      else if (IDWT_type=="DAUB20") set_DWT_type(DAUB20);
      else REPORT_ERROR(1,"DWT::read: Unknown DWT type");
   }
   
   void show() {
      Prog_parameters::show();
      cout << "IDWT type: " << IDWT_type << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "    [-type <str=\"DAUB12\"]    : DWT type. Valid types are:\n"
           << "                                   DAUB4, DAUB12, DAUB20\n";
   }
};

void process_img(ImageXmipp &img, const Prog_parameters *prm) {
   IDWT_parameters *eprm=(IDWT_parameters *) prm;
   IDWT(img(),img());
}

void process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   IDWT_parameters *eprm=(IDWT_parameters *) prm;
   IDWT(vol(),vol());
}

int main (int argc, char **argv) {
   IDWT_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

