/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
 *              Carlos Oscar Sánchez Sorzano (coss@cnb.uam.es)
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

#ifdef _HAVE_VTK
#include "../Prog_adjust_CTF.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippHistograms.hh>
#include <vector.h>   // vector from STL library

/* Number of CTF parameters */
#define CTF_PARAMETERS             23
#define PARAMETRIC_CTF_PARAMETERS  13
#define BACKGROUND_CTF_PARAMETERS  10
#define SQRT_CTF_PARAMETERS         3
#define GAUSSIAN_CTF_PARAMETERS     7

// Global variables -------------------------------------------------------- */
Adjust_CTF_Parameters *global_Adjust_CTF_parameters;
matrix2D<double> *global_ctf_ampl;
matrix2D<double> *global_ctf_ampl2;
matrix2D<int>    *global_ctf_mask;
matrix1D<double> *global_adjust;
vtkImageData     *global_ctf;

histogram1D      *global_hist;
double            global_value_th;

vector<int>       global_minima_X,global_minima_Y;

double            global_min_val;
double            global_max_val;
double            global_min_freq;
double            global_max_freq;
double            global_gamma;
double            global_Tm;
bool              global_penalize;
double            global_penalty;
double            global_current_penalty;
double            global_central_weight;
double            global_current_central_weight;
int               global_action; // 0=adjust background
      	             	      	 // 1=adjust CTF
				 // 2=adjust both
                                 // 3=adjust CTF with few parameters

int               global_show;   // Force CTF_fitness to show debugging info
bool              global_compute_FFT_distance;

// Penalization for forbidden values of the parameters
double		  global_heavy_penalization; 

// Additional debug parameters
int		  global_evaluation_reduction;

// CTF model and noise model
XmippCTF          global_ctfmodel;
XmippCTF         *global_initial_ctfmodel;

/* Assign ctfmodel from a vector and viceversa ----------------------------- */
void assign_CTF_from_parameters(double *p, XmippCTF &ctfmodel,
   int ia, int l, bool astigmatic_noise) {
   ctfmodel.Tm=global_Tm; 
   if (ia<= 0 && l>0) {ctfmodel.K              =p[ 0]; l--;}
   if (ia<= 1 && l>0) {ctfmodel.kV             =p[ 1]; l--;}
   if (ia<= 2 && l>0) {ctfmodel.DeltafU        =p[ 2]; l--;}
   if (ia<= 3 && l>0) {ctfmodel.DeltafV        =p[ 3]; l--;}
   if (ia<= 4 && l>0) {ctfmodel.azimuthal_angle=p[ 4]; l--;}
   if (ia<= 5 && l>0) {ctfmodel.Cs             =p[ 5]; l--;}
   if (ia<= 6 && l>0) {ctfmodel.Ca             =p[ 6]; l--;}
   if (ia<= 7 && l>0) {ctfmodel.espr           =p[ 7]; l--;}
   if (ia<= 8 && l>0) {ctfmodel.ispr           =p[ 8]; l--;}
   if (ia<= 9 && l>0) {ctfmodel.alpha          =p[ 9]; l--;}
   if (ia<=10 && l>0) {ctfmodel.DeltaF         =p[10]; l--;}
   if (ia<=11 && l>0) {ctfmodel.DeltaR         =p[11]; l--;}
   if (ia<=12 && l>0) {ctfmodel.Q0             =p[12]; l--;}
   if (ia<=13 && l>0) {ctfmodel.gaussian_K     =p[13]; l--;}
   if (ia<=14 && l>0) {ctfmodel.base_line      =p[14]; l--;}
   if (ia<=15 && l>0) {ctfmodel.sigmaU         =p[15]; l--;}
   if (ia<=16 && l>0)
      if (astigmatic_noise) {ctfmodel.sigmaV=p[16]; l--;}
      else                  {ctfmodel.sigmaV=p[15]; l--;}
   if (ia<=17 && l>0)
      if (astigmatic_noise) {ctfmodel.gaussian_angle=p[17]; l--;}
      else                  {ctfmodel.gaussian_angle=0;     l--;}
   if (ia<=18 && l>0) {ctfmodel.sqU            =p[18]; l--;}
   if (ia<=19 && l>0)
      if (astigmatic_noise) {ctfmodel.sqV           =p[19]; l--;}
      else                  {ctfmodel.sqV           =p[18]; l--;}
   if (ia<=20 && l>0) {ctfmodel.sqrt_K         =p[20]; l--;}
   if (ia<=21 && l>0) {ctfmodel.cU             =p[21]; l--;}
   if (ia<=22 && l>0)
      if (astigmatic_noise) {ctfmodel.cV            =p[22]; l--;}
      else                  {ctfmodel.cV            =p[21]; l--;}
//   cout << "Estoy 1:" << ctfmodel;
}

void assign_parameters_from_CTF(XmippCTF &ctfmodel, double *p,
   int ia, int l, bool astigmatic_noise) {
   if (ia<= 0 && l>0) {p[ 0]=ctfmodel.K; l--;}
   if (ia<= 1 && l>0) {p[ 1]=ctfmodel.kV; l--;}
   if (ia<= 2 && l>0) {p[ 2]=ctfmodel.DeltafU; l--;}
   if (ia<= 3 && l>0) {p[ 3]=ctfmodel.DeltafV; l--;}
   if (ia<= 4 && l>0) {p[ 4]=ctfmodel.azimuthal_angle; l--;}
   if (ia<= 5 && l>0) {p[ 5]=ctfmodel.Cs; l--;}
   if (ia<= 6 && l>0) {p[ 6]=ctfmodel.Ca; l--;}
   if (ia<= 7 && l>0) {p[ 7]=ctfmodel.espr; l--;}
   if (ia<= 8 && l>0) {p[ 8]=ctfmodel.ispr; l--;}
   if (ia<= 9 && l>0) {p[ 9]=ctfmodel.alpha; l--;}
   if (ia<=10 && l>0) {p[10]=ctfmodel.DeltaF; l--;}
   if (ia<=11 && l>0) {p[11]=ctfmodel.DeltaR; l--;}
   if (ia<=12 && l>0) {p[12]=ctfmodel.Q0; l--;}
   if (ia<=13 && l>0) {p[13]=ctfmodel.gaussian_K; l--;}
   if (ia<=14 && l>0) {p[14]=ctfmodel.base_line; l--;}
   if (ia<=15 && l>0) {p[15]=ctfmodel.sigmaU; l--;}
   if (ia<=16 && l>0)
      if (astigmatic_noise) {p[16]=ctfmodel.sigmaV; l--;}
      else                  {p[16]=0;               l--;}
   if (ia<=17 && l>0)
      if (astigmatic_noise) {p[17]=ctfmodel.gaussian_angle; l--;}
      else                  {p[17]=0;                       l--;}
   if (ia<=18 && l>0) {p[18]=ctfmodel.sqU; l--;}
   if (ia<=19 && l>0)
      if (astigmatic_noise) {p[19]=ctfmodel.sqV; l--;}
      else                  {p[19]=0;            l--;}
   if (ia<=20 && l>0) {p[20]=ctfmodel.sqrt_K; l--;}
   if (ia<=21 && l>0) {p[21]=ctfmodel.cU; l--;}
   if (ia<=22 && l>0)
      if (astigmatic_noise) {p[22]=ctfmodel.cV; l--;}
      else                  {p[22]=0;           l--;}
//   cout << "Estoy 2:" << ctfmodel;
}

/* Read parameters --------------------------------------------------------- */
void Adjust_CTF_Parameters::read(const FileName &fn_param) _THROW {
   adjust.resize(CTF_PARAMETERS);
   ctfmodel.enable_CTF=ctfmodel.enable_CTFnoise=TRUE; 
   ctfmodel.read(fn_param,FALSE);
   Tm=ctfmodel.Tm; 
   assign_parameters_from_CTF(ctfmodel,VEC_ARRAY(adjust),0,CTF_PARAMETERS,
      TRUE); 
   
   FILE *fh_param; 
   if ((fh_param = fopen(fn_param.c_str(), "r")) == NULL) 
   	 REPORT_ERROR(1,(string)"Prog_Adjust_CTF::read: There is a problem " 
            "opening the file "+fn_param); 
   
   fn_ctf=get_param(fh_param,"ctf",0,"");
   fn_outroot=get_param(fh_param,"output_root",0,"");
      if (fn_outroot=="") fn_outroot=fn_ctf.without_extension();
   fn_out_CTF_parameters=get_param(fh_param,"out_CTF_parameters",0,"");
      if (fn_out_CTF_parameters=="")
         fn_out_CTF_parameters=fn_outroot+"_model.param";
   
   show_optimization=check_param(fh_param,"show_optimization"); 
   do_not_optimize=check_param(fh_param,"do_not_optimize"); 
   do_not_fine_search=check_param(fh_param,"do_not_fine_search"); 
   penalty=AtoF(get_param(fh_param,"penalty",0,"128")); 
   central_weight=AtoF(get_param(fh_param,"central_weight",0,"1")); 
   min_freq=AtoF(get_param(fh_param,"min_freq",0,"0.05")); 
   max_freq=AtoF(get_param(fh_param,"max_freq",0,"0.45")); 
   gamma=AtoF(get_param(fh_param,"gamma",0,"3"));
   accuracy=AtoF(get_param(fh_param,"accuracy",0,"1"));
      accuracy/=100; 
   evaluation_reduction=AtoI(get_param(fh_param,"eval_red",0,"4")); 
   astigmatic_noise=check_param(fh_param,"astigmatic_noise");
   
   if (check_param(fh_param,"steps")) 
      steps=get_vector_param(fh_param,"steps",PARAMETRIC_CTF_PARAMETERS); 
   else { 
      steps.resize(PARAMETRIC_CTF_PARAMETERS); 
      steps.init_constant(1); 
      steps(1)=steps(5)=steps(6)=steps(7)=steps(8)=steps(10)=steps(11)=0; 
   } 
   
   value_th=AtoF(get_param(fh_param,"value_th",0,"80")); 
}

/* Write to a file --------------------------------------------------------- */
void Adjust_CTF_Parameters::write(const FileName &fn_prm, bool rewrite)
   _THROW {
   ofstream fh_param;
   if (!rewrite) fh_param.open(fn_prm.c_str(),ios::app);
   else          fh_param.open(fn_prm.c_str(),ios::out);
   if (!fh_param)
      REPORT_ERROR(1,(string)"Adjust_CTF_Parameters::write: There is a problem "
            "opening the file "+fn_prm);
   fh_param << "# Adjust CTF parameters\n";
   if (fn_ctf!="")
      fh_param << "ctf="               << fn_ctf                  << endl;
   if (fn_outroot!="")
      fh_param << "output_root="       << fn_outroot              << endl;
   if (fn_out_CTF_parameters!="" && fn_out_CTF_parameters!="_model.param")
      fh_param << "out_CTF_parameters="<< fn_out_CTF_parameters   << endl;
   fh_param << "penalty="              << penalty                 << endl
            << "central_weight="       << central_weight          << endl
            << "min_freq="             << min_freq                << endl
            << "max_freq="             << max_freq                << endl
            << "gamma="                << gamma                   << endl
            << "accuracy="             << 100*accuracy            << endl
            << "eval_red="             << evaluation_reduction    << endl
            << "value_th="             << value_th                << endl
   ;
   if (show_optimization)  fh_param << "show_optimization=yes\n";
   if (do_not_optimize)    fh_param << "do_not_optimize=yes\n";
   if (do_not_fine_search) fh_param << "do_not_fine_search=yes\n";
   if (astigmatic_noise)   fh_param << "astigmatic_noise=yes\n";

   fh_param << "steps=[" << steps(0);
   for (int i=1; i<PARAMETRIC_CTF_PARAMETERS; i++)
      fh_param << ", " << steps(i);
   fh_param << "]\n";

   fh_param << endl;

   fh_param << ctfmodel << endl;
   fh_param.close();
}


/* Show -------------------------------------------------------------------- */
void Adjust_CTF_Parameters::show() {
   cout << "CTF file:        " << fn_ctf << endl
	<< "Penalty:         " << penalty << endl
	<< "Central weight:  " << central_weight << endl
	<< "Min Freq.:       " << min_freq << endl
	<< "Max Freq.:       " << max_freq << endl
        << "Gamma:           " << gamma << endl
	<< "Accuracy:        " << accuracy << endl
	<< "Eval reduction:  " << evaluation_reduction << endl
	<< "Sampling:        " << Tm << endl
	<< "Threshold Sqrt:  " << value_th << endl
        << "Astigmatic noise:" << astigmatic_noise << endl
	<< "Starting:\n"       << ctfmodel << endl
	<< "Steps:	     " << steps.transpose() << endl
    ;
}

/* Usage ------------------------------------------------------------------- */
void Adjust_CTF_Parameters::Usage() {
   cerr << "This program tries to adjust a parametric model to a CTF file.\n"
        << "Usage: adjust_ctf -i <parameters file>\n"
        << "   Where the parameters file may contain the description of a\n"
        << "      CTF and CTFnoise plus any of the following parameters\n"
	<< "   [ctf=<Fourier Xmipp Image>] : CTF file\n"
        << "   [output_root=<root name>]   : Output files are generated using\n"
        << "                                 this root name. By default, the root\n"
        << "                                 name of the ctf= option is used\n"
        << "   [out_CTF_parameters=<fn_out_CTF>]: The adjusted CTF parameters\n"
        << "                                 are stored in this file. By default,\n"
        << "                                 the ctf= option rootname is used to generate\n"
        << "                                 the input file.\n"
	<< "   [steps=<vector>]            : An Xmipp vector of 13 elements saying what parameters \n"
	<< "                                 of the CTF to adjust (put 1) or not (put 0).\n"  
	<< "   [penalty=<f=128>]           : Penalty applied in background adjust to points whose\n"
	<< "                                 value is greater than the CTF\n"
	<< "   [central_weight=<w=1>]      : Weight for the first two CTF lobes\n"
	<< "   [eval_red=<n=4>]            : Use values of images taken every n pixels. As images are\n"
	<< "                                 usually huge. A value between 4 and 8 is recommended.\n"
        << "   [show_optimization=yes]     : Show optimization process\n"
	<< "   [min_freq=<f=0.05>]         : Minimum digital frequency to use in adjust. Its value\n"
	<< "                                 should be a little lower than the dig. freq. of the first \n"  
	<< "                                 CTF zero.\n"  
 	<< "   [max_freq=<f=0.45>]         : Maximum digital frequency to use in adjust.\n"
	<< "                                 It should be higher than the last zero of the CTF.\n"  
        << "   [gamma=<gamma=3>]           : Gamma correction used when adjusting\n"
        << "                                 the initial defoci without any hint\n"
	<< "   [accuracy=<acc=1>]          : Stop criteria for adjustment iterations\n"
	<< "   [do_not_optimize=yes]       : Only shows adjustment of initial point. Mainly to use\n"
	<< "                                 in combination with -save_images.\n"
	<< "   [do_not_fine_search=yes]    : Do not perform final fine search\n"
	<< "   [value_th=<th0=80>]         : values above this threshold do not contribute\n"
        << "   [astigmatic_noise=yes|no]   : By default, noise is not astigmatic\n"
   ;
}

/* Produce side information ------------------------------------------------ */
void Adjust_CTF_Parameters::produce_side_info() {    
   ctftomodel.FilterBand=ctftomodel.FilterShape=FROM_FILE;
   ctftomodel.fn_mask=fn_ctf;
   ctftomodel.generate_mask(NULL);
}

/* Generate model so far ---------------------------------------------------- */
/* The model is taken from global_adjust and global_ctfmodel is modified */
void generate_model_so_far(ImageXmipp &I, bool apply_log=FALSE) {
   matrix1D<int>    idx(2);  // Indexes for Fourier plane
   matrix1D<double> freq(2); // Frequencies for Fourier plane
   assign_CTF_from_parameters(VEC_ARRAY(*global_adjust),global_ctfmodel,
      0,CTF_PARAMETERS,global_Adjust_CTF_parameters->astigmatic_noise);
   global_ctfmodel.Produce_Side_Info();

   I().resize(*global_ctf_ampl);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(I()) {
      XX(idx)=j; YY(idx)=i;
      FFT_idx2digfreq(global_ctf, idx, freq);
      digfreq2contfreq(freq, freq, global_Tm);
      
      // Decide what to save
      double param;
      if (global_action==0)
         I(i,j)=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
      else {
         I(i,j)=global_ctfmodel.CTF_at(XX(freq),YY(freq));
         I(i,j)*=I(i,j);
      }
      if (apply_log) I(i,j)=log10(I(i,j));
   }
}


/* Save intermidiate results ----------------------------------------------- */
/* First call to generate model so far and then save the image, and a couple
   of cuts along X and Y */
void save_intermidiate_results(const FileName &fn_root) {
   ofstream plotX, plotY;
   ImageXmipp save;
   matrix1D<int>    idx(2);  // Indexes for Fourier plane
   matrix1D<double> freq(2); // Frequencies for Fourier plane
   generate_model_so_far(save,TRUE);
   plotX.open((fn_root+"X.txt").c_str());
   plotY.open((fn_root+"Y.txt").c_str());
   if (!plotX)
      REPORT_ERROR(1,"save_intermidiate_results::Cannot open plot file for writing\n");
   if (!plotY)
      REPORT_ERROR(1,"save_intermidiate_results::Cannot open plot file for writing\n");
   plotX << "# freq gaus ctf2\n";
   plotY << "# freq gaus ctf2\n";

   double w;
   for (int i=STARTINGY(save()); i<=FINISHINGY(save())/2; i++) {
      XX(idx)=0; YY(idx)=i;
      FFT_idx2digfreq(global_ctf, idx, freq); w=freq.module();
      if (w<global_min_freq || w>global_max_freq) continue;
      if (YY(freq)<0) continue;
      digfreq2contfreq(freq, freq, global_Tm);
      plotY << YY(freq) << " " << pow(10,save(i,0)) << " "
	    << (*global_ctf_ampl2)(i,0) << endl;
   }

   for (int j=STARTINGX(save()); j<=FINISHINGX(save())/2; j++) {
      XX(idx)=j; YY(idx)=0;
      FFT_idx2digfreq(global_ctf, idx, freq); w=freq.module();
      if (w<global_min_freq || w>global_max_freq) continue;
      if (XX(freq)<0) continue;
      digfreq2contfreq(freq, freq, global_Tm);
      plotX << XX(freq) << " " << pow(10,save(0,j)) << " "
	    << (*global_ctf_ampl2)(0,j) << endl;
   }
   save.write((fn_root+"_log10.xmp").c_str());
   plotX.close();
   plotY.close();
}

// Estimate sqrt parameters ------------------------------------------------
// Results are written in global_ctfmodel
void estimate_background_sqrt_parameters() {
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   double           w;                   // module of the frequency
   matrix2D<double> *f=global_ctf_ampl2; // A pointer to simplify expresions

   // Estimate the base line taking the value of the CTF for the maximum X and Y frequencies
   double base_line=0; int N=0;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(*f) {
      VECTOR_R2(idx,j,i);
      FFT_idx2digfreq(global_ctf, idx, freq); w=freq.module();
      if (w>0.4) {N++; base_line+=(*f)(i,j);}
   }
   global_ctfmodel.base_line=base_line/N;
   
   // Find the linear least squares solution for the sqrt part
   matrix2D<double> A(2,2); A.init_zeros();
   matrix1D<double> b(2);   b.init_zeros();
   FOR_ALL_ELEMENTS_IN_MATRIX2D(*f) {
      // Disregard this point?
      if ((*f)(i,j)>=global_value_th) continue;
      VECTOR_R2(idx,j,i);
      FFT_idx2digfreq(global_ctf, idx, freq); w=freq.module();
      if (w>=global_max_freq || w<=global_min_freq) continue;
      if (YY(freq)<0) continue;
      digfreq2contfreq(freq, freq, global_Tm); double wcont=freq.module();
      
      // Compute weight for this point
      int n;
      double weight;
      global_hist->val2index(w,n);
      if ((*global_hist)(n)!=0) weight=1/(*global_hist)(n);
      else continue;
      weight*=global_max_freq-w;

      // Compute error
      double explained=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
      double unexplained=(*f)(i,j)-explained;
      if (unexplained<=0) continue;
      unexplained=log(unexplained);

      double X=-sqrt(wcont);
      A(0,0)+=weight*X*X;
      A(0,1)+=weight*X;
      A(1,1)+=weight*1;
      b(0)  +=X*weight*unexplained;
      b(1)  +=  weight*unexplained;
   }
   A(1,0)=A(0,1);
   b=A.inv()*b;
   global_ctfmodel.sqU=global_ctfmodel.sqV=b(0);
   global_ctfmodel.sqrt_K=exp(b(1));
}

// Estimate gaussian parameters --------------------------------------------
void estimate_background_gauss_parameters() {
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   double           w;                   // module of the frequency
   matrix2D<double> *f=global_ctf_ampl2; // A pointer to simplify expresions

   // Estimate the gaussian centers
   VECTOR_R2(freq,0,global_min_freq+0.2*(global_max_freq-global_min_freq));
   digfreq2contfreq(freq, freq, global_Tm);
   global_ctfmodel.cV=YY(freq);
   global_ctfmodel.cV=0;

   VECTOR_R2(freq,global_min_freq+0.2*(global_max_freq-global_min_freq),0);
   digfreq2contfreq(freq, freq, global_Tm);
   global_ctfmodel.cU=XX(freq);
   global_ctfmodel.cU=0;

   // Find the linear least squares solution for the gauss part
   matrix2D<double> A(3,3); A.init_zeros();
   matrix1D<double> b(3);   b.init_zeros();
   FOR_ALL_ELEMENTS_IN_MATRIX2D(*f) {
      // Disregard this point?
      if ((*f)(i,j)>=global_value_th) continue;
      VECTOR_R2(idx,j,i);
      FFT_idx2digfreq(global_ctf, idx, freq); w=freq.module();
      if (w>=global_max_freq || w<=global_min_freq) continue;
      if (YY(freq)<0) continue;
      digfreq2contfreq(freq, freq, global_Tm);

      // Compute weight for this point
      int n;
      double weight;
      global_hist->val2index(w,n);
      if ((*global_hist)(n)!=0) weight=1/(*global_hist)(n);
      else continue;
      weight*=global_max_freq-w;

      // Compute error
      double explained=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
      double unexplained=(*f)(i,j)-explained;
      if (unexplained<=0) continue;
      unexplained=log(unexplained);
      double X=-XX(freq)*XX(freq);
      double Y=-YY(freq)*YY(freq);
      A(0,0)+=weight*X*X;
      A(0,1)+=weight*X*Y;
      A(0,2)+=weight*X;
      A(1,1)+=weight*Y*Y;
      A(1,2)+=weight*Y;
      A(2,2)+=weight*1;
      b(0)  +=X*weight*unexplained;
      b(1)  +=Y*weight*unexplained;
      b(2)  +=  weight*unexplained;
   }
   A(1,0)=A(0,1); A(2,0)=A(0,2); A(2,1)=A(1,2);
   b=A.inv()*b;
   global_ctfmodel.sigmaU=b(0);
   global_ctfmodel.sigmaV=b(1);
   global_ctfmodel.gaussian_K=exp(b(2));
}

/* Compute central region -------------------------------------------------- */
void compute_central_region(double &w1, double &w2, double ang) {
   w1=global_max_freq; w2=global_min_freq; 
   matrix1D<double> freq(2), dir(2);
   
   // Compute first and third zero in the given direction
   VECTOR_R2(dir,COSD(ang),SIND(ang));

   // Detect first zero
   global_ctfmodel.zero(1,dir,freq);
   if (XX(freq)==-1 && YY(freq)==-1) w1=global_min_freq;
   else {
      contfreq2digfreq(freq,freq,global_Tm);
      double w;
      if (XX(dir)>0.1) w=XX(freq)/XX(dir);
      else             w=YY(freq)/YY(dir);
      w1=MAX(global_min_freq,MIN(w1,w));
   }

   // Detect third zero
   global_ctfmodel.zero(4,dir,freq);
   if (XX(freq)==-1 && YY(freq)==-1) w2=global_max_freq;
   else {
      double w;
      contfreq2digfreq(freq,freq,global_Tm);
      if (XX(dir)>0.1) w=XX(freq)/XX(dir);
      else             w=YY(freq)/YY(dir);
      w2=MIN(global_max_freq,MAX(w2,w));
   }
}

/* CTF fitness ------------------------------------------------------------- */
/* All variables are taken as global since there is no other way
   of communicating through Powell function */
double CTF_fitness(double *p) {
   double retval=0;
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   ImageXmipp save;
   matrix2D<double> considered_ctf_ampl2, considered_ctfmodel_ampl2;
   
   // Show input variables for debugging
   if (global_show>=1) {
      int imax;
      switch (global_action) {
         case 0: imax=BACKGROUND_CTF_PARAMETERS; break;
         case 3:
         case 1: imax=PARAMETRIC_CTF_PARAMETERS; break;
         case 2: imax=CTF_PARAMETERS; break;
      }
      cout << "Input vector: ";
      for (int i=0; i<imax; i++) cout << p[i+1] << " ";
      cout << endl;
      cout << "Global vector: " << global_adjust->transpose() << endl;
   }

   if (global_show>=2) {
      save().resize(*global_ctf_ampl2);
   }

   // Do not allow big changes (>10%)
   #ifdef NEVER_DEFINED
   if (global_action==2) {
      for (int i=0; i<CTF_PARAMETERS; i++) {
         if ((*global_adjust)(i)!=0) {
            if (ABS((*global_adjust)(i)-p[i+1])/ABS((*global_adjust)(i))>0.10)
               return global_heavy_penalization;
         }
      }
   }
   #endif

   // Generate CTF model
   switch (global_action) {
      case 0: assign_CTF_from_parameters(p-12,global_ctfmodel,
                13,BACKGROUND_CTF_PARAMETERS,
                global_Adjust_CTF_parameters->astigmatic_noise); break;
      case 3:
      case 1: assign_CTF_from_parameters(p+1,global_ctfmodel,
                 0,PARAMETRIC_CTF_PARAMETERS,
                 global_Adjust_CTF_parameters->astigmatic_noise); break;
      case 2: assign_CTF_from_parameters(p+1,global_ctfmodel,
                 0,CTF_PARAMETERS,
                 global_Adjust_CTF_parameters->astigmatic_noise); break;
   }
   global_ctfmodel.enable_CTF=(global_action!=0);
   global_ctfmodel.enable_CTFnoise=TRUE;
   if (global_show>=1) cout << "Model:\n" << global_ctfmodel << endl;
   global_ctfmodel.Produce_Side_Info();
   if(!global_ctfmodel.physical_meaning()) return global_heavy_penalization;
   
   // Compute central region
   double w1U, w2U, w1V, w2V;
   if (global_current_central_weight!=1) {
      compute_central_region(w1U,w2U,global_ctfmodel.azimuthal_angle);
      compute_central_region(w1V,w2V,global_ctfmodel.azimuthal_angle+90);
   }
   else {w1U=w1V=global_min_freq; w2U=w2V=global_max_freq;}
   
   if (global_action==3) {w1U=w1V=global_min_freq;}
   if (global_show>=1)
      cout << "U=[" << w1U << "," << w2U << "]\n"
           << "V=[" << w1V << "," << w2V << "]\n";

   double N=0;
   double wctf2_sum=0, wparam_sum=0, wctf2_param_sum=0, w_sum=0, w_sum2=0,
          wctf2_sum2=0, wparam_sum2=0;
//   int step=(global_action!=3) ? global_evaluation_reduction:2;
   int step=global_evaluation_reduction;
   
   if (global_action!=0 && global_compute_FFT_distance) {
      considered_ctfmodel_ampl2.init_zeros(CEIL((double)YSIZE(*global_ctf_ampl2)/step),
         CEIL((double)XSIZE(*global_ctf_ampl2)/step));
      considered_ctf_ampl2.init_zeros(considered_ctfmodel_ampl2);
   }
   for(int i=STARTINGY(*global_ctf_ampl2), ii=0;i<=FINISHINGY(*global_ctf_ampl2);i+=step, ii++)
       for(int j=STARTINGX(*global_ctf_ampl2), jj=0;j<=FINISHINGX(*global_ctf_ampl2);j+=step, jj++)
       {
          // Disregard this point?
    	  double ctf2=(*global_ctf_ampl2)(i,j);
      	  if (ctf2>=global_value_th) continue;
    	  XX(idx)=j; YY(idx)=i;
    	  FFT_idx2digfreq(global_ctf, idx, freq);
          if (YY(freq)<0) continue;
    	  double w=freq.module();
	  if (w<=global_min_freq || w>=global_max_freq) continue;
    	  digfreq2contfreq(freq, freq, global_Tm);
	  
	  // Compute each component
          double gaus=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
     	  double param;
	  if (global_action!=0)
             param=global_ctfmodel.CTFpure_at(XX(freq),YY(freq));
	  else param=0;
    	  param=gaus+param*param;

          // Compute weight
      	  double weight;
          if (global_current_central_weight!=1) {
             double ellipsoid_ang=
                atan2(YY(freq),XX(freq))-global_ctfmodel.rad_azimuth;
             double w1Up=w1U*cos(ellipsoid_ang);
             double w1Vp=w1V*sin(ellipsoid_ang);
             double w2Up=w2U*cos(ellipsoid_ang);
             double w2Vp=w2V*sin(ellipsoid_ang);
             double w1=sqrt(w1Up*w1Up+w1Vp*w1Vp);
             double w2=sqrt(w2Up*w2Up+w2Vp*w2Vp);
             if (w<w1)      weight=1.0/global_current_central_weight-1+
                                   exp(log(2.0-1.0/global_current_central_weight)
                                   /(w1-global_min_freq)*(w-global_min_freq));
             else if (w<w2) weight=1;
             else           weight=exp(-log(1.0/global_current_central_weight)
                                   /(w2-global_max_freq)*(w-w2));
          } else
             weight=1;
          // Give more weight to low frequencies
          if (global_gamma!=1)
               weight*=pow(global_max_freq*global_max_freq-w*w,global_gamma);
          else weight*=global_max_freq*global_max_freq-w*w;
          if (global_show>=2) save(i,j)=weight;

      	  // Compute gaussian distance
          double dist_gaus=ABS(ctf2-gaus);
	  if (global_action!=0) dist_gaus/=ctf2;
          if (gaus>ctf2 && global_penalize)
	     dist_gaus*=weight*global_current_penalty*(global_max_freq-w)/
        	(global_max_freq-global_min_freq);

      	  // Compute parametric distance
	  double dist_param;
	  if (global_action==0) dist_param=0;
	  else                  dist_param=ABS(ctf2-param)/ctf2;

      	  // Compute distance
	  double dist=weight*(dist_gaus+dist_param);
	  retval+=dist;
          N++;
	  
	  // Compute correlation contribution
          // Write these values in the considered image
          // for posterior Fourier distance
      	  if (global_action!=0) {
             wctf2_sum2     +=weight*ctf2*ctf2;
             wparam_sum2    +=weight*param*param;
	     wctf2_param_sum+=weight*ctf2*param;
             wctf2_sum      +=weight*ctf2;
             wparam_sum     +=weight*param;
             w_sum          +=weight;

             if (global_compute_FFT_distance) {
                considered_ctfmodel_ampl2(ii,jj)=weight*param;
                considered_ctf_ampl2(ii,jj)=weight*ctf2;
             }
      	  }
       }
   retval/=N;
   if (global_show>=2) cout << "Euclidean error " << retval << endl;

   // Compute correlation distance
   double wcorr, wcovariance, wctf2_var, wparam_var, wctf2_avg, wparam_avg;
   if (global_action!=0) {
      wctf2_avg =wctf2_sum/w_sum;
      wparam_avg=wparam_sum/w_sum;
      wctf2_var =(wctf2_sum2 -wctf2_avg *wctf2_avg *w_sum)/(w_sum-1);
      wparam_var=(wparam_sum2-wparam_avg*wparam_avg*w_sum)/(w_sum-1);
      wcovariance=(wctf2_param_sum-wctf2_avg*wparam_sum-
         wparam_avg*wctf2_sum+wctf2_avg*wparam_avg*w_sum)/(w_sum-1);
      wcorr=wcovariance/sqrt(wctf2_var*wparam_var);
      if (wcorr<0) {
         if (global_show>=1)
            cout << "Negative correlation: " << wcorr << endl;
         return global_heavy_penalization;
      }

      retval+=1-wcorr;
   if (global_show>=2) cout << "Correlation error " << 1-wcorr << endl;
   }

   // Penalize astigmatism
//   if (global_action!=0)
//      retval*=1+0.5*
//         ABS((*global_adjust)(3)-(*global_adjust)(2))/
//         MAX(ABS((*global_adjust)(3)),ABS((*global_adjust)(2)));

   // Compute distance in Fourier Space
   double fft_error=0;
   if (global_action!=0 && global_compute_FFT_distance) {
      // Compute the magnitude of both CTF images
      vtkImageData *fft=NULL;
      ImageXmipp fft_ampl, fftmodel_ampl;
      fft_ampl()=considered_ctf_ampl2;
      considered_ctf_ampl2-=considered_ctf_ampl2.compute_avg();
      FFT_VTK(considered_ctf_ampl2,fft,TRUE);
      FFT_magnitude(fft,fft_ampl());
      fft_ampl().self_log10();
      fft_ampl().set_Xmipp_origin();
 
      considered_ctfmodel_ampl2-=considered_ctfmodel_ampl2.compute_avg();
      FFT_VTK(considered_ctfmodel_ampl2,fft,TRUE);
      FFT_magnitude(fft,fftmodel_ampl());
      fftmodel_ampl().self_log10();
      fftmodel_ampl().set_Xmipp_origin();
      fft->Delete();

      // Compute the differences
      int imax=YSIZE(fftmodel_ampl())/4;
      int jmax=XSIZE(fftmodel_ampl())/4;
      N=0;
      for (int i=0; i<=imax; i++)
          for (int j=-jmax; j<=jmax; j++, N++)
             if (ABS(fft_ampl(i,j))>1e-2)
                fft_error+=ABS(fftmodel_ampl(i,j)-fft_ampl(i,j))/
                   ABS(fft_ampl(i,j));
      retval+=fft_error/N;
      if (global_show>=2) {
         cout << "FFT error " << fft_error/N << endl;
         cout << "Model: ";fftmodel_ampl().print_stats(); cout << endl;
         cout << "Org: "; fft_ampl().print_stats(); cout << endl;
      }
   }

   // Show some debugging information?
   if (global_show>=2) {
      save.write("weight.xmp");
      cout << "Score: " << retval << endl;
      if (global_action!=0)
         cout << "Correlation index: " << wcorr << endl
              << wctf2_avg << " " << wparam_avg << " " << wcovariance << " "
              << wctf2_var << " " << wparam_var << endl;
      matrix1D<double> aux; aux=*global_adjust;
      if (global_action==0) {
        for (int i=1; i<=BACKGROUND_CTF_PARAMETERS; i++)
	   (*global_adjust)(12+i)=p[i];
      } else if (global_action==1 || global_action==3 || global_action==4) {
        for (int i=1; i<=PARAMETRIC_CTF_PARAMETERS; i++)
	   (*global_adjust)(i-1)=p[i];
      } else {
        for (int i=1; i<=CTF_PARAMETERS; i++)
	   (*global_adjust)(i-1)=p[i];
      }
      save_intermidiate_results("inter");
      *global_adjust=aux;
      cout << "Press any key\n"; char c; cin >> c;
   }

   return retval;
}

// Compute CTF stats -------------------------------------------------------
/* Compute image minimum and maximum and the number of points per frequency */
void CTF_stats() {
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   matrix2D<double> *f=global_ctf_ampl2; // Pointer for simplifying expressions
   bool first=TRUE;
   
   // Compute minimum and maximum
   for(int i=STARTINGY(*f);i<=FINISHINGY(*f);i+=global_evaluation_reduction)
       for(int j=STARTINGX(*f);j<=FINISHINGX(*f);j+=global_evaluation_reduction)
       {
    	  XX(idx)=j; YY(idx)=i;
    	  FFT_idx2digfreq(global_ctf, idx, freq);
    	  double w=freq.module();
	  if (w>=global_min_freq && w<=global_max_freq) {
	     if (first) {global_min_val=global_max_val=(*f)(i,j); first=FALSE;}
	     else {
	        global_min_val=MIN(global_min_val,(*f)(i,j));
	        global_max_val=MAX(global_max_val,(*f)(i,j));
	     }
	  }
       }

   // Compute histogram
   global_hist->clear();
   global_hist->init(global_min_freq,global_max_freq,25);
   for(int i=STARTINGY(*f);i<=FINISHINGY(*f);i+=global_evaluation_reduction)
       for(int j=STARTINGX(*f);j<=FINISHINGX(*f);j+=global_evaluation_reduction)
       {
    	  XX(idx)=j; YY(idx)=i;
    	  FFT_idx2digfreq(global_ctf, idx, freq);
    	  double w=freq.module();
	  if (w<global_min_freq || w>global_max_freq) continue;
	  if (YY(freq)<0) continue;
	  global_hist->insert_value(w);
       }
   (*global_hist)/=global_hist->compute_max();
}

/* Estimate best paramtric CTF when the defoci are fixed ------------------- */
void estimate_best_parametric_CTF_for_defoci() {
   matrix1D<int>    idx(2);
   matrix1D<double> freq(2);
   matrix2D<double> *f=global_ctf_ampl2;
   
   // Prepare a parametric CTF
   assign_CTF_from_parameters(VEC_ARRAY(*global_adjust),global_ctfmodel,
      0,CTF_PARAMETERS,global_Adjust_CTF_parameters->astigmatic_noise);
   global_ctfmodel.Produce_Side_Info();
   global_ctfmodel.K=1;

   matrix2D<double> A(2,2); A.init_zeros();
   matrix1D<double> b(2);   b.init_zeros();
   for(int i=STARTINGY(*f);i<=FINISHINGY(*f);i+=global_evaluation_reduction)
       for(int j=STARTINGX(*f);j<=FINISHINGX(*f);j+=global_evaluation_reduction) {
      // Disregard this point?
      if ((*f)(i,j)>=global_value_th) continue;
      VECTOR_R2(idx,j,i);
      FFT_idx2digfreq(global_ctf, idx, freq); double w=freq.module();
      if (w>=global_max_freq || w<=global_min_freq) continue;
      if (YY(freq)<0) continue;
      digfreq2contfreq(freq, freq, global_Tm); double wcont=freq.module();
      
      // Compute weight for this point
      int n;
      double weight;
      global_hist->val2index(w,n);
      if ((*global_hist)(n)!=0) weight=1/(*global_hist)(n);
      else continue;
      weight*=global_max_freq-w;

      // Compute error
      double explained=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
      double unexplained=(*f)(i,j)-explained;
      double plain_ctf=global_ctfmodel.CTFpure_at(XX(freq),YY(freq));
      plain_ctf*=plain_ctf;
      if (unexplained<=0) continue;
      if (plain_ctf<0.2) continue;
      
      unexplained=0.5*log(unexplained/plain_ctf);
      double wcont3=wcont*wcont*wcont;
      double X=1;
      double Y=PI*(global_ctfmodel.Cs*global_ctfmodel.lambda*
                   global_ctfmodel.lambda*wcont3+
                   global_ctfmodel.Deltaf(XX(freq),YY(freq))*wcont);
      Y=-Y*Y;
      A(0,0)+=weight*X*X;
      A(0,1)+=weight*X*Y;
      A(1,1)+=weight*Y*Y;
      b(0)  +=X*weight*unexplained;
      b(1)  +=Y*weight*unexplained;
   }
   A(1,0)=A(0,1);
   b=A.inv()*b;

   (*global_adjust)( 0)=exp(b(0));                 // CTF gain
   if (b(1)>0) (*global_adjust)(9)=sqrt(b(1))*1e3; // alpha (mrad)
   else        (*global_adjust)(9)=0;
}

/* Estimate parametric CTF ------------------------------------------------- */
void estimate_best_parametric_CTF(double (*fitness)(double *p)) {
    double best_defocusU, best_defocusV, best_error;
    bool first=TRUE;
    int i,j;
    double defocusV, defocusU;
    double *p=global_adjust->adapt_for_numerical_recipes();

    double defocusV0=-1e3, defocusU0=-1e3;
    double defocusVF=-100e3, defocusUF=-100e3;
    matrix2D<double> error;

    for (double defocusStep=16e3; defocusStep>=1e3; defocusStep/=2) {
       error.resize(CEIL((defocusV0-defocusVF)/defocusStep),
                    CEIL((defocusU0-defocusUF)/defocusStep));
       if (global_show>=1)
          cout << "V=["<<defocusV0 << "," << defocusVF << "]\n"
               << "U=["<<defocusU0 << "," << defocusUF << "]\n" << defocusStep << endl;
       for (defocusV=defocusV0,i=0; defocusV>defocusVF; defocusV-=defocusStep, i++) {
          for (defocusU=defocusU0,j=0; defocusU>defocusUF; defocusU-=defocusStep, j++) {
	      // Determine the distance between the parametric CTF
	      // (with the defocus and the K values) and the CTF file
              // Select all parameters except defoci
              assign_parameters_from_CTF(*global_initial_ctfmodel,
                 VEC_ARRAY(*global_adjust),0,PARAMETRIC_CTF_PARAMETERS,
                 global_Adjust_CTF_parameters->astigmatic_noise);
	      (*global_adjust)(2)=defocusU;
	      (*global_adjust)(3)=defocusV;

              estimate_best_parametric_CTF_for_defoci();
              error(i,j)=(*fitness)(p);

	      if(error(i,j)<best_error || first) {
	         best_error=error(i,j);
                 best_defocusU=defocusU; best_defocusV=defocusV;
	         first=FALSE;
	         cout << "    (DefocusU,DefocusV)=(" << defocusU << ","
                      << defocusV << ") --> " << error(i,j) << endl;
	      }
          }
       }

       // Find those defoci which are within a 10% of the maximum
       double errmin=error(0,0), errmax=error(0,0);
       for (int ii=STARTINGY(error); ii<=FINISHINGY(error); ii++)
          for (int jj=STARTINGX(error); jj<=FINISHINGX(error); jj++) {
             if (error(ii,jj)!=global_heavy_penalization)
                if (error(ii,jj)<errmin) errmin=error(ii,jj);
                else if (errmax==global_heavy_penalization) errmax=error(ii,jj);
                else if (error(ii,jj)>errmax) errmax=error(ii,jj);
          }

       if (global_show>=1) cout << "Range=" << errmax -errmin << endl;
       double best_defocusVmin=best_defocusV, best_defocusVmax=best_defocusV;
       double best_defocusUmin=best_defocusU, best_defocusUmax=best_defocusU;
       for (defocusV=defocusV0,i=0; defocusV>defocusVF; defocusV-=defocusStep, i++) {
          for (defocusU=defocusU0,j=0; defocusU>defocusUF; defocusU-=defocusStep, j++) {
            if (global_show>=1)
                cout << i << "," << j << " " << error(i,j) << " " << defocusU << " " << defocusV << endl
                     << best_defocusUmin << " " << best_defocusUmax << endl
                     << best_defocusVmin << " " << best_defocusVmax << endl;
             if ((error(i,j)-errmin)/(errmax-errmin)<=0.1) {
                if (defocusV<best_defocusVmin) best_defocusVmin=defocusV;
                if (defocusU<best_defocusUmin) best_defocusUmin=defocusU;
                if (defocusV>best_defocusVmax) best_defocusVmax=defocusV;
                if (defocusU>best_defocusUmax) best_defocusUmax=defocusU;
             }
          }
       }

       defocusVF=MAX(-100e3,best_defocusVmin-defocusStep);
       defocusV0=MIN(-1e3,best_defocusVmax+defocusStep);
       defocusUF=MAX(-100e3,best_defocusUmin-defocusStep);
       defocusU0=MIN(-1e3,best_defocusUmax+defocusStep);
       i=j=0;
       if (global_show>=1) {
          ImageXmipp save; save()=error; save.write("error.xmp");
          cout << "Press any key: Error saved\n";
          char c; cin >> c;
       }
    }

    global_adjust->kill_adaptation_for_numerical_recipes(p);
    (*global_adjust)(2)=best_defocusU;
    (*global_adjust)(3)=best_defocusV;
    estimate_best_parametric_CTF_for_defoci();
}

/* Main routine ------------------------------------------------------------ */
//#define DEBUG
void ROUT_Adjust_CTF(Adjust_CTF_Parameters &prm) {
   matrix2D<double> ctf_ampl,ctf_ampl2; // To store the amplitude images
   ImageXmipp save, gauss;  // To save and keep the gaussian images
   global_Adjust_CTF_parameters=&prm;

   // Show starting parameters
   prm.show();
   global_show=0;
   
   // Read the CTF file
   FourierImageXmipp I;
   I.read(prm.fn_ctf);
   xmippFFT2VTK(I,prm.ctftomodel.mask);

   // Get the amplitude of the cft
   FFT_magnitude(prm.ctftomodel.mask, ctf_ampl);

   // Get the squared amplitude
   ctf_ampl2.resize(ctf_ampl);   
   FOR_ALL_ELEMENTS_IN_MATRIX2D(ctf_ampl)
      ctf_ampl2(i,j)=ctf_ampl(i,j)*ctf_ampl(i,j);
 
   // If images are to be saved, save the amplitude (log)  
   if (prm.show_optimization) {
      save().resize(ctf_ampl);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(save())
         save(i,j)=log10(ctf_ampl2(i,j));
      save.write(prm.fn_ctf.get_root()+"_log10.xmp");
   }
   
   // Prepare global variables
   global_ctf_ampl=&ctf_ampl;
   global_ctf_ampl2=&ctf_ampl2;   
   global_ctf=prm.ctftomodel.mask;
   global_ctfmodel=prm.ctfmodel;
   global_initial_ctfmodel=&prm.ctfmodel;
   global_Tm=prm.Tm;
   global_adjust=&prm.adjust;
   global_penalize=TRUE;
   global_min_freq=prm.min_freq;
   global_max_freq=prm.max_freq;
   global_penalty=prm.penalty;
   global_current_penalty=0;
   global_central_weight=prm.central_weight;
   global_current_central_weight=1;
   global_evaluation_reduction=prm.evaluation_reduction;
   global_heavy_penalization=
      ctf_ampl2(0,0)*ctf_ampl2.RowNo()*ctf_ampl2.ColNo();
   histogram1D hist;
   global_hist=&hist;
   CTF_stats();
   global_value_th=global_min_val+(global_max_val-global_min_val)*
      prm.value_th/100;

   /************************************************************************
     STEP 1:  Find background which best fits the CTF
   /************************************************************************/
   // Some initialization
   double fitness=0;
   int    iter=0;
   global_action=0;
   global_gamma=1;
   
   // If initial parameters weren´t supplied for the gaussian curve,
   // estimate them from the CTF file
   if (prm.adjust(20)==0) {
       cerr << "Computing first sqrt background ...\n";
       estimate_background_sqrt_parameters();
       assign_parameters_from_CTF(global_ctfmodel,
          VEC_ARRAY(prm.adjust),0,CTF_PARAMETERS,
          global_Adjust_CTF_parameters->astigmatic_noise);
       if (prm.show_optimization) save_intermidiate_results("first_sqrt_fit");
   
      // Now optimize 
      if (!prm.do_not_optimize) {
         // Use every parameter to first gausssian optimization
         matrix1D<double> steps(BACKGROUND_CTF_PARAMETERS);
         steps.init_constant(0);
         steps(1)=steps(5)=steps(7)=1;
         if (prm.astigmatic_noise) {steps(4)=steps(6)=1;}
         cerr << "Looking for best fitting sqrt ...\n";
         global_penalize=FALSE;
         Powell_optimizer(prm.adjust, 14, BACKGROUND_CTF_PARAMETERS,
            &CTF_fitness, 2*prm.accuracy, fitness, iter, steps,
            prm.show_optimization);
         if (prm.show_optimization)
            save_intermidiate_results("best_sqrt_fit");
         global_penalize=TRUE;
         cerr << "   Penalizing best fitting sqrt ...\n";
         global_current_penalty=2;
         int imax=CEIL(log(global_penalty)/log(2.0));
         for (int i=1; i<=imax; i++) {
            cerr << "     Iteration " << i
	         << " penalty=" << global_current_penalty << endl;
	    Powell_optimizer(prm.adjust, 14, BACKGROUND_CTF_PARAMETERS, &CTF_fitness,
               2*prm.accuracy, fitness, iter, steps, prm.show_optimization);
	    global_current_penalty*=2;
	    global_current_penalty=MIN(global_current_penalty, global_penalty);
         }
         if (prm.show_optimization)
            save_intermidiate_results("best_penalized_sqrt_fit");
      }
   }

   // If initial parameters weren´t supplied for the gaussian curve,
   // estimate them from the CTF file
   if (prm.adjust(13)==0) {
       cerr << "Computing first background parameters ...\n";
       estimate_background_gauss_parameters();
       assign_parameters_from_CTF(global_ctfmodel,
          VEC_ARRAY(prm.adjust),13,BACKGROUND_CTF_PARAMETERS,
          global_Adjust_CTF_parameters->astigmatic_noise);
       if (prm.show_optimization)
          save_intermidiate_results("first_background_fit");
   
       // Now optimize 
       if (!prm.do_not_optimize) {
          // Use every parameter to first gausssian optimization
          matrix1D<double> steps(BACKGROUND_CTF_PARAMETERS);
          steps.init_constant(1);
          if (!prm.astigmatic_noise) {
             steps(3)=steps(4)=steps(6)=steps(9)=0;
          }
          global_penalize=TRUE;
          cerr << "Looking best fitting background ...\n";
          global_current_penalty=2;
          int imax=CEIL(log(global_penalty)/log(2.0));
          for (int i=1; i<imax; i++) {
             cerr << "     Iteration " << i
 	         << " penalty=" << global_current_penalty << endl;
             Powell_optimizer(prm.adjust, 14, BACKGROUND_CTF_PARAMETERS,
                &CTF_fitness, 2*prm.accuracy, fitness, iter, steps,
                prm.show_optimization);
 	     global_current_penalty*=2;
 	     global_current_penalty=MIN(global_current_penalty, global_penalty);
          }
       }
      if (prm.show_optimization)
         save_intermidiate_results("best_penalized_background_fit");
   }

   /************************************************************************
     STEP 2:  Adjust best parametric CTF with only a few parameters 
              (defocus, basically)
    ************************************************************************/
   // Some initialization
   fitness=0;
   iter=0;
   global_penalize=TRUE;
   global_ctfmodel.enable_CTF=TRUE;
   global_gamma=prm.gamma;
   global_compute_FFT_distance=TRUE;
   
   if (prm.adjust(0)==0) {
      // Find the best initial defocus and K estimation
      if (!prm.do_not_optimize) {
         cerr << "Looking for best fitting parametric CTF ...\n";
         // First perform an exhaustive search
         global_action=3;
         estimate_best_parametric_CTF(&CTF_fitness);
         assign_CTF_from_parameters(VEC_ARRAY(prm.adjust),
            global_ctfmodel,0,CTF_PARAMETERS,
            global_Adjust_CTF_parameters->astigmatic_noise);
         if (prm.show_optimization)
            save_intermidiate_results("first_parametric_fit");

         // Now optimize
         cerr << "Improving parametric fit ...\n";
         global_action=1;
         matrix1D<double> steps(PARAMETRIC_CTF_PARAMETERS);
         for (int i=0; i<PARAMETRIC_CTF_PARAMETERS; i++) 
             steps(i)=prm.steps(i);
         Powell_optimizer(prm.adjust, 1, PARAMETRIC_CTF_PARAMETERS,
            &CTF_fitness, 3*prm.accuracy, fitness, iter, steps,
            prm.show_optimization);
      }
      if (prm.show_optimization)
         save_intermidiate_results("best_parametric_fit");
   }
  
   /************************************************************************
     STEP 3:  Adjust best total CTF with all parameters
    ************************************************************************/
   fitness=0;
   iter=0;
   global_action=2;
//   global_gamma=1;
   
   // Now optimize
   if (!prm.do_not_optimize) {
       matrix1D<double> steps(CTF_PARAMETERS);
       steps.init_constant(1);
       for (int i=0; i<PARAMETRIC_CTF_PARAMETERS; i++) 
           steps(i)=prm.steps(i);
       if (!prm.astigmatic_noise) {
          steps(PARAMETRIC_CTF_PARAMETERS+3)=steps(PARAMETRIC_CTF_PARAMETERS+4)=
             steps(PARAMETRIC_CTF_PARAMETERS+6)=steps(PARAMETRIC_CTF_PARAMETERS+9)=0;
       }
       cerr << "Looking for best penalized fitting CTF ...\n";
       Powell_optimizer(prm.adjust, 1, CTF_PARAMETERS, &CTF_fitness,
          2*prm.accuracy, fitness, iter, steps,
          prm.show_optimization);
   }

   // Save total CTF if images are required
   if (prm.show_optimization)
      save_intermidiate_results("best_penalized_fit");
  
   /************************************************************************
     STEP 4:  Give more weight to central part
    ************************************************************************/
//   global_compute_FFT_distance=FALSE;
   if (!prm.do_not_optimize) {
       matrix1D<double> steps(CTF_PARAMETERS);
       steps.init_constant(1);
       for (int i=0; i<PARAMETRIC_CTF_PARAMETERS; i++) 
           steps(i)=prm.steps(i);
       cerr << "Looking for best central fitting CTF ...\n";
       if (!prm.astigmatic_noise) {
          steps(PARAMETRIC_CTF_PARAMETERS+3)=steps(PARAMETRIC_CTF_PARAMETERS+4)=
             steps(PARAMETRIC_CTF_PARAMETERS+6)=steps(PARAMETRIC_CTF_PARAMETERS+9)=0;
       }
       global_current_central_weight=2;
       
       // Adjust without penalization
       double imax=CEIL(log(prm.central_weight)/log(2.0));
       for (int i=1; i<=imax; i+=2) {
          cerr << "   Central weight= "
               << global_current_central_weight << endl;
          Powell_optimizer(prm.adjust, 1, CTF_PARAMETERS, &CTF_fitness,
                  2*prm.accuracy, fitness, iter, steps,
                  prm.show_optimization);
          #ifdef DEBUG
          if (prm.show_optimization) save_intermidiate_results(
             (string)"best_fit_weight_"+ItoA(i));
          #endif
          global_current_central_weight*=2;
          global_current_central_weight=MIN(prm.central_weight,
             global_current_central_weight);
       }
       if (prm.show_optimization)
          save_intermidiate_results("best_fit_weight");
   }

   /************************************************************************
     STEP 5:  Remove penalization
    ************************************************************************/
   if (!prm.do_not_optimize) {
       matrix1D<double> steps(CTF_PARAMETERS);
       steps.init_constant(1);
       for (int i=0; i<PARAMETRIC_CTF_PARAMETERS; i++) 
           steps(i)=prm.steps(i);
       if (!prm.astigmatic_noise) {
          steps(PARAMETRIC_CTF_PARAMETERS+3)=steps(PARAMETRIC_CTF_PARAMETERS+4)=
             steps(PARAMETRIC_CTF_PARAMETERS+6)=steps(PARAMETRIC_CTF_PARAMETERS+9)=0;
       }
       cerr << "Looking for best not penalized fitting CTF ...\n";
       
       // Adjust with less penalization
       double p=log(global_penalty)/log(2.0);
       int imax;
       if ((int)p==p) imax=(int)p-1;
       else imax=FLOOR(p);
       for (int i=imax; i>=0; i-=2) {
          global_current_penalty=pow(2.0,i);
          if (global_current_penalty==1) global_penalize=FALSE;
          cerr << "     Iteration " << i
 	      << " penalty=" << global_current_penalty << endl;
          Powell_optimizer(prm.adjust, 1, CTF_PARAMETERS, &CTF_fitness,
             2*prm.accuracy, fitness, iter, steps, prm.show_optimization);
          #ifdef DEBUG
          if (prm.show_optimization) save_intermidiate_results(
             (string)"best_fit_less_"+ItoA(imax-i));
          #endif
       }
       if (prm.show_optimization)
          save_intermidiate_results("best_fit_less");
       global_current_penalty=1;
    }

   /************************************************************************
     STEP 6:  Give more weight to central part
    ************************************************************************/
   if (!prm.do_not_optimize) {
       matrix1D<double> steps(CTF_PARAMETERS);
       steps.init_constant(1);
       for (int i=0; i<PARAMETRIC_CTF_PARAMETERS; i++) 
           steps(i)=prm.steps(i);
       if (!prm.astigmatic_noise) {
          steps(PARAMETRIC_CTF_PARAMETERS+3)=steps(PARAMETRIC_CTF_PARAMETERS+4)=
             steps(PARAMETRIC_CTF_PARAMETERS+6)=steps(PARAMETRIC_CTF_PARAMETERS+9)=0;
       }
       cerr << "Looking again for best central fitting CTF ...\n";
       global_current_central_weight=2;
       
       // Adjust without penalization
       double imax=CEIL(log(prm.central_weight)/log(2.0));
       for (int i=1; i<=imax; i+=2) {
          cerr << "   Central weight= "
               << global_current_central_weight << endl;
          Powell_optimizer(prm.adjust, 1, CTF_PARAMETERS, &CTF_fitness,
                  2*prm.accuracy, fitness, iter, steps,
                  prm.show_optimization);
          #ifdef DEBUG
          if (prm.show_optimization) save_intermidiate_results(
             (string)"best_fit_weight_again_"+ItoA(i));
          #endif
          global_current_central_weight*=2;
          global_current_central_weight=MIN(prm.central_weight,
             global_current_central_weight);
       }
       if (prm.show_optimization)
          save_intermidiate_results("best_fit_weight_again");
   }

   /************************************************************************
     STEP 7:  Fine search
    ************************************************************************/
   if (!prm.do_not_optimize && !prm.do_not_fine_search) {
       matrix1D<double> steps(CTF_PARAMETERS);
       steps.init_constant(1);
       for (int i=0; i<PARAMETRIC_CTF_PARAMETERS; i++) 
           steps(i)=prm.steps(i);
       if (!prm.astigmatic_noise) {
          steps(PARAMETRIC_CTF_PARAMETERS+3)=steps(PARAMETRIC_CTF_PARAMETERS+4)=
             steps(PARAMETRIC_CTF_PARAMETERS+6)=steps(PARAMETRIC_CTF_PARAMETERS+9)=0;
       }
       cerr << "Looking for best fitting CTF ...\n";
       global_evaluation_reduction=1;
       
       // Adjust without penalization
       Powell_optimizer(prm.adjust, 1, CTF_PARAMETERS, &CTF_fitness,
          prm.accuracy, fitness, iter, steps,
          prm.show_optimization);
   }

   // Save total CTF if images are required
   if (prm.show_optimization) save_intermidiate_results("best_fit");

   /************************************************************************
     STEP 8:  Write final CTF file
    ************************************************************************/
   FourierImageXmipp Final_CTF;
   generate_model_so_far(save);
   Final_CTF().resize(save());
   FOR_ALL_ELEMENTS_IN_MATRIX2D(Final_CTF()) Final_CTF(i,j)=save(i,j);
   Final_CTF.write(prm.fn_outroot+"_model.fft");

   ofstream output_file((prm.fn_out_CTF_parameters).c_str());
   global_ctfmodel.Produce_Main_Info();
   output_file << global_ctfmodel << endl;
   output_file.close();

   save_intermidiate_results(prm.fn_outroot+"_model");
}
#endif
