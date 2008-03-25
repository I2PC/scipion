/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.uam.es (2004)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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

#include "new_projmatch.h"

// Read arguments ==========================================================
void Prog_new_projection_matching_prm::read(int argc, char **argv)  {

  // Read command line
  if (checkParameter(argc,argv,"-more_options")) { usage(); extendedUsage();}
  fn_ref = getParameter(argc,argv,"-ref","ref");
  fn_exp=getParameter(argc,argv,"-i");
  fn_root=getParameter(argc,argv,"-o","out");
  max_shift=textToFloat(getParameter(argc,argv,"-max_shift","5"));

  // Additional commands
  Ri=textToInteger(getParameter(argc,argv,"-Ri","-1"));
  Ro=textToInteger(getParameter(argc,argv,"-Ro","-1"));

  // Hidden stuff
  verb=textToInteger(getParameter(argc,argv,"-verb","1"));

}

// Show ====================================================================
void Prog_new_projection_matching_prm::show() {

  if (verb>0) {
    std::cerr << "  Input images            : "<< fn_exp << std::endl;
    std::cerr << "  Output rootname         : "<< fn_root << std::endl;
    if (Ri>0)
    std::cerr << "  Inner radius rot-search : "<<Ri<<std::endl;
    if (Ro>0)
    std::cerr << "  Outer radius rot-search : "<<Ro<<std::endl;
    std::cerr << "  -> Limit search of origin offsets to  +/- "<<max_shift<<" pixels"<<std::endl;
    std::cerr << " ================================================================="<<std::endl;
  }
}

// Usage ===================================================================
void Prog_new_projection_matching_prm::usage() {
  std::cerr << "Usage:  projection_matching [options] "<<std::endl;
  std::cerr << "   -i <selfile>                : Selfile with input images \n"
       << "   -vol <volume>               : Reference volume \n"
       << " [ -o <rootname=\"out\"> ]       : Output rootname \n"
       << " [ -max_shift <float=5> ]      : Maximum change in origin offset (+/- pixels) \n"
       << " [ -more_options ]             : Show all program options\n";
}

// Extended usage ===================================================================
void Prog_new_projection_matching_prm::extendedUsage() {
  std::cerr << "Additional options: \n"
       << " [ -Ri <float=1> ]             : Inner radius to limit rotational search \n"
	    << " [ -Ro <float=dim/2 - 1> ]     : Outer radius to limit rotational search \n";
  exit(1);
}

// Side info stuff ===================================================================
void Prog_new_projection_matching_prm::produceSideInfo() {

    ImageXmipp       img,empty;
    Projection       proj;
    DocFile          DF;
    SelFile          SFr,emptySF;
    SymList          SL;
    FileName         fn_tmp, fn_refs;
    double           mean,stddev,psi=0.;
    Matrix1D<double> dataline(3);
    int              nl;

    // Read Selfile and get dimensions
    DFexp.read(fn_exp);

    // Set ring defaults
    if (Ri<1) Ri=1;
    if (Ro<0) Ro=(dim/2)-1;

    // Set up angular sampling
    mysampling.read_sampling_file(fn_ref);

}

void Prog_new_projection_matching_prm::rotationallyAlignOneImage(const Matrix2D<double> &img,
								 int &opt_samplenr,
								 double &opt_psi, 
								 double &opt_flip, 
								 double &maxcorr)
{
    Matrix2D<double>         Maux;
    Polar<double>            P;
    Polar<std::complex <double> > fP,fPm;
    Matrix1D<double>         ang,corr;
    int                      max_index;
    double                   mean,stddev_img;

    // Calculate polar coordinates using gridding
    produceReverseGriddingMatrix2D(img,Maux,kb);
    P.getPolarFromCartesianGridding(Maux,kb,Ri,Ro);
    mean = P.computeSum(true);
    stddev_img = P.computeSum2(true);
    stddev_img = sqrt(stddev_img - mean * mean);
    P -= mean; // for normalized cross-correlation coefficient
    fP = P.fourierTransformRings(false);
    fPm = P.fourierTransformRings(true);

    // Loop over all relevant "neighbors" (i.e. directions within the search range)
    maxcorr = -99.e99;

    for (int i = 0; i < id_ref.size(); i++)
    {
	// A. Check straight image
	rotationalCorrelation(fP,fP_ref[i],ang,corr);
	corr /= stddev_ref[i] * stddev_img; // for normalized ccf
	for (int k = 0; k < XSIZE(corr); k++)
	{
	    if (corr(k)> maxcorr)
	    {
		maxcorr = corr(k);
		opt_psi = ang(k);
		opt_samplenr = i;
		opt_flip = 0.;
	    }
	}

	// B. Check mirrored image
	rotationalCorrelation(fPm,fP_ref[i],ang,corr);
	corr /= stddev_ref[i] * stddev_img; // for normalized ccf
	for (int k = 0; k < XSIZE(corr); k++)
	{
	    if (corr(k)> maxcorr)
	    {
		maxcorr = corr(k);
		opt_psi = ang(k);
		opt_samplenr = i;
		opt_flip = 1.;
	    }
	}
    }

}

void Prog_new_projection_matching_prm::translationallyAlignOneImage(const Matrix2D<double> &img,
								    const int &samplenr, const double &psi, 
								    const double &opt_flip, double &opt_xoff, 
								    double &opt_yoff, double &maxcorr)
{
    Projection       proj;
    Matrix2D<double> Mtrans,Mimg,Mref;

    // Rotate stored reference projection by phi degrees
    proj_ref[samplenr].rotateBSpline(3,-psi,Mref,DONT_WRAP);

    if (opt_flip > 0.)
    {
	// Flip experimental image
	Matrix2D<double> A(3,3);
	A.initIdentity();
	A(0, 0) *= -1.;
	A(0, 1) *= -1.;
	applyGeometry(Mimg,A, img, IS_INV, DONT_WRAP);
    }
    else
	Mimg = img;

    // Perform the actual search for the optimal shift
    if (max_shift>0) 
	best_shift(Mref,Mimg,opt_xoff,opt_yoff);
    else 
	opt_xoff = opt_yoff = 0.;
    if (opt_xoff * opt_xoff + opt_yoff * opt_yoff > max_shift * max_shift) 
	opt_xoff = opt_yoff = 0.;

    // Calculate standard cross-correlation coefficient
    Mimg.translate(vectorR2(opt_xoff,opt_yoff),Mtrans,true);
    maxcorr = correlation_index(Mref,Mtrans);

    // Correct X-shift for mirrored images
    if (opt_flip>0.)
	opt_xoff *= -1.;	

}

void Prog_new_projection_matching_prm::getCurrentImage(int imgno, ImageXmipp img)
{

    FileName         fn_img;
    DocLine          DL;
    Matrix2D<double> A;

    // jump to line imgno in DFexp, get data and filename
    DFexp.locate(imgno);
    DL = DFexp.get_current_line();
    DFexp.previous();
    if (DFexp.get_current_line().Is_comment()) 
    {
	fn_img = (DFexp.get_current_line()).get_text();
    }
    else
    {
	REPORT_ERROR(1,"BUG: no comment in DFexp where expected....");
    }

    // Read actual image
    img.read(fn_img);
    img().setXmippOrigin();

    // Store translation in header and apply it to the actual image
    img.Xoff() = DL[3];
    img.Yoff() = DL[4];
    img.rot()  = 0.;
    img.tilt() = 0.;
    img.psi()  = 0.;
    img.flip() = 0.;

    A = img.get_transformation_matrix(true);
    if (!A.isIdentity())
	img().selfApplyGeometryBSpline(A, 3, IS_INV, WRAP);

}

void Prog_new_projection_matching_prm::getCurrentReferences(int imgno)
{

    FileName                      fnt;
    ImageXmipp                    img;
    int                           refno;
    double                        mean,stddev;
    Matrix2D<double>              Maux;
    Polar<double>                 P;
    Polar<std::complex <double> > fP;

    //TODO: KEEP REFS THAT I STILL WANT TO USE INSTEAD OF READING ALL
    //AGAIN FROM DISC!!!
    fP_ref.clear();
    stddev_ref.clear();
    proj_ref.clear();

    // Loop over all relevant references
    for (int i = 0; i < mysampling.my_neighbors[imgno].size(); i++)
    {
	// all my neighboyurs are:
	refno = mysampling.my_neighbors[imgno][i];
	fnt.compose(fn_ref,refno,"xmp");
	img.read(fnt);
	  
	// Calculate FTs of polar rings (store its complex conjugated!)
	produceReverseGriddingMatrix2D(img(),Maux,kb);
	P.getPolarFromCartesianGridding(Maux,kb,Ri,Ro);
	mean = P.computeSum(true);
	stddev = P.computeSum2(true);
	stddev = sqrt(stddev - mean * mean);
	P -= mean;
	fP = P.fourierTransformRings(true);
	fP_ref.push_back(fP);
	stddev_ref.push_back(stddev);
	// Also store reference projections in cartesian coordinates 
	proj_ref.push_back(img());
	id_ref.push_back(refno);
    }

}



void Prog_new_projection_matching_prm::processSomeImages(int * my_images, double * my_output) 
{


  ImageXmipp img;
  double opt_rot, opt_tilt, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr;
  int opt_samplenr;

  int nr_images = my_images[0];
  // Prepare my_output
  my_output[0] = nr_images * MY_OUPUT_SIZE;

  // Loop over all images
  for (int imgno = 0; imgno < nr_images; imgno++)
  {

      int this_image = my_images[imgno];

      // read experimental image and all corresponding references in memory
      getCurrentImage(this_image,img);

      // Read all references and set vectors with FTs of their polar coordinates
      getCurrentReferences(this_image);
      
      // Align the image (for now 3D+2D search only)
      rotationallyAlignOneImage(img(), opt_samplenr, opt_psi, opt_flip, maxcorr);
      opt_rot  = XX(mysampling.no_redundant_sampling_points_angles[id_ref[opt_samplenr]]);
      opt_tilt = YY(mysampling.no_redundant_sampling_points_angles[id_ref[opt_samplenr]]);
      opt_psi += ZZ(mysampling.no_redundant_sampling_points_angles[id_ref[opt_samplenr]]);
      translationallyAlignOneImage(img(), opt_samplenr, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr);

      // Add previously applied translation to the newly found one
      opt_xoff += img.Xoff();
      opt_yoff += img.Yoff();

      // Output
      my_output[imgno * MY_OUPUT_SIZE + 1] = imgno;
      my_output[imgno * MY_OUPUT_SIZE + 2] = opt_rot;
      my_output[imgno * MY_OUPUT_SIZE + 3] = opt_tilt;
      my_output[imgno * MY_OUPUT_SIZE + 4] = opt_psi;
      my_output[imgno * MY_OUPUT_SIZE + 5] = opt_xoff;
      my_output[imgno * MY_OUPUT_SIZE + 6] = opt_yoff;
      my_output[imgno * MY_OUPUT_SIZE + 7] = id_ref[opt_samplenr];
      my_output[imgno * MY_OUPUT_SIZE + 8] = opt_flip;
      my_output[imgno * MY_OUPUT_SIZE + 9] = maxcorr;

  }

}

