/***************************************************************************
 *
 * Authors: Carlos Oscar (coss@cnb.uam.es)
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

#include "../xmippMicrograph.hh"
#include "../xmippArgs.hh"
#include "../xmippSelFiles.hh"
#include "../xmippMasks.hh"
#include "../xmippGeometry.hh"
#include <fstream>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#ifdef LINUX
   #include <unistd.h>
#endif


/* Clear ------------------------------------------------------------------- */
void Micrograph::clear() {
   single_particle.clear();
   coords.clear();
   fn_coords=fn_micrograph="";
   X_window_size=Y_window_size=-1;
   fh_micrograph=-1;
   Xdim=Ydim=-1;
   __depth=-1;
   m8=NULL;
   m16=NULL;
   um16=NULL;
   m32=NULL;
   compute_transmitance=false;
   compute_inverse=false;
   __scaling_valid=false;
   /* __in_core=FALSE;*/
}

/* Open micrograph --------------------------------------------------------- */
void Micrograph::open_micrograph(const FileName &_fn_micrograph,
   /*bool in_core,*/ bool reversed) {
   struct stat info;

   // Micrograph name
   fn_micrograph = _fn_micrograph;

   // Look for micrograph dimensions
   // check if the file format is spider
   if(Is_ImageXmipp(fn_micrograph)){
      headerXmipp     header;
      header.read(fn_micrograph);
      float fXdim,fYdim;
      header.get_dimension(fYdim,fXdim);
      Xdim = (int) fXdim; Ydim= (int)fYdim;
      __offset=header.get_header_size();
      __depth=32;
      reversed=header.reversed();
   }
   else {
      fn_inf=fn_micrograph.add_extension("inf");
      FILE *fh_inf=fopen(fn_inf.c_str(),"r");
      if (!fh_inf)
	 REPORT_ERROR(1,(string)"Micrograph::open_micrograph: Cannot find "+
            fn_inf);
      Xdim=AtoI(get_param(fh_inf,"Xdim"));
      Ydim=AtoI(get_param(fh_inf,"Ydim"));
      __depth=AtoI(get_param(fh_inf,"bitspersample"));
      if(check_param(fh_inf,"offset"))
	 __offset=AtoI(get_param(fh_inf,"offset"));
      else
	 __offset=0;   
      if (check_param(fh_inf,"is_signed"))
	 __is_signed=(get_param(fh_inf,"is_signed")=="true" ||
                      get_param(fh_inf,"is_signed")=="TRUE");
      else __is_signed=false;
      fclose(fh_inf);
   }
   // Open micrograph and map
   fh_micrograph=open(fn_micrograph.c_str(),O_RDWR,S_IREAD|S_IWRITE);
   if (fh_micrograph==-1)
      REPORT_ERROR(1,(string)"Micrograph::open_micrograph: There is a "
         "problem opening " +fn_micrograph);
   char *aux_ptr;
   switch (__depth) {
      case 8:
         /* if (!in_core) { */
            m8=(unsigned char *) mmap(0,(__depth/8)*Ydim*Xdim+__offset,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
            if (m8==MAP_FAILED)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
	     m8+=__offset;
         /* } else {
            m8=new unsigned char (Ydim*Xdim*__depth/8);
            int length=Ydim*Xdim*__depth/8;
            int read_length=read(fh_micrograph,m8,length);
            cout << SSIZE_MAX << endl;
            cout << read_length << endl;
            if (read_length!=length)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot read "+
                  _fn_micrograph+" in memory");
         } */
         break;
      case 16:
         /* if (!in_core) { */
	 if(__is_signed){
            m16=(short int *) mmap(0,(__depth/8)*Ydim*Xdim+__offset,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
            if (m16==MAP_FAILED) {
               /*
	       switch (errno) {
	          case EACCES:    cout << "EACCES:   \n"; break;
      	          case EAGAIN:    cout << "EAGAIN:   \n"; break;
		  case EBADF:     cout << "EBADF:    \n"; break;
		  case EINVAL:    cout << "EINVAL:   \n"; break;
		  case EMFILE:    cout << "EMFILE:   \n"; break;
		  case ENODEV:    cout << "ENODEV:   \n"; break;
		  case ENOMEM:    cout << "ENOMEM:   \n"; break;
		  case ENOTSUP:   cout << "ENOTSUP:  \n"; break;
		  case ENXIO:     cout << "ENXIO:    \n"; break;
		  case EOVERFLOW: cout << "EOVERFLOW:\n"; break;
	       }
	       */
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
	    }
          aux_ptr=(char *)m16;

	  aux_ptr+=__offset;
	  m16=(short int *) aux_ptr;
	  }
	  else {
            um16=(unsigned short int *) mmap(0,(__depth/8)*Ydim*Xdim+__offset,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
            if (um16==MAP_FAILED) {
               /*
	       switch (errno) {
	          case EACCES:    cout << "EACCES:   \n"; break;
      	          case EAGAIN:    cout << "EAGAIN:   \n"; break;
		  case EBADF:     cout << "EBADF:    \n"; break;
		  case EINVAL:    cout << "EINVAL:   \n"; break;
		  case EMFILE:    cout << "EMFILE:   \n"; break;
		  case ENODEV:    cout << "ENODEV:   \n"; break;
		  case ENOMEM:    cout << "ENOMEM:   \n"; break;
		  case ENOTSUP:   cout << "ENOTSUP:  \n"; break;
		  case ENXIO:     cout << "ENXIO:    \n"; break;
		  case EOVERFLOW: cout << "EOVERFLOW:\n"; break;
	       }
	       */
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
	    }
          aux_ptr=(char *)um16;
	  aux_ptr+=__offset;
	  um16=(unsigned short int *) aux_ptr;
	  } 
         /* } else {
            m16=new short int (Ydim*Xdim*__depth/8);
            int length=Ydim*Xdim*__depth/8;
            if (read(fh_micrograph,m16,length)!=length)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot read "+
                  _fn_micrograph+" in memory");
         } */
         break;
      case 32:
      	    // Map file in memory
            m32=(float*) mmap(0,(__depth/8)*Ydim*Xdim+__offset,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph,0);
            if (m32==MAP_FAILED)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
            aux_ptr=(char *)m32;
	    aux_ptr+=__offset;
	    m32=(float *) aux_ptr;
	    break;
      default:
         REPORT_ERROR(1,"Micrograph::open_micrograph: depth is not 8, 16 nor 32");
   }
   /*__in_core=in_core; */
   __reversed=reversed;
}

/* Close micrograph -------------------------------------------------------- */
void Micrograph::close_micrograph() {
   if (fh_micrograph!=-1) {
      close(fh_micrograph);
      /* if (!__in_core) { */
         if      (__depth== 8) {
	    m8-=__offset;
	    munmap((char *)m8,Ydim*Xdim*__depth/8+__offset);
         } else if (__depth==16) {
	    if(__is_signed){
		char *aux_ptr=(char *)m16;
		aux_ptr-=__offset;
		m16=(short int *)aux_ptr;
		munmap((char *)m16,Ydim*Xdim*__depth/8+__offset);
            }
	    else{
	        char *aux_ptr=(char *)um16;
	        aux_ptr-=__offset;
	        um16=(unsigned short int *)aux_ptr;
	        munmap((char *)um16,Ydim*Xdim*__depth/8+__offset);
	    }	
	 } else if (__depth==32) {
	    char *aux_ptr=(char *)m32;
	    aux_ptr-=__offset;
	    m32=(float *)aux_ptr;
	    munmap((char *)m32,Ydim*Xdim*__depth/8+__offset);
	 }
      /* } else {
         if      (__depth== 8) delete m8;
         else if (__depth==16) delete m16;
      } */
   }
}

/* Compute 8 bit scaling --------------------------------------------------- */
void Micrograph::compute_8_bit_scaling() {
   cerr << "Computing 8 bit scaling ...\n";

   // Compute minimum and maximum value
   float minval, maxval;
   minval=maxval=(*this)(0,0);
   for (int i=0; i<Ydim; i++) {
      for (int j=0; j<Xdim; j++) {
         float tmp=(*this)(j,i);
         if      (tmp<minval) minval=tmp;
         else if (tmp>maxval) maxval=tmp;
/*
if(maxval > 32000)
  cout << "(i,j) max min valuefloat value" << i << " " << j 
     << " " << maxval << " " << minval << " "<< tmp 
     << " " << (*this)(j,i) << endl;
*/
      }
   }

   // Compute output range
   float minF, maxF;
   if      (minval<0)         {minF=0; maxF=MIN(255,maxval+minval);}
   else if (maxval>255)       {minF=MAX(0,minval-(maxval-255)); maxF=255;}
   else if (maxval-minval<32) {minF=0; maxF=255;}
   else                       {minF=minval; maxF=maxval;}

   // Compute scaling
   __a=(maxF-minF)/(maxval-minval);
   __b=minF-__a*minval;
   __scaling_valid=true;
}

/* Save coordinates to disk ------------------------------------------------ */
void Micrograph::write_coordinates(int label, const FileName &_fn_coords) {
   ofstream fh;
   if (_fn_coords!="") fn_coords=_fn_coords;
   fh.open(fn_coords.c_str(), ios::out);
   if (!fh)
      REPORT_ERROR(1,(string)"Micrograph::write: File "+fn_coords+
         " cannot be openned for output");
   int imax=coords.size();
   fh << "# <X position> <Y position>\n";
   for (int i=0; i<imax; i++)
      if (coords[i].valid && coords[i].label==label)
         fh << coords[i].X << " " << coords[i].Y << endl;
   fh.close();
}

/* Read coordinates from disk ---------------------------------------------- */
void Micrograph::read_coordinates(int label, const FileName &_fn_coords) {
   ifstream  fh;
   int       line_no=0;
   string    line;
   
   fn_coords=_fn_coords;
   fh.open(fn_coords.c_str(), ios::in);
   if (!fh)
      REPORT_ERROR(1,(string)"Micrograph::read: File "+fn_coords+" not found");

   // Count the number of lines
   fh.peek();
   while (!fh.eof()) {
      getline(fh,line);
      if (line.length()>0 && line[0]!='#' && line[0]!=';') line_no++;
      fh.peek();
   }
   fh.close();
   fh.clear();

   // Resize coordinate list and read
   fh.open(fn_coords.c_str(), ios::in);
   coords.reserve(line_no);
   struct Particle_coords aux;
   aux.valid=true;
   aux.label=label;
   fh.peek();
   while (!fh.eof()) {
      getline(fh,line);
      if (line.length()>0 && line[0]!='#' && line[0]!=';') {
      	 int converted_elements=sscanf(line.c_str(),"%d %d",
	    &aux.X, &aux.Y);
	 if (converted_elements!=2)
	      cerr << "Ignoring line: " << line << endl;
         else coords.push_back(aux);
      }
      fh.peek();
   }
   fh.close();
}

/* Scissor ----------------------------------------------------------------- */
int Micrograph::scissor(const Particle_coords &P, Image &result,
   double Dmin, double Dmax, double scaleX, double scaleY, 
   bool only_check) {
   if (X_window_size==-1 || Y_window_size==-1)
      REPORT_ERROR(1,"Micrograph::scissor: window size not set");

   result().resize(Y_window_size, X_window_size);
   int i0=ROUND(scaleY*P.Y)+FIRST_XMIPP_INDEX(Y_window_size);
   int iF=ROUND(scaleY*P.Y)+LAST_XMIPP_INDEX(Y_window_size);
   int j0=ROUND(scaleX*P.X)+FIRST_XMIPP_INDEX(X_window_size);
   int jF=ROUND(scaleX*P.X)+LAST_XMIPP_INDEX(X_window_size);
   int retval=1;
   double range,temp;
   range=Dmax-Dmin;
   if (i0<0 || iF>=Ydim || j0<0 || jF>=Xdim) {
      result().init_zeros();
      retval=0;
   } else
      if (!only_check){
	 for (int i=i0; i<=iF; i++)
            for (int j=j0; j<=jF; j++){
	       if(compute_transmitance){
		  if((*this)(j,i)<1)
		     temp = (*this)(j,i);
		  else   
		     temp = log10((double)(*this)(j,i));
		  if(compute_inverse)
		     result(i-i0,j-j0)= (Dmax-temp)/range;
		  else  
		     result(i-i0,j-j0)= (temp-Dmin)/range;
		  }   
	        else{
		   if(compute_inverse)
	                result(i-i0,j-j0)= (Dmax-(*this)(j,i))/range;
	           else
	                result(i-i0,j-j0)= (*this)(j,i);
	        }
	     }
      }	     	  
   return retval;
}

/* Produce all images ------------------------------------------------------ */
void Micrograph::produce_all_images(int label, const FileName &fn_root,
   int starting_index, const FileName &fn_image, double ang, double tilt,
   double psi) {
   SelFile SF;
   FileName fn_out;
   ImageXmipp I;
   Micrograph *M;

   // Set Source image
   if (fn_image=="") M=this;
   else {
      M=new Micrograph;
      M->open_micrograph(fn_image,__reversed);
      M->set_window_size(X_window_size, Y_window_size);
      M->set_transmitance_flag(compute_transmitance);
      M->set_inverse_flag(compute_inverse);
   }
   
   // Set scale for particles
   int MXdim, MYdim, thisXdim, thisYdim;
   M->size(MXdim, MYdim);
   this->size(thisXdim, thisYdim);
   double scaleX=(double)MXdim/thisXdim;
   double scaleY=(double)MYdim/thisYdim;
   
   // Compute max and minimum if compute_transmitance 
   // or compute_inverse flags are ON
   double Dmax, Dmin;
   if(compute_transmitance || compute_inverse){
      (*this).compute_double_minmax(Dmin,Dmax);
      //#define DEBUG66
      #ifdef DEBUG66
        cout << "Min= " << Dmin << " Dmax" << Dmax << endl;
      #endif
      #undef DEBUG66
      if(compute_transmitance)
         {
	 if( Dmin > 1) Dmin = log10(Dmin);
         if( Dmax > 1) Dmax = log10(Dmax);
	 }
   }
   // Scissor all particles
   if (ang!=0)
      cout << "Angle from Y axis to tilt axis " << ang << endl
	   << "   applying apropriate rotation\n";
   int i=starting_index;
   int nmax=ParticleNo();
   for (int n=0; n<nmax; n++) 
      if (coords[n].valid && coords[n].label==label) {
         fn_out.compose(fn_root,i++,"xmp");
         if (!M->scissor(coords[n],(Image &) I, Dmin, Dmax, scaleX, scaleY)) {
            cout << "Particle " << fn_out << " is very near the border, "
                 << "corresponding image is set to blank\n";
            SF.insert(fn_out,SelLine::DISCARDED);
         } else SF.insert(fn_out);
//	 if (ang!=0) I().rotate(-ang);
         I.rot()=(float)ang;
         I.tilt()=(float)tilt;
         I.psi()=(float)psi;
         I.write(fn_out);
      }
   if (labels[label]!="") {
      SF.write(fn_micrograph.remove_directories()+"."+labels[label]+".sel");
      write_coordinates(label,fn_micrograph+"."+labels[label]+".pos");
   } else {
      SF.write(fn_micrograph.remove_directories()+".sel");
      write_coordinates(label,fn_micrograph+".pos");
   }
   
   // Free source image??
   if (fn_image!="") {
      M->close_micrograph();
      delete M;
   }
}

/* Search coordinate near a position --------------------------------------- */
int Micrograph::search_coord_near(int x, int y, int prec) const {
   int imax=coords.size();
   int prec2=prec*prec;
   for (int i=0; i<imax; i++)
      if ((coords[i].X-x)*(coords[i].X-x)+(coords[i].Y-y)*(coords[i].Y-y)<prec2
          && coords[i].valid) return i;
   return -1;
}

/* Invalidate a coordinate ------------------------------------------------- */
void Micrograph::invalidate_coord(int n) {
   if (n<0 || n>=ParticleNo())
      REPORT_ERROR(1,"Micrograph::invalidate_coord: Index out of range");
   coords[n].valid=false;   
}

/* Add coordinate ---------------------------------------------------------- */
int Micrograph::add_coord(int x, int y, int label) {
   struct Particle_coords aux;
   aux.valid=true;
   aux.X=x;
   aux.Y=y;
   aux.label=label;
   coords.push_back(aux);
   return coords.size()-1;
}

/* Move last coordinate ---------------------------------------------------- */
void Micrograph::move_last_coord_to(int x, int y) {
   if (coords.size()>0) {
      coords.back().X=x;
      coords.back().Y=y;
   }
}

/* Downsample -------------------------------------------------------------- */
void downsample(const Micrograph &M, int Xstep, int Ystep,
   const matrix2D<double> &kernel, Micrograph &Mp) {
   // Find first and last indexes in each direction
   // it is assumed that (0,0) is within the kernel
   int Ydim, Xdim, Ypdim, Xpdim;
   M.size(Xdim,Ydim);
   Mp.size(Xpdim,Ypdim);
   int x0=0;
   int y0=0;
   int xF=Xdim;
   int yF=Ydim;
   
   double pixval;
   int ii, y, jj, x;
   time_config();

   // Look for input/output ranges
   double a=1;
   double b=0;
   double scale=1;

   if (Mp.depth()!=32) {
      double imin, imax;
      double omin, omax;
      bool ifirst=true, ofirst=true;

      if (M.depth()!=32)
         scale=(pow(2.0,Mp.depth())-1.0)/(pow(2.0,M.depth())-1.0);
      else if (M.depth()==32) scale=1;

      init_progress_bar(yF/Ystep);
      for (ii=0, y=y0; y<yF; y+=Ystep, ii++) {
	  for (jj=0, x=x0; x<xF; x+=Xstep, jj++) {
             pixval=0;
             FOR_ALL_ELEMENTS_IN_MATRIX2D(kernel) {
	        int j2=intWRAP(j+x,0,xF-1);
	        int i2=intWRAP(i+y,0,yF-1);
        	if (ifirst) {imin=imax=M(j2,i2); ifirst=false;}
        	else {imin=MIN(imin,M(j2,i2)); imax=MAX(imax,M(j2,i2));}
        	pixval += kernel(i,j)*M(j2,i2);
             }
             pixval *= scale;
             if (ii<Ypdim && jj<Xpdim) {
        	if (ofirst) {omin=omax=pixval; ofirst=false;}
        	else {omin=MIN(omin,pixval); omax=MAX(omax,pixval);}
             }
	  }
	  if (ii%50==0) progress_bar(ii);
      }
      progress_bar(yF/Ystep);

      // Compute range transformation
      double irange=imax-imin;
      double orange=omax-omin;

      if (M.depth()!=32) {
         a=scale*irange/orange;
         b=-omin;
      } else {
         a=(pow(2.0,Mp.depth())-1.0)/orange;
         scale=1;
         b=-omin;
      }
   }

   // Really downsample
   init_progress_bar(yF/Ystep);
   for (ii=0, y=y0; y<yF; y+=Ystep, ii++) {
       for (jj=0, x=x0; x<xF; x+=Xstep, jj++) {
          pixval=0;
          FOR_ALL_ELEMENTS_IN_MATRIX2D(kernel) {
	     int j2=intWRAP(j+x,0,xF-1);
	     int i2=intWRAP(i+y,0,yF-1);
             pixval += kernel(i,j)*M(j2,i2);
	  }
          if (ii<Ypdim && jj<Xpdim) 
	     if (Mp.depth()!=32) Mp.set_val(jj,ii,FLOOR(a*(pixval*scale+b)));
	     else Mp.set_val(jj,ii,pixval);
       }
       if (ii%50==0) progress_bar(ii);
   }
   progress_bar(yF/Ystep);
}

/* Normalizations ---------------------------------------------------------- */
void normalize_OldXmipp(Image *I) {
   double avg, stddev, min, max;
   (*I)().compute_stats(avg, stddev, min, max);
   (*I)() -= avg;
   (*I)() /= stddev;
}

void normalize_Near_OldXmipp(Image *I, const matrix2D<int> &bg_mask) {
   double avg, stddev, min, max;
   double avgbg, stddevbg, minbg, maxbg;
   (*I)().compute_stats(avg, stddev, min, max);
   compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
      stddevbg);
   (*I)() -= avg;
   (*I)() /= stddevbg;
}

void normalize_Oldmipp_decomposition(Image *I, const matrix2D<int> &bg_mask,
   const matrix2D<double> *mask) {
   double avgbg, stddevbg, minbg, maxbg;
   compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
      stddevbg);
   (*I)() -= avgbg;
   (*I)() /= stddevbg;
   if (mask!=NULL) (*I)() *= *mask;
   normalize_OldXmipp(I);
}

void normalize_Michael(Image *I, const matrix2D<int> &bg_mask) {
   double avg, stddev, min, max;
   double avgbg, stddevbg, minbg, maxbg;
   (*I)().compute_stats(avg, stddev, min, max);
   compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
      stddevbg);
   if (avgbg>0) {
      (*I)() -= avgbg;
      (*I)() /= avgbg;
   } else { // To avoid the contrast inversion
      (*I)() -= (avgbg-min);
      (*I)() /= (avgbg-min);
   }
}

void normalize_NewXmipp(Image *I, const matrix2D<int> &bg_mask) {
   double avgbg, stddevbg, minbg, maxbg;
   compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
      stddevbg);
   (*I)() -= avgbg;
   (*I)() /= stddevbg;
}

void normalize_NewXmipp2(Image *I, const matrix2D<int> &bg_mask) {
   double avg, stddev, min, max;
   double avgbg, stddevbg, minbg, maxbg;
   (*I)().compute_stats(avg, stddev, min, max);
   compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
      stddevbg);
   (*I)() -= avgbg;
   (*I)() /= ABS(avg-avgbg);
}

void normalize_ramp(Image *I, const matrix2D<int> &bg_mask) {
  fit_point          onepoint;
  vector<fit_point>  allpoints;
  double             pA,pB,pC;
  double             avgbg, stddevbg, minbg, maxbg;

  // Fit a least squares plane through the background pixels
  allpoints.clear();
  (*I)().set_Xmipp_origin();
  FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)()) {
    if (MAT_ELEM(bg_mask,i,j)) {
      onepoint.x=j;
      onepoint.y=i;
      onepoint.z=MAT_ELEM((*I)(),i,j);
      onepoint.w=1.;
      allpoints.push_back(onepoint);
    }    
  }
  least_squares_plane_fit(allpoints,pA,pB,pC);
  // Substract the plane from the image
  FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)()) {
    MAT_ELEM((*I)(),i,j)-=pA*j+pB*i+pC;
  }  
  // Divide by the remaining std.dev. in the background region
  compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
				   stddevbg);
  (*I)() /= stddevbg;

}

#ifdef NEVER_DEFINED
// This version doesn't work because of the high variance of avg-avg_bg
void normalize_NewXmipp(Image *I, const matrix2D<int> &bg_mask) {
   double avg, stddev, min, max;
   double avgbg, stddevbg, minbg, maxbg;
   (*I)().compute_stats(avg, stddev, min, max);
   compute_stats_within_binary_mask(bg_mask, (*I)(), minbg, maxbg, avgbg,
      stddevbg);
   (*I)() -= avgbg;
   (*I)() /= avg-avgbg;
}
#endif
