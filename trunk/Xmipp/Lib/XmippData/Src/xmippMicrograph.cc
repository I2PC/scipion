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
#include <fstream.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
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
   m32=NULL;
   __scaling_valid=FALSE;
   /* __in_core=FALSE;*/
}

/* Open micrograph --------------------------------------------------------- */
void Micrograph::open_micrograph(const FileName &_fn_micrograph,
   /*bool in_core,*/ bool reversed) _THROW {
   struct stat info;

   // Micrograph name
   fn_micrograph = _fn_micrograph;

   // Look for micrograph dimensions
   fn_inf=fn_micrograph.add_extension("inf");
   FILE *fh_inf=fopen(fn_inf.c_str(),"r");
   if (!fh_inf)
      REPORT_ERROR(1,(string)"Micrograph::open_micrograph: Cannot find "+
         fn_inf);
   Xdim=AtoI(get_param(fh_inf,"Xdim"));
   Ydim=AtoI(get_param(fh_inf,"Ydim"));
   __depth=AtoI(get_param(fh_inf,"bitspersample"));
   fclose(fh_inf);

   // Open micrograph and map
   fh_micrograph=open(fn_micrograph.c_str(),O_RDWR,S_IREAD|S_IWRITE);
   if (fh_micrograph==-1)
      REPORT_ERROR(1,(string)"Micrograph::open_micrograph: There is a "
         "problem opening "+fn_micrograph);
   switch (__depth) {
      case 8:
         /* if (!in_core) { */
            m8=(unsigned char *) mmap(0,Ydim*Xdim*__depth/8,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
            if (m8==MAP_FAILED)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
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
            m16=(short int *) mmap(0,Ydim*Xdim*__depth/8,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph, 0);
            if (m16==MAP_FAILED)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
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
            m32=(float*) mmap(0,Ydim*Xdim*__depth/8,
               PROT_READ|PROT_WRITE, MAP_SHARED, fh_micrograph,0);
            if (m32==MAP_FAILED)
               REPORT_ERROR(1,(string)"Micrograph::open_micrograph: cannot map "+
                  _fn_micrograph+" in memory");
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
         if      (__depth== 8) munmap((char *)m8,Ydim*Xdim*__depth/8);
         else if (__depth==16) munmap((char *)m16,Ydim*Xdim*__depth/8);
         else if (__depth==32) munmap((char *)m32,Ydim*Xdim*__depth/8);
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
   for (int j=0; j<Xdim; j++)
      for (int i=0; i<Ydim; i++) {
         if ((*this)(j,i)<minval) minval=(*this)(j,i);
         if ((*this)(j,i)>maxval) maxval=(*this)(j,i);
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
   __scaling_valid=TRUE;
}

/* Save coordinates to disk ------------------------------------------------ */
void Micrograph::write_coordinates(int label, const FileName &_fn_coords)
   _THROW {
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
void Micrograph::read_coordinates(int label, const FileName &_fn_coords)
   _THROW {
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
   aux.valid=TRUE;
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
   double scaleX, double scaleY) _THROW {
   if (X_window_size==-1 || Y_window_size==-1)
      REPORT_ERROR(1,"Micrograph::scissor: window size not set");

   result().resize(Y_window_size, X_window_size);
   int i0=ROUND(scaleY*P.Y)+FIRST_XMIPP_INDEX(Y_window_size);
   int iF=ROUND(scaleY*P.Y)+LAST_XMIPP_INDEX(Y_window_size);
   int j0=ROUND(scaleX*P.X)+FIRST_XMIPP_INDEX(X_window_size);
   int jF=ROUND(scaleX*P.X)+LAST_XMIPP_INDEX(X_window_size);
   int retval=1;
   if (i0<0 || iF>=Ydim || j0<0 || jF>=Xdim) {
      result().init_zeros();
      retval=0;
   } else
      for (int i=i0; i<=iF; i++)
         for (int j=j0; j<=jF; j++)
            result(i-i0,j-j0)=(*this)(j,i);
   return retval;
}

/* Produce all images ------------------------------------------------------ */
void Micrograph::produce_all_images(int label, const FileName &fn_root,
   int starting_index, const FileName &fn_image, double ang) _THROW {
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
   }
   
   // Set scale for particles
   int MXdim, MYdim, thisXdim, thisYdim;
   M->size(MXdim, MYdim);
   this->size(thisXdim, thisYdim);
   double scaleX=(double)MXdim/thisXdim;
   double scaleY=(double)MYdim/thisYdim;

   // Scissor all particles
   if (ang!=0)
      cout << "Angle from Y axis to tilt axis " << ang << endl
	   << "   applying apropriate rotation\n";
   int i=starting_index;
   int nmax=ParticleNo();
   for (int n=0; n<nmax; n++)
      if (coords[n].valid && coords[n].label==label) {
         fn_out.compose(fn_root,i++,"xmp");
         if (!M->scissor(coords[n],(Image &) I, scaleX, scaleY)) {
            cout << "Particle " << fn_out << " is very near the border, "
                 << "corresponding image is set to blank\n";
            SF.insert(fn_out,SelLine::DISCARDED);
         } else SF.insert(fn_out);
	 if (ang!=0) I().rotate(-ang);
         I.write(fn_out);
      }
   if (labels[label]!="") {
      SF.write(fn_micrograph+"."+labels[label]+".sel");
      write_coordinates(label,fn_micrograph+"."+labels[label]+".pos");
   } else {
      SF.write(fn_micrograph+".sel");
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
void Micrograph::invalidate_coord(int n) _THROW {
   if (n<0 || n>=ParticleNo())
      REPORT_ERROR(1,"Micrograph::invalidate_coord: Index out of range");
   coords[n].valid=FALSE;   
}

/* Add coordinate ---------------------------------------------------------- */
void Micrograph::add_coord(int x, int y, int label) {
   struct Particle_coords aux;
   aux.valid=TRUE;
   aux.X=x;
   aux.Y=y;
   aux.label=label;
   coords.push_back(aux);
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
      bool ifirst=TRUE, ofirst=TRUE;

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
        	if (ifirst) {imin=imax=M(j2,i2); ifirst=FALSE;}
        	else {imin=MIN(imin,M(j2,i2)); imax=MAX(imax,M(j2,i2));}
        	pixval += kernel(i,j)*M(j2,i2);
             }
             pixval *= scale;
             if (ii<Ypdim && jj<Xpdim) {
        	if (ofirst) {omin=omax=pixval; ofirst=FALSE;}
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
