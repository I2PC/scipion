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

#include "../xmippFilters.hh"
#include <slist.h>

/* Substract background ---------------------------------------------------- */
void substract_background_plane(Image *I) {
   matrix2D<double> A(3,3);
   matrix1D<double> x(3), b(3);
   
   // Solve the plane 'x'
   A.init_zeros(); b.init_zeros();
   FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(*I)) {
      A(0,0) += j*j;    A(0,1) +=j*i;    A(0,2) +=j;
      A(1,0)  = A(0,1); A(1,1) +=i*i;    A(1,2) +=i;
      A(2,0)  = A(0,2); A(2,1)  =A(1,2); A(2,2) +=1;
      b(0)   += j*IMGPIXEL(*I,i,j);   
      b(1)   += i*IMGPIXEL(*I,i,j);
      b(2) +=IMGPIXEL(*I,i,j);
   }
   solve(A,b,x);

   // Now substract the plane
   FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(*I))
      IMGPIXEL(*I,i,j) -= x(0)*i+x(1)*j+x(2);
}

/* Contranst enhancement --------------------------------------------------- */
void contrast_enhancement(Image *I) {
   (*I)().range_adjust(0,255);
}

/* Region growing for images ----------------------------------------------- */
 void region_growing(const Image *I_in, Image *I_out, int i, int j,
    float stop_colour, float filling_colour, bool less)
{
   slist<int> iNeighbours;       /* A list for neighbour pixels */
   int iCurrenti, iCurrentj;     /* Coordinates of the current pixel considered */
   
   /* First task is copying the input image into the output one */
   *I_out=*I_in;
   
   /**** Then the region growing is done ****/
   /* Insert at the beginning of the list the seed coordinates */
   iNeighbours.push_front(j);
   iNeighbours.push_front(i);

   /* Fill the seed coordinates */
   IMGPIXEL(*I_out,i,j)=filling_colour;
   
   while(!iNeighbours.empty())
   {
      matrix1D<double> r(2);

      /* Take the current pixel to explore */
      iCurrenti=iNeighbours.front();
      iNeighbours.pop_front();
      iCurrentj=iNeighbours.front();
      iNeighbours.pop_front();

      #define CHECK_POINT(i,j) \
         XX(r)=j; YY(r)=i; \
         if (!(*I_out)().outside(r))  { \
            if (IMGPIXEL(*I_out,i,j)!=filling_colour) \
               if ((less && IMGPIXEL (*I_out,i,j) < stop_colour) || \
                   (!less && IMGPIXEL (*I_out,i,j) > stop_colour)) { \
                  IMGPIXEL (*I_out,i,j)=filling_colour; \
                  iNeighbours.push_front(j); \
                  iNeighbours.push_front(i); \
               } \
         }

      /* Make the exploration of the pixel´s neighbours */
      CHECK_POINT(iCurrenti  ,iCurrentj-1);
      CHECK_POINT(iCurrenti  ,iCurrentj+1);
      CHECK_POINT(iCurrenti-1,iCurrentj  );
      CHECK_POINT(iCurrenti-1,iCurrentj-1);
      CHECK_POINT(iCurrenti-1,iCurrentj+1);
      CHECK_POINT(iCurrenti+1,iCurrentj  );
      CHECK_POINT(iCurrenti+1,iCurrentj-1);
      CHECK_POINT(iCurrenti+1,iCurrentj+1);
   }             
}

/* Region growing for volumes ----------------------------------------------- */
void region_growing(const Volume *V_in, Volume *V_out, int k, int i, int j,
   float stop_colour, float filling_colour, bool less)
{
   slist<int> iNeighbours;       /* A list for neighbour voxels */
   int iCurrentk, iCurrenti, iCurrentj;     /* Coordinates of the current voxel considered */
   
   /* First task is copying the input volume into the output one */
   *V_out=*V_in;
   
   /**** Then the region growing is done in output volume ****/
   /* Insert at the beginning of the list the seed coordinates */
   iNeighbours.push_front(j);
   iNeighbours.push_front(i);
   iNeighbours.push_front(k);
   
   /* Fill the seed coordinates */
   VOLVOXEL(*V_out,k,i,j)=filling_colour;
   
   while(!iNeighbours.empty())
   {
      matrix1D<double> r(3);

      /* Take the current pixel to explore */
      iCurrentk=iNeighbours.front();
      iNeighbours.pop_front();     
      iCurrenti=iNeighbours.front();
      iNeighbours.pop_front();
      iCurrentj=iNeighbours.front();
      iNeighbours.pop_front();

      /* a macro for doing exploration of a voxel. If the voxel has a value 
      lower than stop_colour, its filled with filling colour and added to the
      list for exploring its neighbours */
      #define CHECK_POINT_3D(k,i,j) \
         XX(r)=j; YY(r)=i; ZZ(r)=k; \
         if (!(*V_out)().outside(r))  { \
            if (VOLVOXEL(*V_out,k,i,j)!=filling_colour) \
               if ((less && VOLVOXEL (*V_out,k,i,j) < stop_colour)|| \
                   (!less &&VOLVOXEL (*V_out,k,i,j) > stop_colour)) { \
                  VOLVOXEL (*V_out,k,i,j)=filling_colour; \
                  iNeighbours.push_front(j); \
                  iNeighbours.push_front(i); \
                  iNeighbours.push_front(k); \
               } \
         }

      /* Make the exploration of the pixel´s neighbours */
      CHECK_POINT_3D(iCurrentk  ,iCurrenti  ,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk  ,iCurrenti  ,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk  ,iCurrenti-1,iCurrentj  );
      CHECK_POINT_3D(iCurrentk  ,iCurrenti-1,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk  ,iCurrenti-1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk  ,iCurrenti+1,iCurrentj  );
      CHECK_POINT_3D(iCurrentk  ,iCurrenti+1,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk  ,iCurrenti+1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk-1,iCurrenti  ,iCurrentj  );
      CHECK_POINT_3D(iCurrentk-1,iCurrenti  ,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk-1,iCurrenti  ,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk-1,iCurrenti-1,iCurrentj  );
      CHECK_POINT_3D(iCurrentk-1,iCurrenti-1,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk-1,iCurrenti-1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk-1,iCurrenti+1,iCurrentj  );
      CHECK_POINT_3D(iCurrentk-1,iCurrenti+1,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk-1,iCurrenti+1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti  ,iCurrentj  );
      CHECK_POINT_3D(iCurrentk+1,iCurrenti  ,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti  ,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti-1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti-1,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti-1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti+1,iCurrentj+1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti+1,iCurrentj-1);
      CHECK_POINT_3D(iCurrentk+1,iCurrenti+1,iCurrentj+1);
   }             
}

/* Label image ------------------------------------------------------------ */
int label_image(const Image *I, Image *label) {
   (*label)()=(*I)();
   int colour=32000;
   bool found;
   FOR_ALL_ELEMENTS_IN_MATRIX2D((*label)()) {
      if ((*label)(i,j)!=1) continue;
      region_growing(label,label,i,j,0,colour,FALSE);
      colour++;
   }
   FOR_ALL_ELEMENTS_IN_MATRIX2D((*label)())
      if ((*label)(i,j)!=0) (*label)(i,j)=(*label)(i,j)-31999;
   return colour-32000;
}

/* Label volume ------------------------------------------------------------ */
int label_volume(const Volume *V, Volume *label) {
   (*label)()=(*V)();
   int colour=32000;
   bool found;
   FOR_ALL_ELEMENTS_IN_MATRIX3D((*label)()) {
      if ((*label)(k,i,j)!=1) continue;
      region_growing(label,label,k,i,j,0,colour,FALSE);
      colour++;
   }
   FOR_ALL_ELEMENTS_IN_MATRIX3D((*label)())
      if ((*label)(k,i,j)!=0) (*label)(k,i,j)=(*label)(k,i,j)-31999;
   return colour-32000;
}

/* Correlation ------------------------------------------------------------- */
template <class T> 
double correlation(matrix1D<T> &x,matrix1D<T> &y,const matrix1D<int> *mask,int l)
{
   /* Note: l index is for rows and m index for columns */
   
   SPEED_UP_temps;
   double retval=0;        // returned value
   int i,ip;               // indexes 
   int Rows;               // of the matrices
   
   Rows=x.get_dim();
   for(i=0;i<Rows;i++)
   {   
      ip=i-l;
      if(ip>=0 && ip<Rows)
	  {
		 if (mask!=NULL)
            if (!(*mask)(i)) continue;
		 retval+=DIRECT_VEC_ELEM(x,i)*DIRECT_VEC_ELEM(y,ip);
	  }
   }
   
//     return retval/N;
     return retval/(Rows);

}

template <class T>
double correlation(matrix2D<T> &x,matrix2D<T> &y,
                    const matrix2D<int> *mask, int l, int m)
{
   /* Note: l index is for rows and m index for columns */
   
   SPEED_UP_temps;
   double retval=0;         // returned value
   int i,j,ip,jp;           // indexes 
   int Rows, Cols;          // of the matrices
   // do the computation
   Cols=x.ColNo();
   Rows=x.RowNo();
   
   for(i=0;i<Rows;i++)
      for(j=0;j<Cols;j++)
      {
      	 ip=i-l;
	     jp=j-m;
      	 if(ip>=0 && ip<Rows && jp>=0 && jp<Cols)
	     {
		    if (mask!=NULL)
               if (!(*mask)(i,j)) continue;
            retval+=DIRECT_MAT_ELEM(x,i,j)*DIRECT_MAT_ELEM(y,ip,jp);
	     }
      }
    
//    return retval/N;
   return retval/(Cols*Rows);
}

template <class T>
double correlation(matrix3D<T> &x, matrix3D<T> &y,
          const matrix3D<int> *mask,int l, int m, int q) 
{
   /* Note: l index is for rows and m index for columns */
   
   SPEED_UP_temps;
   double retval=0;                 // returned value
   int i,j,k,ip,jp,kp;             // indexes 
   int Rows, Cols, Slices;          // of the volumes
   
    // do the computation
   Cols=x.ColNo(); Rows=x.RowNo(); Slices=x.SliNo();
   long N=0;   
   for(k=0;k<Slices;k++)
      for(i=0;i<Rows;i++)
         for(j=0;j<Cols;j++)
         {
      	    ip=i-l;
	        jp=j-m;
			kp=k-q;
      	    if(ip>=0 && ip<Rows && jp>=0 && jp<Cols && kp>=0 && kp<Slices)
	        {
		       if (mask!=NULL)
                  if (!(*mask)(k,i,j)) continue;
			   retval+=VOL_ELEM(x,k,i,j)*VOL_ELEM(y,kp,ip,jp);
	        }
         }
   
 
//   return retval/N;
   return retval/(Slices*Rows*Cols);
}

/* correlation_index ------------------------------------------------------------ */
template <class T>
   double correlation_index(const matrix1D<T> &x, const matrix1D<T> &y)
{
   SPEED_UP_temps;
   double retval=0;
   double mean_x, mean_y;
   double stddev_x, stddev_y;
   double aux;
   long n=0;
   
  
   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x,y)
   {
         retval+=(VEC_ELEM(x,i)-mean_x)*(VEC_ELEM(y,i)-mean_y);
         n++;
   }

   if (n!=0) return retval/((stddev_x*stddev_y)*n); else return 0;
}

template <class T>
   double correlation_index(const matrix2D<T> &x, const matrix2D<T> &y,
      const matrix2D<int> *mask,matrix2D<double> *Contributions) {
   SPEED_UP_temps;
   double retval=0,aux;
   double mean_x, mean_y;
   double stddev_x, stddev_y;
   T dummy;  
   long n=0;

   if(mask==NULL)
   {
       x.compute_stats(mean_x, stddev_x, dummy, dummy);
       y.compute_stats(mean_y, stddev_y, dummy, dummy);
   
   }
   else
   {
       compute_stats_within_binary_mask( *mask, x, dummy, dummy, mean_x, stddev_x);
       compute_stats_within_binary_mask( *mask, y, dummy, dummy, mean_y, stddev_y);
   }   

   // If contributions are desired. Please, be careful interpreting individual 
   // contributions to the covariance! One pixel value afect others.
   if(Contributions!=NULL)
   {

      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(i,j)) continue;
		 aux=(MAT_ELEM(x,i,j)-mean_x)*(MAT_ELEM(y,i,j)-mean_y);
		 MAT_ELEM(*Contributions,i,j)=aux;
		 retval+=aux;
		 n++;
      }

      FOR_ALL_ELEMENTS_IN_MATRIX2D(*Contributions)
	  {
		 MAT_ELEM(*Contributions,i,j)/=((stddev_x*stddev_y)*n);
      }
	  
   }
   // In other case, normal process (less computationaly expensive)   
   else
   {

      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(i,j)) continue;
		 retval+=(MAT_ELEM(x,i,j)-mean_x)*(MAT_ELEM(y,i,j)-mean_y);;
		 n++;
      }

   }
   
   if (n!=0) return retval/((stddev_x*stddev_y)*n); else return 0;
}

template <class T>
   double correlation_index(const matrix3D<T> &x, const matrix3D<T> &y,
      const matrix3D<int> *mask,matrix3D<double> *Contributions) {
   SPEED_UP_temps;
   double retval=0,aux;
   T dummy;  
   double mean_x, mean_y;
   double stddev_x, stddev_y;
	  
   long n=0;
   
   if(mask==NULL)
   {
       x.compute_stats(mean_x, stddev_x, dummy, dummy);
       y.compute_stats(mean_y, stddev_y, dummy, dummy);
   
   }
   else
   {
       compute_stats_within_binary_mask( *mask, x, dummy, dummy, mean_x, stddev_x);
       compute_stats_within_binary_mask( *mask, y, dummy, dummy, mean_y, stddev_y);
   }
   
   // If contributions are desired. Please, be careful interpreting individual 
   // contributions to the covariance! One pixel value afect others.
   if(Contributions!=NULL)
   {
       FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x,y)
       {
          if (mask!=NULL)
             if (!(*mask)(k,i,j)) continue;
          aux=(VOL_ELEM(x,k,i,j)-mean_x)*(VOL_ELEM(y,k,i,j)-mean_y);;
		  VOL_ELEM(*Contributions,k,i,j)=aux;
          retval+=aux;
		  n++;
       }

      FOR_ALL_ELEMENTS_IN_MATRIX3D(*Contributions)
	  {
		 VOL_ELEM(*Contributions,k,i,j)/=((stddev_x*stddev_y)*n);
      }
	   
   }
   else
   {
       FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x,y)
       {
          if (mask!=NULL)
             if (!(*mask)(k,i,j)) continue;
          retval+=(VOL_ELEM(x,k,i,j)-mean_x)*(VOL_ELEM(y,k,i,j)-mean_y);
          n++;
       }
   }
   
   if (n!=0) return retval/((stddev_x*stddev_y)*n); else return 0;
}

/* RMS --------------------------------------------------------------------- */
template <class T>
   double rms(const matrix1D<T> &x, const matrix1D<T> &y,
      const matrix1D<int> *mask,matrix1D<double> *Contributions) {
   SPEED_UP_temps;
   double retval=0,aux;
   int n=0;
 
   // If contributions are desired.
   if(Contributions!=NULL)
   {
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(i)) continue;
		 aux=(VEC_ELEM(x,i)-VEC_ELEM(y,i))*(VEC_ELEM(x,i)-VEC_ELEM(y,i));
		 VEC_ELEM(*Contributions,i)=aux;
		 retval+=aux;
		 n++;
      }

      FOR_ALL_ELEMENTS_IN_MATRIX1D(*Contributions)
	  {
		 VEC_ELEM(*Contributions,i)=sqrt(VEC_ELEM(*Contributions,i)/n);
      }
	  
   
   }
   // In other case, normal process (less computationaly expensive)
   else
   {
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(i)) continue;
         retval+=(VEC_ELEM(x,i)-VEC_ELEM(y,i))*(VEC_ELEM(x,i)-VEC_ELEM(y,i));
         n++;
      }
   
   }
   if (n!=0) return sqrt(retval/n); else return 0;
}

template <class T>
   double rms(const matrix2D<T> &x, const matrix2D<T> &y,
      const matrix2D<int> *mask,matrix2D<double> *Contributions) {
   SPEED_UP_temps;
   double retval=0,aux;
   int n=0;
   
      // If contributions are desired.
   if(Contributions!=NULL)
   {
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(i,j)) continue;
 
 		 aux=(MAT_ELEM(x,i,j)-MAT_ELEM(y,i,j))*(MAT_ELEM(x,i,j)-MAT_ELEM(y,i,j));
		 MAT_ELEM(*Contributions,i,j)=aux;
		 retval+=aux;
		 n++;
      }

      FOR_ALL_ELEMENTS_IN_MATRIX2D(*Contributions)
	  {
		 MAT_ELEM(*Contributions,i,j)=sqrt(MAT_ELEM(*Contributions,i,j)/n);
      }
	  
   
   }
   // In other case, normal process (less computationaly expensive)
   else
   {
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(i,j)) continue;
         retval+=(MAT_ELEM(x,i,j)-MAT_ELEM(y,i,j))*
            (MAT_ELEM(x,i,j)-MAT_ELEM(y,i,j));
         n++;
      }
   }

   if (n!=0) return sqrt(retval/n); else return 0;
}

template <class T>
   double rms(const matrix3D<T> &x, const matrix3D<T> &y,
      const matrix3D<int> *mask,matrix3D<double> *Contributions) {
   SPEED_UP_temps;
   double retval=0;
   double aux;
   int n=0;

   // If contributions are desired.
   if(Contributions!=NULL)
   {
   
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x,y)
	  {
         if (mask!=NULL)
            if (!(*mask)(k,i,j)) continue;
		 aux=(VOL_ELEM(x,k,i,j)-VOL_ELEM(y,k,i,j))*(VOL_ELEM(x,k,i,j)-VOL_ELEM(y,k,i,j));
		 VOL_ELEM(*Contributions,k,i,j)=aux;
		 retval+=aux;
        n++;
      }
	  
      FOR_ALL_ELEMENTS_IN_MATRIX3D(*Contributions)
	  {
		 VOL_ELEM(*Contributions,k,i,j)=sqrt(VOL_ELEM(*Contributions,k,i,j)/n);;
      }

   }
   else
   {
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x,y) {
         if (mask!=NULL)
            if (!(*mask)(k,i,j)) continue;
         retval+=(VOL_ELEM(x,k,i,j)-VOL_ELEM(y,k,i,j))*
            (VOL_ELEM(x,k,i,j)-VOL_ELEM(y,k,i,j));
         n++;
      }   
   }
   if (n!=0) return sqrt(retval/n); else return 0;
}

/* Fourier-Bessel decomposition. ------------------------------------------- */
void Fourier_Bessel_decomposition(const matrix2D<double> &img_in,
   matrix2D<double> &m_out, double r1, double r2, int k1, int k2) {
   for (int k=k1; k<=k2; k++) {
      int k_1=k-1;

      // Compute h and a,b coefficients
      double coefca=0,coefcb=0,coefsa=0,coefsb=0;
      double h=0, my5=0;
      if (k_1!=0) {
          double my=1+PI*r2/2/k_1;
          double my2=2*my;
          double my4=my*k_1;
          double my5=my4-1;
          double ntot=4*my4;
                 h=2*PI/ntot;
          double hdpi=h/PI;
          double th=k_1*h;
          double ys=sin (th);
          double zs=cos (th);
          double ys2 = sin (2*th);
          double b1=2/(th*th)*(1 + zs*zs - ys2/th);
          double g1=4/(th*th)*(ys/th-zs);
          double d1=2*th/45*ys2;
          double e1=d1*ys*2;
                 coefca=(b1+e1)*hdpi;
                 coefcb=(g1-d1)*hdpi;
                 coefsa=(b1-e1)*hdpi;
                 coefsb=(g1+d1)*hdpi;
       } else {
          double my=1+PI*r2/2;
          double my2=2*my;
          double my4=my;
          double my5=my4-1;
          double ntot=4*my4;
                 h=2*PI/ntot;
                 coefca=h/PI/2.;
     }
     
     matrix1D<double> sine(CEIL(my5));
     FOR_ALL_ELEMENTS_IN_MATRIX1D(sine) sine(i)=sin((i+1)*h);
     
   }
}

/* Harmonic decomposition. ------------------------------------------------- */
void harmonic_decomposition(const matrix2D<double> &img_in,
   matrix1D<double> &v_out) {
}

/* Median filter -----------------------------------------------------------*/
template <class T>
void sort(T a,T b,T c,matrix1D<T> &v) 
{                
		if(a<b)                       
		  if(b<c)                     
		  {  v(0)=a;v(1)=b;v(2)=c; }  
		  else if(a<c)                
		  {  v(0)=a;v(1)=c;v(2)=b; }  
		  else                        
		  {  v(0)=c;v(1)=a;v(2)=b; }  
		else                          
		  if(a<c)                     
		  {  v(0)=b;v(1)=a;v(2)=c; }  
		  else if(b<c)                
		  {  v(0)=b;v(1)=c;v(2)=a; }  
		  else                        
		  {  v(0)=c;v(1)=b;v(2)=a; } 
}		  

template <class T>
void merge_sort(matrix1D<T> &v1,matrix1D<T> &v2,matrix1D<T> &v)  
{
		  int i1=0,i2=0,i=0; 
		  while((i1 < 3) && (i2 < 3)) 
		  { 
    		if(v1(i1) < v2(i2)) 
    		  v(i++) = v1(i1++); 
    		else 
    		  v(i++) = v2(i2++); 
		  } 
		  while(i1 < 3) 
    		 v(i++) = v1(i1++); 
		  while(i2 < 3) 
    		 v(i++) = v2(i2++); 
}
			 
// This ugly function performs a fast merge sort for the case of vectors of 3
// elements. This way is guaranteed a minimum number of comparisons (maximum
// number of comparisons to perform the sort, 5)
template <class T> 
void fast_merge_sort(matrix1D<T> &x,matrix1D<T> &y,matrix1D<T> &v)
{                          
if(x(0)<y(0))                                           
{														
   v(0)=x(0);											
   if(x(1)<y(0))										
   {													
      v(1)=x(1);										
	  if(x(2)<y(0)) 									
	  { v(2)=x(2);v(3)=y(0);v(4)=y(1);v(5)=y(2); }		
	  else												
	  { 												
	    v(2)=y(0);										
		if(x(2)<y(1))									
		{ v(3)=x(2);v(4)=y(1);v(5)=y(2); }				
		else											
		{												
		  v(3)=y(1);                                    
		  if(x(2)<y(2)) 								
		  { v(4)=x(2);v(5)=y(2); }						
		  else											
		  { v(4)=y(2);v(5)=x(2); }						
		}												
	  } 												
   }													
   else 												
   {													
      v(1)=y(0);										
	  if(x(1)<y(1)) 									
	  { 												
	     v(2)=x(1); 									
		 if(x(2)<y(1))									
		 { v(3)=x(2);v(4)=y(1);v(5)=y(2); } 			
		 else                                           
		 {												
		   v(3)=y(1);									
		   if(x(2)<y(2))								
		   { v(4)=x(2);v(5)=y(2); } 					
		   else 										
		   { v(4)=y(2);v(5)=x(2); } 					
		 }												
	  } 												
	  else												
	  { 												
	     v(2)=y(1); 									
		 if(x(1)<y(2))									
		 {												
		    v(3)=x(1);									
			if(x(2)<y(2))								
		    { v(4)=x(2);v(5)=y(2); }                    
		    else										
		    { v(4)=y(2);v(5)=x(2); }					
		 }												
		 else											
		 { v(3)=y(2); v(4)=x(1); v(5)=x(2); }		 	
	  } 												
   }													
}														
else													
{														
   v(0)=y(0);											
   if(x(0)<y(1))										
   {													
      v(1)=x(0);										
	  if(x(1)<y(1)) 									
	  {                                                 
	     v(2)=x(1);                                     
		 if(x(2)<y(1))									
		 { v(3)=x(2);v(4)=y(1);v(5)=y(2); } 			
		 else											
		 {												
		   v(3)=y(1);									
		   if(x(2)<y(2))								
		   { v(4)=x(2);v(5)=y(2); } 					
		   else 										
		   { v(4)=y(2);v(5)=x(2); } 					
		 }												
	  } 												
	  else												
	  { 												
	     v(2)=y(1); 									
		 if(x(1)<y(2))                                  
		 {                                              
		    v(3)=x(1);                                  
			if(x(2)<y(2))								
		    { v(4)=x(2);v(5)=y(2); }					
		    else										
		    { v(4)=y(2);v(5)=x(2); }					
		 }												
		 else											
		 { v(3)=y(2);v(4)=x(1);v(5)=x(2); } 			
	  } 												
   }													
   else 												
   {													
      v(1)=y(1);										
	  if(x(0)<y(2)) 									
	  {                                                 
	     v(2)=x(0); 									
		 if(x(1)<y(2))									
		 {												
		    v(3)=x(1);									
			if(x(2)<y(2))								
		    { v(4)=x(2);v(5)=y(2); }					
		    else										
		    { v(4)=y(2);v(5)=x(2); }					
		 }												
		 else											
		 { v(3)=y(2);v(4)=x(1);v(5)=x(2); } 			
	  } 												
	  else												
	  { v(2)=y(2);v(3)=x(0);v(4)=x(1);v(5)=x(2); }		
   }													
}
}	
			 
template <class T> 
void median(matrix1D<T> &x,matrix1D<T> &y,T &m) 
{                       
if(x(0)<y(1))								 
   if(x(1)<y(1))							 
      if(x(2)<y(1))  m=y(1);				 
	  else           m=x(2);				 
   else 									 
      if(x(1)<y(2))  m=y(2);				 
      else           m=x(1);				 
else										 
   if(x(0)<y(2))							 
      if(x(1)<y(2))  m=y(2);				 
	  else           m=x(1);				 
   else 									 
      if(x(0)<y(3))  m=y(3);				 
	  else           						 
	     if(x(0)<y(4)) m=x(0);				 
		 else          m=y(4);				 
}


template <class T> 
void median_filter3x3(matrix2D<T> &m,matrix2D<T> &out)
{
    matrix1D<T> v1(3),v2(3),v3(3),v4(3);
	matrix1D<T> v(6); 
	
    // Set the output matrix size
	out.resize(m);
	// Set the initial and final matrix indices to explore
	int initialY=1, initialX=1;
	int finalY=m.RowNo()-2;
	int finalX=m.ColNo()-2;
	// For every row 
	for(int i=initialY;i<=finalY;i++)
	{
	    // For every pair of pixels (mean is computed obtaining
		// two means at the same time using an efficient method)
    	for(int j=initialX;j<=finalX;j+=2)
		{
		   // If we are in the border case
		   if(j==1)
		   {
		      // Order the first and second vectors of 3 elements
			  sort(DIRECT_MAT_ELEM(m,i-1,j-1),DIRECT_MAT_ELEM(m,i,j-1),
			       DIRECT_MAT_ELEM(m,i+1,j-1),v1);
			  sort(DIRECT_MAT_ELEM(m,i-1,j),DIRECT_MAT_ELEM(m,i,j),
			       DIRECT_MAT_ELEM(m,i+1,j),v2);
		   }
		   else
		   {
		      // Simply take ordered vectors from previous
			  v1=v3;
			  v2=v4;
		   }
		   
	       // As we are computing 2 medians at the same time, if the matrix has
	       // an odd number of columns, the last column isn't calculated. Its done here
		   if(j==finalX)
		   { 
		      v1=v3;
			  v2=v4; 
			  sort(DIRECT_MAT_ELEM(m,i-1,j+1),DIRECT_MAT_ELEM(m,i,j+1),
			       DIRECT_MAT_ELEM(m,i+1,j+1),v3);
              fast_merge_sort(v2,v3,v);
		      median(v1,v,out(i,j));
		   }
		   else
		   {
		       // Order the third and fourth vectors of 3 elements
		          sort(DIRECT_MAT_ELEM(m,i-1,j+1),DIRECT_MAT_ELEM(m,i,j+1),
		        	   DIRECT_MAT_ELEM(m,i+1,j+1),v3);
		          sort(DIRECT_MAT_ELEM(m,i-1,j+2),DIRECT_MAT_ELEM(m,i,j+2),
		        	   DIRECT_MAT_ELEM(m,i+1,j+2),v4);
		       // Merge sort the second and third vectors
		       fast_merge_sort(v2,v3,v);
		       // Find the first median and assign it to the output
		       median(v1,v,out(i,j));
		       // Find the second median and assign it to the output
		       median(v4,v,out(i,j+1));
		   }
		}
	}
		
}

/* Instantiation ----------------------------------------------------------- */
template <class T>
   void instantiate_filtersT(T t) {
   matrix1D<T> x1,y1;
   matrix2D<T> x2,y2;
   matrix3D<T> x3,y3;

   matrix1D<int> mask1;
   matrix2D<int> mask2;
   matrix3D<int> mask3;
   matrix1D<double> cont1;
   matrix2D<double> cont2;
   matrix3D<double> cont3;

   double mean;
   int position;

   correlation(x3,y3,&mask3,position,position,position);
   correlation(x2,y2,&mask2,position,position);
   correlation(x1,y1,&mask1,position);
   correlation(x1,y1,&mask1);
   correlation(x2,y2,&mask2);
   correlation(x3,y3,&mask3);
   correlation(x1,y1);
   correlation(x2,y2);
   correlation(x3,y3);

   correlation_index(x1,y1);
   correlation_index(x2,y2);
   correlation_index(x3,y3);
   correlation_index(x2,y2,&mask2,&cont2);
   correlation_index(x3,y3,&mask3,&cont3);
   
   rms(x1,y1);
   rms(x2,y2);
   rms(x3,y3);
   rms(x1,y1,&mask1,&cont1);
   rms(x2,y2,&mask2,&cont2);
   rms(x3,y3,&mask3,&cont3);
   
   T  m;
    sort(m,m,m,x1);
	merge_sort(x1,y1,x1);
	fast_merge_sort(x1,y1,y1);  
    median(x1,x1,m);
   median_filter3x3(x2,y2);
}

void instantiate_filters() {
    int    i; instantiate_filtersT(i);
    double d; instantiate_filtersT(d);
}

