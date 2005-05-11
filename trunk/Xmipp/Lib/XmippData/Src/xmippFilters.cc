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
#include "../xmippFFT.hh"
#include<list>

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
 void region_growing(const matrix2D<double> &I_in, matrix2D<double> &I_out,
    int i, int j,
    float stop_colour, float filling_colour, bool less, int neighbourhood)
{
   list<int> iNeighbours;       /* A list for neighbour pixels */
   int iCurrenti, iCurrentj;     /* Coordinates of the current pixel considered */
   
   /* First task is copying the input image into the output one */
   I_out=I_in;
   
   /**** Then the region growing is done ****/
   /* Insert at the beginning of the list the seed coordinates */
   iNeighbours.push_front(j);
   iNeighbours.push_front(i);

   /* Fill the seed coordinates */
   MAT_ELEM(I_out,i,j)=filling_colour;
   
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
         if (!I_out.outside(r))  { \
            if (MAT_ELEM(I_out,i,j)!=filling_colour) \
               if ((less && MAT_ELEM (I_out,i,j) < stop_colour) || \
                   (!less && MAT_ELEM (I_out,i,j) > stop_colour)) { \
                  MAT_ELEM (I_out,i,j)=filling_colour; \
                  iNeighbours.push_front(j); \
                  iNeighbours.push_front(i); \
               } \
         }

      /* Make the exploration of the pixel´s neighbours */
      CHECK_POINT(iCurrenti  ,iCurrentj-1);
      CHECK_POINT(iCurrenti  ,iCurrentj+1);
      CHECK_POINT(iCurrenti-1,iCurrentj  );
      CHECK_POINT(iCurrenti+1,iCurrentj  );
      if (neighbourhood==8) {
	 CHECK_POINT(iCurrenti-1,iCurrentj-1);
	 CHECK_POINT(iCurrenti-1,iCurrentj+1);
	 CHECK_POINT(iCurrenti+1,iCurrentj-1);
	 CHECK_POINT(iCurrenti+1,iCurrentj+1);
      }
   }             
}

/* Region growing for volumes ----------------------------------------------- */
void region_growing(const matrix3D<double> &V_in, matrix3D<double> &V_out,
   int k, int i, int j,
   float stop_colour, float filling_colour, bool less)
{
   list<int> iNeighbours;       /* A list for neighbour voxels */
   int iCurrentk, iCurrenti, iCurrentj;     /* Coordinates of the current voxel considered */
   
   /* First task is copying the input volume into the output one */
   V_out=V_in;
   
   /**** Then the region growing is done in output volume ****/
   /* Insert at the beginning of the list the seed coordinates */
   iNeighbours.push_front(j);
   iNeighbours.push_front(i);
   iNeighbours.push_front(k);
   
   /* Fill the seed coordinates */
   VOL_ELEM(V_out,k,i,j)=filling_colour;
   
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
         if (!V_out.outside(r))  { \
            if (VOL_ELEM(V_out,k,i,j)!=filling_colour) \
               if ((less && VOL_ELEM (V_out,k,i,j) < stop_colour)|| \
                   (!less &&VOL_ELEM (V_out,k,i,j) > stop_colour)) { \
                  VOL_ELEM (V_out,k,i,j)=filling_colour; \
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
int label_image(const matrix2D<double> &I, matrix2D<double> &label,
   int neighbourhood) {
   label=I;
   int colour=32000;
   bool found;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label) {
      if (label(i,j)!=1) continue;
      region_growing(label,label,i,j,0,colour,FALSE,neighbourhood);
      colour++;
   }
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
      if (label(i,j)!=0) label(i,j)=label(i,j)-31999;
   return colour-32000;
}

/* Label volume ------------------------------------------------------------ */
int label_volume(const matrix3D<double> &V, matrix3D<double> &label) {
   label=V;
   int colour=32000;
   bool found;
   FOR_ALL_ELEMENTS_IN_MATRIX3D(label) {
      if (label(k,i,j)!=1) continue;
      region_growing(label,label,k,i,j,0,colour,FALSE);
      colour++;
   }
   FOR_ALL_ELEMENTS_IN_MATRIX3D(label)
      if (label(k,i,j)!=0) label(k,i,j)=label(k,i,j)-31999;
   return colour-32000;
}

/* Remove small components ------------------------------------------------- */
void remove_small_components(matrix2D<double> &I, int size,
   int neighbourhood) {
   matrix2D<double> label;
   int imax=label_image(I,label,neighbourhood);
   matrix1D<int> nlabel(imax+1);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label) nlabel((int)(label(i,j)))++;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
      if (nlabel((int)(label(i,j)))<size) I(i,j)=0;
}

/* Keep biggest component -------------------------------------------------- */
void keep_biggest_component(matrix2D<double> &I, double percentage,
   int neighbourhood) {
   matrix2D<double> label;
   int imax=label_image(I,label,neighbourhood);
   matrix1D<int> nlabel(imax+1);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label) {
      int idx=(int)(label(i,j));
      if (idx==0) continue;
      nlabel(idx)++;
   }
   matrix1D<int> best=nlabel.index_sort();
   best-=1;
   int nbest=XSIZE(best)-1;
   double total=nlabel.sum();
   double explained=nlabel(best(nbest));
   while (explained<percentage*total) {
      nbest--;
      explained+=nlabel(best(nbest));
   }
   
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label) {
      bool among_the_best=false;
      for (int n=nbest; n<imax+1; n++)
         if (label(i,j)==best(n)) {among_the_best=true; break;}
      if (!among_the_best) I(i,j)=0;
   }
}

/* Fill object ------------------------------------------------------------- */
void fill_binary_object(matrix2D<double> &I, int neighbourhood) {
   matrix2D<double> label;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(I) I(i,j)=1-I(i,j);
   int imax=label_image(I,label,neighbourhood);
   double l0=label(STARTINGY(I),STARTINGX(I));
   FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
      if (label(i,j)==l0) I(i,j)=0;
      else                I(i,j)=1;
}

/* Correlation ------------------------------------------------------------- */
template <class T> 
double correlation(const matrix1D<T> &x,const matrix1D<T> &y,const matrix1D<int> *mask,int l)
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
		   //            if (!(*mask)(i)) continue;
               if (!DIRECT_VEC_ELEM((*mask),i)) continue;
		 retval+=DIRECT_VEC_ELEM(x,i)*DIRECT_VEC_ELEM(y,ip);
	  }
   }
   
//     return retval/N;
     return retval/(Rows);

}

template <class T>
double correlation(const matrix2D<T> &x,const matrix2D<T> &y,
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
	//  Sjors: 
	//if (!(*mask)(i,j)) continue;
               if (!DIRECT_MAT_ELEM((*mask),i,j)) continue;
		    retval+=DIRECT_MAT_ELEM(x,i,j)*DIRECT_MAT_ELEM(y,ip,jp);
	     }
      }
    
//    return retval/N;
   return retval/(Cols*Rows);
}

template <class T>
double correlation(const matrix3D<T> &x, const matrix3D<T> &y,
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
			 //  Sjors: 
			 //if (!(*mask)(k,i,j)) continue;
			 //retval+=VOL_ELEM(x,k,i,j)*VOL_ELEM(y,kp,ip,jp);
			 if (!DIRECT_VOL_ELEM((*mask),k,i,j)) continue;
		       retval+=DIRECT_VOL_ELEM(x,k,i,j)*DIRECT_VOL_ELEM(y,kp,ip,jp);
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
   T dummy;
   long n=0;
   
   x.compute_stats(mean_x, stddev_x, dummy, dummy);
   y.compute_stats(mean_y, stddev_y, dummy, dummy);

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

/* Best shift -------------------------------------------------------------- */
void best_shift(const matrix2D<double> &I1, const matrix2D<double> &I2,
   double &shiftX, double &shiftY, const matrix2D<int> *mask) {
   int              imax,jmax,i_actual,j_actual;
   double           max,xmax,ymax,sumcorr,avecorr,stdcorr,dummy;
   float            xshift,yshift,shift;
   int              n_max=-1;
   bool             neighbourhood=TRUE;
   matrix2D<double> Mcorr;

   correlation_matrix(I1,I2,Mcorr);
 
   // Adjust statistics within shiftmask to average 0 and stddev 1
   if (mask!=NULL) {
      compute_stats_within_binary_mask(*mask,Mcorr,dummy,dummy,avecorr,stdcorr);
      FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr)
        if (MAT_ELEM(*mask,i,j)) 
	   MAT_ELEM(Mcorr,i,j)=(MAT_ELEM(Mcorr,i,j)-avecorr)/stdcorr;
        else MAT_ELEM(Mcorr,i,j)=0.;
   } else
      Mcorr.statistics_adjust(0,1);
   Mcorr.max_index(imax,jmax);
   max=MAT_ELEM(Mcorr,imax,jmax);

   while (neighbourhood) {
      n_max ++;
      for (int i=-n_max; i <= n_max; i++)
          for (int j=-n_max; j <= n_max; j++) {
              i_actual = i+imax;
              j_actual = j+jmax;
              if (i_actual < STARTINGY(Mcorr)  || j_actual < STARTINGX(Mcorr) ||
        	  i_actual > FINISHINGY(Mcorr) || j_actual > FINISHINGX(Mcorr) )
        	 neighbourhood=FALSE;
              else if (max/1.414 > MAT_ELEM(Mcorr,i_actual,j_actual))
        	  neighbourhood=FALSE;
	  }
   }

   // We have the neighbourhood => looking for the gravity centre
   xmax = ymax = sumcorr = 0.;
   for (int i=-n_max; i <= n_max; i++)
       for (int j=-n_max; j <= n_max; j++) {
	 i_actual = i+imax;
	 j_actual = j+jmax;
         if (i_actual >= STARTINGY(Mcorr)  && j_actual >= STARTINGX(Mcorr) &&
             i_actual <= FINISHINGY(Mcorr) && j_actual <= FINISHINGX(Mcorr)) {
            ymax += i_actual*MAT_ELEM(Mcorr,i_actual,j_actual);
            xmax += j_actual*MAT_ELEM(Mcorr,i_actual,j_actual);
            sumcorr += MAT_ELEM(Mcorr,i_actual,j_actual);
	 }
   }
   shiftX=xmax /sumcorr;
   shiftY=ymax /sumcorr;
}

/* Euclidian distance ------------------------------------------------------ */
template <class T>
   double euclidian_distance(const matrix1D<T> &x, const matrix1D<T> &y) {
   SPEED_UP_temps;
   double retval=0;
   long n=0;
   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x,y) {
     retval+=(VEC_ELEM(x,i)-VEC_ELEM(y,i))*(VEC_ELEM(x,i)-VEC_ELEM(y,i));
     n++;
   }
   if (n!=0) return sqrt(retval); else return 0;
    
}

/* Euclidian distance ------------------------------------------------------ */
template <class T>
   double euclidian_distance(const matrix2D<T> &x, const matrix2D<T> &y,
      const matrix2D<int> *mask) {
   SPEED_UP_temps;
   double retval=0;
   long n=0;

   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x,y) {
     if (mask!=NULL)
       if (!(*mask)(i,j)) continue;
         retval+=(MAT_ELEM(x,i,j)-MAT_ELEM(y,i,j))*(MAT_ELEM(x,i,j)-MAT_ELEM(y,i,j));
         n++;
   }
   if (n!=0) return sqrt(retval); else return 0;
}

template <class T>
   double euclidian_distance(const matrix3D<T> &x, const matrix3D<T> &y,
      const matrix3D<int> *mask) {
   SPEED_UP_temps;
   double retval=0;
   long n=0;

   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x,y) {
     if (mask!=NULL)
       if (!(*mask)(k,i,j)) continue;
         retval+=(VOL_ELEM(x,k,i,j)-VOL_ELEM(y,k,i,j))*
                 (VOL_ELEM(x,k,i,j)-VOL_ELEM(y,k,i,j));
         n++;
   }
   if (n!=0) return sqrt(retval); else return 0;
}

/* Mutual information --------------------------------------------------- */
template <class T>
   double mutual_information(const matrix1D<T> &x, const matrix1D<T> &y,
      int nx, int ny) {
   SPEED_UP_temps;

   long n=0;
   histogram1D       histx,histy;
   histogram2D       histxy;
   matrix1D<T>       aux_x,aux_y;
   matrix1D<double>  mx,my;
   matrix2D<double>  mxy;
   int xdim,ydim;
   double MI=0.;
   double HAB=0.;
   double retval=0;

   aux_x.resize(x);
   aux_y.resize(y);

   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x,y) {
     aux_x(n)= VEC_ELEM(x,i);
     aux_y(n)= VEC_ELEM(y,i);
     n++;
   }
   aux_x.resize(n);
   aux_y.resize(n);

   if (n!=0) {

     if (nx==0) nx=(int)((log((double)n)/LOG2)+1); //Assume Gaussian distribution
     if (ny==0) ny=(int)((log((double)n)/LOG2)+1); //Assume Gaussian distribution
     compute_hist(aux_x,histx,nx);
     compute_hist(aux_y,histy,ny);
     compute_hist(aux_x,aux_y,histxy,nx,ny);

     mx=histx;
     my=histy;
     mxy=histxy;
     for(int i=0;i<nx;i++) {
       double histxi=(histx(i))/n;
       for(int j=0;j<ny;j++) {
	 double histyj=(histy(j))/n;
	 double histxyij=(histxy(i,j))/n;
	 if (histxyij>0) retval +=histxyij*log(histxyij/(histxi*histyj))/LOG2;
       }
     }
     return retval;
   } else return 0;
    
}

/* Mutual information ------------------------------------------------------ */
template <class T>
   double mutual_information(const matrix2D<T> &x, const matrix2D<T> &y,
      int nx, int ny, const matrix2D<int> *mask) {
   SPEED_UP_temps;
   long n=0;
   histogram1D       histx,histy;
   histogram2D       histxy;
   matrix1D<T>       aux_x,aux_y;
   matrix1D<double>  mx,my;
   matrix2D<double>  mxy;
   int xdim,ydim;
   double retval=0.;

   x.get_dim(xdim,ydim);
   aux_x.resize(xdim*ydim);
   y.get_dim(xdim,ydim);
   aux_y.resize(xdim*ydim);

   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x,y) {
     if (mask!=NULL)
       if (!(*mask)(i,j)) continue;
	 aux_x(n)= MAT_ELEM(x,i,j);
	 aux_y(n)= MAT_ELEM(y,i,j);
         n++;
   }

   aux_x.resize(n);
   aux_y.resize(n);

   if (n!=0) {

     if (nx==0) nx=(int)((log((double)n)/LOG2)+1); //Assume Gaussian distribution
     if (ny==0) ny=(int)((log((double)n)/LOG2)+1); //Assume Gaussian distribution
     compute_hist(aux_x,histx,nx);
     compute_hist(aux_y,histy,ny);
     compute_hist(aux_x,aux_y,histxy,nx,ny);

     mx=histx;
     my=histy;
     mxy=histxy;
     for(int i=0;i<nx;i++) {
       double histxi=(histx(i))/n;
       for(int j=0;j<ny;j++) {
	 double histyj=(histy(j))/n;
	 double histxyij=(histxy(i,j))/n;
	 if (histxyij>0) retval +=histxyij*log(histxyij/(histxi*histyj))/LOG2;
       }
     }
     return retval;
   } else return 0;
}

template <class T>
   double mutual_information(const matrix3D<T> &x, const matrix3D<T> &y,
      int nx, int ny, const matrix3D<int> *mask) {

   SPEED_UP_temps;
   long n=0;
   histogram1D       histx,histy;
   histogram2D       histxy;
   matrix1D<T>       aux_x,aux_y;
   matrix1D<double>  mx,my;
   matrix2D<double>  mxy;
   int xdim,ydim,zdim;
   double retval=0.;

   x.get_dim(ydim,xdim,zdim);
   aux_x.resize(xdim*ydim*zdim);
   y.get_dim(ydim,xdim,zdim);
   aux_y.resize(xdim*ydim*zdim);

   FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x,y) {
     if (mask!=NULL)
       if (!(*mask)(k,i,j)) continue;
	 aux_x(n)= VOL_ELEM(x,k,i,j);
	 aux_y(n)= VOL_ELEM(y,k,i,j);
         n++;
   }
   aux_x.resize(n);
   aux_y.resize(n);

   if (n!=0) {

     if (nx==0) nx=(int)((log((double)n)/LOG2)+1); //Assume Gaussian distribution
     if (ny==0) ny=(int)((log((double)n)/LOG2)+1); //Assume Gaussian distribution
     compute_hist(aux_x,histx,nx);
     compute_hist(aux_y,histy,ny);
     compute_hist(aux_x,aux_y,histxy,nx,ny);

     mx=histx;
     my=histy;
     mxy=histxy;
     for(int i=0;i<nx;i++) {
       double histxi=(histx(i))/n;
       for(int j=0;j<ny;j++) {
	 double histyj=(histy(j))/n;
	 double histxyij=(histxy(i,j))/n;
	 if (histxyij>0) retval +=histxyij*log(histxyij/(histxi*histyj))/LOG2;
       }
     }
     return retval;
   } else return 0;
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

/* Shah energy ------------------------------------------------------------- */
/* This function computes the current functional energy */
double Shah_energy(const matrix2D<double> &img, 
   const matrix2D<double> &surface_strength,
   const matrix2D<double> &edge_strength,
   double K, const matrix1D<double> &W) {
   int Ydim1 = YSIZE(img) - 1;
   int Xdim1 = XSIZE(img) - 1;

   double Kinv = 1.0/K;
   
   /* Calculate surface energy */
   double E1 = 0.0, E2 = 0.0, E3 = 0.0, E4 = 0.0;
   for (int i=1; i<Ydim1; i++)
      for (int j=1; j<Xdim1; j++) {
	    /* Calculate data matching terms */
	    double D = dMij(img,i,j);
	    double F = dMij(surface_strength,i,j);
	    double S = dMij(edge_strength,i,j);
	    E1 += W(0) * (D - F) * (D - F);
	    E3 += W(2) * K * S * S;

	    /* Calculate first derivative terms */
	    double Fx = (dMij(surface_strength,i,j+1) - dMij(surface_strength,i,j-1)) / 2;
	    double Fy = (dMij(surface_strength,i+1,j) - dMij(surface_strength,i-1,j)) / 2;
	    double Sx = (dMij(edge_strength,i,j+1)    - dMij(edge_strength,i,j-1)) / 2;
	    double Sy = (dMij(edge_strength,i+1,j)    - dMij(edge_strength,i-1,j)) / 2;
	    E2 += W(1) * (1 - S) * (1 - S) * (Fx * Fx + Fy * Fy);
	    E4 += W(3) * Kinv * (Sx * Sx + Sy * Sy);
      }

   return E1 + E2 + E3 + E4; // Total energy
}

/* Update Surface Shah ----------------------------------------------------- */
/* This routine performs one update to the edge estimate based     
   on a finite differences solution to the following equation:     
       0 = dE/df = dF/df - d(dF/dfx)/dx - d(dF/dfy)/dy                 
           + dd(dF/dfxx)/dxx + dd(dF/dfxy)/dxy + dd(dF/dfyy)/dyy */
double Update_surface_Shah(matrix2D<double> &img, 
   matrix2D<double> &surface_strength, matrix2D<double> &edge_strength,
   const matrix1D<double> &W) {
   double Diff = 0.0, Norm=0.0;

   int Ydim1 = YSIZE(img) - 1;
   int Xdim1 = XSIZE(img) - 1;

   /* Update surface estimate */
   for (int i=1; i<Ydim1; i++)
      for (int j=1; j<Xdim1; j++) {
	 /* Calculate edge partial derivative terms */
	 double S  =  dMij(edge_strength,i,j);
	 double Sx = (dMij(edge_strength,i,j+1)    - dMij(edge_strength,i,j-1)) / 2;
	 double Sy = (dMij(edge_strength,i+1,j)    - dMij(edge_strength,i-1,j)) / 2;

	 double nS  = 1-S;
	 double nS2 = nS*nS;

	 /* Calculate surface partial derivative terms (excluding central pixel) */
	 double F,D; F=D=dMij(img,i,j);
	 double Fx = (dMij(surface_strength,i,j+1) - dMij(surface_strength,i,j-1)) / 2;
	 double Fy = (dMij(surface_strength,i+1,j) - dMij(surface_strength,i-1,j)) / 2;
	 double Fxx=  dMij(surface_strength,i,j+1) + dMij(surface_strength,i,j-1);
	 double Fyy=  dMij(surface_strength,i+1,j) + dMij(surface_strength,i-1,j);

	 /* Calculate surface partial derivative weights */
	 double wFx = 4 * W(1) * nS * Sx;
	 double wFy = 4 * W(1) * nS * Sy;
	 double wFxx = -2 * W(1) * nS2;
	 double wFyy = -2 * W(1) * nS2;

	 /* Calculate new surface value */
	 double Constant = -2 * W(0) * D;
	 double Central  = -2 * W(0) + 2 * wFxx + 2 * wFyy;
	 double Neighbors = wFx * Fx + wFy * Fy + wFxx * Fxx + wFyy * Fyy;

	 if (ABS(Central)>XMIPP_EQUAL_ACCURACY)
	    F=(Constant + Neighbors) / Central;
	 F=CLIP(F,0,1);

	 // Compute the difference. 
	 Diff += ABS(dMij(surface_strength,i,j) - F);
	 Norm += ABS(dMij(surface_strength,i,j));

         // Update the new value.
	 dMij(surface_strength,i,j) = F;
      }
   return Diff/Norm; // Return the relative difference.
}

/* Update Edge Shah -------------------------------------------------------- */
/* This routine performs one update to the edge estimate based	 
   on a finite differences solution to the following equation:
   0 = dE/ds = dF/ds - d(dF/dsx)/dx - d(dF/dsy)/dy */
double Update_edge_Shah(matrix2D<double> &img, 
   matrix2D<double> &surface_strength, matrix2D<double> &edge_strength,
   double K, const matrix1D<double> &W) {
   double Diff = 0.0, Norm=0.0;

   int Ydim1 = YSIZE(img) - 1;
   int Xdim1 = XSIZE(img) - 1;
   double Kinv = 1.0 / K;

   /* Update edge estimate */
   for (int i=1; i<Ydim1; i++)
      for (int j=1; j<Xdim1; j++) {
	 /* Calculate first and second derivative terms */
	 double Fx = (dMij(surface_strength,i,j+1) - dMij(surface_strength,i,j-1)) / 2;
	 double Fy = (dMij(surface_strength,i+1,j) - dMij(surface_strength,i-1,j)) / 2;
	 double Constant = W(1) * (Fx * Fx + Fy * Fy);

	 /* Calculate weights for central pixel and neighbors */
	 double Central   = W(2) * K + W(3) * Kinv * 4;
	 double Neighbors = W(3) * Kinv * (
	     dMij(edge_strength,i-1,j)+dMij(edge_strength,i+1,j)
	    +dMij(edge_strength,i,j-1)+dMij(edge_strength,i,j+1));

	 /* Calculate new S value */
	 double Old_edge_strength = dMij(edge_strength,i,j);
	 double S = (Constant + Neighbors) / (Constant + Central);
	 if	 (S < 0) dMij(edge_strength,i,j)/=2;
	 else if (S > 1) dMij(edge_strength,i,j)=(dMij(edge_strength,i,j)+1)/2;
	 else		 dMij(edge_strength,i,j)=S;

	 // Compute the difference. 
	 Diff += ABS(dMij(edge_strength,i,j)-Old_edge_strength);
	 Norm += ABS(Old_edge_strength);
      }
   return Diff/Norm; // Return the relative difference.
}

/* Smoothing Shah ---------------------------------------------------------- */
#define SHAH_CONVERGENCE_THRESHOLD  0.0001
void Smoothing_Shah(matrix2D<double> &img, 
   matrix2D<double> &surface_strength, matrix2D<double> &edge_strength,
   const matrix1D<double> &W, int OuterLoops, int InnerLoops,
   int RefinementLoops, bool adjust_range) {

   type_cast(img,surface_strength);
   if (adjust_range) surface_strength.range_adjust(0,1);
   edge_strength.resize(img);

   for(int k=1; k<=RefinementLoops; k++) {
       // Initialize Edge Image.
       edge_strength.init_zeros();

       double diffsurface = MAXFLOAT; // Reset surface difference
       for (int i=0; ((i<OuterLoops) && OuterLoops) ||
               ((diffsurface>SHAH_CONVERGENCE_THRESHOLD) && !OuterLoops); i++) {

      	  /* cout << "Iteration ..." << i+1;*/
          /* Iteratively update surface estimate */
          for (int j=0; j<InnerLoops; j++)
             diffsurface=
	        Update_surface_Shah(img,surface_strength,edge_strength,W);
	     
          /* Iteratively update edge estimate */
	  double diffedge;
          for (int j=0; j<InnerLoops; j++)
             diffedge=
	        Update_edge_Shah(img,surface_strength,edge_strength,k,W);
	     	     
          /* Calculate new functional energy */
          double energy=Shah_energy(img,surface_strength,edge_strength,k,W);
	  /* cout << " Energy " << energy
	       << " ... Relative Diff " << diffsurface
	       << " ... Edge diff " << diffedge
	       << endl; */
       }
    }
}

/* Rotational invariant moments -------------------------------------------- */
void rotational_invariant_moments(const matrix2D<double> &img,
   const matrix2D<int> *mask,
   matrix1D<double> &v_out) {
   // Prepare some variables
   double m_11=0, m_02=0, m_20=0, m_12=0, m_21=0, m_03=0, m_30=0; //, m_00=0;
   double normalize_x=2.0/XSIZE(img);
   double normalize_y=2.0/YSIZE(img);

   // Compute low-level moments
   FOR_ALL_ELEMENTS_IN_MATRIX2D(img) {
      if (mask!=NULL)
         if (!(*mask)(i,j)) continue;
      double val=img(i,j);
      double dx=j*normalize_x;
      double dy=i*normalize_y;
      double dx2=dx*dx, dx3=dx2*dx;
      double dy2=dy*dy, dy3=dy2*dy;
      // m_00+=val;
      m_11+=val*dx*dy;  m_02+=val*dy2;    m_20+=val*dx2;
      m_12+=val*dx*dy2; m_21+=val*dx2*dy; m_03+=val*dy3; m_30+=val*dx3;
   }
   //m_11/=m_00; m_02/=m_00; m_20/=m_00;
   //m_12/=m_00; m_21/=m_00; m_03/=m_00; m_30/=m_00;
   
   // Compute high-level rotational invariant moments
   v_out.resize(5);
   v_out( 0)= m_20+m_02;
   v_out( 1)=(m_20-m_02)*(m_20-m_02)+4*m_11*m_11;
   v_out( 2)=(m_30-3*m_12)*(m_30-3*m_12)+
             (3*m_21-m_03)*(3*m_21-m_03);
   v_out( 3)=(m_30+m_12)*(m_30+m_12)+
             (m_21+m_03)*(m_21+m_03);
   v_out( 4)=(m_30-3*m_12)*(m_30+m_12)*(
                   (m_30+m_12)*(m_30+m_12)
	        -3*(m_21+m_03)*(m_21+m_03))+
	     (3*m_21-m_03)*(m_21+m_03)*(
	         3*(m_30+m_12)*(m_30+m_12)
		  -(m_21+m_03)*(m_21+m_03));
/*
   v_out( 5)=(m_20+m_02)*(
             	  (m_30+m_12)*(m_30+m_12)
	       -3*(m_21+m_03)*(m_21+m_03))
	     +4*m_11*(m_30+m_12)*(m_03+m_21);
   v_out( 6)=(3*m_21-m_03)*(m_12+m_30)*(
             	    (m_30+m_12)*(m_30+m_12)
	     	 -3*(m_21+m_03)*(m_21+m_03))
	     -(m_30-3*m_12)*(m_12+m_03)*(
	     	  3*(m_30+m_12)*(m_30+m_12)
		   -(m_21+m_03)*(m_21+m_03));
*/
}   

/* Inertia moments --------------------------------------------------------- */
void inertia_moments(const matrix2D<double> &img,
   const matrix2D<int> *mask,
   matrix1D<double> &v_out) {
   // Prepare some variables
   double m_11=0, m_02=0, m_20=0;
   double normalize_x=2.0/XSIZE(img);
   double normalize_y=2.0/YSIZE(img);

   // Compute low-level moments
   FOR_ALL_ELEMENTS_IN_MATRIX2D(img) {
      if (mask!=NULL)
         if (!(*mask)(i,j)) continue;
      double val=img(i,j);
      double dx=j*normalize_x;
      double dy=i*normalize_y;
      double dx2=dx*dx;
      double dy2=dy*dy;
      m_11+=val*dx*dy;  m_02+=val*dy2;    m_20+=val*dx2;
   }
   
   // Compute the eigen values of the inertia matrix
   // [m_02 m_11
   //  m_11 m_20]
   matrix2D<double> A(2,2);
   A(0,0)=m_02;
   A(0,1)=A(1,0)=m_11;
   A(1,1)=m_20;
   matrix2D<double> u,v;
   svdcmp(A,u,v_out,v);
   v_out=v_out.sort();
}   

/* Fill triangle ----------------------------------------------------------- */
void fill_triangle(matrix2D<double> &img, int *tx, int *ty, double color) {
    /*
     * Order in y values
     */
    int y1 = 0;
    while (!y1) {
       y1 = 1;
       for (int y2 = 0; y2 < 2; y2++) {
          if (ty[y2] > ty[y2 + 1] ||
             (ty[y2] == ty[y2 + 1] && tx[y2] < tx[y2 + 1])) {
	     int x1 = ty[y2];
	     ty[y2] = ty[y2 + 1];
	     ty[y2 + 1] = x1;
	     x1 = tx[y2];
	     tx[y2] = tx[y2 + 1];
	     tx[y2 + 1] = x1;
	     y1 = 0;
          }
       }
    }

    int dx1 = tx[1] - tx[0];
    int dx2 = tx[2] - tx[0];
    int dy1 = ty[1] - ty[0];
    int dy2 = ty[2] - ty[0];

    int sx1 = SGN0(dx1);
    int sx2 = SGN0(dx2);
    int sy1 = SGN0(dy1);

    int ix1 = ABS(dx1);
    int ix2 = ABS(dx2);
    int iy1 = ABS(dy1);
    int iy2 = ABS(dy2);

    int inc1 = MAX(ix1, iy1);
    int inc2 = MAX(ix2, iy2);

    int x1, x2, y2, xl, xr;
    x1 = x2 = y1 = y2 = 0;
    xl = xr = tx[0];
    int y = ty[0];

    while (y != ty[1]) {
       int doit1 = 0;
       int doit2 = 0;

       while (!doit1) {   /* Wait until y changes */
	   x1 += ix1;
	   y1 += iy1;
	   if (x1 > inc1) {
	      x1 -= inc1;
	      xl += sx1;
	   }
	   if (y1 > inc1) {
	      y1 -= inc1;
	      y += sy1;
	      doit1 = 1;
	   }
       }

       while (!doit2) {   /* Wait until y changes */
	   x2 += ix2;
	   y2 += iy2;
	   if (x2 > inc2) {
	      x2 -= inc2;
	      xr += sx2;
	   }
	   if (y2 > inc2) {
	      y2 -= inc2;
	      doit2 = 1;
	   }
       }

       for (int myx=xl; myx<=xr; myx++)
          img(y,myx)=color;
    }

    dx1 = tx[2] - tx[1];
    dy1 = ty[2] - ty[1];

    sx1 = SGN0(dx1);
    sy1 = SGN0(dy1);

    ix1 = ABS(dx1);
    iy1 = ABS(dy1);

    inc1 = MAX(ix1, iy1);
    xl = tx[1];
    x1 = 0;

    while (y != ty[2]) {
       int doit1 = 0;
       int doit2 = 0;

       while (!doit1) {   /* Wait until y changes */
	   x1 += ix1;
	   y1 += iy1;
	   if (x1 > inc1) {
	      x1 -= inc1;
	      xl += sx1;
	   }
	   if (y1 > inc1) {
	      y1 -= inc1;
	      y += sy1;
	      doit1 = 1;
	   }
       }

       while (!doit2) {   /* Wait until y changes */
	   x2 += ix2;
	   y2 += iy2;
	   if (x2 > inc2) {
	      x2 -= inc2;
	      xr += sx2;
	   }
	   if (y2 > inc2) {
	      y2 -= inc2;
	      doit2 = 1;
	   }
       }

       for (int myx=xl; myx<=xr; myx++)
          img(y,myx)=color;
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

   euclidian_distance(x1,y1);
   euclidian_distance(x2,y2);
   euclidian_distance(x3,y3);
   euclidian_distance(x2,y2,&mask2);
   euclidian_distance(x3,y3,&mask3);
   
   mutual_information(x1,y1,position,position);
   mutual_information(x2,y2,position,position);
   mutual_information(x3,y3,position,position);
   mutual_information(x2,y2,position,position,&mask2);
   mutual_information(x3,y3,position,position,&mask3);
   
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
    float  f; instantiate_filtersT(f);
}

