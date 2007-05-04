/***************************************************************************
 *
 * Authors:    Roberto Marabini
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
#include "sampling.h"

/* Default Constructor */
XmippSampling::XmippSampling()
{  

   vertices aux;
   aux.rot =   -PI/5.; aux.tilt = PI/2.-cte_w  ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot =    PI/5.; aux.tilt = PI/2.-cte_w  ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot = 3.*PI/5.; aux.tilt = PI/2.-cte_w  ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot = 5.*PI/5.; aux.tilt = PI/2.-cte_w  ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot =-3.*PI/5.; aux.tilt = PI/2.-cte_w  ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot =    0./5.; aux.tilt = -cte_w+PI/2. ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot = 2.*PI/5.; aux.tilt = -cte_w+PI/2. ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot = 4.*PI/5.; aux.tilt = -cte_w+PI/2. ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot =-4.*PI/5.; aux.tilt = -cte_w+PI/2. ; aux.psi =  0.; vertices_angles.push_back(aux);
   aux.rot =-2.*PI/5.; aux.tilt = -cte_w+PI/2. ; aux.psi =  0.; vertices_angles.push_back(aux);

   vertices_vectors.push_back(vector_R3(0.,0.,1.));
   vertices_vectors.push_back(vector_R3(0.723606900230461,-0.525731185781806,0.447213343087301));
   vertices_vectors.push_back(vector_R3(0.723606900230461,0.525731185781806,0.447213343087301));
   vertices_vectors.push_back(vector_R3(-0.276393239417711,0.850650928976665,0.447213343087301));
   vertices_vectors.push_back(vector_R3(-0.8944273172062,0.,0.447213343087301));
   vertices_vectors.push_back(vector_R3(-0.276393239417711,-0.850650928976665,0.447213343087301));
   vertices_vectors.push_back(vector_R3(0.8944273172062,0.,-0.447213343087301));
   vertices_vectors.push_back(vector_R3(0.276393242471372,0.850650927984471,-0.447213343087301));
   vertices_vectors.push_back(vector_R3(-0.723606898343194,0.525731188379405,-0.447213343087301));
   vertices_vectors.push_back(vector_R3(-0.723606898343194,-0.525731188379405,-0.447213343087301));
   vertices_vectors.push_back(vector_R3(0.276393242471372,-0.850650927984471,-0.447213343087301));
   vertices_vectors.push_back(vector_R3(0.,0.,-1.));
   //#define DEBUG1
   #ifdef  DEBUG1  
    for(    int i = 0; 
            i < vertices_vectors.size(); 
            i++)
       cout  <<  vertices_vectors[].transpose()  << endl;
   #endif
   #undef DEBUG1
}

void XmippSampling::SetSampling(double sampling)
{
   sampling_rate_rad = DEG2RAD(sampling);
   number_of_samples = ROUND(cte_w/sampling_rate_rad);
}

/* Compute edge sampling points using Baumgardner  1995 */

void XmippSampling::Compute_sampling_points(bool only_half_sphere)
 {
    /** vector to decimate the triangles */  
   vector <matrix1D<double> >edge_vector_start;
   /** vector to decimate the triangles */  
   vector <matrix1D<double> >edge_vector_end;
  // I need 10 auxiliary vector for edges
   matrix1D<double> starting_point, ending_point;
   //01a   
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[1];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[6]; ending_point=vertices_vectors[1];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //01b
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[2];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[6]; ending_point=vertices_vectors[2];
   fill_edge(starting_point,ending_point,edge_vector_end,true);
   //02a   
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[2];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[7]; ending_point=vertices_vectors[2];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //02b
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[3];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[7]; ending_point=vertices_vectors[3];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //03a   
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[3];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[8]; ending_point=vertices_vectors[3];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //03b
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[4];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[8]; ending_point=vertices_vectors[4];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //04a   
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[4];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[9]; ending_point=vertices_vectors[4];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //04b
   starting_point=vertices_vectors[0]; ending_point=vertices_vectors[5];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[9]; ending_point=vertices_vectors[5];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //05a   
   starting_point=vertices_vectors[0];  ending_point=vertices_vectors[5];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[10]; ending_point=vertices_vectors[5];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //05b
   starting_point=vertices_vectors[0];  ending_point=vertices_vectors[1];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[10]; ending_point=vertices_vectors[1];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //06a   
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[10];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[5];  ending_point=vertices_vectors[10];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //06b
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[9];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[5];  ending_point=vertices_vectors[9];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //07a   
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[9];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[4];  ending_point=vertices_vectors[9];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //07b
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[8];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[4];  ending_point=vertices_vectors[8];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //08a   
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[8];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[3];  ending_point=vertices_vectors[8];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //08b
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[7];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[3];  ending_point=vertices_vectors[7];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //09a   
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[7];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[2];  ending_point=vertices_vectors[7];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //09b
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[6];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[2];  ending_point=vertices_vectors[6];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //10a   
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[6];
   fill_edge(starting_point,ending_point,edge_vector_start,false);
   starting_point=vertices_vectors[1]; ending_point=vertices_vectors[6];
   fill_edge(starting_point,ending_point,edge_vector_start,true);
   //10b
   starting_point=vertices_vectors[11]; ending_point=vertices_vectors[10];
   fill_edge(starting_point,ending_point,edge_vector_end,false);
   starting_point=vertices_vectors[1];  ending_point=vertices_vectors[10];
   fill_edge(starting_point,ending_point,edge_vector_end,true);

   //#define DEBUG2
   #ifdef  DEBUG2
    for(    int i = 0; 
            i < edge_vector_start.size(); 
            i++)
       {     
       cout  <<  edge_vector_start[i].transpose()  << " 1 1 " << endl;
       cout  <<  edge_vector_end[i].transpose()  << " 1 2 " << endl;
       }
   //cout  <<  ending_point.transpose()    << " 1.1 1.5 " << endl;
   #endif
   #undef DEBUG2
   // add  main corners
       for(    int i = 0; 
            i < vertices_vectors.size(); 
            i++)
           { 
           if (only_half_sphere && ZZ(vertices_vectors[i])<0.0)
              continue;
           else   
	      sampling_points_vector.push_back(vertices_vectors[i]);
           }   
   // add edges    
       for(    int i = 0; 
            i < edge_vector_start.size(); 
            i++)
            {
            if(i<number_of_samples * 10 -15)
               { 
               if (only_half_sphere && ZZ(edge_vector_start[i])<0.0)
                  continue;
               else   
	          sampling_points_vector.push_back(edge_vector_start[i]);
               }   
            else
               { 
               if (only_half_sphere && ZZ(edge_vector_end[i])<0.0)
                  continue;
               else   
	          sampling_points_vector.push_back(edge_vector_end[i]);
               }   
            }
    // add in between points
    int j=0;
    bool j_flag=false;
    for(    int i = 0; 
            i < edge_vector_start.size(); 
           i++)
       {     
           if((j % (number_of_samples-1))==0 && j!=0) {j=0; j_flag=true;}
           if((j % (number_of_samples-2))==0 && j!=0  && j_flag==true) 
               {j=0; j_flag=false;}
           fill_distance(edge_vector_start[i],
                          edge_vector_end[i],
                          sampling_points_vector,
                          (j+1)%number_of_samples,
                          only_half_sphere);
            j++;              
       }
   //#define DEBUG3
   #ifdef  DEBUG3
    for(    int i = 0; 
            i < sampling_points_vector.size(); 
            i++)
       {     
       cout  <<  sampling_points_vector[i].transpose()  << " 1 1 " << endl;
       }
   #endif
   #undef DEBUG3
   matrix1D<double> aux(3), aux1(3);
   ZZ(aux)=0.;
   double rot, tilt,psi;
    for(    int i = 0; 
            i < sampling_points_vector.size(); 
            i++)
         {
         XX(aux)=atan2(YY(sampling_points_vector[i]),
                       XX(sampling_points_vector[i]));
         YY(aux)=acos(ZZ(sampling_points_vector[i])); 
         if (YY(aux)<0.) YY(aux) += PI;            
         sampling_points_angles.push_back(aux*180./PI);  
         }
#ifdef NEVERDEFINED
  //sort points
  int k;
  int bound = sampling_points_angles.size()-1;
  int t;
  int last_swap;
  matrix1D<double> aux_angle(3);
  matrix1D<double> aux_vector(3);

  while (bound) {
    last_swap = 0;
    for ( k=0; k<bound; k++ ){
       aux = sampling_points_angles[k]; /* aux is a maximum of A[0]..A[k] */
       aux1 = sampling_points_vector[k]; /* aux is a maximum of A[0]..A[k] */
       if ( sort_func(aux, sampling_points_angles[k+1])) {
         sampling_points_angles[k] = sampling_points_angles[k+1];
         sampling_points_angles[k] = sampling_points_angles[k+1];
         sampling_points_angles[k+1] = aux; /*swap*/
         sampling_points_vector[k] = sampling_points_vector[k+1];
         sampling_points_vector[k] = sampling_points_vector[k+1];
         sampling_points_vector[k+1] = aux1; /*swap*/
         last_swap = k; /* mark the last swap position */
       }//if
    }//for
    bound=last_swap; /*elements after bound already sorted */       
  }//while       
#endif
   //#define DEBUG3
   #ifdef  DEBUG3
    for(    int i = 0; 
            i < sampling_points_angles.size(); 
            i++)
       {     
       cout  <<  sampling_points_angles[i].transpose()  << " 1 1 " << endl;
       }
   #endif
   #undef DEBUG3

 }  
 
 // return 1 if a should go first 0 is equal -1 if before
int XmippSampling::sort_func(matrix1D<double> &t, matrix1D<double> &a)
 {
 if (YY(t) - 0.000001> YY(a)) { return 1;}
 else if (YY(t) + 0.000001 < YY(a)) { return 0;};
 if (XX(t)  > XX(a)) { return 1;}
 else  { return 0;};
 }
    
void XmippSampling::fill_edge(matrix1D<double> starting_point,
                                      matrix1D<double> ending_point,
                                      vector <matrix1D<double> > & edge_vector,
                                      bool END_FLAG
                                      )
   {
    matrix1D<double> v_aux(3);
   
    double alpha;
    double beta;
    double gamma;
   // skip first corener, already computed;
   double upsilon = acos(dot_product(starting_point,ending_point));
   for (int i1=1; i1< number_of_samples; i1++)
      { 
	gamma  = (double)i1/(number_of_samples-1);
        alpha  = sin((1.-gamma)*upsilon)/(sin(upsilon));
        beta   = sin(gamma*upsilon)/sin(upsilon);
	v_aux = alpha*starting_point + beta*ending_point;
	v_aux=v_aux.normalize();
        if(beta>0.9999 && END_FLAG) continue;
	edge_vector.push_back(v_aux);
      }
 }
void XmippSampling::fill_distance(matrix1D<double> starting_point,
                                          matrix1D<double> ending_point,
                                          vector <matrix1D<double> > &
                                          sampling_points_vector,
                                          int my_number_of_samples,
                                          bool only_half_sphere
                                          )
   {
    matrix1D<double> v_aux(3);
   
    double alpha;
    double beta;
    double gamma;
   // skip first corener, already computed;
   double upsilon = acos(dot_product(starting_point,ending_point));
   for (int i1=1; i1< my_number_of_samples; i1++)
      { 
	gamma  = (double)i1/ (my_number_of_samples);
        alpha  = sin((1.-gamma)*upsilon)/(sin(upsilon));
        beta   = sin(gamma*upsilon)/sin(upsilon);
	v_aux = alpha*starting_point + beta*ending_point;
	v_aux=v_aux.normalize();
        if (only_half_sphere && ZZ(v_aux)<0.0)
           continue;
        else   
	   sampling_points_vector.push_back(v_aux);
      }
 }

void XmippSampling::remove_redundant_points(string symmetry,
                                            int sym_order )
{
  matrix2D<double>  L(4,4), R(4,4);
  matrix2D<double>  aux(3,3);
  matrix1D<double>  row1(3),row2(3);
  bool valid=true;
  double rot, tilt, psi=0;
  double rotp,tiltp, psip =0.;
  no_redundant_sampling_points_vector.clear();
  double aux1,aux2;
  bool match=false;
 
  //int j_end=0;
  matrix1D<double>  row(3);


  no_redundant_sampling_points_vector.clear();     
  no_redundant_sampling_points_angles.clear();     
  double my_dot_product;
  if (symmetry.compare("cn")==0)
     {//OK
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= (-180./sym_order) &&
               XX(sampling_points_angles[i])<=  (180./sym_order) )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("ci")==0 ||
           symmetry.compare("cs")==0 )
     {//OK
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (YY(sampling_points_angles[i])<=90)
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("cnv")==0)
     {//OK
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>=    0./sym_order &&
               XX(sampling_points_angles[i])<=  180./sym_order )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("cnh")==0)
     {//OK
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= -180./sym_order &&
               XX(sampling_points_angles[i])<=  180./sym_order &&
               YY(sampling_points_angles[i])<=    90.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("sn")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= -180./(2*sym_order) &&
               XX(sampling_points_angles[i])<=  180./(2*sym_order) &&
               YY(sampling_points_angles[i])<=    90.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("dn")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= -180./(sym_order) &&
               XX(sampling_points_angles[i])<=  180./(sym_order) &&
               YY(sampling_points_angles[i])<=    90.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("dnv")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>=    0. &&
               XX(sampling_points_angles[i])<=  180./(sym_order) &&
               YY(sampling_points_angles[i])<=    90.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("dnh")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= 0. &&
               XX(sampling_points_angles[i])<=  180./(sym_order) &&
               YY(sampling_points_angles[i])<=   90.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("t")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= -180./(3) &&
               XX(sampling_points_angles[i])<=  180./(3) &&
               YY(sampling_points_angles[i])>=    0.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("td")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= 0. &&
               XX(sampling_points_angles[i])<=  180./(3) &&
               YY(sampling_points_angles[i])>=    0.
               )
              {
              no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
              no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
              }
        }// for i  
     }
  else if (symmetry.compare("th")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
        int flag = false;
           if (XX(sampling_points_angles[i])>= 0. &&
               XX(sampling_points_angles[i])<=  180./(2) 
               )
             {
              if(XX(sampling_points_angles[i])>= 0 && 
                 XX(sampling_points_angles[i])<= 45 &&
                 YY(sampling_points_angles[i])<= XX(sampling_points_angles[i]) 
                 )  
                 flag=true;
              else if(XX(sampling_points_angles[i])> 45 && 
                 XX(sampling_points_angles[i])<= 90 &&
                 YY(sampling_points_angles[i])<= 90.-XX(sampling_points_angles[i]) 
                 )  
                 flag=true;
               
             if (flag)
                {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
             }
        }// for i  
     }
  else if (symmetry.compare("o")==0)
     {
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
        int flag = false;
           if (XX(sampling_points_angles[i])>= 0. &&
               XX(sampling_points_angles[i])<=  180./(2) &&
               YY(sampling_points_angles[i])>=45
               )
                {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
        }// for i  
     }
  else if (symmetry.compare("oh")==0)
     {
        /** coordinates of the icosahedron P5 axis
            +-31.7174745559,90 
            coordinates of the icosahedron P3 axis
            0,69.094843368*/
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
        int flag = false;
           if (XX(sampling_points_angles[i])>= 90./2. &&
               XX(sampling_points_angles[i])<=  180./(2.) &&
               YY(sampling_points_angles[i])>=45
               )
                {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
         }// for i  
     }
  else if (symmetry.compare("i")==0)
     {//?
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])>= -31.7174745559 &&
               XX(sampling_points_angles[i])<=  0 &&
               YY(sampling_points_angles[i])<=90 &&
               YY(sampling_points_angles[i])>= (69.094843368 - 
                                                XX(sampling_points_angles[i])*
                                                (90.-69.094843368)/31.7174745559
                                               )
               )
                {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
            else if (   XX(sampling_points_angles[i]) <= 31.7174745559 &&
                        XX(sampling_points_angles[i]) >   0 &&
                        YY(sampling_points_angles[i]) <= 90 &&
                        YY(sampling_points_angles[i]) >= (69.094843368 -RAD2DEG(sampling_rate_rad) + 
                                                          XX(sampling_points_angles[i])*
                                                          (90.-69.094843368)/31.7174745559
                                                         )
               )
                {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
         }// for i  
    }
  else if (symmetry.compare("ih")==0)
      for (int i=0; i< sampling_points_angles.size(); i++)
        {
           if (XX(sampling_points_angles[i])<= 31.7174745559 &&
               XX(sampling_points_angles[i])>=  0 &&
               YY(sampling_points_angles[i])<=90 &&
               YY(sampling_points_angles[i])>= (69.094843368 + 
                                              XX(sampling_points_angles[i])*
                                              (90.-69.094843368)/31.7174745559
                                              )
               )
                {
                no_redundant_sampling_points_angles.push_back(sampling_points_angles[i]);
                no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
                }
         }// for i  
  else  
     {
     cerr << "ERROR: Symmetry " << symmetry  << "is not known" << endl;
     exit(0);
     }


#ifdef NEVERDEFINED  


/*  
  for (int j=1; j<SL.v_rotation_axis.size(); j++)
    Euler_direction2angles(SL.v_rotation_axis[j], rot, tilt,psi);  
    for (int isym=0; isym<=SL.SymsNo(); isym++) 
      {
      if(isym!=SL.SymsNo())
           {
           SL.get_matrices(isym,L,R);
           L.resize(3,3); // Erase last row and column
           R.resize(3,3); // as only the relative orientation
           Euler_apply_transf(L,R,rot,tilt,psi,rotp,tiltp,psip);
           Euler_direction(rotp, tiltp, psip,row);
           aux1=ABS(dot_product(row1,SL.v_rotation_axis[0]));
           //row = R * no_redundant_sampling_points_vector[j];
           }
      }
   }  
*/
//  cout << SL.v_rotation_axis[0];
//  cout << SL.v_ang_incr[0];
#ifdef NEVERDEFINED  
  int j_end=SL.v_rotation_axis.size();
  int i_end=sampling_points_vector.size();
  for (int i=0; i< i_end; i++)
  { 
    valid =true;
    for (int j=0; j<j_end; j++)
        {
         Euler_direction2angles(SL.v_rotation_axis[j], rot, tilt,psi);  
         Euler_angles2matrix(rot, tilt,psi,aux);
         aux.getRow(0,row1);
         row1=row1.transpose();
         Euler_angles2matrix(rot, tilt,psi+SL.v_ang_incr[j],aux);
         aux.getRow(0,row2);
         row2=row2.transpose();
         row1=vector_product(row1,sampling_points_vector[i]);
         row2=vector_product(row2,sampling_points_vector[i]);
         
         aux1=dot_product(row1,SL.v_rotation_axis[j]);
         aux2=dot_product(row2,SL.v_rotation_axis[j]);
         //cerr << i << " " << j << " " << aux1 << " " << aux2 <<endl;
         
         if( (aux1<0 && aux2>0) || (aux1==0 && aux2==0))
            continue;
         else
           valid =false;   
         }
    if(valid)     
       no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
  }      
#endif
  //#define TEST
  #ifdef TEST
  matrix1D<double>  two(3)=vector_R3(0.,0.,1.).normalize();
  matrix1D<double>  five(3)=vector_R3(1.,0.,1.618033989).normalize();
  matrix1D<double>  three(3)=vector_R3(0.,0.53934467,1.4120227).normalize();
  #endif
  no_redundant_sampling_points_vector.clear();     
  int j_end=0;
  int i_end=sampling_points_vector.size();
  double my_dot_product;
  for (int i=0; i< i_end; i++)
    {cerr << ".";
    for (int j=0; j<no_redundant_sampling_points_vector.size(); j++) 
       {
       Euler_direction2angles(no_redundant_sampling_points_vector[j], rot, tilt,psi);  
       for (int isym=0; isym<=SL.SymsNo(); isym++) 
          {
          if(isym!=SL.SymsNo())
               {
               SL.get_matrices(isym,L,R);
               L.resize(3,3); // Erase last row and column
               R.resize(3,3); // as only the relative orientation
               Euler_apply_transf(L,R,rot,tilt,psi,rotp,tiltp,psip);
               Euler_direction(rotp, tiltp, psip,row);
               //row = R * no_redundant_sampling_points_vector[j];
               }
          else
               row =   no_redundant_sampling_points_vector[j];   
//          if (ZZ(row)<0) continue;
//          cerr << i << " " << j << " " << isym << " " << 
//               dot_product(row,sampling_points_vector[i]) << " " <<
//               acos(dot_product(row,sampling_points_vector[i]))
//               << " " << sampling_rate_rad
//               << " " << row.module()
//               << endl;
           my_dot_product= dot_product(row,sampling_points_vector[i]);
           if (my_dot_product >=1.) my_dot_product=0.999999;
           
           if(acos(my_dot_product) < (sampling_rate_rad*.75))
              {
              match=true;
              break;
              }//if
           }//for isym   
      if(match) break;
       }//for j
    if (!match)
       {
       no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
       //cout << "non_redu " << no_redundant_sampling_points_vector.size()<< endl;;
       j_end++;
       }
       match=false;
    }// for i  
#ifdef NEVERDEFINED  
  for (int i=0; i< i_end; i++)
    {
    for (int j=0; j<j_end; j++) 
       {
       cerr << "j="     << j     << " " << endl;
       cerr << "j_end=" << j_end << " " << endl;
       Euler_direction2angles(no_redundant_sampling_points_vector[j], rot, tilt,psi);  
       for (int isym=0; isym<SL.SymsNo(); isym++) 
          {
          SL.get_matrices(isym,L,R);

          L.resize(3,3); // Erase last row and column
          R.resize(3,3); // as only the relative orientation
          Euler_apply_transf(L,R,rot,tilt,psi,rotp,tiltp,psip);
          Euler_direction(rotp, tiltp, psip,row);
          cerr << isym << " " ;
          cerr << 
                  row.transpose() << "   " << 
                  sampling_points_vector[i].transpose() <<
                  " " << ABS(acos(dot_product(row,sampling_points_vector[i]) ) )<<
                  " " << sampling_rate_rad*0.7 << "\n";
          cerr.flush();       

          if(ABS(acos(dot_product(row,sampling_points_vector[i]) ) ) < 
                 (sampling_rate_rad*0.7))
                 {
                 match=true;
                 break;
                 }

          }
       if(match) break;
       }    
    if (!match)
       {
       no_redundant_sampling_points_vector.push_back(sampling_points_vector[i]);
       cout << sampling_points_vector[i].transpose() << " 1 1 " << endl;
       match=false;
       j_end++;
       } 
   }
#endif
#endif

}     
/* Create symmetry file----------------------------------------------------- */
void XmippSampling::create_sym_file(string symmetry, int sym_order) {
   symmetry_file=symmetry+".sys";
   ofstream SymFile;
   SymFile.open(symmetry_file.c_str(), ios::out);
       if (symmetry.compare("cn")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          }
       else if (symmetry.compare("ci")==0)
          {
          SymFile << "invertion ";
          }
       else if (symmetry.compare("cs")==0)
          {
          SymFile << "mirror_plane 0 0 1";
          }
       else if (symmetry.compare("cnv")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          SymFile << endl;
          SymFile << "mirror_plane 1 0 0";
          }
       else if (symmetry.compare("cnh")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          SymFile << endl;
          SymFile << "mirror_plane 0 0 1";
          }
       else if (symmetry.compare("sn")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          SymFile << endl;
          SymFile << "invertion ";
          }
       else if (symmetry.compare("dn")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          SymFile << endl;
          SymFile << "rot_axis " << "2" << " 1 0 0";
          }
       else if (symmetry.compare("dnv")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          SymFile << endl;
          SymFile << "rot_axis " << "2" << " 1 0 0";
          SymFile << endl;
          SymFile << "mirror_plane 1 0 0";
          }
       else if (symmetry.compare("dnh")==0)
          {
          SymFile << "rot_axis " << sym_order << " 0 0 1";
          SymFile << endl;
          SymFile << "rot_axis " << "2" << " 1 0 0";
          SymFile << endl;
          SymFile << "mirror_plane 0 0 1";
          }
       else if (symmetry.compare("t")==0)
          {
          SymFile << "rot_axis " << "3" << "  .57735026918962576451  .57735026918962576451 .57735026918962576451";
          SymFile << endl;
          SymFile << "rot_axis " << "2" << " 0 0 1";
          }
       else if (symmetry.compare("td")==0)
          {
          SymFile << "rot_axis " << "3" << "  .57735026918962576451  .57735026918962576451 .57735026918962576451";
          SymFile << endl;
          SymFile << "rot_axis " << "2" << " 0 0 1";
          SymFile << endl;
          SymFile << "mirror_plane 0 1 1";
          }
       else if (symmetry.compare("th")==0)
          {
          SymFile << "rot_axis " << "3" << "  .57735026918962576451  .57735026918962576451 .57735026918962576451";
          SymFile << endl;
          SymFile << "rot_axis " << "2" << " 0 0 1";
          SymFile << endl;
          SymFile << "inversion";
          }
       else if (symmetry.compare("o")==0)
          {
          SymFile << "rot_axis " << "3" << "  .57735026918962576451  .57735026918962576451 .57735026918962576451";
          SymFile << endl;
          SymFile << "rot_axis " << "4" << " 0 0 1";
          SymFile << endl;
          SymFile << "inversion";
          }
       else if (symmetry.compare("oh")==0)
          {
          SymFile << "rot_axis " << "3" << "  .57735026918962576451  .57735026918962576451 .57735026918962576451";
          SymFile << endl;
          SymFile << "rot_axis " << "4" << " 0 0 1";
          SymFile << endl;
          SymFile << "mirror_plane 0 1 1";
          }
       else if (symmetry.compare("i")==0)
          {
          SymFile << "rot_axis 2  0             0          1";
          SymFile << endl;
          SymFile << "rot_axis 5 -1.618033989  -1           0";
          SymFile << endl;
          SymFile << "rot_axis 3 -0.53934467   -1.4120227   0";
          }
       else if (symmetry.compare("ih")==0)
          {
          SymFile << "rot_axis 2  0             0          1";
          SymFile << endl;
          SymFile << "rot_axis 5 -1.618033989  -1           0";
          SymFile << endl;
          SymFile << "rot_axis 3 -0.53934467   -1.4120227   0";
          SymFile << endl;
          SymFile << "mirror_plane 1 0 0";
          }
       else  
          {
          cerr << "ERROR: Symmetry " << symmetry  << "is not known" << endl;
          exit(0);
          }
       SymFile.close();
       SL.read_sym_file(symmetry_file);

}
