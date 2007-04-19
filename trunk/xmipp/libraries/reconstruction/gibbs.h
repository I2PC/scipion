/***************************************************************************
 *
 * Authors:     R. Marabini
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
/* ------------------------------------------------------------------------- */
/* FCC_GRIDS                                                                 */
/* ------------------------------------------------------------------------- */

#ifndef _GIBBS_HH
#   define _GIBBS_HH

#include <vector>

#include <data/funcs.h>
#include <data/volume.h>
#include <data/geometry.h>
#include <data/args.h>
#include <interface/opendx.h>
#include <interface/vrml.h>

#include "grids.h"

#include <stdlib.h>

/**@name FCC_Grids
    This class contains the basic bricks to work with Gibbs prios in a FCC grid.
    This class is based in the Gridvolume and Grid classes and whenever possible
    I have expanded these two classes instead the FCC-Grids class. But those
    tricks that only work for integer volumes in a FCC grid and that can
    not be easily generalizated are here
*/
//A few defines
#define FCC_no_VOLUMES    1//This may be helpfull if we define the fcc
                           //as four sc
#define FCC_NEIGHBORS    12
#define NEIGHBORS_IN     13// # neighbors including itself
#define RANDOM_NUMBER_POOL_SIZE 10000000 //200000
#define DIF_CONF_EXT 0x2000// number of different clicks=2^13
#define DIF_CONF 0x4000// 2*number of different clicks
                       // the extra ones are for filling valid positions with
                       // invalid clicks, that is spell in the periphery

#define MASKBIT        0x100000
#define PASSMASKBITS   0x0fffff
#define REMASKBIT        0x200000
#define PASSREMASKBITS   0x1fffff
#define PRECISION 100.//Maximum number of decimal in user defined
                         // constant iCo_Border, iCo_homo, iBeta;
/** FCC-Grids class.
*/
template <class T> class FCC_Gibs {
// By the moment only int will be instanciated.
public:
// Structure ---------------------------------------------------------------
    //volume
   GridVolumeT<T> gridvolume;
   //volume with click values in each point
   GridVolumeT<T> gridvolumewithClicks;
   // aux grid volume to save the volume during the process
   GridVolumeT<T> gridvolume_debug;
   // this pointers will help to speed up calculations
   VolumeT<T> * _FCC, * _aux_FCC;
   // these point to their respectibe volumes, in the a
   VolumeT<T> *  _FCC_Click;
   // The easiest way to generate a gridvolume is through a grid.
   Grid FCC_grid;
   //inicializate the volume with...
   typedef enum { ZERO = 0, ONE    = 1,
                  TWO  = 2, RANDOM = -1,
		  FROM_FILE = -2} InitMode;
   // this guys point towards the neighbors, the first 16 are the right
   // ones for a point in the first SC grid and so on.

//   Marsaglia<short> short_Random_pool; // indexes
   Marsaglia<unsigned int> Int_Random_pool;//indexes
   Marsaglia<float> float_Random_pool;//statistics

   typedef enum { OUTSIDE= -2,
                  NON_VALID = -1,
                  BORDER  = 1,
                  HOMOGENEOUS  = 2} ClickMode;
   typedef enum { MY_TRUE= 1,
                  MY_FALSE = 0} MY_SWITCH;

   ClickMode Click_Table[DIF_CONF];
   int EnergyChange[DIF_CONF][NEIGHBORS_IN];
   int ClickChange[DIF_CONF][NEIGHBORS_IN];
   long int iNu_iter; //number of iterations
   //random vector with valid coordinates
   //I know this should be private but I am triying to speed up things
   int * valid_coordinates;
private:
   //history file
    ofstream fh_out;
   // vector that point to neighbours in grid indexes
   matrix1D<int> FCC_Vectors[FCC_NEIGHBORS+1];
   // vector that point to neighbours in real Space
   matrix1D<double> FCC_Vectors_RS[FCC_NEIGHBORS+1];
   // Number of valid points (points insside radius R)
   int iNumber_of_points_inside_sphere;
//   // Number of  points with valid cliques
//   int iNumber_of_valid_points_with_valid_cliques;
//    ^This is iValid_clicks. we do not need it in advance
   // Number of homogeneous clicks
   int iHomo_clicks;
   // Number of border clicks
   int iBorder_clicks;
   // Number of border clicks
   int iValid_clicks;
   // Number of grid points
   int iGridPointNumber;
   //input parameters normalized to integres
   int iCte_Border, iCte_homo, iBeta;
   //Total Energy
   long int Total_Energy;
   //Output file_name seed
   FileName       fn_out_seed;
   //file with initial gridvolume
   FileName       fn_init_grid;
public:
// Fuctions ---------------------------------------------------------------
//  default constructor
   FCC_Gibs()
   {;}
//  default destructor
   ~FCC_Gibs()
   {fh_out.close();}

/** Initialice the class (create a grid, a volume and initialize
    vectors related with the position of the neighbors */

   void Init(double relative_size, double R,
                              FileName Random_Number_File,
                              InitMode initmode) {
          GridVolumeT<T> temporal_gridvolume;
          if(initmode == FCC_Gibs<T>::FROM_FILE)//read init values from old gridvolume
            {
            temporal_gridvolume.read(fn_init_grid);
            R=temporal_gridvolume.grid(0).get_interest_radius();
            relative_size=temporal_gridvolume.grid(0).get_relative_size();

            }
          //Create Grid
          FCC_Gibs<T>::FCC_grid=Create_FCC_grid (relative_size, R);
          //create Gridvolume and aux gridvolume with Clicks
          FCC_Gibs<T>::gridvolume.adapt_to_grid(FCC_grid);
          FCC_Gibs<T>::gridvolume_debug.adapt_to_grid(FCC_grid);
          FCC_Gibs<T>::gridvolumewithClicks.adapt_to_grid(FCC_grid);
          //Initialize pointer to individual SC volumes to speed up things
           _FCC   = &gridvolume(0);
           _FCC_Click   = &gridvolumewithClicks(0);

          //Open history file
          FCC_Gibs<T>::fh_out.open((fn_out_seed+".hist").c_str() , ios::out);
          //save some info
          FCC_Gibs<T>::fh_out
            << "iCte_homo: "          << FCC_Gibs<T>::iCte_homo << endl
            << "iCte_Border: "        << FCC_Gibs<T>::iCte_Border << endl
            << "iBeta: "              << FCC_Gibs<T>::iBeta << endl
            << "iNu_iter: "           << FCC_Gibs<T>::iNu_iter << endl
            << "fn_out_seed: "        << FCC_Gibs<T>::fn_out_seed << endl
            << "fn_init_grid: "        << FCC_Gibs<T>::fn_init_grid << endl
            << "relative_size: "      << relative_size << endl
            << "R: "                  << R << endl
            << "Random_Number_File: " << Random_Number_File << endl
            << "InitMode: "           << initmode << endl
            << "Precision: "          << PRECISION<< endl<< endl
            ;

          FCC_Gibs<T>::fh_out << "GRID DESCRIPTION\n " << FCC_Gibs<T>::FCC_grid;

          // Init random number generator. By the moment I am going to fix the size
          try {
              //select dinamic range for the different random pools
              //Spell positions
              //As soon as we know the number of points inside the
              //sphere Int_Random_pool will be reinizializated.
              FCC_Gibs<T>::Int_Random_pool.Init(Random_Number_File,
                                             RANDOM_NUMBER_POOL_SIZE);

              FCC_Gibs<T>::float_Random_pool.Init(Random_Number_File,
                                          RANDOM_NUMBER_POOL_SIZE);
                           float_Random_pool.Marsaglia_log();
                           float_Random_pool.mul((float)(PRECISION*PRECISION));


          }
          catch (Xmipp_error &XE) {
             cout << XE << endl;
          }//upper program should handle this

          //fill the volume with the initial value
          if(initmode == FCC_Gibs<T>::RANDOM)
           {//init with integer random numbers between 0 an 1
                     FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(* _FCC))
                        {
                        VOLVOXEL((*_FCC),k,i,j) =
                                      (Int_Random_pool.Get_One_Number() & 1) ;
                        }
           }// if end
          else if(initmode == FCC_Gibs<T>::FROM_FILE)
               {
	       VolumeT<T> * temporal_FCC;
	       temporal_FCC = &temporal_gridvolume(0);
                     FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(* _FCC))
                        {
                        VOLVOXEL((*_FCC),k,i,j) =
                                   (VOLVOXEL((*temporal_FCC),k,i,j) & PASSMASKBITS) ;
                        }

	       }
          else{
             T constant;
             if(initmode==FCC_Gibs<T>::ZERO)    // inizialization though gridvolume
                 constant=(T)0;                 // does not work,  I do not know why.
             else if(initmode==FCC_Gibs<T>::ONE)
                 constant=(T)1;
             else if(initmode==FCC_Gibs<T>::TWO)
                 constant=(T)2;
             else
                 REPORT_ERROR(1,(string)"FCC_grid::Init:Wrong Initmode code");
             FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(* FCC_Gibs<T>::_FCC))
                VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j) = constant;
           }//else end

      //#define SINGLEPOINT
      //#define PLOT_SINGLEPOINT
      #ifdef SINGLEPOINT
      {
      matrix1D<double> XYZ=vector_R3(0.,0.,0.);
      //               z y x
      cout << "SINGLEPOINT_DEBUG"<<endl;
      VOLVOXEL((*_FCC),0,0,2) = 0x1;
      cout << "(z,y,x)=0,0,2" <<endl;
      FCC_grid(0).grid2universe(vector_R3(0.,0.,2.),XYZ);
      cout << XYZ.transpose() << endl;
      }
      #ifdef PLOT_SINGLEPOINT
      {
      //print neigh positions and create and openddx file
         matrix1D<double> XYZ=vector_R3(0.,0.,0.);
          cout << "\n\tPLOT_SINGLEPOINT enabled\n";
         openDX DX_r,DX_g,DX_b;
         DX_r.openDXFile( (string) "SINGLEPOINT_DEBUG_red.dx" );
         DX_b.openDXFile( (string) "SINGLEPOINT_DEBUG_blue.dx" );

         //print grid in blue
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
      //cout <<  VOLVOXEL((*_FCC),k,i,j) << endl;

              if( VOLVOXEL((*_FCC),k,i,j) == 0)
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 DX_r.Add_Item(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D

         FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( VOLVOXEL((*_FCC),k,i,j) == 1)
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 DX_b.Add_Item(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D
      }
      #endif
      #endif
      #undef SINGLEPOINT
      #undef PLOT_SINGLEPOINT
      //#define TWOPLANES
      #ifdef TWOPLANES
      /**/
      //VOLVOXEL((*FCC_Gibs<T>::_FCC0),0,0,0) = 0;
      { matrix1D<double> aux_matrix, aux_matrix2;
        double _dot_product;

              aux_matrix.resize(3);aux_matrix2.resize(3);
              aux_matrix2=vector_R3(0.,1.,1.);//one vector
      //            _FCC  = &gridvolume(0);
                     FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
                     {
                        gridvolume.grid(jj).grid2universe(vector_R3( (double)j,
                                                               (double)i,
                                                               (double)k),
                                                                aux_matrix);
                        _dot_product=
                                   dot_product(aux_matrix2,aux_matrix);
                        if(_dot_product> (+0.1))// should be integer except
                                            // for rounding errors
                                  VOLVOXEL((*_FCC),k,i,j) = 0x0;
      //cout << "(x,y,z)= " << j << " " << i << " " << k << " "  << _dot_product << endl;
      //cout << aux_matrix.transpose() << aux_matrix2.transpose()<<endl <<endl;
                     }
      }
      /**/
      #endif
      #undef TWOPLANES


          //alloc space and fill table with the auxiliar vectors
          //these vectors connect each grid point with their neighbours
          InitNeighVector();
      //#define BORDERCLICK_DEBGU
      #ifdef BORDERCLICK_DEBGU
      {
      matrix1D<int> XYZ=vector_R3(0,0,0);
      matrix1D<int> point=vector_R3(1,0,0);
      cout << "BORDERCLICK_DEBGU"<<endl;

      XYZ = point+FCC_Vectors[0x0];
      VOLVOXEL((*_FCC),ZZ(XYZ),YY(XYZ),XX(XYZ)) = 0x1;
      XYZ = point+FCC_Vectors[0x2];
      VOLVOXEL((*_FCC),ZZ(XYZ),YY(XYZ),XX(XYZ)) = 0x1;
      XYZ = point+FCC_Vectors[0x4];
      VOLVOXEL((*_FCC),ZZ(XYZ),YY(XYZ),XX(XYZ)) = 0x1;
      XYZ = point+FCC_Vectors[0x6];
      VOLVOXEL((*_FCC),ZZ(XYZ),YY(XYZ),XX(XYZ)) = 0x1;
      XYZ = point+FCC_Vectors[0x8];
      VOLVOXEL((*_FCC),ZZ(XYZ),YY(XYZ),XX(XYZ)) = 0x1;

      }
      #endif// BORDERCLICK_DEBGU
      #undef BORDERCLICK_DEBGU

          // fill table with valid clicks
          GenerateValidClicks();
          // create mask. Only points inside mask will be modified
          CreateMask();
          //normalize  FCC_Gibs<T>::Int_Random_pool to a maximum
          Int_Random_pool.M_max(Random_Number_File,
                                (T)(iNumber_of_points_inside_sphere-1));
          //only points inside the second mask will have meaningfull clicks
          Create_Second_Mask();
          // alloc memory for vector with random coordinates
          // iNumber_of_points_inside_sphere is calculated in CreateMask
          Alloc_and_Fill_valid_coordinates_vector();
          //calculate configuration energy
          Calculate_Total_Energy();
          //Fill gridvolumewithClicks with Click values
          FillauxGridVolumewithClicks();
          //Klick change due to spell change in klick k
          ClickChangeDueSpellChange();
          //Energy change due to spell change in klick k
          EnergyChangeDueSpellChange();
          //Write initial state in history file
          fh_out << "BEFORE FIRST ITERATION ";
          FCC_Gibs<T>::Print_System_State(FCC_Gibs<T>::fh_out);

      }
/**  Initialice a table containing all valid clicks.
      A Click is a valid one when:
      \begin{itemize}
      \item Is all 0's or all 1's
      \item Given a p, if we draw a vector from each point to the neighborhood
      then all those points which doc product positive are 1, all those
       points with product negative are 0 and at least one of the points with
        dot product zero is equal to the point at the origin.
       \end{itemize}
       The total number of valid click is 72+2
      */
   void GenerateValidClicks(void) {
      double _dot_product;
      int valid_hh_10, valid_hh_1;

       for(int hh=DIF_CONF_EXT; hh < DIF_CONF ; hh++)
          Click_Table[hh]=FCC_Gibs<T>::NON_VALID;// fake clicks for outside of the
                                                 // mask elements
       int center_hh;
       for(int hh=0; hh < DIF_CONF_EXT ; hh++){
          center_hh= hh & (1 << FCC_NEIGHBORS);
          Click_Table[hh]=FCC_Gibs<T>::NON_VALID;
          for(int ii=0; ii < FCC_NEIGHBORS; ii++)
             {
             valid_hh_10=0;// counter for points with dot product different from zero
             valid_hh_1=0; // counter for points with dot product equal to zero
             for(int jj=0; jj < FCC_NEIGHBORS; jj++)
                {
                _dot_product= dot_product(FCC_Vectors_RS[ii],FCC_Vectors_RS[jj]);

                if(_dot_product > XMIPP_EQUAL_ACCURACY )
                  {
                  if( ((hh & (1 << jj)) == 0 && (hh & (1 << ii)) ==0) ||
                      ((hh & (1 << jj)) != 0 && (hh & (1 << ii)) !=0)    )
                         valid_hh_10++;
                  }
                else if (_dot_product > -XMIPP_EQUAL_ACCURACY )
                  {
                  if( (hh & (1 << jj)) == 0 && (center_hh ==0) ||
                      (hh & (1 << jj)) != 0 && (center_hh !=0)    )
                      valid_hh_1=1;
                  }
                else
                {
                  if( ((hh & (1 << jj)) == 0 && (hh & (1 << ii)) !=0) ||
                      ((hh & (1 << jj)) != 0 && (hh & (1 << ii)) ==0)    )
                         valid_hh_10++;
                 }

           if (valid_hh_10==10 && valid_hh_1>=1)
              {
              Click_Table[hh]=FCC_Gibs<T>::BORDER;
              ii=jj=FCC_NEIGHBORS;
              }//valid_hh_10 end

          }/*for jj end*/
        }/*for ii end*/


        }// hh end
        Click_Table[0x0]=FCC_Gibs<T>::HOMOGENEOUS;
        Click_Table[0x1fff]=FCC_Gibs<T>::HOMOGENEOUS;

      //-----------------------
      //Only DEBUG code bellow
      //-----------------------

      //#define   GenerateValidClicks_DEBUG
      #ifdef GenerateValidClicks_DEBUG
          {//begin
          cout << "\nGenerateValidClicks_DEBUG enabled\n";
          for(int hh=0; hh < DIF_CONF ; hh++)
             if(Click_Table[hh]!=FCC_Gibs<T>::NON_VALID)
                cout << hh << endl;

          char fh_FileName[32];
          matrix1D<double>  type_cast_matrix;
          matrix1D<double> RGB=vector_R3(1.,0.,0.);
          matrix1D<double> XYZ;
          int ii;

          for(int hh=0; hh < DIF_CONF ; hh++)
              {
              ii=-1;
              if(Click_Table[hh]==FCC_Gibs<T>::NON_VALID)   continue;
              sprintf(fh_FileName,"klick_%04x.wrl",hh);
              VrmlFile  *_VRML = new VrmlFile((string)fh_FileName);
              RGB=vector_R3(1.,0.,0.);
              XYZ=vector_R3(0.,0.,0.);
              _VRML->Sphere( XYZ,RGB,0.05);
              for(int jj=0; jj < NEIGHBORS_IN ; jj++)
                  {
                  if( (hh & (1<<jj))!=0)
                      {
                      type_cast(FCC_Vectors_RS[jj],type_cast_matrix);
                      _VRML->Add_sphere(type_cast_matrix);
                      }
                  else ii=jj;
                  }
              RGB=vector_R3(0.,0.,1.);
              if(ii != -1)
                {
                  type_cast(FCC_Vectors[ii],XYZ);
                 _VRML->Sphere( XYZ,RGB,0.0525);
                 for(int jj=0; jj < NEIGHBORS_IN ; jj++)
                     {
                     if( (hh & (1<<jj))==0)
                         {
                         type_cast(FCC_Vectors_RS[jj],type_cast_matrix);
                         _VRML->Add_sphere(type_cast_matrix);
                         }
                     }
                }
              //mark center
              RGB=vector_R3(0.,1.,0.);
              XYZ=vector_R3(0.,0.,0.08);
              _VRML->Sphere( XYZ,RGB,0.02);
              delete _VRML;
              }
          }//end
      #endif
      #undef GenerateValidClicks_DEBUG
      }

/**   I do not have time to check if a point is inside or outside the volume
      so I alloc memory for a volume bigger than needed so all the access to
       memory position are always valid. This extra points are marked with
       a 1 in the bit 15*/
   void CreateMask(void) {
      double R;
      matrix1D<double>  distance_vector;
      int flag;


      iNumber_of_points_inside_sphere=0;
      R= gridvolume.grid(0).get_interest_radius()+XMIPP_EQUAL_ACCURACY;

      FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
          {
                                              //x,y,z
          FCC_grid(0).grid2universe(vector_R3((double)j,
                                              (double)i,
                                              (double)k),
                                              distance_vector);

           if( distance_vector.module()>R )
                            //z,y,x
              VOLVOXEL((*_FCC),k,i,j) = (MASKBIT +
                                (VOLVOXEL((*_FCC),k,i,j)&PASSMASKBITS));
            else
              {//check that all the neigh are inside the matrix
               //we will visit them without checking in the future
              flag=0;//for counting neigh
              for (int ii=0; ii<FCC_NEIGHBORS+1; ii++)
                 {
                  if ( VOLMATRIX(*_FCC).outside(k+ZZ(FCC_Vectors[ii]),
                                                i+YY(FCC_Vectors[ii]),
                                                j+XX(FCC_Vectors[ii]) ) )
                         {
                         VOLVOXEL((*_FCC),k,i,j) = (MASKBIT +
                                  (VOLVOXEL((*_FCC),k,i,j)&PASSMASKBITS));
      //cout << "OUT: " << j+XX(FCC_Vectors[ii]) << " "
      //                << i+YY(FCC_Vectors[ii]) << " "
      //                << k+ZZ(FCC_Vectors[ii]) << endl;
                         flag=1;
                         break;
                         }

                 }
               if(flag!=1)
                  iNumber_of_points_inside_sphere++;
              }
          }// FOR ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))

      //----
      //Only DEBUG code Bellow
      //-----------------------
      //#define CreateMask_DEBUG
      #ifdef CreateMask_DEBUG
      {
      //print neigh positions and create and wrl
          cout << "\n\tCreateMask_DEBUG enabled\n";
         matrix1D<double> RGB=vector_R3(0.,0.,1.);
         matrix1D<double> XYZ=vector_R3(0.,0.,0.);

         VrmlFile _VRML( (string) "CreateMask_DEBUG.wrl" );
         //print grid in blue
         _VRML.Sphere( XYZ,RGB,0.05);
         FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==0)
                 continue;
                                                  //x,y,z
              FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
      //       cout << XYZ.transpose() << " " <<vector_R3(j,i,k).transpose() << endl;
              _VRML.Add_sphere(XYZ);
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D
          RGB=vector_R3(1.,0.,0.);
          XYZ=vector_R3(0.,0.,0.);//reset
         _VRML.Sphere( XYZ,RGB,0.05);
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==0)
                 {
                                                     //x,y,z
              FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
              _VRML.Add_sphere(XYZ);
                 }

              }//FOR_ALL_ELEMENTS_IN_MATRIX3D
         //mark center with grey sphere
      //   RGB=vector_R3(.2,.2,.2);
      //   XYZ=vector_R3(0.0,0.,0.);
      //   _VRML.Trans_Sphere( XYZ,RGB, 2.0);
          _VRML.Axis(2.,0.01);

      }
      #endif
      #undef CreateMask_DEBUG

      }


/**   Only point with a complete neighbourh will be used to compute the energy
      This mask marks all this points that do NOT have a complete neighbour.
      The mask makes the bit 16 of the volume equal to 1 */
   void Create_Second_Mask(void) {
      double R;
      matrix1D<double>  distance_vector;

      //iNumber_of_valid_points_with_valid_cliques=0;

      R= gridvolume.grid(0).get_interest_radius()+XMIPP_EQUAL_ACCURACY;
      FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
          {
          for (int ii=0; ii<FCC_NEIGHBORS+1; ii++)
             {
             FCC_grid(0).grid2universe(vector_R3((double)j,
                                                 (double)i,
                                                 (double)k),
                                                 distance_vector);

             distance_vector = distance_vector+FCC_Vectors_RS[ii];
             if( distance_vector.module()>R ||
                  (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==MASKBIT )
                {
                  VOLVOXEL((*_FCC),k,i,j) = (REMASKBIT +
                                    (VOLVOXEL((*_FCC),k,i,j)&PASSREMASKBITS));

                  break;
                  }//if end
      //        if(ii==FCC_NEIGHBORS)
      //           iNumber_of_valid_points_with_valid_cliques++;
               }//for ii
            }// FOR ALL

      //----
      //Only DEBUG code Bellow
      //-----------------------
      //#define SecondMask_DEBUG
      #ifdef SecondMask_DEBUG
      {
      //print neigh positions and create and wrl
          cout << "\n\tSecondMask_DEBUG enabled\n";
         VrmlFile _VRML( (string) "SecondMask_DEBUG.wrl" );
         matrix1D<double> RGB=vector_R3(0.,0.,1.);
         matrix1D<double> XYZ=vector_R3(0.,0.,0.);

         //print grid in blue
         _VRML.Sphere( XYZ,RGB,0.05);
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==MASKBIT)
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 _VRML.Add_sphere(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D
          RGB=vector_R3(0.,1.,0.);
          XYZ=vector_R3(0.,0.,0.);//reset
         _VRML.Sphere( XYZ,RGB,0.05);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & REMASKBIT) ==REMASKBIT
               && (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) !=MASKBIT )
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 _VRML.Add_sphere(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D
          RGB=vector_R3(1.,0.,0.);
          XYZ=vector_R3(0.,0.,0.);//reset
         _VRML.Sphere( XYZ,RGB,0.05);
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==0
                &&(VOLVOXEL((*_FCC),k,i,j) & REMASKBIT) ==0)
                 {
                                                     //x,y,z
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                     (double)i,
                                                     (double)k),
                                                     XYZ);
                 _VRML.Add_sphere(XYZ);
                 }

              }//FOR_ALL_ELEMENTS_IN_MATRIX3D
         //mark center with grey sphere
          _VRML.Axis(2.,0.01);

      }
      #endif
      #undef SecondMask_DEBUG
      //#define SecondMask_DX_DEBUG
      #ifdef SecondMask_DX_DEBUG
      {
      //print neigh positions and create and openddx file
         matrix1D<double> XYZ=vector_R3(0.,0.,0.);
          cout << "\n\tSecondMask_DX_DEBUG enabled\n";
         openDX DX_r,DX_g,DX_b;
         DX_r.openDXFile( (string) "SecondMask_DX_DEBUG_red.dx" );
         DX_g.openDXFile( (string) "SecondMask_DX_DEBUG_green.dx" );
         DX_b.openDXFile( (string) "SecondMask_DX_DEBUG_blue.dx" );

         //print grid in blue
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==MASKBIT)
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 DX_b.Add_Item(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D



        FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & REMASKBIT) ==REMASKBIT
               && (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) !=MASKBIT )
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 DX_g.Add_Item(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              {
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) ==0
                &&(VOLVOXEL((*_FCC),k,i,j) & REMASKBIT) ==0)
                 {
                                                     //x,y,z
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                     (double)i,
                                                     (double)k),
                                                     XYZ);
                 DX_r.Add_Item(XYZ);
                 }

              }//FOR_ALL_ELEMENTS_IN_MATRIX3D
         //mark center with grey sphere
      }
      #endif
      #undef SecondMask_DX_DEBUG


      }


/**    To speed up the process the coordinates of the valid points are stored            (all the valid coordinates, even without valid click go here) */
   void Alloc_and_Fill_valid_coordinates_vector(void) {
      //alloc memory
      if( (valid_coordinates=(int *)
                 malloc(iNumber_of_points_inside_sphere*3*sizeof(int)))==NULL)
         {
         cerr << "ERROR: allocating memory for valid_coordinates\n"
              << "Volume to big?" <<endl; exit(0);
         }
      //fill it
      int icounter=0;
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
              if( (VOLVOXEL((*_FCC),k,i,j) & MASKBIT) !=MASKBIT)
                   {
                   valid_coordinates[3*icounter+0]=k;
                   valid_coordinates[3*icounter+1]=i;
                   valid_coordinates[3*icounter+2]=j;
                   icounter++;
                   }//if end
      if(icounter!=(iNumber_of_points_inside_sphere))
         {
         cerr << "\nERROR: this should "
              << "\nnever happend(Alloc_and_Fill_valid_coordinates_vector)\n"
              << "\nicounter="<< icounter <<endl
              << "\niNumber_of_points_inside_sphere=" <<
                  iNumber_of_points_inside_sphere;
         }

      //#define Alloc_and_Fill_valid_coordinates_vector_before
      #ifdef Alloc_and_Fill_valid_coordinates_vector_before
      {
      cout << "\niNumber_of_points_inside_sphere-before\n";
      for(int ii=0; ii<iNumber_of_points_inside_sphere; ii++)
          {
          cout << valid_coordinates[3*ii+0] << " "
               << valid_coordinates[3*ii+1] << " "
               << valid_coordinates[3*ii+2] << endl;
          }//if end
      }
      #endif// Alloc_and_Fill_valid_coordinates_vector_before
      #undef Alloc_and_Fill_valid_coordinates_vector_before

      //randomize it
      //find next power_of_two
      int random1, random2;
      int XX,YY,ZZ;

      for(int ii=0; ii<iNumber_of_points_inside_sphere; ii++)
         {
              random1=Int_Random_pool.Get_One_Number();
              random2=Int_Random_pool.Get_One_Number();
              ZZ=valid_coordinates[3*random1+0];
              YY=valid_coordinates[3*random1+1];
              XX=valid_coordinates[3*random1+2];
              valid_coordinates[3*random1+0]=valid_coordinates[3*random2+0];
              valid_coordinates[3*random1+1]=valid_coordinates[3*random2+1];
              valid_coordinates[3*random1+2]=valid_coordinates[3*random2+2];
              valid_coordinates[3*random2+0]=ZZ;
              valid_coordinates[3*random2+1]=YY;
              valid_coordinates[3*random2+2]=XX;
         }//randomize end (for ii=0
      //#define Alloc_and_Fill_valid_coordinates_vector_after
      #ifdef  Alloc_and_Fill_valid_coordinates_vector_after
      {
      cout << "\niNumber_of_points_inside_sphere-after\n";
      for(int ii=0; ii<iNumber_of_points_inside_sphere; ii++)
          {
          cout << valid_coordinates[3*ii+0] << " "
               << valid_coordinates[3*ii+1] << " "
               << valid_coordinates[3*ii+2] << endl;
          }//if end
      }
      #endif// Alloc_and_Fill_valid_coordinates_vector_after
      #undef Alloc_and_Fill_valid_coordinates_vector_after
      //#define Alloc_and_Fill_valid_coordinates_vector_after_DX
      #ifdef Alloc_and_Fill_valid_coordinates_vector_after_DX
      {
      //print neigh positions and create and openddx file
         matrix1D<double> XYZ=vector_R3(0.,0.,0.);
          cout << "\nAlloc_and_Fill_valid_coordinates_vector_after_DX enabled\n";
         openDX DX_r,DX_g,DX_b;
         DX_r.openDXFile( (string) "alloc_and_fill.dx" );

         //print grid in red
      for(int ii=0; ii<iNumber_of_points_inside_sphere; ii++)
              {
                                                  //x,y,z
               FCC_grid(0).grid2universe(vector_R3(
                                                  (double)valid_coordinates[3*ii+2],
                                                  (double)valid_coordinates[3*ii+1],
                                                  (double)valid_coordinates[3*ii+0]),
                                                  XYZ);
                 DX_r.Add_Item(XYZ);
              }
      }
      #endif// Alloc_and_Fill_valid_coordinates_vector_after_DX
      #undef Alloc_and_Fill_valid_coordinates_vector_after_DX

      }


/**   Check if the point is the center of a valid Click
      It the point do not have 12 neighbours it returns OUTSIDE,
      If the point is the center of a valid click it returns either
      HOMOGENEOUS, BORDER or NON_VALID (that is neither HOMOGENEOUS nor
      BORDER*/
   int IsThisValidClick(int j, int i, int k) //x y z
      {

      int ii=WhatClickIsThis(  k, i, j);
      if(ii == FCC_Gibs<T>::OUTSIDE)
         return (FCC_Gibs<T>::OUTSIDE);
      else
         return( Click_Table[ii] );//either border, homogeneous or nonvalid

      }//end IsThisValidClick

/**   Check if the point is the center of a valid Click
      It the point do not have 12 neighbours it returns OUTSIDE,
      If the point is the center of a valid click it returns the
      click code*/
   int WhatClickIsThis(int j, int i, int k) //x y z
      {
       int neigh_value;

          if(  (VOLVOXEL( *FCC_Gibs<T>::_FCC,k,i,j) & REMASKBIT)==REMASKBIT )
             return (FCC_Gibs<T>::OUTSIDE);    //z,y,x
          neigh_value =
      //this should match FCC_Vectors
              ((VOLVOXEL((*_FCC),k+0,i+0,j+0) & 0x1) << 0xc) +

              ((VOLVOXEL((*_FCC),k+1,i+0,j+0) & 0x1) << 0x0) +
              ((VOLVOXEL((*_FCC),k-1,i+0,j+0) & 0x1) << 0x1) +
              ((VOLVOXEL((*_FCC),k+0,i-1,j+1) & 0x1) << 0x2) +
              ((VOLVOXEL((*_FCC),k+0,i+1,j-1) & 0x1) << 0x3) +

              ((VOLVOXEL((*_FCC),k+1,i+0,j-1) & 0x1) << 0x4) +
              ((VOLVOXEL((*_FCC),k-1,i+0,j+1) & 0x1) << 0x5) +
              ((VOLVOXEL((*_FCC),k+0,i-1,j+0) & 0x1) << 0x6) +
              ((VOLVOXEL((*_FCC),k+0,i+1,j+0) & 0x1) << 0x7) +

              ((VOLVOXEL((*_FCC),k+1,i-1,j+0) & 0x1) << 0x8) +
              ((VOLVOXEL((*_FCC),k-1,i+1,j+0) & 0x1) << 0x9) +
              ((VOLVOXEL((*_FCC),k+0,i+0,j-1) & 0x1) << 0xa) +
              ((VOLVOXEL((*_FCC),k+0,i+0,j+1) & 0x1) << 0xb) ;

      return(neigh_value);

      }

/** Count clicks of each type, So long homogeneous or border */
   void Count_clicks_of_each_type(void) {
        int click_code;
        iHomo_clicks=iBorder_clicks=iValid_clicks=0;

        try {
      //          _FCC   = &gridvolume(0);
                FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC)){
                   click_code=IsThisValidClick(k,i,j);
                   if(click_code!=FCC_Gibs<T>::OUTSIDE) iValid_clicks++;
                   if(click_code==FCC_Gibs<T>::HOMOGENEOUS){
                       iHomo_clicks++;
                       }
                   else if(click_code==FCC_Gibs<T>::BORDER)
                       iBorder_clicks++;
               }//FOR_ALL
          }//try
          catch (Xmipp_error &XE) {
          cout << XE << endl;
       }//upper program should handle this

      }

/** Calculate total Energy */
   void Calculate_Total_Energy(void){
   Count_clicks_of_each_type();
   Total_Energy  = (long int)iCte_Border * (long int)iBorder_clicks +
                   (long int)iCte_homo   * (long int)iHomo_clicks;
   Total_Energy *= (long int) iBeta;
   }
/** Fill auxiliar Grid volume with the value of the clicks centered on
    each point (of the principal grid volume */
   void FillauxGridVolumewithClicks(void) {
        try {
                FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC_Click)){

                   VOLVOXEL((*_FCC_Click),k,i,j) =
                                WhatClickIsThis(k,i,j);
                   if (VOLVOXEL((*_FCC_Click),k,i,j)==FCC_Gibs<T>::OUTSIDE)
                       VOLVOXEL((*_FCC_Click),k,i,j)=0x2000;
                                 // 0x2000 is a click that is not changed
                                 // by the update click mechanism
                                 // in a valid click

      //#define ONE_POINT_DEBUG
      #ifdef ONE_POINT_DEBUG
      if(VOLVOXEL((*_FCC_Click),k,i,j)!=0x2000 && VOLVOXEL((*_FCC_Click),k,i,j)!=0)
      {
      cout << k << " " << i << " "<< j <<"\t\t" ;
      printb(cout,VOLVOXEL((*_FCC_Click),k,i,j)); cout << endl;
      //remember that you should get opposite vectors
      }
      #endif
      #undef ONE_POINT_DEBUG
               }//FOR_ALL
          }//try
          catch (Xmipp_error &XE) {
          cout << XE << endl;
       }//upper program should handle this
      //#define FillauxGridVolumewithClicks_DEBUG
      #ifdef FillauxGridVolumewithClicks_DEBUG
      {
      //print neigh positions and create and openddx file
         matrix1D<double> XYZ=vector_R3(0.,0.,0.);
          cout << "\n\tFillauxGridVolumewithClicks_DEBUG enabled\n";
         openDX DX_r,DX_g,DX_b;
         DX_r.openDXFile( (string) "FillauxGridVolumewithClicks_DEBUG_red.dx" );
         DX_b.openDXFile( (string) "FillauxGridVolumewithClicks_DEBUG_blue.dx" );

         //print grid in blue
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC_Click))
              {
              if( VOLVOXEL((*_FCC_Click),k,i,j) !=0x2000)
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 DX_r.Add_Item(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D



        FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC_Click))
              {
              if( VOLVOXEL((*_FCC_Click),k,i,j) ==0x2000)
                                                  //x,y,z
                 {
                 FCC_grid(0).grid2universe(vector_R3((double)j,
                                                  (double)i,
                                                  (double)k),
                                                  XYZ);
                 DX_b.Add_Item(XYZ);
                 }
          }//FOR_ALL_ELEMENTS_IN_MATRIX3D
      }
      #endif
      #undef FillauxGridVolumewithClicks_DEBUG

      }

/** Store the user provided parameters as integers */
   void Set_User_Parameters(int homo, int border, int beta,
                            long int nu_iter, FileName fh_out_name,
			                       FileName fn_initfile)
   {
   iCte_homo   = homo;
   iCte_Border = border;
   iBeta       = beta;
   iNu_iter    = nu_iter;
   fn_out_seed = fh_out_name;
   fn_init_grid = fn_initfile;
   }
/** Create a table (of size  DIF_CONF x NEIGHBORS_IN) with the value of the
    new click when in the klick k the voxel v is changed*/
   void ClickChangeDueSpellChange(void){
       for(int hh=DIF_CONF_EXT; hh < DIF_CONF ; hh++)
          for(int jj=0; jj < NEIGHBORS_IN ; jj++)
             ClickChange[hh][jj]=0x2000;   // > 0x1fff
       for(int hh=0; hh < DIF_CONF_EXT ; hh++){
          ClickChange[hh][0x0]=hh^0x0001; //01
          ClickChange[hh][0x1]=hh^0x0002; //02
          ClickChange[hh][0x2]=hh^0x0004; //03
          ClickChange[hh][0x3]=hh^0x0008; //04
          ClickChange[hh][0x4]=hh^0x0010; //05
          ClickChange[hh][0x5]=hh^0x0020; //06
          ClickChange[hh][0x6]=hh^0x0040; //07
          ClickChange[hh][0x7]=hh^0x0080; //08
          ClickChange[hh][0x8]=hh^0x0100; //09
          ClickChange[hh][0x9]=hh^0x0200; //10
          ClickChange[hh][0xa]=hh^0x0400; //11
          ClickChange[hh][0xb]=hh^0x0800; //12
          ClickChange[hh][0xc]=hh^0x1000; //13-4096
         } // end for(int hh=0; hh < DIF_CONF ; hh++){
      //#define PRINT_CLICKCHANGETABLE_DEBUG
      #ifdef PRINT_CLICKCHANGETABLE_DEBUG
      cout << "PRINT_CLICKCHANGETABLE_DEBUG" << endl;
       for(int hh=0; hh < DIF_CONF_EXT ; hh++)
          {
      //    printb(cout, hh); cout << " ";
           cout << hh << "\t-> ";
          for(int ii=0; ii < NEIGHBORS_IN ; ii++)
               {
      //         printb(cout,ClickChange[hh][ii]); cout << " ";
                 cout << ClickChange[hh][ii] << " ";
               }
          cout << endl;
          }
      #endif// PRINT_CLICKCHANGETABLE_DEBUG
      #undef PRINT_CLICKCHANGETABLE_DEBUG

      }

/** Create a table (of size  DIF_CONF x NEIGHBORS_IN) with the energy change
    when in the click k the voxel v is changed */
   void EnergyChangeDueSpellChange(void) {
       int iCte_Border_by_iBeta;
       int iCte_homo_by_iBeta;
       iCte_Border_by_iBeta = iCte_Border * iBeta;
       iCte_homo_by_iBeta   = iCte_homo * iBeta;

       for(int hh=DIF_CONF_EXT; hh < DIF_CONF ; hh++)
          for(int jj=0; jj < NEIGHBORS_IN ; jj++)
             EnergyChange[hh][jj]=0x0;   // > 0x1fff
       for(int hh=0; hh < DIF_CONF_EXT ; hh++){
      // cout << endl << hh << "\t";
           for(int jj=0; jj < NEIGHBORS_IN ; jj++){
               if(Click_Table[hh]==FCC_Gibs<T>::NON_VALID)
                                            //either border, homogeneous or nonvalid
                  EnergyChange[hh][jj]=0;
               else{//b
                  if(Click_Table[hh]==FCC_Gibs<T>::HOMOGENEOUS)
                       EnergyChange[hh][jj] = -iCte_homo_by_iBeta;
                  else if (Click_Table[hh] == FCC_Gibs<T>::BORDER)
                       EnergyChange[hh][jj] = -iCte_Border_by_iBeta;
                  else//a
                     {
                     cout << "ERROR1: If you see this line something is going"
                          << "wrong in the routine EnergyChangeDueSpellChange" <<endl;
                     exit(1);
                     }//end else a
               }//end else b
               if(Click_Table[ClickChange[hh][jj]]==FCC_Gibs<T>::NON_VALID)
                  ;
               else{//b
                  if(Click_Table[ClickChange[hh][jj]]==FCC_Gibs<T>::HOMOGENEOUS)
                       EnergyChange[hh][jj] += iCte_homo_by_iBeta;
                  else if (Click_Table[ClickChange[hh][jj]] == FCC_Gibs<T>::BORDER)
                       EnergyChange[hh][jj] += iCte_Border_by_iBeta;
                  else//a
                     {
                     cout << "ERROR2: If you see this line something is going"
                          << "wrong in the routine EnergyChangeDueSpellChange" <<endl;
                     exit(1);
                     }//end else a
               }//end else b
      // cout << EnergyChange[hh][jj] << "\t";
             } // end for(int jj=0; jj < DIF_CONF ; jj++)
         } // end for(int hh=0; hh < DIF_CONF ; hh++){

      //#define ENERGY_CHANGE_DUE_TO_SPELL_CHANGE_DEBUG
      #ifdef ENERGY_CHANGE_DUE_TO_SPELL_CHANGE_DEBUG
      cout << "ENERGY_CHANGE_DUE_TO_SPELL_CHANGE_DEBUG"<<endl;
       for(int hh=0; hh < DIF_CONF_EXT ; hh++)
          {
      //    printb(cout, hh); cout << " ";
           cout << hh << "\t-> ";
          for(int ii=0; ii < NEIGHBORS_IN ; ii++)
               {
      //         printb(cout,ClickChange[hh][ii]); cout << " ";
                 cout << EnergyChange[hh][ii] << " ";
               }
          cout << endl;
          }
      #endif// ENERGY_CHANGE_DUE_TO_SPELL_CHANGE_DEBUG
      #undef ENERGY_CHANGE_DUE_TO_SPELL_CHANGE_DEBUG

      }

/** Calculate variation of energy due to spell flip */
   int Get_Energy_Change( int k, int i, int j) //z y x
      {

      int energy_change;
          energy_change=
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+0,j+0)][0xc] +

            EnergyChange[VOLVOXEL(*_FCC_Click,k+1,i+0,j+0)][0x1] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k-1,i+0,j+0)][0x0] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i-1,j+1)][0x3] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+1,j-1)][0x2] +

            EnergyChange[VOLVOXEL(*_FCC_Click,k+1,i+0,j-1)][0x5] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k-1,i+0,j+1)][0x4] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i-1,j+0)][0x7] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+1,j+0)][0x6] +

            EnergyChange[VOLVOXEL(*_FCC_Click,k+1,i-1,j+0)][0x9] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k-1,i+1,j+0)][0x8] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+0,j-1)][0xb] +
            EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+0,j+1)][0xa] ;

      //#define GET_ENERGY_CHANGE
      #ifdef GET_ENERGY_CHANGE
      //WARNING: if this debug activated all changes are acepted!!!!!!!!!!!!!!!!
      //compare energy_change with manual computing
      if (energy_change>=0)
      {
      cout << "\nGET_ENERGY_CHANGE_AFTER"<<endl;
      //---------

      cout <<  "0:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+0,j+0)][0xc] <<endl;

      cout <<  "1:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+1,i+0,j+0)][0x1] <<endl;
      cout <<  "2:" << EnergyChange[VOLVOXEL(*_FCC_Click,k-1,i+0,j+0)][0x0] <<endl;
      cout <<  "3:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i-1,j+1)][0x3] <<endl;
      cout <<  "4:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+1,j-1)][0x2] <<endl;

      cout <<  "5:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+1,i+0,j-1)][0x5] <<endl;
      cout <<  "6:" << EnergyChange[VOLVOXEL(*_FCC_Click,k-1,i+0,j+1)][0x4] <<endl;
      cout <<  "7:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i-1,j+0)][0x7] <<endl;
      cout <<  "8:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+1,j+0)][0x6] <<endl;

      cout <<  "9:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+1,i-1,j+0)][0x9] <<endl;
      cout << "10:" << EnergyChange[VOLVOXEL(*_FCC_Click,k-1,i+1,j+0)][0x8] <<endl;
      cout << "11:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+0,j-1)][0xb] <<endl;
      cout << "12:" << EnergyChange[VOLVOXEL(*_FCC_Click,k+0,i+0,j+1)][0xa] <<endl;

      //--------
      //    Set_Click_Change(k,i,j);
      //    Set_Spell_Change(k,i,j);

          cout << k << " " << i << " "<< j <<": "<< energy_change;
      //    FCC_Gibs<T>::Calculate_Total_Energy();// Count_clicks_of_each_type   called from
                                        // Calculate_Total_Energy
      //    FCC_Gibs<T>::Print_System_State(cout);
      }
      #endif// GET_ENERGY_CHANGE
      #undef GET_ENERGY_CHANGE
          return(energy_change);

      }

/** Get iNumber_of_points_inside_sphere */
    int Get_iNumber_of_points_inside_sphere(void)
     { return(iNumber_of_points_inside_sphere);}
/** Update volume grid with clicks */
   void Set_Click_Change(int k, int i, int j) //z y x
      {

            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i+0,j+0) ^= 0x1000;

            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+1,i+0,j+0) ^= 0x0002;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k-1,i+0,j+0) ^= 0x0001;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i-1,j+1) ^= 0x0008;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i+1,j-1) ^= 0x0004;

            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+1,i+0,j-1) ^= 0x0020;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k-1,i+0,j+1) ^= 0x0010;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i-1,j+0) ^= 0x0080;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i+1,j+0) ^= 0x0040;

            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+1,i-1,j+0) ^= 0x0200;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k-1,i+1,j+0) ^= 0x0100;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i+0,j-1) ^= 0x0800;
            VOLVOXEL(*FCC_Gibs<T>::_FCC_Click,k+0,i+0,j+1) ^= 0x0400;
      }

/** Update volume grid with voxels  */
   void Set_Spell_Change(int k, int i, int j) //z y x
      {
                VOLVOXEL( *(_FCC),z,y,x) ^= 1;
      }

/** Add text to the history file */
   void Add_text_to_the_history_file(char * line)
   {
   if (!fh_out.is_open ())
      {
       cout << "Error: history file is not open"
            << endl;
       exit (1);
      }
   FCC_Gibs<T>::fh_out << line ;
   }
/** Count Number of spells with value n */
   int Count_Number_Spells_Value_n (T n)
      {
        int n_counter=0;

        FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
           if ( (VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j)&PASSMASKBITS)==n
                && (VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j)&MASKBIT)!=MASKBIT)
                 n_counter++;
        return(n_counter);

      }

/** Remove Mask and save gridvolume */
    void Save_gridvolume(void)
    {
    fh_out << "AFTER LAST ITERATION ";
//    fh_out << Get_Total_energy() << endl;
    Calculate_Total_Energy();
    Print_System_State(FCC_Gibs<T>::fh_out);
    Remove_Mask();
    gridvolume.write((fn_out_seed+".fcc").c_str());
    }
/** Copy volume to auxiliar volume, Remove Mask and save aux gridvolume */
    void Save_debug_gridvolume(void)
    {
    static int priv_counter=0;
//    Copy_to_Gridvolume_debug_and_Remove_Mask();
    Copy_to_Gridvolume_debug();
    gridvolume_debug.write((fn_out_seed+ItoA(priv_counter++,4)+".fcc").c_str());
    }

/** Get Total Energy */
   long int Get_Total_energy(void) {return (Total_Energy);}
/** Update Total Energy */
   void Update_Total_energy(int energy_change) {Total_Energy += energy_change;}  /** Remove gridvolume mask */
   void Remove_Mask(void)
      {

          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
             VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j) &= 0x1;

      }

/** auxiliar function to save intermediate results. */
   void Copy_to_Gridvolume_debug(void)
      {

        _aux_FCC = &gridvolume_debug(0);
             FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
                VOLVOXEL((*FCC_Gibs<T>::_aux_FCC),k,i,j) =
                VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j) & 0x1;

      }

/** Get Beta value */
   int Get_Beta(void) { return(iBeta);}
/** Get Homogeneous clique  constant */
   int Get_Cte_homo(void) { return(iCte_homo);}
#ifdef NEVERDEFINED//a
/** Count Number of spells */
   void Count_Number_Spells(void);

//DEBUG FUNCTIONS

#endif//if NEVERDEFINED//a
void Print_System_State(ostream &o){
 o << "\nSYSTEM STATE"
   << "\n  iHomo_clicks= " << iHomo_clicks
   << "\n  iBorder_clicks=   " << iBorder_clicks
   << "\n  iValid_clicks=    " << iValid_clicks
   << "\n  iNumber_of_points_inside_sphere= "
                              << iNumber_of_points_inside_sphere
//   << "\n iNumber_of_valid_points_with_valid_cliques= "
//                              << iNumber_of_valid_points_with_valid_cliques
   << "\n  Total_Energy= " << Total_Energy  << endl;
   }
void Print_All_Clicks(void){
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(* _FCC_Click))
             if(  (VOLVOXEL( *_FCC,k,i,j) & REMASKBIT)!=REMASKBIT )
                cout << k << "\t" << i << "\t" << j
                     << "\t" << VOLVOXEL((*_FCC_Click),k,i,j) << endl;

   }

private:

/**  Alloc space and fill table with the auxiliar vectors called
     FCC_Vectors. These vectors connect each grid point with their
     neighbours */

     void InitNeighVector(void) {
          matrix1D<double>  type_cast_matrix;
          //alloc space for auxiliar vector
          //these vectors conect each grid point with their neighbours
          for(int ii=0; ii< FCC_NEIGHBORS+1; ii++)
              {
              FCC_Gibs<T>::FCC_Vectors[ii].resize(3);
              FCC_Gibs<T>::FCC_Vectors_RS[ii].resize(3);
              }

      //----neigh in INDEX space
      // The order of this vectors is crucial
      // The order of this vectors is crucial and should not be changed
      // DO NOT CHANGE IT!!!!!!!!!
      // just in case I did not make myself clear  DO NOT CHANGE THE ORDER!!!!!!!!!
      // rmember order here is x,y,z but in volumes is z,y,x
                                            //x,y,z
          FCC_Vectors[0x0]  = vector_R3( +0, +0, +1);
          FCC_Vectors[0x1]  = vector_R3( +0, +0, -1);
          FCC_Vectors[0x2]  = vector_R3( +1, -1, +0);
          FCC_Vectors[0x3]  = vector_R3( -1, +1, +0);

          FCC_Vectors[0x4]  = vector_R3( -1, +0, +1);
          FCC_Vectors[0x5]  = vector_R3( +1, +0, -1);
          FCC_Vectors[0x6]  = vector_R3( +0, -1, +0);
          FCC_Vectors[0x7]  = vector_R3( +0, +1, +0);

          FCC_Vectors[0x8]  = vector_R3( +0, -1, +1);
          FCC_Vectors[0x9]  = vector_R3( +0, +1, -1);
          FCC_Vectors[0xa]  = vector_R3( -1, +0, +0);
          FCC_Vectors[0xb]  = vector_R3( +1, +0, +0);

          FCC_Vectors[0xc]  = vector_R3( +0, +0, +0);

          for(int ii=0; ii< FCC_NEIGHBORS + 1; ii++)
              {
              type_cast(FCC_Vectors[ii],type_cast_matrix);
              FCC_grid(0).grid2universe(type_cast_matrix,FCC_Vectors_RS[ii]);
              }
      //-----------------------
      //Only DEBUG code bellow
      //-----------------------
      //#define InitNeighVector_DEBUG1
      #ifdef InitNeighVector_DEBUG1
      {
      //print neigh positions and create and wrl
      //the printed data is in index not real space coordinates
          matrix1D<double>  matrix_aux;
          type_cast(FCC_Vectors[0],type_cast_matrix);
          FCC_grid(0).grid2universe(type_cast_matrix, matrix_aux);
          matrix1D<double> RGB=vector_R3(0.,0.,1.);
          matrix1D<double> XYZ=matrix_aux;
          cout << "\n\tInitNeighVector_DEBUG1 enabled\n";

          for(int ii=0; ii< FCC_NEIGHBORS+1; ii++)
              cout << FCC_Gibs<T>::FCC_Vectors[ii].transpose() << endl;

          VrmlFile _VRML( (string) "InitNeighVector_DEBUG1.wrl" );
          _VRML.Sphere( XYZ,RGB,0.05);
          for(int ii=0; ii< FCC_NEIGHBORS+1; ii++)
              {
              type_cast(FCC_Vectors[ii],type_cast_matrix);
              FCC_grid(0).grid2universe(type_cast_matrix, matrix_aux);
              _VRML.Add_sphere(matrix_aux);
              }//for(int ii=0; ii< FCC_NEIGHBORS; ii++)
         //mark center
          RGB=vector_R3(1.,0.,0.);
          XYZ=vector_R3(0.,0.,0.08);
          _VRML.Sphere( XYZ,RGB,0.02);
      }
      #endif
      #undef InitNeighVector_DEBUG1

      //#define InitNeighVector_DEBUG2
      #ifdef InitNeighVector_DEBUG2
      {
      //print neigh positions using FCC_NEIGHBORS_RS and create and wrl
          cout << "\n\tInitNeighVector_DEBUG2 enabled\n";
          for(int ii=0; ii< FCC_NEIGHBORS+1; ii++)
              cout << FCC_Gibs<T>::FCC_Vectors_RS[ii].transpose() << endl;

          matrix1D<double>  matrix_aux;
          matrix1D<double> RGB=vector_R3(1.,0.,0.);
          matrix1D<double> XYZ=matrix_aux;
          VrmlFile _VRML( (string) "InitNeighVector_DEBUG2.wrl" );
          RGB=vector_R3(1.,0.,0.);
          XYZ=vector_R3(0.,0.,0.);
          _VRML.Sphere( XYZ,RGB,0.05);
          for(int ii=0; ii< FCC_NEIGHBORS+1; ii++)
              {
              _VRML.Add_sphere(FCC_Vectors_RS[ii]);
               }//for(int ii=0; ii< FCC_NEIGHBORS_RS; ii++)
         //mark center
          RGB=vector_R3(0.,0.,1.);
          XYZ=vector_R3(0.,0.,0.08);
          _VRML.Sphere( XYZ,RGB,0.02);

      }
      #endif
      #undef InitNeighVector_DEBUG2


      }
};
#endif
