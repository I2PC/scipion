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

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippGeometry.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippInterface/xmippOpenDX.hh>
#include <XmippInterface/xmippVrml.hh>
#include "grids.hh"
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
   long long iNu_iter; //number of iterations
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
   long long Total_Energy;
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
                              InitMode initmode) /*_THROW*/;
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
   void GenerateValidClicks(void); // called from init

/**   I do not have time to check if a point is inside or outside the volume
      so I alloc memory for a volume bigger than needed so all the access to
       memory position are always valid. This extra points are marked with 
       a 1 in the bit 15*/
   void CreateMask(void);
/**   Only point with a complete neighbourh will be used to compute the energy
      This mask marks all this points that do NOT have a complete neighbour.
      The mask makes the bit 16 of the volume equal to 1 */
   void Create_Second_Mask(void);
/**    To speed up the process the coordinates of the valid points are stored            (all the valid coordinates, even without valid click go here) */
   void Alloc_and_Fill_valid_coordinates_vector(void);

/**   Check if the point is the center of a valid Click 
      It the point do not have 12 neighbours it returns OUTSIDE,
      If the point is the center of a valid click it returns either
      HOMOGENEOUS, BORDER or NON_VALID (that is neither HOMOGENEOUS nor 
      BORDER*/
   int IsThisValidClick(int j, int i, int k);//x y z
/**   Check if the point is the center of a valid Click 
      It the point do not have 12 neighbours it returns OUTSIDE,
      If the point is the center of a valid click it returns the 
      click code*/
   int WhatClickIsThis(int j, int i, int k);//x y z

/** Count clicks of each type, So long homogeneous or border */
   void Count_clicks_of_each_type(void);

/** Calculate total Energy */
   void Calculate_Total_Energy(void){
   Count_clicks_of_each_type();
   Total_Energy  = (long long)iCte_Border * (long long)iBorder_clicks +
                   (long long)iCte_homo   * (long long)iHomo_clicks;
   Total_Energy *= (long long) iBeta;
   }
/** Fill auxiliar Grid volume with the value of the clicks centered on 
    each point (of the principal grid volume */
   void FillauxGridVolumewithClicks(void);
/** Store the user provided parameters as integers */
   void Set_User_Parameters(int homo, int border, int beta, 
                            long long nu_iter, FileName fh_out_name,
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
   void ClickChangeDueSpellChange(void);
/** Create a table (of size  DIF_CONF x NEIGHBORS_IN) with the energy change
    when in the click k the voxel v is changed */
   void EnergyChangeDueSpellChange(void);
/** Calculate variation of energy due to spell flip */
   int Get_Energy_Change( int k, int i, int j);//z y x
/** Get iNumber_of_points_inside_sphere */
    int Get_iNumber_of_points_inside_sphere(void)
     { return(iNumber_of_points_inside_sphere);}
/** Update volume grid with clicks */
   void Set_Click_Change(int k, int i, int j);//z y x   
/** Update volume grid with voxels  */
   void Set_Spell_Change(int k, int i, int j);//z y x   
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
   int Count_Number_Spells_Value_n (T n); 
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
   long long Get_Total_energy(void) {return (Total_Energy);}   
/** Update Total Energy */
   void Update_Total_energy(int energy_change) {Total_Energy += energy_change;}  /** Remove gridvolume mask */
   void Remove_Mask(void); 
/** auxiliar function to save intermediate results. */
   void Copy_to_Gridvolume_debug(void); 
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
     
     void InitNeighVector(void);   
};
#endif
