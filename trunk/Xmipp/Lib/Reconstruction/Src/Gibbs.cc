/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

#include "../Gibbs.hh"

template <class T> void
FCC_Gibs<T>::Init(double relative_size, double R, 
                                     FileName Random_Number_File, 
                                    InitMode initmode) /*_THROW*/
{
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

}// end FCC_Grids::Init()
//-----------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::InitNeighVector(void)

{
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
   

}// end InitNeighVector
//--------------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::GenerateValidClicks(void)
{
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
}//GenerateValidClicks

//-------------------------------------------------------------------------
//mask out all those points OUTSIDE an sphere of radius R
//------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::CreateMask(void)
{
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

//--------------------------------------------------------
//Mask out all those points that do not have a full neighbourhood
//------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::Create_Second_Mask(void)//z y x
{
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


}//end SecondMask

//--------------------------------------------------------
//MAX of randompool should be iNumber_of_points_inside_sphere
template <class T>
void  FCC_Gibs<T>::Alloc_and_Fill_valid_coordinates_vector(void)
{
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

}//Alloc_and_Fill_valid_coordinates_vector

//--------------------------------------------------------
//--------------------------------------------------------
template <class T>
int FCC_Gibs<T>::WhatClickIsThis( int k, int i, int j)//z y x
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

}//end WhatClickIsThis

//--------------------------------------------------------
template <class T>
int FCC_Gibs<T>::IsThisValidClick( int k, int i, int j)//z y x
{
    
int ii=WhatClickIsThis(  k, i, j);
if(ii == FCC_Gibs<T>::OUTSIDE)
   return (FCC_Gibs<T>::OUTSIDE);
else   
   return( Click_Table[ii] );//either border, homogeneous or nonvalid

}//end IsThisValidClick

//-----------------------------------------------------------------------
// Fill auxiliar Grid volume with the value of the clicks centered on 
//    each point (of the principal grid volume 
//-----------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::FillauxGridVolumewithClicks(void){

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

}// end FillauxGridVolumewithClicks

//-----------------------------------------------------------------------
// count clicks of each type 
//-----------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::Count_clicks_of_each_type(void){
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

}//end Count_clicks_of_each_type
 
//-----------------------------------------------------------------------
// Create a table (of size  DIF_CONF x NEIGHBORS_IN) with the value of the
//    new click when in the klick k the voxel v is changed
//-----------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::ClickChangeDueSpellChange(void){
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

}// end ClickChangeDueSpellChange
//-----------------------------------------------------------------------
// Create a table (of size  DIF_CONF x NEIGHBORS_IN) with the value of the
//    new click when in the klick k the voxel v is changed
//-----------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::EnergyChangeDueSpellChange(void){
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

}// end EnergyChangeDueSpellChange


//--------------------------------------------------------
//if I change the spell s in the volume v what is the energy change
//------------------------------------------------------------
template <class T>
int FCC_Gibs<T>::Get_Energy_Change(int k, int i, int j)//z y x
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

}//end Get_Energy_Change

//--------------------------------------------------------
//if I change the spell s in the volume v what is the energy change
//------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::Set_Click_Change(int k, int i, int j)//z y x
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
}//end Set_Click_Change


template <class T>
void FCC_Gibs<T>::Set_Spell_Change(int z, int y, int x)//z y x
{
          VOLVOXEL( *(_FCC),z,y,x) ^= 1;
}//end Set_Spell_Change
//this application should maximice energy

//-----------------------------------------------------------------------
//Count numer of Spells with value n
//-----------------------------------------------------------------------
// count clicks of each type
template <class T>
int FCC_Gibs<T>::Count_Number_Spells_Value_n(T n){
  int n_counter=0;
  
  FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
     if ( (VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j)&PASSMASKBITS)==n
          && (VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j)&MASKBIT)!=MASKBIT)
           n_counter++;
  return(n_counter);           

}//end Update_Number_Spells_Value_n

//-----------------------------------------------------------------------
//Remove Mask 
//-----------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::Remove_Mask(void){
  
    FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
       VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j) &= 0x1;

}//end Remove_Mask
//-----------------------------------------------------------------------
//Copy to gridvolume_debug and Remove Mask 
//-----------------------------------------------------------------------
//my apologies,  instead of this function a should have define
// the = function in gridvolume.

template <class T>
void FCC_Gibs<T>::Copy_to_Gridvolume_debug(void){
  
  _aux_FCC = &gridvolume_debug(0);
       FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
          VOLVOXEL((*FCC_Gibs<T>::_aux_FCC),k,i,j) =
          VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j) & 0x1;

}//end Copy_to_Gridvolume_debug_and_Remove_Mask
#ifdef NEVERDEFINED
//-----------------------------------------------------------------------
//Count number of Spells 
//-----------------------------------------------------------------------
template <class T>
void FCC_Gibs<T>::Count_Number_Spells(void){

  iGridPointNumber=0;
          FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*_FCC))
             if ( (VOLVOXEL((*FCC_Gibs<T>::_FCC),k,i,j)&MASKBIT)!=MASKBIT)
                   iGridPointNumber++;

}//end Count_Number_Spells
 #endif
//------------------------------------------------------------------- 


void Instanciate_FCC_Gibs()
{
double relative_size;
double R=1.;

FCC_Gibs<int> fcc_gibs;
   fcc_gibs.Init(1., R, "m",(FCC_Gibs<int>::InitMode) 1);
   fcc_gibs.IsThisValidClick(0,0,0);
   fcc_gibs.WhatClickIsThis(0,0,0);
   fcc_gibs.Count_clicks_of_each_type();
   fcc_gibs.Get_Energy_Change(0,0,0);//z y x
   fcc_gibs.Set_Click_Change(0,0,0);//z y x
   fcc_gibs.Set_Spell_Change(0,0,0);//z y x
   fcc_gibs.Count_Number_Spells_Value_n(1);
   fcc_gibs.Remove_Mask();
   fcc_gibs.Copy_to_Gridvolume_debug();
   fcc_gibs.Save_debug_gridvolume();
}




