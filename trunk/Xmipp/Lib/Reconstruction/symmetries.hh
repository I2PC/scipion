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
/* ------------------------------------------------------------------------- */
/* SYMMETRIES                                                                */
/* ------------------------------------------------------------------------- */
#ifndef _SYMMETRIES_HH
#  define _SYMMETRIES_HH

//define the  crystallographic groups
//symmetry matices from http://cci.lbl.gov/asu_gallery/
#define sym_undefined    -1
#define sym_P1            1
#define sym_P2            3//c-axis, b-axis
#define sym_P2_1          4
#define sym_C2           12
#define sym_P222         16
#define sym_P222_1       17
//#define sym_P22_12       17
#define sym_P22_12_1     18
#define sym_P4           75
#define sym_P422         89
#define sym_P42_12       90
#define sym_P3          143
#define sym_P312        150
#define sym_P6          168
#define sym_P622        177
// This constant is obsolete and should be revised with the use of
// matrix euality accuracy
#define SYM_ACCURACY 1e-6

#include <XmippData/xmippMatrices2D.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippArgs.hh>
#include <Reconstruction/grids.hh>
//#include <Reconstruction/Programs/Prog_art_crystal.hh>

/**@name Symmetry lists
    The symmetry lists are, simply, lists of 2D matrices. It's the way of
    taking symmetry into account in the reconstruction programs. The
    symmetry list must contain matrices which express equivalent views to
    the actual one due to the underlying volume symmetry. The identity matrix
    is not within the list. You know that symmetry matrices should form
    a subgroup, when reading a file the subgroup is automatically computed
    and, when you add or remove a new matrix, the subgroup must be
    manually computed.
*/
//@{
/** Number of an image in the reconstruction list.
    This macro returns the index of a symmetry image (after the symmetry matrix
    number sym_no) within a list where the first images are true images and the
    last ones, the symmetrized copies (all copies of a same image are
    together). The total number of real images is numIMG, and i is the index
    within this first numIMG images of the image we want to symmetrize The
    first image in the list is the number 0 */
#define SYMINDEX(SL, sym_no, i, numIMG) \
   numIMG+SL.__L.ydim/4*i+sym_no

/** Symmetry List class.
    Internally the symmetry list class is implemented as a single 2D matrix,
    where every 4 rows (remember that in 3D the geometrical transformation
    matrices are 4x4) comprise a symmetry matrix. Access, and ways to modify
    the symmetry list are supplied. Remind that any symmetry is expressed
    in terms of two matrices L and R, so that any Euler matrix must be
    transformed by L*Euler*R resulting into a new perspective of the volume
    which is equivalent to the original one.
    
    The typical use of the symmetry lists is to read the symmetry file, and
    do nothing else but reading matrices from it.
    
    The symmetry file format is
    \begin{verbatim}
    #This is a comment
    # The following line is a 6-fold rotational symmetry axis along Z-axis.
    # The fold is the number of times that the volume can be rotated along
    # the symmetry axis giving the same view from different view points.
    # the structure for the rotational axis is
    # rot_axis      <fold> <X0> <Y0> <Z0>
    # mirror_plane         <X0> <Y0> <Z0>
    rot_axis      6 0 0 1
    mirror_plane    0 0 1
    \end{verbatim}
*/
class SymList {
private:
   // Crystallographic space group. This is only a guess based on angles.
   // No check on vectors magnitude is made
   int		space_group;
public:
   // L and R matrices
   matrix2D<double> __L, __R;
   matrix2D<double> __shift;  // It is used for crystallographic symmetries
   matrix1D<int>    __chain_length;

   // As the symmetry elements form a subgroup, this is the number of
   // true symmetry elements belonging to the list, the rest of
   // the list are simply the elements to fill the subgroup
   int true_symNo;

   // Number of Axis, mirrors, ...
   int              __sym_elements;
   
public:
   /** Create an empty list.
       The 2D matrices are 0x0.
       \\ Ex: SymList SL; */
   SymList() {__sym_elements=true_symNo=space_group=0;}
   
   /** Create Symmetry List from a Symmetry file.
       All the subgroup elements are computed automatically.
       \\ Ex: SymList SL("sym.txt"); */
   SymList(FileName fn_sym, double accuracy=SYM_ACCURACY)
       {read_sym_file(fn_sym,accuracy);}

   /** Get matrices from the symmetry list.
       The number of matrices inside the list is given by \Ref{SymsNo}.
       This function return the 4x4 transformation matrices associated to
       the one in the list which occupies the position 'i'. The matrix
       numbering within the list starts at 0. The output transformation
       matrices is given as a pointer to gain speed.
       \\ Ex:
       \begin{verbatim}
          for (i=0; i<SL.SymsNo; i++) {
              SL.get_matrices(i,L,R);
              ...
          }
       \end{verbatim} */
   void get_matrices(int i, matrix2D<double> &L, matrix2D<double> &R) const;

   /** Set a couple of matrices in the symmetry list.
       The number of matrices inside the list is given by \Ref{SymsNo}.
       This function sets the 4x4 transformation matrices associated to
       the one in the list which occupies the position 'i'. The matrix
       numbering within the list starts at 0.
       \\ Ex:
       \begin{verbatim}
          for (i=0; i<SL.SymsNo; i++) {
              SL.set_matrix(i,L,R);
              ...
          }
       \end{verbatim} */
   void set_matrices(int i, const matrix2D<double> &L,
      const matrix2D<double> &R);

   /** Get shift.
       Returns the shift associated to a certain symmetry. */
   void get_shift(int i, matrix1D<double> &shift) const;

   /** Set shift.
       Set the shift associated to a certain symmetry. */
   void set_shift(int i, const matrix1D<double> &shift) _THROW;

   /** Add shift.
       Add a shift vector to the shift matrix. An exception is thrown if
       the input vector is not a 3x1 vector.*/
   void add_shift(const matrix1D<double> &shift) _THROW;

   /** Read a symmetry file into a symmetry list.
       The former symmetry list is overwritten with the new one. All the
       subgroup members are added to the list. If the accuracy is negative
       then the subgroup is not generated.
       \\ Ex: SL.read_sym_file("sym.txt");*/
   void read_sym_file(FileName fn_sym, double accuracy=SYM_ACCURACY);

   /** Add symmetry matrices to the symmetry list.
       The given matrix must specify a point of view equivalent to the
       actual point of view. The matrices are added to the subgroup generator
       but the subgroup is not updated, you must do it manually using
       \Ref{compute_subgroup}. What is more, the subgroup after the insertion
       is corrupted.
       
       The chain length is the number of single matrices multiplication of
       which the inserted one is compound.*/
   void add_matrices(const matrix2D<double> &L, const matrix2D<double> &R,
       int chain_length) _THROW;

   /** Compute subgroup for this structure.
       After adding or setting a matrix, the subgroup information
       is lost, you must recalculate it using this function. The different
       matrices are multiplied until no more different matrices are produced.
       The accuracy is used in order to compare when two matrix elements are
       the same.
       
       So far, all the shifts associated to generated matrices are set to 0*/
   void compute_subgroup(double accuracy=SYM_ACCURACY);

   /** Number of symmetry matrices inside the structure.
       This is the number of all the matrices inside the subgroup.
       \\ Ex:
       \begin{verbatim}
          for (i=0; i<SL.SymsNo; i++) {
              SL.get_matrix(i,A);
              ...
          }
       \end{verbatim} */
   int SymsNo() const {return __L.ydim/4;}

   /** Number of symmetry matrices which generated the structure.
       This is the number of the matrices which generated the structure,
       notice that it should be always less or equal to the total number
       of matrices in the subgroup. */
   int TrueSymsNo() const {return true_symNo;}
   
   /** Guess Crystallographic space group.  
       Return the  \URL[space group]{
       http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen} number. So
       far it has only been implemented for P1 (1), P2212 (17), P4 (75), 
       P4212 (90) and P6 (168).
       
       Mag_a and Mag_b are the crystal vector magnitude. */
       
   int  crystallographic_space_group (double mag_a, 
                                           double mag_b,
					   double ang_a2b_deg) const;

}; 

/** Applies to the crystal vectors de n-th symmetry  matrix, It also
   inizializates the shift vector. The crystal vectors and the basis must be
   the same  except for a constant!!
   A note: Please realize that we are not repeating code here.
   The class SymList deals with symmetries when expressed in 
   Cartesian space, that is the basis is orthonormal. Here
   we describe symmetries in the crystallographic way
   that is, the basis and the crystal vectors are the same.
   For same symmetries both representations are almost the same
   but in general they are rather different.
 */


    void symmetrize_crystal_vectors(matrix1D<double> &aint, 
			          matrix1D<double> &bint,
				  matrix1D<double> &shift,
				  int space_group,
				  int sym_no,
				  const matrix1D<double> &eprm_aint,
				  const matrix1D<double> &eprm_bint);   

/** Symmetrizes a crystal volume.
 */

    void symmetrize_crystal_volume(GridVolume &vol,
                              const matrix1D<double> &eprm_aint, 
			      const matrix1D<double> &eprm_bint,
			      int eprm_space_group,const matrix2D<int> &mask,
			      int grid_type );

/** Symmetrizes a simple grid with P2221  symmetry
*/
     void symmetry_P222_1(Volume &vol, const SimpleGrid &grid,
                                const matrix1D<double> &eprm_aint, 
			        const matrix1D<double> &eprm_bint,
				const matrix2D<int> &mask, int volume_no,
				int grid_type);
/** Symmetrizes a simple grid with P4  symmetry
*/
     void symmetry_P4(Volume &vol, const SimpleGrid &grid,
                                const matrix1D<double> &eprm_aint, 
			        const matrix1D<double> &eprm_bint,
				const matrix2D<int> &mask, int volume_no,
				int grid_type);
/** Symmetrizes a simple grid with P4212 symmetry
*/
     void symmetry_P42_12(Volume &vol, const SimpleGrid &grid, 
                                const matrix1D<double> &eprm_aint, 
			        const matrix1D<double> &eprm_bint,
				const matrix2D<int> &mask, int volume_no,
				int grid_type);
/** Symmetrizes a simple grid with P6 symmetry
*/
     void symmetry_P6(Volume &vol, const SimpleGrid &grid,
                                const matrix1D<double> &eprm_aint, 
			        const matrix1D<double> &eprm_bint,
				const matrix2D<int> &mask, int volume_no,
				int grid_type);

//@}

#endif
