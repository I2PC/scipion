/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano   (coss@cnb.csic.es)
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

#include "kSVD.h"
#include <vector>

/* ------------------------------------------------------------------------- */
/* Orthogonal Matching Pursuit                                               */
/* ------------------------------------------------------------------------- */
double orthogonalMatchingPursuit(const Matrix1D<double> &x,
    const Matrix2D<double> &D, int S, Matrix1D<double> &alpha)
{
    // Compute approximation error
    double approximationError=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(x)
        approximationError+=x(i)*x(i);
    double errorLimit=approximationError*1e-10;

    // Initialize the number of available atoms
    int N=YSIZE(D); // Dimension of each atom
    int K=XSIZE(D); // Number of atoms in the dictionary
    Matrix1D<int> availableAtoms;
    availableAtoms.initZeros(K);
    availableAtoms.initConstant(1);
    
    // Compute the projection of x onto each of the atoms and look for the
    // maximum
    Matrix1D<double> c, w, e, u;
    c.initZeros(K);
    e.initZeros(K); e.initConstant(1);
    u.initZeros(K); u.initConstant(1);

    // Compute dot product between x and di
    int km=0;
    double ckm=-1;
    for (int k=0; k<K; k++)
    {
        for (int j=0; j<N; j++)
            DIRECT_VEC_ELEM(c,k)+=DIRECT_MAT_ELEM(D,j,k)*DIRECT_VEC_ELEM(x,j);
        
        // Check if this is the maximum
        if (ABS(c(k))>ckm)
        {
            km=k;
            ckm=ABS(c(k));
        }
    }
    
    // Delete Js from the available indexes and update the approximation error
    availableAtoms(km)=0;
    approximationError-=ckm*ckm;
    
    // Annotate Jk
    std::vector<int> J;
    int s=0;
    J.push_back(km);

    // Build the R matrix
    Matrix2D<double> R(S,K);
    R(0,km)=u(km);

    while (s<S-1 && approximationError>errorLimit)
    {
        // Update all the rest of R except Js
        ckm=0;
        int nextKm=-1;
        for (int k=0; k<K; k++)
        {
            // If already used, skip it
            if (!availableAtoms(k)) continue;
        
            // Compute the product between the atoms km and k
            double aux=0;
            for (int j=0; j<N; j++)
                aux+=DIRECT_MAT_ELEM(D,j,km)*DIRECT_MAT_ELEM(D,j,k);
            DIRECT_MAT_ELEM(R,s,k)=aux;
            
            // Readjust R
            for (int n=0; n<s; n++)
                DIRECT_MAT_ELEM(R,s,k)-=DIRECT_MAT_ELEM(R,n,km)*
                    DIRECT_MAT_ELEM(R,n,k);
            if (DIRECT_VEC_ELEM(u,km)>XMIPP_EQUAL_ACCURACY)
                DIRECT_MAT_ELEM(R,s,k)/=DIRECT_VEC_ELEM(u,km);
            
            // Update the rest of variables
            c(k)=c(k)*u(k)-c(km)*R(s,k);
            e(k)-=R(s,k)*R(s,k);
            u(k)=sqrt(ABS(e(k)));
            if (u(k)!=0) c(k)/=u(k);

            // Check if this is the best c(k)
            if (ABS(c(k))>ckm)
            {
                nextKm=k;
                ckm=ABS(c(k));
            }
        }
        
        // Delete km from the available indexes and update the approximation error
        km=nextKm;
        J.push_back(km);
        availableAtoms(km)=0;
        s++;
        R(s,km)=u(km);
        approximationError-=ckm*ckm;
    }

    // Perform backsubstitution
    alpha.initZeros(K);
    for (int k=s; k>=0; k--)
    {
        int Jk=J[k];
        for (int n=s; n>=k+1; n--)
            c(Jk)-=R(k,J[n])*c(J[n]);
        if (R(k,Jk)!=0) c(Jk)/=R(k,Jk);
        alpha(Jk)=c(Jk);
    }
    
    return sqrt(ABS(approximationError));
}
