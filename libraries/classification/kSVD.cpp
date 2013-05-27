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
 *  e-mail address 'xmipp@cnb.csic.es'                                  
 ***************************************************************************/

#include "kSVD.h"

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
    int N=D.Ydim(); // Dimension of each atom
    int K=D.Xdim(); // Number of atoms in the dictionary
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
            c(k)+=D(j,k)*x(j);
        
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
                aux+=D(j,km)*D(j,k);
            R(s,k)=aux;
            
            // Readjust R
            for (int n=0; n<s; n++)
                R(s,k)-=R(n,km)*R(n,k);
            if (u(km)>XMIPP_EQUAL_ACCURACY)
                R(s,k)/=u(km);
            
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

/* ------------------------------------------------------------------------- */
/* Lasso (shooting)                                                          */
/* ------------------------------------------------------------------------- */
double lasso(const Matrix1D<double> &x,
    const Matrix2D<double> &D, 
    const Matrix2D<double> &DtD,
    const Matrix2D<double> &DtDlambdaInv,
    double lambda, Matrix1D<double> &alpha,
    const int maxIter, const double tol)
{
    int K=D.Xdim();

    // Compute the ridge least squares solution
    // Compute D^t*x
    Matrix1D<double> Dtx(K);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(Dtx)
        // Compute the dot product of the i-th column of D and x
        for (size_t k=0; k<D.Ydim(); k++)
            Dtx(i)+= D(k,i)*x(k);

    // Now multiply DtdLambdaInv * D^t * x
    alpha.initZeros(K);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(DtDlambdaInv)
        alpha(i)+=DtDlambdaInv(i,j)*Dtx(j);

    // Iterate over alpha
    Matrix1D<double> alphaOld;
    int iter=0;
    bool finished=false;
    while (!finished)
    {
        // Keep the current alpha
        alphaOld=alpha;
        
        // Update alpha
        FOR_ALL_ELEMENTS_IN_MATRIX1D(alpha)
        {
            double S=-Dtx(i);
            for (size_t j=0; j<alpha.size(); j++)
                if (i!=j) S+=DtD(i,j)*alpha(j);
            if (fabs(S)<lambda)
                alpha(i)=0;
            else if (S>lambda)
                alpha(i)=(lambda-S)/DtD(i,i);
            else
                alpha(i)=(-lambda-S)/DtD(i,i);
        }
        
        // Prepare for next iteration
        iter++;
        double normDiff=0, normAlpha=0;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(alpha)
        {
            double diff=alpha(i)-alphaOld(i);
            normDiff+=diff*diff;
            normAlpha+=alpha(i)*alpha(i);
        }
        finished=(iter>=maxIter) || normDiff/normAlpha<tol;
    }
    
    // Compute the approximation error
    Matrix1D<double> xp(x.size());
    for (size_t j=0; j<D.Xdim(); j++)
        if (alpha(j)!=0)
            for (size_t i=0; i<D.Ydim(); i++)
                xp(i)+=D(i,j)*alpha(j);
    double approximationError=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(xp)
    {
        double diff=xp(i)-x(i);
        approximationError+=diff*diff;
    }
    
    return sqrt(ABS(approximationError));
}

/* ------------------------------------------------------------------------- */
/* kSVD                                                                      */
/* ------------------------------------------------------------------------- */
//#define DEBUG
double kSVD(const std::vector< Matrix1D<double> > &X, int S,
    Matrix2D<double> &D, std::vector< Matrix1D<double> > &Alpha,
    bool keepFirstColumn, int maxIter, double minChange,
    int projectionMethod, double lambda)
{
	size_t Nvectors=X.size();
	size_t K=D.Xdim(); // Number of atoms
	size_t N=D.Ydim(); // Dimension of the atoms

    // Ask for memory for Alpha if necessary
    if (Alpha.size()!=Nvectors)
    {
        int Nalpha=Alpha.size();
        Matrix1D<double> dummy;
        dummy.initZeros(K);
        for (size_t i=0; i<Nvectors-Nalpha; i++)
            Alpha.push_back(dummy);
    }

    // Ask for memory for which vectors use which atoms
    std::vector< std::vector<int> > listUsers;
    for (size_t k=0; k<K; k++)
    {
        std::vector<int> dummy;
        listUsers.push_back(dummy);
    }

    // Compute the power of the input vectors
    Matrix1D<double> power;
    power.initZeros(Nvectors);
    for (size_t n=0; n<Nvectors; n++)
        FOR_ALL_ELEMENTS_IN_MATRIX1D(X[n])
            power(n)+=X[n](i)*X[n](i);
    double avgPower=power.sum(true);

    // Perform kSVD
    int iterations=0;
    Matrix1D<double> error;
    error.initZeros(Nvectors);
    double previousError=0, currentError=0;
    bool limitAchieved=false;
    do {
        // Sparse coding step
        Matrix2D<double> DtD, DtDlambda, DtDlambdaInv;
        if (projectionMethod==LASSO_PROJECTION)
        {
            DtD.initZeros(K,K);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(DtD)
                // Compute the dot product between the columns i and j of D
                for (size_t k=0; k<D.Ydim(); k++)
                    DtD(i,j)+=D(k,i)*D(k,j);
            DtDlambda=DtD;
            for (size_t i=0; i<K; i++)
                DtDlambda(i,i)+=lambda;
            DtDlambda.inv(DtDlambdaInv);
        }
        std::cout << "Sparse coding\n";
        init_progress_bar(Nvectors);
        for (size_t n=0; n<Nvectors; n++)
        {
            // Compute alpha
            if (projectionMethod==LASSO_PROJECTION)
                error(n)=lasso(X[n],D,DtD,DtDlambdaInv,lambda,Alpha[n]);
            else
                error(n)=orthogonalMatchingPursuit(X[n],D,S,Alpha[n]);
            
            // Check which are the atoms this vector is using
            for (size_t k=0; k<K; k++)
                if (Alpha[n](k)!=0)
                    listUsers[k].push_back(n);
            if (n%100==0) progress_bar(n);
        }
        progress_bar(Nvectors);
        double Noise=error.sum(true);
        std::cout << "kSVD error=" << Noise << " "
                  << " (" << 100*Noise/avgPower << "%)" << std::endl;
        if (projectionMethod==LASSO_PROJECTION)
        {
            int countNonZeros=0;
            for (size_t n=0; n<Nvectors; n++)
                FOR_ALL_ELEMENTS_IN_MATRIX1D(Alpha[n])
                    if (Alpha[n](i)!=0) countNonZeros++;
            std::cout << "Average sparsity = " << (double)countNonZeros/
                Nvectors << std::endl;
        }
        
        // Codebook update
        Matrix2D<double> EkR, U, V;
        Matrix1D<double> xp, W;
        int firstKtoUpdate=0;
        if (keepFirstColumn) firstKtoUpdate=1;
        std::cout << "Updating codebook\n";
        init_progress_bar(K);
        for (size_t k=firstKtoUpdate; k<K; k++)
        {
            // Compute the error that would be commited if the
            // atom k were not used
        	size_t Nk=listUsers[k].size();
            // std::cout << "Atom k=" << k << " is used by " << Nk << " vectors\n";
            if (Nk>1)
            {
                EkR.initZeros(N,Nk);
                for (size_t nk=0; nk<Nk; nk++)
                {
                    // Select vector nk
                    int n=listUsers[k][nk];

                    // Compute the represented vector if atom is not used
                    xp.initZeros(N);
                    for (size_t kp=0; kp<K; kp++)
                    {
                        double w=Alpha[n](kp);
                        if (kp!=k && w!=0)
                            for (size_t j=0; j<N; j++)
                                xp(j)+=w*D(j,kp);
                    }

                    // Fill the error matrix EkR
                    const Matrix1D<double> &x=X[n];
                    for (size_t j=0; j<N; j++)
                        EkR(j,nk)=x(j)-xp(j);
                }

                // Compute the SVD decomposition of the EkR matrix
                svdcmp(EkR,U,W,V); // EkR=U * diag(W) * V^t

                // Update the dictionary with the first column of U
                double normU0=0;
                for (size_t j=0; j<N; j++)
                {
                    double uj0=U(j,0);
                    normU0+=uj0*uj0;
                }
                normU0=sqrt(normU0);
                double inormU0=1/normU0;
                for (size_t j=0; j<N; j++)
                    D(j,k)=U(j,0)*inormU0;

                // Update the coefficients of the users of atom k
                double W0=W(0)*normU0;
                for (size_t nk=0; nk<Nk; nk++)
                {
                    // Select vector nk
                    int n=listUsers[k][nk];
                    Alpha[n](k)=V(nk,0)*W0;
                }

                #ifdef DEBUG
                    for (int n=0; n<Nvectors; n++)
                    {
                        Matrix1D<double> xp=D*Alpha[n];
                        std::cout << "x= " << X[n].transpose() << std::endl
                                  << "xp=" << xp.transpose() << std::endl
                                  << "diff norm=" << (X[n]-xp).module() << std::endl;
                    }
                #endif
            }
            else if (Nk==0)
            {
                // Pick the vector that is worse represented
                int iworse;
                error.maxIndex(iworse);
                error(iworse)=0;
                for (size_t j=0; j<N; j++)
                    D(j,k)=X[iworse](j);
            }
            progress_bar(k);
        }
        progress_bar(K);
        
        // Prepare for next iteration
        iterations++;
        currentError=error.sum(true);
        if (previousError==0) previousError=currentError;
        else {
            limitAchieved=
                ABS((previousError-currentError)/previousError)<minChange;
            std::cout << "Change=" << ((previousError-currentError)/previousError) << std::endl;
            previousError=currentError;
        }
    } while (iterations<maxIter && !limitAchieved);

    return currentError;
}
#undef DEBUG
