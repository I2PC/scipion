/***************************************************************************
 *
 * Authors:     Jorge García de la Nava Ruiz (gdl@ac.uma.es)
 *              Carlos Oscar Sanchez Sorzano
 *              Alberto Pascual Montano (pascual@cnb.uam.es)
 *
 * Departamento de Arquitectura de Computadores, Universidad de Málaga
 *
 * Copyright (c) 2001 , CSIC/UMA.
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'pascual@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef XMIPPPC_H
#define XMIPPPC_H

using namespace std;

#include "base_algorithm.h"
#include "data_types.h"
#include "training_vector.h"

#include <data/funcs.h>
#include <data/matrix2d.h>

/**@name PCA classes*/
//@{
/** Basic PCA class */
class xmippPC {
public:

	/**
	* Make an empty xmippPC
	*/
	xmippPC(void){}

	/**
	* Construct a xmippPC object with eigenvectors & eigenvalues
	* from ts.
	* @param ts The vectors.
	*/
	xmippPC(xmippCTVectors const &ts) {reset(ts);}

	/**
	* Construct a xmippPC object with eigenvectors & eigenvalues
	* from ts, using only the items given in idx.
	* @param ts The vectors.
	* @param idx The indexes of the vectors to use
	*/
	xmippPC(xmippCTVectors const &ts, vector<unsigned> const & idx) {reset(ts,idx);}

	/**
	* Calculate the eigenval/vecs
	* @param ts The vectors.
	*/
	void reset(xmippCTVectors const &ts);

	/**
	* Calculate the eigenval/vecs
	* @param ts The vectors.
	* @param idx The indexes of the vectors to use
	*/
	void reset(xmippCTVectors const &ts, vector<unsigned> const & idx);

	/**
	* The eigenvectors
	*/
	vector<xmippVector> eigenvec;

	/**
	* The eigenvalues
	*/
	xmippVector eigenval;

	/**
	* Mean of the input training set
	*/
        xmippVector mean;

        /** Number of relevant eigenvectors */
	int D;

        /** Product <mean,mean> */
	double prod_mean_mean;
	
	/** Products <ei,mean> */
	vector<double> prod_ei_mean;
	
	/** Products <ei,ei> */
	vector<double> prod_ei_ei;

      	/** Average of the mean vector */
	double avg_mean;
	
	/** Average of the ei vectors */
	vector<double> avg_ei;

	/**Set identity matrix as eigenvector matrix
		@param n The number of eigenvectors*/
	void setIdentity(int n);
  	
        /** Clear.
	    Clean the eigenvector, eigenvalues and D */
        void clear();
	
      	/** Prepare for correlation.
	    This function computes the inner products emong the ei,
	    and the mean vector, among the ei, and the average of all
	    vectors*/
        void prepare_for_correlation();

        /** Set the number of relevant eigenvalues. */
	void set_Dimension(int _D) {D=_D;}

        /** Get the number of relevant eigenvalues. */
	int get_Dimension() const {return D;}

      	/** Get the dimension of the eigenvectors. */
	int get_eigenDimension() const {
	   if (eigenvec.size()>1) return eigenvec[0].size();
	   else return 0;
	}

      	/** Number of components for a given accuracy explanation.
	    This function returns the number of components to be taken
	    if th_var% of the variance should be explained */
	int Dimension_for_variance(double th_var);

        /** Project onto PCA space.
	    Given a vector of the same size as the PCA vectors this function
	    returns the projection of size D onto the first D eigenvectors.
	    D is set via the set_Dimension function
	
	    An exception is thrown if the input vectors are not of the same size
	    as the PCA ones.*/
        void Project(xmippVector &input, xmippVector &output);

	/** Defines Listener class
  	*/
  	void setListener(xmippBaseListener* _listener) { listener = _listener; };

	/** Show relevant eigenvectors and eigenvalues */
	friend ostream& operator << (ostream &out, const xmippPC &PC);

	/** Read a set of PCA just as shown */
	friend istream& operator >> (istream &in, xmippPC &PC);

private:

  	xmippBaseListener* listener;   // Listener class


};


/** Set of PCA classes */
class PCA_set {
public:
   /** Set of PCA analysis. */
   vector<xmippPC *> PCA;

public:
   /** Destructor */
   ~PCA_set();

   /** Create empty PCA.
       Creates space for n new PCA and returns the index of the first one */
   int create_empty_PCA(int n=1);

   /** Returns the number of PCA analysis.*/
   int PCANo() const {return PCA.size();}

   /** Returns a pointer to PCA number i*/
   xmippPC * operator ()(int i) const {return PCA[i];}

   /** Show all PCA */
   friend ostream& operator << (ostream &out, const PCA_set &PS);

   /** Read a set of PCA just as shown */
   friend istream& operator >> (istream &in, PCA_set &PS);
};

/** Running PCA.
    Running PCA is an algorithm that estimates iteratively the
    principal components of a dataset that is provided to the algorithm
    as vectors are available.

    See J. Weng, Y. Zhang, W.S. Hwang. Candid covariance-free incremental
    principal component analysis. IEEE Trans. On Pattern Analysis and
    Machine Intelligence, 25(8): 1034-1040 (2003).
*/
class Running_PCA {
public:
   /// Total number of eigenvectors to be computed.
   int J;

   /// Dimension of the sample vectors
   int d;

   /// Current estimate of the population mean
   matrix1D<double> current_sample_mean;

   /// Current number of samples seen
   long n;

   /** Current estimate of the eigenvectors.
       Each column is an eigenvector. */
   matrix2D<double> eigenvectors;

   /** Constructor.
       J is the number of eigenvectors to compute. d is the
       dimension of the sample vectors. */
   Running_PCA(int _J, int _d);

   /** Update estimates with a new sample. */
   void new_sample(const matrix1D<double> &sample);

   /** Project a sample vector on the PCA space. */
   void project(const matrix1D<double> &input, matrix1D<double> &output) const;

   /// Get a certain eigenvector.
   void get_eigenvector(int j, matrix1D<double> &result) const {
      eigenvectors.getCol(j,result);
   }

   /// Get the variance associated to a certain eigenvector.
   double get_eigenvector_variance(int j) const {
      if (n<=j) return 0.0;
      double mean_proj=sum_proj(j)/n;
      return sum_proj2(j)/n-mean_proj*mean_proj;
   }
public:
   // Sum of all samples so far
   matrix1D<double> sum_all_samples;

   // Sum of all projections so far
   matrix1D<double> sum_proj;

   // Sum of all projections squared so far
   matrix1D<double> sum_proj2;
};
//@}
#endif

