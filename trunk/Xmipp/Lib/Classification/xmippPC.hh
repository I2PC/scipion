/***************************************************************************
 *
 * Authors:     Jorge García de la Nava Ruiz (gdl@ac.uma.es)
 *              Carlos Oscar Sanchez Sorzano
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

#include <Classification/xmippBaseAlgo.hh>
#include <Classification/xmippCDataTypes.hh>
#include <Classification/xmippCTVectors.hh>
#include <XmippData/xmippFuncs.hh>

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
	void reset(xmippCTVectors const &ts, vector<unsigned> const & idx) _THROW;

	/**
	* The eigenvectors
	*/
	vector<xmippVector> eigenvec;

	/**
	* The eigenvalues
	*/
	xmippVector eigenval;

	/**Set identity matrix as eigenvector matrix
		@param n The number of eigenvectors*/
	void setIdentity(int n);
  	
	/** Defines Listener class
  	*/
  	void setListener(xmippBaseListener* _listener) { listener = _listener; };

private:

  	xmippBaseListener* listener;   // Listener class   


};
#endif

