#ifndef CRISBEND_H
#define CRISBEND_H

///////////////////////////// GENERAL LIBRARIES /////////////////////////
// I/O
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <XmippData/xmippArgs.hh>

//Data 
#include <vector>
#include <string>

///////////////////////////// SPECIFIC LIBRARIES /////////////////////////
//#include "DelTriang.h"
#include <XmippData/xmippDelTriang.hh>
#include <XmippInterface/xmippCCLattice_IO.hh>

/////////////////////////////// DATA TYPES //////////////////////////////

typedef struct{

	int dim[2];
	vector <float> O;
	vector <float> a;
	vector <float> b;
	int Na[2];
	int  Nb[2];
}LatParam;

// Struct to store corelation peak related information
typedef struct{
	//Lattice indexes
	int i,j;
	//Point Position
	float x,y;
	//Deviation from Theoric Position
	float Incrx,Incry;
	//Interpolation flag
	bool Interp;
}LatPoint;


////////////////////////////////  FUNCTIONS //////////////////////////

//////////////////////////  Peaks Correspondance

// Peaks Correspondance 
void PeaksCorresp(CCLattice_IO & , vector <LatPoint> &, double );


/////////////////////////// Interpolation

//Displacement Interpolation
void BendingInterp(vector <LatPoint> &);

/**Displacement Interpolation from Triangulation of Irregular Grid*/
void LinInterp(vector <LatPoint> &, LatPoint &, vector <ITRIANGLE> & );

////////////////////////////  Triangulation

// Lattice Triangulation (using DelTriang.hh)
void LatTriang(vector <LatPoint> &, vector <ITRIANGLE> & );

// Returns index of triangle in vector LatTri
int FindNearest(vector <LatPoint> &, LatPoint &, vector <ITRIANGLE> & );


//////////////////////////// I/O functions

// Read *.cor file. Currently using CCLattice_IO class ReadMRCCord function
// void ReadMRCCord(const char *, LatParam &,vector <float> &, vector <float> &);
//Save Lat. displacements
void  SaveLatBending(const char * ,CCLattice_IO &, vector <LatPoint> &);

#endif
