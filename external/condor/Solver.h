/*

CONDOR 1.06 - COnstrained, Non-linear, Direct, parallel Optimization 
              using trust Region method for high-computing load, 
              noisy functions
Copyright (C) 2004 Frank Vanden Berghen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation version 2
of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

If you want to include this tools in any commercial product, 
you can contact the author at fvandenb@iridia.ulb.ac.be

*/

#ifndef _INCLUDE_SOLVER_H
#define _INCLUDE_SOLVER_H

#include "IntPoly.h"
#include "Vector.h"

Vector L2NormMinimizer(Polynomial q, double delta, 
                       int *infoOut=NULL, int maxIter=1000, double *lambda1=NULL);
Vector L2NormMinimizer(Polynomial q, Vector pointXk, double delta, 
                       int *infoOut=NULL, int maxIter=1000, double *lambda1=NULL);
Vector L2NormMinimizer(Polynomial q, Vector pointXk, double delta, 
                       int *infoOut, int maxIter, double *lambda1, Vector minusG, Matrix H);

Vector LAGMAXModified(Polynomial q, double rho,double &VMAX);
Vector LAGMAXModified(Polynomial q, Vector pointXk, double rho,double &VMAX);
Vector LAGMAXModified(Vector G, Matrix H, double rho,double &VMAX);

void CONDOR(double rhoStart, double rhoEnd, int niter, 
            ObjectiveFunction *of, int nnode=0);
#endif

