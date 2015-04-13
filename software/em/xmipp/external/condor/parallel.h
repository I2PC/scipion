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

#ifndef _PARALLEL_H_
#define _PARALLEL_H_

void parallelImprove(InterPolynomial *p, int *_k, double _rho, double *_valueFk, Vector _Base);
void parallelInit(int _nnode, int _dim, ObjectiveFunction *_of);
void parallelFinish();
void startParallelThread();

#endif

