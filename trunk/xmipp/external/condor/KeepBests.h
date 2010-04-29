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

#ifndef __INCLUDE_KEEPBEST__
#define __INCLUDE_KEEPBEST__

#define INF 1.7E+308 

typedef struct cell_tag
    {
        double K;
        double value;
        double *optValue;
        struct cell_tag *prev;
    } cell;

class KeepBests
{
public:
    KeepBests(int n);
    KeepBests(int n, int optionalN);
    void setOptionalN(int optinalN);
    ~KeepBests();
    void reset();
    void add(double key, double value);
    void add(double key, double value, double optionalValue);
    void add(double key, double value, double *optionalValue);
    void add(double key, double value, double *optionalValue, int nn);
    double getKey(int i);
    double getValue(int i);
    double getOptValue(int i, int n);
    double* getOptValue(int i);
    int sz() {return n;};
private:
    void init();
    cell *ctable,*end,*_local_getOptValueC;
    int n,optionalN,_local_getOptValueI;
};

#endif
