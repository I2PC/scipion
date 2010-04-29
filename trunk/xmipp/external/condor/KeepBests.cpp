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

#include "KeepBests.h"
#include <stdlib.h>
#include <string.h>

KeepBests::KeepBests(int _n): n(_n), optionalN(0)
{
    init();
}

KeepBests::KeepBests(int _n, int _optionalN): n(_n), optionalN(_optionalN)
{
    init();
}

void KeepBests::init()
{
    int i;
    double *t;
    ctable=(cell*)malloc(n*sizeof(cell));
    if (optionalN) t=(double*)malloc(optionalN*n*sizeof(double));
    for (i=0; i<n; i++)
    {
        if (optionalN)
        {
            ctable[i].optValue=t;
            t+=optionalN;
        }
        ctable[i].K=INF;
        ctable[i].prev=ctable+(i-1);
    }
    ctable[0].prev=NULL;
    end=ctable+(n-1);
    _local_getOptValueI=-1;
}

void KeepBests::setOptionalN(int _optionalN)
{
    int i;
    double *t;
    if (optionalN) t=(double*)realloc(ctable[0].optValue,_optionalN*n*sizeof(double));
    else t=(double*)malloc(_optionalN*n*sizeof(double));
    for (i=0; i<n; i++) 
    {
        ctable[i].optValue=t;
        t+=_optionalN;
    }
    optionalN=_optionalN;
}

KeepBests::~KeepBests()
{
    if (optionalN) free(ctable[0].optValue);
    free(ctable);
}

void KeepBests::reset()
{
    int i;
    for (i=0; i<n; i++) ctable[i].K=INF;
//    if (optionalN) memset(ctable[0].optValue,0,optionalN*n*sizeof(double));
}

void KeepBests::add(double key, double value)
{
    add(key,value,NULL,0);
}
void KeepBests::add(double key, double value, double optionalValue)
{
    add(key,value,&optionalValue,1);
}

void KeepBests::add(double key, double value, double *optionalValue)
{
    add(key,value,optionalValue,optionalN);
}

void KeepBests::add(double key, double value, double *optionalValue, int nn)
{
    cell *t=end, *prev, *t_next=NULL;
    while ((t)&&(t->K>key)) { t_next=t; t=t->prev; };
    if (t_next)
    {
        if (t_next==end)
        {
            end->K=key;
            end->value=value;
            if ((optionalN)&&(optionalValue))
            {
                memcpy(end->optValue, optionalValue, nn*sizeof(double));
                if (optionalN-nn>0) 
                    memset(end->optValue+nn,0,(optionalN-nn)*sizeof(double));
            }
        } else
        {
            prev=end->prev;
            end->prev=t;
            t_next->prev=end;
    
            end->K=key;
            end->value=value;
            if ((optionalN)&&(optionalValue))
            {
                memcpy(end->optValue, optionalValue, nn*sizeof(double));
                if (optionalN-nn) 
                    memset(end->optValue+nn,0,(optionalN-nn)*sizeof(double));
            }
            end=prev;
        };
    };
}

double KeepBests::getValue(int i)
{
    cell *t=end;
    i=n-i-1;
    while (i) { t=t->prev; i--; }
    return t->value;
}

double KeepBests::getKey(int i)
{
    cell *t=end;
    i=n-i-1;
    while (i) { t=t->prev; i--; }
    return t->K;
}

double KeepBests::getOptValue(int i, int no)
{
    if (i==_local_getOptValueI) return _local_getOptValueC->optValue[no];
    _local_getOptValueI=i;
    cell *t=end;
    i=n-i-1;
    while (i) { t=t->prev; i--; }
    _local_getOptValueC=t;
    return t->optValue[no];
}

double *KeepBests::getOptValue(int i)
{
    cell *t=end;
    i=n-i-1;
    while (i) { t=t->prev; i--; }
    return t->optValue;
}
