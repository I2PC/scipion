/***************************************************************************
 *
 * Author:     Monica Chagoyen          monica@b12sg1.cnb.uam.es
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
/****************************************************************************/
/* Program for finding the center of an image. The original in Fortran was  */
/* written by A. Santisteban. For any information about the program, contact*/
/* him. Translated into C by Juan P. Secilla (MSC)  Jun/86		    */
/****************************************************************************/
/****************************************************************************/
/* Created a wrapper to fit NewXmipp 					    */
/* Alberto Pascual, October 2001 					    */
/* pascual@cnb.uam.es 					    		    */
/****************************************************************************/

/********************** Include's and macro definitions *********************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))

/************************* Global variables *********************************/

extern float coseno [1281];
extern unsigned char *imagen[513];
extern float r1, r2, r3, cxp1, cxp2, cxp3, cxm1, cxm2, cxm3, rh, xc0, yc0, del,
      rbajo, ralto, zzz0;
extern int ir, m, mu, mu1, mu4, ntot, ncic, mt, idz, ncic2, ntot4, in, ni, largo,
      lancho, indmul;
/*
extern float conv1x (double, double);
void busca (), suaviza (), lectur (), ergrot (double, double, float *);
*/
extern void busca (), suaviza (), ergrot (double, double, float *);


/****************************************************************************/

void busca ()

{
    static float f[22][22];
    float pi = 3.141592653;
    short n1, in1, mu5, i, mu3, ind, ind1, indice, iflag, j, ix, iy, k;
    float h, th, costh, sinth, fi, anc, xk0, yk0, hpi, x, y, z, xx, yy;
    float b, g, d1, e1, ext;

    suaviza();
/*    printf("      x         y        funci¢n    delta  radio M x.\n");*/
    in=min (in, 10);
    mu1=mu-1;
    m=2*mu;
    mu4=mu*ncic;
    mt=2*m;
    ntot=mt*ncic;
    ntot4=ntot+mu4+ncic;
    zzz0=ntot;
    ncic2=2*ncic;
    in1=in+1;
    h=2.*pi/ntot;
    idz=2-ntot;
    th=ncic*h;
    costh = cos (th);
    sinth = sin (th);
    b=2./(th*th)*(1.+costh*costh-sin(2.*th)/th);
    g=4./(th*th)*(sinth/th-costh);
    d1=2.*th/45.;
    e1=d1*sinth*2.;
    d1*=sin(2.*th);
    mu5=mu4-1;
    for (i = 1; i <= mu5; i++)
    {	fi=i*h;
	coseno[i]=sin(fi);
    }
    coseno[mu4]=1.;
    mu3=2*mu4;
    for (i = 1; i<= mu5; i++)
    {	coseno[mu3-i]=coseno[i];
	coseno[mu3+i]=-coseno[i];
	coseno[ntot-i]=-coseno[i];
	coseno[ntot+i]=coseno[i];
    }
    coseno[mu3]=0.;
    coseno[mu3+mu4]=-1.;
    coseno[ntot]=0.;
    coseno[ntot+mu4]=1.;
    ind=2*in+1;
    ind1=ind-1;
    indice=0;
    iflag=0;
e9: if(indice >= ni)
	goto e10;
    if(iflag == 2)
	goto e22;
    anc=(int)(3.+in*del);
    xk0=xc0;
    yk0=yc0;
    rh=min (rh, r2);
    rh=min(rh,xk0-anc-1.);
    rh=min(rh,yk0-anc-1.);
    rh=min(rh, largo-xk0-anc);
    rh=min(rh, lancho-yk0-anc);
    ir= (int) ((rh-r1)/r3+1);
    hpi=h/2./pi/ir;
    cxp1=2.*b*e1*hpi;
    cxp2=-2.*g*d1*hpi;
    cxp3=2.*(g*e1-b*d1)*hpi;
    cxm1=(b*b+e1*e1)*hpi;
    cxm2=(g*g+d1*d1)*hpi;
    cxm3=2.*(b*g-d1*e1)*hpi;
/*    if(iflag == 1) goto 21					       */
    if(iflag == 2)
	goto e22;
    x=xc0-in*del;
    for (i = 1; i <= ind; i++) /* do 5 */
    {	y=yc0-in*del;
	for (j = 1; j <= ind; j++) /*do 4 */
	{   ergrot(x,y,&z);
/*	    printf ("%10.3f%10.3f%10.5f%10.3f%10.3f AA\n", x,y,z,del,rh);*/
            printf(".");
            fflush(stdout);
	    z=100.+indmul*z;
	    f[i][j]=z;
	    y+=del;
	}
	x+=del;
    }
e23:ext=-1000000.;
    for (i = 1; i <= ind; i++) /* do 7 */
	for (j = 1; j <= ind; j++) /* do 7 */
	{   if(f[i][j] > ext)
	    {	ix=i;
		iy=j;
		ext=f[i][j];
	    }
	}
    xc0=xc0+(ix-in-1)*del;
    yc0=yc0+(iy-in-1)*del;
    if(ix == 1 || ix == ind || iy == 1 || iy == ind)
	goto e8;
    del/=in;
    iflag=2;
    indice++;
    goto e9;
e8: iflag=1;
    goto e9;
e22:f[1][1]=f[ix-1][iy-1];
    f[ind][1]=f[ix+1][iy-1];
    f[1][ind]=f[ix-1][iy+1];
    f[ind][ind]=f[ix+1][iy+1];
    f[1][in1]=f[ix-1][iy];
    f[in1][1]=f[ix][iy-1];
    f[ind][in1]=f[ix+1][iy];
    f[in1][ind]=f[ix][iy+1];
    f[in1][in1]=f[ix][iy];
    x=xc0-(in-1)*del;
    for (i = 2; i < ind1; i++) /* do 11 */
    {	y=yc0-(in-1)*del;
	for (j = 2;j <= ind1; j++) /* do 12 */
	{   if(i == in1 && j == in1)
		 goto e13;
	    ergrot(x,y,&z);
/*	    printf ("%10.3f%10.3f%10.5f%10.3f BB \n", x,y,z,del);*/
            printf(".");
            fflush(stdout);
	    z=100.+indmul*z;
	    f[i][j]=z;
e13:	    y+=del;
	}
      x+=del;
    }
    x=xc0-(in-1)*del;
    y=yc0-(in-1)*del;
    for (k = 2; k <= ind1; k++) /* do 16 */
    {	if(k == in1)
	    goto e17;
	xx=xc0-in*del;
	ergrot(xx,y,&z);
/*	printf ("%10.3f%10.3f%10.5f%10.3f CC \n", xx,y,z,del);*/
        printf(".");
        fflush(stdout);
	z=100.+indmul*z;
	f[1][k]=z;
	yy=yc0-in*del;
	ergrot(x,yy,&z);
/*	printf ("%10.3f%10.3f%10.5f%10.3f DD \n", xx,y,z,del);*/
        printf(".");
        fflush(stdout);
	z=100.+indmul*z;
	f[k][1]=z;
	xx=xc0+in*del;
	ergrot(xx,y,&z);
/*	printf ("%10.3f%10.3f%10.5f%10.3f EE \n", xx,y,z,del);*/
        printf(".");
        fflush(stdout);
	z=100.+indmul*z;
	f[ind][k]=z;
	yy=yc0+in*del;
	ergrot(x,yy,&z);
/*	printf ("%10.3f%10.3f%10.5f%10.3f FF\n", x,yy,z,del);*/
        printf(".");
        fflush(stdout);
	z=100.+indmul*z;
	f[k][ind]=z;
e17:	x+=del;
	y+=del;
    }
      goto e23;
e10:  printf ("\nOptimal center coordinates: x =%10.3f, y = %10.3f\n",
	      yc0-1,xc0-1);
      return;
}

/****************************************************************************/
