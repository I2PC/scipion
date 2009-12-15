/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *		      Jean-Noël PIOCHE (jnp95@hotmail.com) 
 *		
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoléculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Châlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 27/05/2009 15:52:56  
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

#include "projection_real_shears.h"

///Parameters reading. Note that all parameters are required.
void Prog_Project_Parameters_2::read(int argc, char **argv)
{
	fn_proj_param = getParameter(argc, argv, "-i");
	fn_sel_file   = getParameter(argc, argv, "-o", "");
	display = !checkParameter(argc, argv, "-quiet");
}

///Description of the projection_real_shears function.
void Prog_Project_Parameters_2::usage()
{
	printf("\nUsage:\n\n");
    	printf("projection_real_shears -i <Parameters File> \n"
               "                      [-o <sel_file>]\n"
	       "                      [-quiet]\n");
	printf("\tWhere:\n"
               "\t<Parameters File>:  File containing projection parameters\n"
	       "\t                    Note that only the top-four parameters lines are read\n"
               "\t                    Check the manual for a description of the parameters\n"
               "\t<sel_file>:         This is a selection file with all the generated\n"
               "\t                    projections.\n");
}

//-----------------------------------------------------------------------------------------------
///Reads the projection parameters of the input file and inserts it into Projection_Parameters fields.\n
///This is an "overloaded" function in order to use translation parameters.
void Projection_Parameters_withShift::read_withShift(const FileName &fn_proj_param)
{
	FILE    *fh_param;
    	char    line[201];
    	int     lineNo = 0;
    	char    *auxstr;

    	if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
        	REPORT_ERROR(3005,
                (std::string)"Projection_Parameters_withShift::read: There is a problem "
                 "opening the file " + fn_proj_param);

    	while (fgets(line, 200, fh_param) != NULL)
    	{
        	if (line[0] == 0)    continue;
        	if (line[0] == '#')  continue;
        	if (line[0] == '\n') continue;
        	switch (lineNo)
        	{
        		case 0: //***** Line 1 *****
				//Volume file
            			fnPhantom = firstWord(line, 3007,
                                    "Projection_Parameters_withShift::read: Phantom name not found");

				if (!exists(fnPhantom))
                			REPORT_ERROR(3007, (std::string)"Projection_Parameters_withShift::read: "
                             			"file " + fnPhantom + " doesn't exist");

            			lineNo = 1;
            			break;
       		 	case 1: //***** Line 2 *****
            			fnProjectionSeed = firstWord(line, 3007,
                           	     "Projection_Parameters_withShift::read: Error in Projection seed");

            			auxstr = nextToken();
            			if (auxstr != NULL) starting =
                    			textToInteger(auxstr, 3007,
                         		"Projection_Parameters_withShift::read: Error in First "
                         		"projection number");

            			fn_projection_extension = nextWord(3007, (std::string)"Projection_Parameters_withShift::read: "
                             					"Error in Projection extension");
				
            			lineNo = 2;
            			break;
        		case 2: //***** Line 3 *****
            			proj_Xdim = textToInteger(firstToken(line), 3007,
                             		"Projection_Parameters_withShift::read: Error in projection dimension");
            			proj_Zdim = proj_Ydim = proj_Xdim ;

            			lineNo = 3;
            			break;
        		case 3: //***** Line 4 *****
            			// Angle file
            			fn_angle = firstWord(line, 3007,
                                    "Projection_Parameters_withShift::read: Angle file name not found");

            			if (!exists(fn_angle))
                			REPORT_ERROR(3007, (std::string)"Projection_Parameters_withShift::read: "
                             			"file " + fn_angle + " doesn't exist");
					
            			lineNo = 4;
            			break;
        		default:
				break;
        	} // switch end
    } // while end

    if (lineNo != 4) //If all parameters was not read
        REPORT_ERROR(3007, (std::string)"Projection_Parameters_withShift::read: I "
                     "couldn't read all parameters from file " + fn_proj_param);
    fclose(fh_param);
}

//-----------------------------------------------------------------------------------------------
///Computes one iteration of the projection. The resulting projection is into the pointer parameter called "Projection".\n
///Returns possible error.
int Projection_real_shears::do_compute_projection	(double *VWOrigin, 
										 long   Ndv, 
										 long Ndw, 
										 double *Identity_orientV, 
										 double *Identity_orientW, 
										 double dv, 
										 double dw, 
										 double *CoefVolume, 
										 double absscale, 
										 double *Binv, 
										 double *BinvCscaled, 
										 long ksimax,  
										 int *arr, 
										 long CoefVolumeNx, 
										 long CoefVolumeNy, 
										 long lmax, 
										 long mmax, 
										 double *Projection){  
 
 
int		Status = !ERROR; 
long	i, n, l, l1, l2, m, m1, m2, ksi, CC1, CC2, CC3, row, column, index; 
double	CVinc[4], CWinc[4], Operhlp[4], X[4], K[4], ToAdd[4], Arg[4], idw, ndv; 
double	Proj, sc, g, h, rows, columns, Coeff, Difference; 
double	gminusl, hminusm; 
 
CC1 = CoefVolumeNx * CoefVolumeNy; 
 
for (i = 0L; i < Ndw; i++){ 
	idw = (double) i * dw; 
	if (VectorScale(Identity_orientW, CWinc, idw, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					 "Error returned by VectorScale");
		return(ERROR); 
	} 
	if (VectorAdd(VWOrigin, CWinc, Operhlp, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					 "Error returned by VectorAdd"); 
		return(ERROR); 
	} 
	for (n = 0L; n < Ndv; n++){		 
		ndv = (double) n * dv; 
		if (VectorScale(Identity_orientV, CVinc, ndv, 4L) == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					      "Error returned by VectorScale"); 
			return(ERROR); 
		} 
		if (VectorAdd(Operhlp, CVinc, X, 4L) == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					      "Error returned by VectorAdd"); 
			return(ERROR); 
		} 
		if (MatrixTimesVector(Binv, X, K, 4L, 4L) == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					      "Error returned by MatrixTimesVector"); 
			return(ERROR); 
		} 
		Proj = 0.0; 
		for (ksi = 0L; ksi < ksimax; ksi++){ 
			CC2 = CC1 * ksi; 
			sc = (double) ksi - K[arr[0]]; 
			if (VectorScale(BinvCscaled, ToAdd, sc, 4L) == ERROR){  
				REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					 		 "Error returned by VectorScale"); 
				return(ERROR); 
			} 
			if (VectorAdd(K, ToAdd, Arg, 4L) == ERROR){  
				REPORT_ERROR(1, "Projection_real_shears::do_compute_projection: "
					 		 "Error returned by VectorAdd"); 
				return(ERROR); 
			} 
			g = Arg[arr[1]]; 
			h = Arg[arr[2]]; 
					 
			l1 = (long) ceil (g - 2.0); 
			l2 = l1 + 3L; 
			m1 = (long) ceil (h - 2.0); 
			m2 = m1 + 3L; 
			columns = 0.0; 
			for (m = m1; m <= m2; m++){ 
				if (m < mmax && m > -1L) { 
					CC3 = CC2 + CoefVolumeNx * m; 
					rows = 0.0; 
					for (l = l1; l <= l2; l++){ 
						if ((l < lmax && l > -1L)) { 
							gminusl = g - (double) l; 
							Coeff = (double) CoefVolume[CC3 + l]; 
							rows += Coeff * Bspline03(gminusl); 
						} 
					} 
					hminusm = h - (double) m; 
					columns +=  rows * Bspline03(hminusm); 
				} 
			} 
			Proj += columns; 
		} 
		Proj *= absscale; 
		 
		*Projection++ = Proj;											 
	} 
} 
 
return(!ERROR); 
}/* End of do_compute_projection */ 

//-----------------------------------------------------------------------------------------------
///Computes projection. The resulting projection is into the pointer parameter called "Projection".\n
///Returns possible error.
int Projection_real_shears::Compute_projection(double *Parameters, 
                                   		  double *Coef_x, 
									  double *Coef_y, 
									  double *Coef_z, 
									  long Nx, 
									  long Ny, 
									  long Nz, 
									  short *Proj_dims, 
									  double *Identity_orientN, 
									  double *Identity_orientV, 
									  double *Identity_orientW, 
									  double *IdentityOrigin, 
									  double *PeriodOfSamplingInVDirection, 
									  double *PeriodOfSamplingInWDirection, 
									  double *RightOperHlp, 
									  double *Ac, 
									  double *Projection, 
									  double *B				) { 
 
int	Status = !ERROR, arr[3]; 
long	Ndv, Ndw; 
long	CoefVolumeNx, CoefVolumeNy, lmax, mmax, ksimax; 
double	dv, dw, psi, theta, phi, Sinphi, Cosphi, Sinpsi, Cospsi, Sintheta, Costheta; 
double	scale, scale_x, scale_y, scale_z, m_x, m_y, m_z, minm; 
double	*hlp, *Rz, *Ry, *Rx, *At; 
double	*Help1, *Help2, *Help3, *Help4, *Binv; 
double	*C1, *C2, *C3, *VWOrigin, *BinvC, *BinvCscaled; 
double	*Coef_xyz, *Pr; 
 
Pr = Projection; 
	 
Ndv = (long) *Proj_dims++; 
Ndw = (long) *Proj_dims; 
	 
dv = *PeriodOfSamplingInVDirection; 
dw = *PeriodOfSamplingInWDirection; 
 
psi = Parameters[0]; 
theta =	Parameters[1]; 
phi = Parameters[2]; 
 
Sinphi = sin(phi); 
Cosphi = cos(phi); 
Sinpsi = sin(psi); 
Cospsi = cos(psi); 
Sintheta = sin(theta); 
Costheta = cos(theta); 
 
At = (double *)malloc((size_t) 16L * sizeof(double)); 
if (At == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for At"); 
	return(ERROR); 
} 
 
if (GetIdentitySquareMatrix(At, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by GetIdentitySquareMatrix"); 
	free(At); 
	return(ERROR); 
} 
	 
hlp = At + (ptrdiff_t)3L; 
*hlp = Parameters[3]; 
hlp += (ptrdiff_t)4L; 
*hlp = Parameters[4]; 
hlp += (ptrdiff_t)4L; 
*hlp = Parameters[5]; 
 
Help1 = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Help1 == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Help1"); 
	free(At); 
	return(ERROR); 
} 
 
if (MatrixMultiply(At, Ac, Help1, 4L, 4L, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by MatrixMultiply"); 
	free(At); 
	free(Help1); 
	return(ERROR); 
} 
 
Rx = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Rx == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Rx"); 
	free(At); 
	free(Help1); 
	return(ERROR); 
} 
	 
if (GetIdentitySquareMatrix(Rx, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by GetIdentitySquareMatrix"); 
	free(At); 
	free(Help1); 
	free(Rx); 
	return(ERROR); 
} 
	 
hlp = Rx; 
hlp += (ptrdiff_t)5L; 
*hlp++ = Cosphi; 
*hlp = - Sinphi; 
hlp += (ptrdiff_t)3L; 
*hlp++ = Sinphi; 
*hlp = Cosphi; 
 
Help2 = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Help2 == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Help2"); 
	free(At); 
	free(Help1); 
	free(Rx); 
	return(ERROR); 
} 
 
if (MatrixMultiply(Help1, Rx, Help2, 4L, 4L, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by MatrixMultiply"); 
	free(At); 
	free(Help1); 
	free(Rx); 
	free(Help2); 
	return(ERROR); 
} 
 
free(Rx); 
 
Ry = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Ry == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Ry"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	return(ERROR); 
} 
 
if (GetIdentitySquareMatrix(Ry, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by GetIdentitySquareMatrix"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	return(ERROR); 
} 
	 
hlp = Ry; 
*hlp = Costheta; 
hlp += (ptrdiff_t)2L; 
*hlp = Sintheta; 
hlp += (ptrdiff_t)6L; 
*hlp = - Sintheta; 
hlp += (ptrdiff_t)2L; 
*hlp = Costheta; 
 
Help3 = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Help3 == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Help3"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	return(ERROR); 
} 
 
if (MatrixMultiply(Help2, Ry, Help3, 4L, 4L, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by MatrixMultiply"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	return(ERROR); 
} 
 
Rz = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Rz == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Rz"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	return(ERROR); 
} 
 
if (GetIdentitySquareMatrix(Rz, 4L) == ERROR){ 
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by GetIdentitySquareMatrix"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	free(Rz); 
	return(ERROR); 
} 
	 
hlp = Rz; 
*hlp++ = Cospsi; 
*hlp = - Sinpsi; 
hlp += (ptrdiff_t)3L; 
*hlp++ = Sinpsi; 
*hlp = Cospsi; 
 
Help4 = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Help4 == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Help4"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	free(Rz); 
	return(ERROR); 
} 
 
if (MatrixMultiply(Help3, Rz, Help4, 4L, 4L, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by MatrixMultiply"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	free(Rz); 
	free(Help4); 
	return(ERROR); 
} 
 
 
if (MatrixMultiply(Help4, RightOperHlp, B, 4L, 4L, 4L) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by MatrixMultiply"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	free(Rz); 
	free(Help4); 
	return(ERROR); 
} 
 
free(Help4); 
 
Binv = (double *)malloc((size_t) 16L * sizeof(double)); 
if (Binv == (double *)NULL){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for Binv"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	free(Rz); 
	return(ERROR); 
} 
 
if (SquareMatrixInvertGauss(B, Binv, 4L, DBL_EPSILON, &Status) == ERROR){  
	REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by SquareMatrixInvertGauss"); 
	free(At); 
	free(Help1); 
	free(Help2); 
	free(Ry); 
	free(Help3); 
	free(Rz); 
	free(Binv); 
	return(ERROR); 
} 
 
free(Rz); 
free(Ry); 
free(At); 
 
 
	C1 = (double *)malloc((size_t) 4L * sizeof(double)); 
	if (C1 == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for C1"); 
		free(Binv); 
		return(ERROR); 
	} 
 
	hlp = C1; 
	*hlp++ = (double) Identity_orientN[0]; 
	*hlp++ = (double) Identity_orientN[1]; 
	*hlp++ = (double) Identity_orientN[2]; 
	*hlp = 0.0; 
	 
	BinvC = (double *)malloc((size_t) 4L * sizeof(double)); 
	if (BinvC == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for BinvC"); 
		free(Binv); 
		free(C1); 
		return(ERROR); 
	} 
	 
	if (MatrixTimesVector(Binv, C1, BinvC, 4L, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by MatrixTimesVector"); 
		free(Binv); 
		free(C1); 
		free(BinvC); 
		return(ERROR); 
	}
 
	if (BinvC[0] != 0.0) 
		scale_x = 1.0 / BinvC[0]; 
	else 
		scale_x = DBL_MAX; 
		 
	if (BinvC[1] != 0.0) 
		scale_y = 1.0 / BinvC[1]; 
	else 
		scale_y = DBL_MAX; 
		 
	if (BinvC[2] != 0.0) 
		scale_z = 1.0 / BinvC[2]; 
	else 
		scale_z = DBL_MAX; 
		 
	m_x = fabs(scale_x); 
	m_y = fabs(scale_y); 
	m_z = fabs(scale_z); 
	 
	BinvCscaled = (double *)malloc((size_t) 4L * sizeof(double)); 
	if (BinvCscaled == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for BinvCscaled"); 
		free(Binv); 
		free(C1); 
		free(BinvC); 
		return(ERROR); 
	} 
 
	minm = m_x; 
	scale = scale_x; 
	if (VectorScale(BinvC, BinvCscaled, scale_x, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by VectorScale"); 
		free(Binv); 
		free(C1); 
		free(BinvC); 
		free(BinvCscaled); 
		return(ERROR); 
	} 
	ksimax = Nx; 
	arr[0] = 0; 
	arr[1] = 1; 
	arr[2] = 2; 
	CoefVolumeNx = Ny; 
	CoefVolumeNy = Nz; 
	lmax = Ny; 
	mmax = Nz; 
	Coef_xyz = Coef_x; 
	if (m_y < minm) { 
		minm = m_y; 
		scale = scale_y; 
		if (VectorScale(BinvC, BinvCscaled, scale_y, 4L) == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by VectorScale"); 
			free(Binv); 
			free(C1); 
			free(BinvC); 
			free(BinvCscaled); 
			return(ERROR); 
		} 
		ksimax = Ny; 
		arr[0] = 1; 
		arr[1] = 0; 
		arr[2] = 2; 
		CoefVolumeNx = Nx; 
		CoefVolumeNy = Nz; 
		lmax = Nx; 
		mmax = Nz; 
		Coef_xyz = Coef_y; 
	} 
	if (m_z < minm) { 
		minm = m_z; 
		scale = scale_z; 
		if (VectorScale(BinvC, BinvCscaled, scale_z, 4L) == ERROR){ 
			free(Binv); 
			free(C1); 
			free(BinvC); 
			free(BinvCscaled);  
			return(ERROR); 
		} 
		ksimax = Nz; 
		arr[0] = 2; 
		arr[1] = 0; 
		arr[2] = 1; 
		CoefVolumeNx = Nx; 
		CoefVolumeNy = Ny; 
		lmax = Nx; 
		mmax = Ny; 
		Coef_xyz = Coef_z; 
	} 
 
	 
	free(BinvC);	 
	free(C1); 
	 
	C2 = (double *)malloc((size_t) 4L * sizeof(double)); 
	if (C2 == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for C2"); 
		free(Binv); 
		free(BinvCscaled); 
		return(ERROR); 
	} 
 
	C3 = (double *)malloc((size_t) 4L * sizeof(double)); 
	if (C3 == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for C3"); 
		free(Binv); 
		free(BinvCscaled); 
		free(C2); 
		return(ERROR); 
	} 
	 
	VWOrigin = (double *)malloc((size_t) 4L * sizeof(double)); 
	if (VWOrigin == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "ERROR - Not enough memory for VWOrigin"); 
		free(Binv); 
		free(BinvCscaled); 
		free(C2); 
		free(C3); 
		return(ERROR); 
	} 
	 
	hlp = C2; 
	*hlp++ = (double) Identity_orientV[0]; 
	*hlp++ = (double) Identity_orientV[1]; 
	*hlp++ = (double) Identity_orientV[2]; 
	*hlp = 0.0; 
	 
	hlp = C3; 
	*hlp++ = (double) Identity_orientW[0]; 
	*hlp++ = (double) Identity_orientW[1]; 
	*hlp++ = (double) Identity_orientW[2]; 
	*hlp = 0.0; 
	 
	hlp = VWOrigin; 
	*hlp++ = (double) IdentityOrigin[0]; 
	*hlp++ = (double) IdentityOrigin[1]; 
	*hlp++ = (double) IdentityOrigin[2]; 
	*hlp = 1.0; 
	 
	if (do_compute_projection(VWOrigin, Ndv, Ndw, C2, C3, 
		dv, dw, Coef_xyz, minm, Binv, BinvCscaled, ksimax, arr, CoefVolumeNx, 
		CoefVolumeNy, lmax, mmax, Pr) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::Compute_projection: "
					 "Error returned by do_compute_projection"); 
		free(Binv); 
		free(BinvCscaled); 
		free(C2); 
		free(C3); 
		free(VWOrigin); 
		return(ERROR); 
	} 
	 
free(BinvCscaled); 
free(C2); 
free(C3); 
free(VWOrigin);	 
free(Binv);
free(Help1);
free(Help2);
free(Help3); 
			 
return(!ERROR); 
}/* End of Compute_projection */

///Main compute function. Returns possible error.
int Projection_real_shears::ROUT_project_execute(VolumeStruct &Data2)
{
	int	Status = !ERROR; 
	long	DesProjSize; 
	long	i, m, n, l, Nx, Ny, Nz; 
	double	lambda; 
	double	*Parameters; 
	double	*hlp, *Av, *As, *Ac, *Acinv, *RightOperHlp, *B;	 
	double	*Projection, *VolumeCoef, *InputVolume, *InputVolumePlane; 
	double	*InputVolumeRow, *Coef_x, *Coef_y, *Coef_z; 
		 
	Nx = Data2.nx_Volume; 
	Ny = Data2.ny_Volume; 
	Nz = Data2.nz_Volume; 
	 
	InputVolume = Data2.Volume;

	AllocateVolumeDouble( &Coef_x, Ny, Nz, Nx, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Coef_x"); 
		return(ERROR); 
	} 
 
	AllocateVolumeDouble( &Coef_y, Nx, Nz, Ny, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Coef_y"); 
		FreeVolumeDouble(&Coef_x); 
		return(ERROR); 
	} 
 
	AllocateVolumeDouble( &Coef_z, Nx, Ny, Nz, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Coef_z"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		return(ERROR); 
	} 
 
	 
	AllocateVolumeDouble( &VolumeCoef, Ny, Nz, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for VolumeCoef"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		return(ERROR); 
	} 
	 
	AllocateVolumeDouble( &InputVolumePlane, Ny, Nz, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for InputVolumePlane"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&VolumeCoef); 
		return(ERROR); 
	} 
	 
	AllocateVolumeDouble( &InputVolumeRow, 1L, Ny, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for InputVolumePlane"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&VolumeCoef); 
		FreeVolumeDouble(&InputVolumePlane); 
		return(ERROR); 
	} 
 
	for (l = 0L; l < Nx; l++){
		 
		for (m = 0L; m < Nz; m++){  
			if (CopyDoubleToDouble(InputVolume, Nx, Ny, Nz, l,	0L,	m,	InputVolumeRow, 1L, Ny, 1L, 0L, 0L, 0L, 1L, Ny, 1L) == ERROR){ 
				REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
				FreeVolumeDouble(&Coef_x); 
				FreeVolumeDouble(&Coef_y); 
				FreeVolumeDouble(&Coef_z); 
				FreeVolumeDouble(&VolumeCoef); 
				FreeVolumeDouble(&InputVolumePlane); 
				FreeVolumeDouble(&InputVolumeRow); 
				return(ERROR); 
			} 
			for (i = 0L; i < Ny; i++){ 
				InputVolumePlane[Ny * m + i] = InputVolumeRow[i]; 
			} 
		} 
		 
		ChangeBasisVolume(InputVolumePlane, VolumeCoef, Ny, Nz, 1L, CardinalSpline, 
					BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status); 
		if (Status == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR"); 
			FreeVolumeDouble(&Coef_x); 
			FreeVolumeDouble(&Coef_y); 
			FreeVolumeDouble(&Coef_z); 
			FreeVolumeDouble(&VolumeCoef); 
			FreeVolumeDouble(&InputVolumePlane); 
			FreeVolumeDouble(&InputVolumeRow); 
			return(ERROR); 
		} 
		 
		if (CopyDoubleToDouble(VolumeCoef, Ny, Nz, 1L, 0L,	0L,	0L,	Coef_x, Ny, Nz, Nx, 0L, 0L, l, Ny, Nz, 1L) == ERROR){ 
			REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
			FreeVolumeDouble(&Coef_x); 
			FreeVolumeDouble(&Coef_y); 
			FreeVolumeDouble(&Coef_z); 
			FreeVolumeDouble(&VolumeCoef); 
			FreeVolumeDouble(&InputVolumePlane); 
			FreeVolumeDouble(&InputVolumeRow); 
			return(ERROR); 
		} 
		 
	} 
	 
	FreeVolumeDouble(&VolumeCoef);  
	FreeVolumeDouble(&InputVolumePlane); 
	FreeVolumeDouble(&InputVolumeRow);
	 
 
	 
	AllocateVolumeDouble( &VolumeCoef, Nx, Nz, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for VolumeCoef"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		return(ERROR); 
	} 
	 
	AllocateVolumeDouble( &InputVolumePlane, Nx, Nz, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for InputVolumePlane"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&VolumeCoef); 
		return(ERROR); 
	} 
	 
	AllocateVolumeDouble( &InputVolumeRow, Nx, 1L, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for InputVolumePlane"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&VolumeCoef); 
		FreeVolumeDouble(&InputVolumePlane); 
		return(ERROR); 
	} 
 
	for (l = 0L; l < Ny; l++){ 
		for (m = 0L; m < Nz; m++){  
			if (CopyDoubleToDouble(InputVolume, Nx, Ny, Nz, 0L, l,	m, InputVolumeRow, Nx, 1L, 1L, 0L, 0L, 0L, Nx, 1L, 1L) == ERROR){ 
				REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
				FreeVolumeDouble(&Coef_x); 
				FreeVolumeDouble(&Coef_y); 
				FreeVolumeDouble(&Coef_z); 
				FreeVolumeDouble(&VolumeCoef); 
				FreeVolumeDouble(&InputVolumePlane); 
				FreeVolumeDouble(&InputVolumeRow); 
				return(ERROR); 
			} 
			 
			if (CopyDoubleToDouble(InputVolumeRow, Nx, 1L, 1L, 0L,	0L,	0L,	InputVolumePlane, Nx, Nz, 1L, 0L, m, 0L, Nx, 1L, 1L) == ERROR){ 
				REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
				FreeVolumeDouble(&Coef_x); 
				FreeVolumeDouble(&Coef_y); 
				FreeVolumeDouble(&Coef_z); 
				FreeVolumeDouble(&VolumeCoef); 
				FreeVolumeDouble(&InputVolumePlane); 
				FreeVolumeDouble(&InputVolumeRow); 
				return(ERROR); 
			} 
			 
		} 
	 
		ChangeBasisVolume((double*)InputVolumePlane, (double*)VolumeCoef, Nx, Nz, 1L, CardinalSpline, 
					BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status); 
		if (Status == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR"); 
			FreeVolumeDouble(&Coef_x); 
			FreeVolumeDouble(&Coef_y); 
			FreeVolumeDouble(&Coef_z); 
			FreeVolumeDouble(&VolumeCoef); 
			FreeVolumeDouble(&InputVolumePlane); 
			FreeVolumeDouble(&InputVolumeRow); 
			return(ERROR); 
		} 
		 
		if (CopyDoubleToDouble(VolumeCoef, Nx, Nz, 1L, 0L,	0L,	0L,	Coef_y, Nx, Nz, Ny, 0L, 0L, l, Nx, Nz, 1L) == ERROR){ 
			REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
			FreeVolumeDouble(&Coef_x); 
			FreeVolumeDouble(&Coef_y); 
			FreeVolumeDouble(&Coef_z); 
			FreeVolumeDouble(&VolumeCoef); 
			FreeVolumeDouble(&InputVolumePlane); 
			FreeVolumeDouble(&InputVolumeRow); 
			return(ERROR); 
		} 
	 
	} 
	 
	FreeVolumeDouble(&VolumeCoef); 
	FreeVolumeDouble(&InputVolumePlane); 
	FreeVolumeDouble(&InputVolumeRow);	 
	 
	 
	AllocateVolumeDouble( &VolumeCoef, Nx, Ny, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for VolumeCoef"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		return(ERROR); 
	} 
	 
	AllocateVolumeDouble( &InputVolumePlane, Nx, Ny, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for InputVolumePlane"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&VolumeCoef); 
		return(ERROR); 
	} 
	 
	AllocateVolumeDouble( &InputVolumeRow, Nx, 1L, 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for InputVolumePlane"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&VolumeCoef); 
		FreeVolumeDouble(&InputVolumePlane); 
		return(ERROR); 
	} 
 
	for (l = 0L; l < Nz; l++){ 
		for (m = 0L; m < Ny; m++){  
			if (CopyDoubleToDouble(InputVolume, Nx, Ny, Nz, 0L, m,	l, InputVolumeRow, Nx, 1L, 1L, 0L, 0L, 0L, Nx, 1L, 1L) == ERROR){ 
				REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
				FreeVolumeDouble(&Coef_x); 
				FreeVolumeDouble(&Coef_y); 
				FreeVolumeDouble(&Coef_z); 
				FreeVolumeDouble(&VolumeCoef); 
				FreeVolumeDouble(&InputVolumePlane); 
				FreeVolumeDouble(&InputVolumeRow); 
				return(ERROR); 
			} 
			 
			if (CopyDoubleToDouble(InputVolumeRow, Nx, 1L, 1L, 0L,	0L,	0L,	InputVolumePlane, Nx, Ny, 1L, 0L, m, 0L, Nx, 1L, 1L) == ERROR){ 
				REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
				FreeVolumeDouble(&Coef_x); 
				FreeVolumeDouble(&Coef_y); 
				FreeVolumeDouble(&Coef_z); 
				FreeVolumeDouble(&VolumeCoef); 
				FreeVolumeDouble(&InputVolumePlane); 
				FreeVolumeDouble(&InputVolumeRow); 
				return(ERROR); 
			} 
			 
		} 
	 
		ChangeBasisVolume(InputVolumePlane, VolumeCoef, Nx, Ny, 1L, CardinalSpline, 
					BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status); 
		if (Status == ERROR){  
			REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR"); 
			FreeVolumeDouble(&Coef_x); 
			FreeVolumeDouble(&Coef_y); 
			FreeVolumeDouble(&Coef_z); 
			FreeVolumeDouble(&VolumeCoef); 
			FreeVolumeDouble(&InputVolumePlane); 
			FreeVolumeDouble(&InputVolumeRow); 
			return(ERROR); 
		} 
		 
		if (CopyDoubleToDouble(VolumeCoef, Nx, Ny, 1L, 0L,	0L,	0L,	Coef_z, Nx, Ny, Nz, 0L, 0L, l, Nx, Ny, 1L) == ERROR){ 
			REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by CopyDoubleToDouble");  
			FreeVolumeDouble(&Coef_x); 
			FreeVolumeDouble(&Coef_y); 
			FreeVolumeDouble(&Coef_z); 
			FreeVolumeDouble(&VolumeCoef); 
			FreeVolumeDouble(&InputVolumePlane); 
			FreeVolumeDouble(&InputVolumeRow); 
			return(ERROR); 
		} 
	 
	} 
	 
	FreeVolumeDouble(&VolumeCoef); 
	FreeVolumeDouble(&InputVolumePlane); 
	FreeVolumeDouble(&InputVolumeRow);	 
 
 
	 
	AllocateVolumeDouble( &Projection, *(Data2.Proj_dims), *(Data2.Proj_dims + (ptrdiff_t) 1L), 1L, &Status); 
	if (Status == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Projection"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		return(ERROR); 
	} 
	 
	Av = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (Av == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Av"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		return(ERROR); 
	} 
		 
	if (GetIdentitySquareMatrix(Av, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by GetIdentitySquareMatrix"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		return(ERROR); 
	} 
 
	hlp = Av + (ptrdiff_t)3L; 
	*hlp = (double) Data2.K123[0]; 
	hlp += (ptrdiff_t)4L; 
	*hlp = (double) Data2.K123[1]; 
	hlp += (ptrdiff_t)4L; 
	*hlp = (double) Data2.K123[2]; 
	 
	As = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (As == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for As"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		return(ERROR); 
	} 
	 
	if (GetIdentitySquareMatrix(As, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by GetIdentitySquareMatrix"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		return(ERROR); 
	} 
	 
	hlp = As; 
	*hlp = (double) Data2.Lambda123[0]; 
	hlp += (ptrdiff_t)5L; 
	*hlp = (double) Data2.Lambda123[1]; 
	hlp += (ptrdiff_t)5L; 
	*hlp = (double) Data2.Lambda123[2]; 
	 
	Ac = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (Ac == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Ac"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		return(ERROR); 
	} 
	 
	if (GetIdentitySquareMatrix(Ac, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by GetIdentitySquareMatrix"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		free(Ac); 
		return(ERROR); 
	} 
	 
	hlp = Ac + (ptrdiff_t)3L; 
	*hlp = (double) Data2.Gama123[0]; 
	hlp += (ptrdiff_t)4L; 
	*hlp = (double) Data2.Gama123[1]; 
	hlp += (ptrdiff_t)4L; 
	*hlp = (double) Data2.Gama123[2]; 
	 
	Acinv = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (Acinv == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Acinv"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		free(Ac); 
		return(ERROR); 
	} 
	 
	if (SquareMatrixInvertGauss(Ac, Acinv, 4L, DBL_EPSILON, &Status) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by SquareMatrixInvertGauss"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		free(Ac); 
		free(Acinv); 
		return(ERROR); 
	} 
	 
	RightOperHlp = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (RightOperHlp == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for RightOperHlp"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		free(Ac); 
		free(Acinv); 
		return(ERROR); 
	} 
		 
	if (multiply_3Matrices(Acinv, As, Av, RightOperHlp, 4L, 4L, 4L, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by multiply_3Matrices"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Av); 
		free(As); 
		free(Ac); 
		free(Acinv); 
		free(RightOperHlp); 
		return(ERROR); 
	} 
	 
	free(Av); 
	free(As); 
	free(Acinv); 
	 
	Parameters = (double *)malloc((size_t) 6L * sizeof(double)); 
	if (Parameters == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for Parameters"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Ac); 
		free(RightOperHlp); 
		return(ERROR); 
	} 
						 
	Parameters[2] = Data2.InitPsiThetaPhi[0];	 
	Parameters[1] = Data2.InitPsiThetaPhi[1]; 
	Parameters[0] = Data2.InitPsiThetaPhi[2]; 
	Parameters[3] = Data2.InitDelta123[0]; 
	Parameters[4] = Data2.InitDelta123[1]; 
	Parameters[5] = Data2.InitDelta123[2]; 
	 
	B = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (B == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "ERROR - Not enough memory for B"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Ac); 
		free(RightOperHlp); 
		free(Parameters); 
		return(ERROR); 
	} 
	 
	if (Compute_projection(Parameters, Coef_x, Coef_y, Coef_z, Nx, Ny, Nz, Data2.Proj_dims,  
							Data2.Identity_orientN, Data2.Identity_orientV, Data2.Identity_orientW,  
							Data2.IdentityOrigin, &Data2.PeriodOfSamplingInVDirection, 
							&Data2.PeriodOfSamplingInWDirection, RightOperHlp, Ac,  
							Projection, B) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_execute: "
					 "Error returned by Compute_projection"); 
		FreeVolumeDouble(&Coef_x); 
		FreeVolumeDouble(&Coef_y); 
		FreeVolumeDouble(&Coef_z); 
		FreeVolumeDouble(&Projection); 
		free(Ac); 
		free(RightOperHlp); 
		free(Parameters); 
		free(B); 
		return(ERROR); 
	}

	//Copy of the result pointer (Projection is a dynamic allocation)
	Data2.Output = Projection;

	free(Parameters);
	free(Ac);
	free(RightOperHlp);
	free(B);
	FreeVolumeDouble(&Coef_x); 
	FreeVolumeDouble(&Coef_y); 
	FreeVolumeDouble(&Coef_z); 

	return(!ERROR);
 }

///Returns the sign of the parameter : +1 if the number is positive or null, -1 else.
inline double sign(double number)
{
	if(number>=0.) return +1.;
	else           return -1.;
}

///Transforms the degree parameter into radian number and returns it.
inline double toRad(double number)
{
	return (number*PI/180.0);
}

///Returns a pointer of the multiplication of the 5 matrices built with identity matrices and the parameters.\n
///It returns NULL if there is an error.
double *Projection_real_shears::MatrixBem(double phi, double theta, double psi, double x0, double y0, double scale_x, double scale_y, double scale_z)
{
	double *pointer;
 
	double ss = sin(phi); 
	double cc = cos(phi);

//------------------------------------------------ 
	/*Rz1 [4][4] = {{ cc,  ss, 0.0, 0.0}, 
			{-ss,  cc, 0.0, 0.0}, 
			{0.0, 0.0, 1.0, 0.0}, 
			{0.0, 0.0, 0.0, 1.0}};*/
//------------------------------------------------

	double *Rz1 = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (Rz1 == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "ERROR - Not enough memory for Rz1"); 
		free(Rz1); 
		return(NULL); 
	} 
 
	if (GetIdentitySquareMatrix(Rz1, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "Error returned by GetIdentitySquareMatrix"); 
		free(Rz1); 
		return(NULL); 
	} 
	 
	pointer = Rz1; 
	*pointer = cc; 
	pointer += (ptrdiff_t)1L; 
	*pointer = ss; 
	pointer += (ptrdiff_t)3L; 
	*pointer = - ss; 
	pointer += (ptrdiff_t)1L; 
	*pointer = cc;
 
 
	ss = sin(theta); 
	cc = cos(theta);

//------------------------------------------------ 
	/*Ry [4][4] = {{ cc, 0.0, -ss, 0.0}, 
		       {0.0, 1.0, 0.0, 0.0}, 
		       { ss, 0.0,  cc, 0.0}, 
		       {0.0, 0.0, 0.0, 1.0}};*/
//------------------------------------------------ 

	double *Ry = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (Ry == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "ERROR - Not enough memory for Ry"); 
		free(Rz1);
		free(Ry); 
		return(NULL); 
	} 
 
	if (GetIdentitySquareMatrix(Ry, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "Error returned by GetIdentitySquareMatrix"); 
		free(Rz1);
		free(Ry); 
		return(NULL); 
	} 
	 
	pointer = Ry; 
	*pointer = cc; 
	pointer += (ptrdiff_t)2L; 
	*pointer = -ss; 
	pointer += (ptrdiff_t)6L; 
	*pointer = ss; 
	pointer += (ptrdiff_t)2L; 
	*pointer = cc;


 
	ss = sin(psi); 
	cc = cos(psi);

//------------------------------------------------ 
	/*Rz2 [4][4] = {{ cc,  ss, 0.0, 0.0}, 
			{-ss,  cc, 0.0, 0.0}, 
			{0.0, 0.0, 1.0, 0.0}, 
			{0.0, 0.0, 0.0, 1.0}};*/
//------------------------------------------------

	double *Rz2 = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (Rz2 == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "ERROR - Not enough memory for Rz2"); 
		free(Rz1);
		free(Ry);
		free(Rz2); 
		return(NULL); 
	} 
 
	if (GetIdentitySquareMatrix(Rz2, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "Error returned by GetIdentitySquareMatrix"); 
		free(Rz1);
		free(Ry);
		free(Rz2); 
		return(NULL); 
	} 
	 
	pointer = Rz2; 
	*pointer = cc; 
	pointer += (ptrdiff_t)1L; 
	*pointer = ss; 
	pointer += (ptrdiff_t)3L; 
	*pointer = -ss; 
	pointer += (ptrdiff_t)1L; 
	*pointer = cc;



	//(Identity Matrix)
//------------------------------------------------ 
	/*At [4][4] = {{1., 0., 0., 0.},
		       {0., 1., 0., 0.},
		       {0., 0., 1., 0.},
		       {x0, y0, 0., 1.}};*/
//------------------------------------------------

	double *At = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (At == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "ERROR - Not enough memory for At"); 
		free(Rz1);
		free(Ry);
		free(Rz2);
		free(At); 
		return(NULL); 
	} 
 
	if (GetIdentitySquareMatrix(At, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "Error returned by GetIdentitySquareMatrix"); 
		free(Rz1);
		free(Ry);
		free(Rz2);
		free(At); 
		return(NULL); 
	} 
	 
	pointer = At; 
	pointer += (ptrdiff_t)8L; 
	*pointer = x0; 
	pointer += (ptrdiff_t)1L; 
	*pointer = y0;




	//(Identity Matrix)
//------------------------------------------------ 
	/*As [4][4] = {{scale_x,      0.,      0., 0.},
		       {     0., scale_y,      0., 0.},
		       {     0.,      0., scale_z, 0.},
		       {     0.,      0.,      0., 1.}};*/
//------------------------------------------------

	double *As = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (As == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "ERROR - Not enough memory for As"); 
		free(Rz1);
		free(Ry);
		free(Rz2);
		free(At);
		free(As); 
		return(NULL); 
	} 
 
	if (GetIdentitySquareMatrix(As, 4L) == ERROR){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "Error returned by GetIdentitySquareMatrix"); 
		free(Rz1);
		free(Ry);
		free(Rz2);
		free(At);
		free(As); 
		return(NULL); 
	} 
	 
	pointer = As;
	*pointer = scale_x; 
	pointer += (ptrdiff_t)5L; 
	*pointer = scale_y; 
	pointer += (ptrdiff_t)5L; 
	*pointer = scale_z;



//------------------------------------------------
	double *matrB = (double *)malloc((size_t) 16L * sizeof(double)); 
	if (matrB == (double *)NULL){  
		REPORT_ERROR(1, "Projection_real_shears::MatrixBem: "
					 "ERROR - Not enough memory for matrB"); 
		free(Rz1);
		free(Ry);
		free(Rz2);
		free(At);
		free(As);
		free(matrB); 
		return(NULL); 
	}


	multiply_5Matrices(At, As, Rz2, Ry, Rz1, matrB, 4L, 4L, 4L, 4L, 4L, 4L);

	free(Rz1);
	free(Ry);
	free(Rz2);
	free(At);
	free(As); 
 
	return matrB;
}

///Transforms angles from (Ry, Rz, Ry) to (Rx, Ry, Rz) system. Returns possible error.
int Projection_real_shears::angles_transcription(double *angles, double *Lambda123)
{	
	double phi = angles[0];
	double theta = angles[1];
	double psi = angles[2];
	double x0 = 0.;
	double y0 = 0.;
	double scale_x = Lambda123[0];
	double scale_y = Lambda123[1];
	double scale_z = Lambda123[2];
	
	double *Bem_1D = MatrixBem(phi, theta, psi, x0, y0, scale_x, scale_y, scale_z);
	if(Bem_1D==NULL) return (ERROR);

	double Bem [4][4];
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			Bem[i][j] = Bem_1D[i*4+j];
	free(Bem_1D);

	double A00 = Bem[0][0];
	double A10 = Bem[1][0];
	double A20 = Bem[2][0];
	double A21 = Bem[2][1];
	double A22 = Bem[2][2];
	double A12 = Bem[1][2];
	double A02 = Bem[0][2];
	

	double abs_cosay = sqrt(A22*A22+A21*A21);

	//Results angles
	double ax, ay, az;

	if(abs_cosay > 0.0)
	{
		double sign_cosay;

		ax = atan2(-A21, A22);
    	 	az = atan2(-A10, A00);

     	if(abs(cos(ax)) == 0.0)
     		sign_cosay = sign(-A21/sin(ax));
     	else
		{
			if (cos(ax) > 0.0) 
				sign_cosay = sign(A22);
			else
				sign_cosay = - sign(A22);
		}

		ay  = atan2(A20, sign_cosay * abs_cosay);
	}
	else
	{
		//Let's consider the matrix as a rotation around Z
		if (sign(A20) > 0.0)
		{
			ax = 0.0;
     		ay  = PI/2.0;
     		az = atan2(A12, -A02);
		}
		else
		{
			ax = 0.0;
     		ay  = -PI/2.0;
     		az = atan2(-A12, A02);
		}
	}

	angles[0] = ax;
	angles[1] = ay;
	angles[2] = az;

	return (!ERROR);
}

///Allocates and fixes some VolumeStruct fields
void Projection_real_shears::allocAndInit_VolumeStruct(VolumeStruct &Data2)
{
	Data2.Volume = (double*) malloc((size_t)(Data2.nx_Volume * Data2.ny_Volume * Data2.nz_Volume) * sizeof(double));

	Data2.Proj_dims = (short*) malloc((size_t)2L * sizeof(short));
	Data2.Proj_dims[0] = Data2.nx_Volume;
	Data2.Proj_dims[1] = Data2.ny_Volume;

	Data2.Identity_orientN = (double*) malloc((size_t)3L * sizeof(double));
	Data2.Identity_orientN[0] = 0.;
	Data2.Identity_orientN[1] = 0.;
	Data2.Identity_orientN[2] = 1.;
 
	Data2.Identity_orientV = (double*) malloc((size_t)3L * sizeof(double));
	Data2.Identity_orientV[0] = 1.;
	Data2.Identity_orientV[1] = 0.;
	Data2.Identity_orientV[2] = 0.;  
	Data2.Identity_orientW = (double*) malloc((size_t)3L * sizeof(double));
	Data2.Identity_orientW[0] = 0.;
	Data2.Identity_orientW[1] = 1.;
	Data2.Identity_orientW[2] = 0.; 
 
	Data2.IdentityOrigin = (double*) malloc((size_t)3L * sizeof(double));
	Data2.IdentityOrigin[0] = 0.;
	Data2.IdentityOrigin[1] = 0.;
	Data2.IdentityOrigin[2] = 0.;

	Data2.K123 = (double*) malloc((size_t)3L * sizeof(double));
	Data2.K123[0] = 0.;
	Data2.K123[1] = 0.;
	Data2.K123[2] = 0.;

	Data2.Lambda123 = (double*) malloc((size_t)3L * sizeof(double));
	Data2.Lambda123[0] = 1.;
	Data2.Lambda123[1] = 1.;
	Data2.Lambda123[2] = 1.;

	Data2.Gama123 = (double*) malloc((size_t)3L * sizeof(double));
	Data2.Gama123[0] = Data2.nx_Volume/2.;
	Data2.Gama123[1] = Data2.Gama123[0];
	Data2.Gama123[2] = Data2.Gama123[0];

	Data2.InitDelta123 = (double*) malloc((size_t)3L * sizeof(double));
	Data2.InitDelta123[2] = 0.;

	Data2.InitPsiThetaPhi = (double*) malloc((size_t)3L * sizeof(double)); 

	Data2.PeriodOfSamplingInVDirection = 1.;
	Data2.PeriodOfSamplingInWDirection = 1.;
}

///Desallocates VolumeStruct fields
void Projection_real_shears::del_VolumeStruct(VolumeStruct &Data2)
{
	free(Data2.Volume); 
	free(Data2.Proj_dims); 
	free(Data2.Identity_orientN); 
	free(Data2.Identity_orientV);
	free(Data2.Identity_orientW); 
	free(Data2.IdentityOrigin); 
	free(Data2.K123); 
	free(Data2.Lambda123); 
	free(Data2.Gama123);
	free(Data2.InitDelta123);
	free(Data2.InitPsiThetaPhi);
}

///Writes the projection file obtained. Returns possibles errors.
int Projection_real_shears::write_projection_file(int numFile)
{
	proj.clear();
	proj.reset(Data.nx_Volume, Data.ny_Volume);

	int indexX_start = -Data.nx_Volume / 2L; 
	int indexY_start = -Data.ny_Volume / 2L;

	//Projection object filling
	for(int i=0, k=indexY_start; i<Data.ny_Volume; i++, k++)
		for(int j=0, l=indexX_start; j<Data.nx_Volume; j++, l++)
			proj.setPixel(k, l, Data.Output[i*Data.nx_Volume+j]);

	//Composition of the projection name
	FileName fn_proj;
     fn_proj.compose(prm.fnProjectionSeed, numFile, prm.fn_projection_extension);

	//Projection save
	proj.write(fn_proj);

     SF.insert(fn_proj, SelLine::ACTIVE);

	return(!ERROR);
}

///Reads a DocLine and fills Data fields. Returns possibles errors.
int Projection_real_shears::read_a_DocLine()
{
	DocLine dl = DF.get_current_line(); // This is a data line of the fn_angle file 

	double rot    ; 
	double tilt   ; 
	double psi    ;
	double shiftX ;
	double shiftY ;

	//If the line is a comment
	if(dl.Is_comment())
	{
		DF.next_data_line();
		return (!ERROR);
	}
	else if(dl.get_no_components()<3) //If there aren't enough parameters
	{
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_real_shears: "
					 "Error returned by do_oneProjection : too few arguments from the angle file");
		return (ERROR);
	}
	else if(dl.get_no_components()==3) //If there are angles parameters but not shifts parameters
	{
		shiftX = 0.;
		shiftY = 0.;
	}
	else
	{
		shiftX = dl[3];
		shiftY = dl[4];
	}

	rot    = dl[0]; 
	tilt   = dl[1]; 
	psi    = dl[2]; 

	Data.InitPsiThetaPhi[2] = -toRad(rot);
	Data.InitPsiThetaPhi[1] = -toRad(tilt);
	Data.InitPsiThetaPhi[0] = -toRad(psi);

	Data.InitDelta123[0] = shiftX;
	Data.InitDelta123[1] = shiftY;

	return(!ERROR);
}


///////////////////////// MAIN INSTRUCTION FOR MPI ////////////////////////////////	
///Execute instructions for one projection
int Projection_real_shears::do_oneProjection(VolumeStruct &Data2)
{		
		//Transcription of the angles
		if(angles_transcription(Data2.InitPsiThetaPhi, Data2.Lambda123)==ERROR) return (ERROR);

		//Calls the compute function
		if(ROUT_project_execute(Data2) == ERROR) return(ERROR);

	return (!ERROR);
}
//////////////////////////////////////////////////////////////////////////////////


///Does start instructions. Returns possibles errors.
int Projection_real_shears::start_to_process()
{
	if(prog_param.display) cout<<endl;

	prm.read_withShift(prog_param.fn_proj_param);

	DF = DocFile(prm.fn_angle);  // Reads the fn_angle file

	//Reads the reference volume
	VolumeXmipp V;
	V.read(prm.fnPhantom);
	V().setXmippOrigin();
	Matrix3D <double> Volume = V();

	int z, y, x;
	Volume.getDimension(z, y, x);
	Data.nx_Volume = x;
	Data.ny_Volume = y;
	Data.nz_Volume = z;
	if(Data.nx_Volume<=0 || Data.ny_Volume<=0 || Data.nz_Volume<=0
		|| Data.nx_Volume!=Data.ny_Volume || Data.ny_Volume!=Data.nz_Volume)
	{
		REPORT_ERROR(1, "Projection_real_shears::ROUT_project_real_shears: "
					 "Error returned by ROUT_project_real_shears : the dimensions of the volume specified are not correct");
		return (ERROR);
	}

	if(Data.nx_Volume!=prm.proj_Xdim || Data.ny_Volume!=prm.proj_Ydim || Data.nz_Volume!=prm.proj_Zdim)
	{
		cout<<"\n\tWarning : the dimension specified in the input file is different to the volume dimension.";
		cout<<"\n\tThe program will only keep the volume dimensions.\n"<<endl;
	}

	//Does dynamics allocations and sets values
	allocAndInit_VolumeStruct(Data);

	//Filling of the Data Volume field :
	double *temp = MULTIDIM_ARRAY(Volume);
	int arraySize = x*y*z;
	for(int i=0; i<arraySize; i++)
		Data.Volume[i] = temp[i]; 


	SF.clear();
	SF.reserve(DF.dataLineNo());

	num_file = prm.starting;
	DF.go_first_data_line();

	return (!ERROR);
}

///Does finish instructions. Returns possibles errors.
int Projection_real_shears::finish_to_process()
{
	//Destruction of the dynamics allocations of Data
	del_VolumeStruct(Data);

	//SelFile save
	if(prog_param.display)
	{
		if (prog_param.fn_sel_file == "") //If the name of the output file is not specified
		{
			prog_param.fn_sel_file = "sel"+prm.fnProjectionSeed+".sel";
			cout<<"\n\tOutput file : "+prog_param.fn_sel_file<<endl;
		}

		cout<<endl;
		cout<<"\t**************\n";
		cout<<"\tSelFile save :"<<endl;
		SF.write(prog_param.fn_sel_file);
		cout<<"\t     OK"<<endl;
		cout<<"\t**************\n"<<endl;
	}
	else
	{
		if (prog_param.fn_sel_file == "") //If the name of the output file is not specified
			prog_param.fn_sel_file = "sel"+prm.fnProjectionSeed+".sel";

		SF.write(prog_param.fn_sel_file);
	}

	return (!ERROR);
}

//-------------------------------------- Main function ----------------------------------------
///Main function of the projection_real_shears program.
int Projection_real_shears::ROUT_project_real_shears()
{
	if(start_to_process() == ERROR)
		return (ERROR);

	while(!DF.eof())
	{
		if(prog_param.display) cout<<"\tProjection "<<num_file<<"..."<<endl;

		if(read_a_DocLine() == ERROR)
			return (ERROR);
		
		if(do_oneProjection(Data) == ERROR)
			return (ERROR);

		if(write_projection_file(num_file) == ERROR)
			return (ERROR);

		FreeVolumeDouble(&Data.Output);

		num_file++;
		DF.next_data_line();
	}

	if(finish_to_process() == ERROR)
		return (ERROR);

	return(!ERROR);	
}


