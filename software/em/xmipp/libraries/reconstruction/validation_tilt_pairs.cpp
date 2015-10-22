/***************************************************************************
 *
 * Authors:    Jose Luis Vilas          (jlvilas@cnb.csic.es)
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
#include "validation_tilt_pairs.h"
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <complex>
//#include <cmath>

/*
 * xmipp_validation_tilt_pairs --tilt /home/vilas/ScipionUserData/projects/rct/Runs/001623_XmippProtValidateTilt/extra/tilted/angles_iter001_00.xmd --untilt /home/vilas/ScipionUserData/projects/rct/Runs/001623_XmippProtValidateTilt/extra/untilted/angles_iter001_00.xmd -o caca
 * */

//Define Program parameters
void ProgValidationTiltPairs::defineParams()
{
    //Usage
    addUsageLine("Takes two coordinates sets and defines the coordinate transformation between them");
	addUsageLine("First set defines the untilted coordinates, second set defines the tilted coordinates");
	addParamsLine(" --tilt <metadata> : Metadata with angular assignment for the tilted images");
	addParamsLine(" --untilt <metadata> : Metadata with angular assignment for the untilted images");
	addParamsLine(" -o <metadata> : Metadata with matrix transformation");
}



//Read params
void ProgValidationTiltPairs::readParams()
{
    fntiltimage_In = getParam("--tilt");  //Set of tilted coordinates
    fnuntiltimage_In = getParam("--untilt");
	fnOut = getParam("-o");  //Output file
}


void ProgValidationTiltPairs::quaternion2Paulibasis(double rot, double tilt, double psi, std::complex<double> (&L)[4])
{
	double cr, ct, cp, sr, st, sp;

	cr = cos(rot/2);
	ct = cos(tilt/2);
	cp = cos(psi/2);
	sr = sin(rot/2);
	st = sin(tilt/2);
	sp = sin(psi/2);

	L[0] = cr*ct*cp - sr*ct*sp;
	L[1] = 1i*(sr*st*cp - cr*st*sp);
	L[2] = 1i*(cr*st*cp + sr*st*sp);
	L[3] = 1i*(sr*ct*cp+cr*ct*sp);
}


void ProgValidationTiltPairs::matrix2Paulibasis(std::complex<double> M[4],
		std::complex<double> (&P)[4])
{
	//M[0] = m11; M[1]=m12; M[2]=m21; M[3]=m22

	std::complex<double> I=1i;
	std::complex<double> aux=0.5;
	P[0]=(M[0]+M[3])*aux;
	P[1]=(M[1]+M[2])*aux;
	P[2]=I*(M[1]-M[2])*aux;
	P[3]=(M[0]-M[3])*aux;
}

void ProgValidationTiltPairs::InversefromPaulibasis(std::complex<double> Original[4],
		std::complex<double> (&Inver)[4])
{
	//It takes a Pauli expression and returns its inverse expressed in Pauli basis

	//TODO Raise an exception is the Original matrix does not belong to SU(2) group

	std::complex<double> Inver_matrix[4], mat[4];
	std::complex<double> NOriginal;
	double aux=0.5;

	Paulibasis2matrix(Original,mat);

	Inver_matrix[0] = mat[3];
	Inver_matrix[1] = -mat[1];
	Inver_matrix[2] = -mat[2];
	Inver_matrix[3] = mat[0];

	matrix2Paulibasis(Inver_matrix,Inver);
}

void ProgValidationTiltPairs::inverse_matrixSU2(std::complex<double> Original[4],
		std::complex<double> (&Inver)[4])
{
	//It takes a matrix and returns its inverse expressed in Pauli basis

	//TODO Raise an exception is the Original matrix does not belong to SU(2) group

	std::complex<double> Inver_matrix[4];

	Inver_matrix[0] = Original[3];
	Inver_matrix[1] = -Original[1];
	Inver_matrix[2] = -Original[2];
	Inver_matrix[3] = Original[0];

	matrix2Paulibasis(Inver_matrix,Inver);
}

void ProgValidationTiltPairs::Paulibasis2matrix(std::complex<double> P[4], std::complex<double> (&M)[4])
{
	std::complex<double> I=1i;
	M[0] = (P[0]+P[3]);
	M[1] = (P[1]-I*P[2]);
	M[2] = (P[1]+I*P[2]);
	M[3] = (P[0]-P[3]);
}

void ProgValidationTiltPairs::Pauliproduct(std::complex<double> A[4], std::complex<double> B[4],
		std::complex<double> (&P)[4])
{
	std::complex<double> A_matrix[4], B_matrix[4], aux[4];

	Paulibasis2matrix(A,A_matrix);
	Paulibasis2matrix(B,B_matrix);

	aux[0] = A_matrix[0]*B_matrix[0] + A_matrix[1]*B_matrix[2];
	aux[2] = A_matrix[2]*B_matrix[0] + A_matrix[3]*B_matrix[2];
	aux[1] = A_matrix[0]*B_matrix[1] + A_matrix[1]*B_matrix[3];
	aux[3] = A_matrix[2]*B_matrix[1] + A_matrix[3]*B_matrix[3];

	matrix2Paulibasis(aux, P);
}

void ProgValidationTiltPairs::extrarotationangles(std::complex<double> R[4], double &alpha_x, double &alpha_y)
{
	std::complex<double> I=1i;
	std::complex<double> aux1 = I*R[1]/R[0],  aux2 = I*R[2]/R[0], alpha_aux_x, alpha_aux_y;

	if ((aux1.imag() == 0) && (aux2.imag() == 0))
	{
		alpha_aux_x = 2*atan(aux1.real());
		alpha_aux_y = 2*atan(aux2.real());
		alpha_x = alpha_aux_x.real()*180/PI;
		alpha_y = alpha_aux_y.real()*180/PI;
	}
	else
	{
		std::cout << "Error, argument of atan is complex" << std::endl;
	}
}


void ProgValidationTiltPairs::angles2tranformation(double untilt_angles[3],
		double tilt_angles[3], double alpha_x, double alpha_y)
{
	double rotu = untilt_angles[0], tiltu = untilt_angles[1], psiu = untilt_angles[2];
	double rott = tilt_angles[0], tiltt = tilt_angles[1], psit = tilt_angles[2];
	std::complex<double> qu[4], M[4], Inv_qu[4], qt[4], R[4];

	quaternion2Paulibasis(rotu, tiltu, psiu, qu);
    Paulibasis2matrix(qu,M);
    inverse_matrixSU2(M, Inv_qu);
    quaternion2Paulibasis(rott, tiltt, psit, qt);
    Pauliproduct(qt, Inv_qu, R);

    extrarotationangles(R, alpha_x, alpha_y);
    std::cout << "alpha_x = " << alpha_x << std::endl;
    std::cout << "alpha_y = " << alpha_y << std::endl;
}

void ProgValidationTiltPairs::run()
{
	MetaData MD_tilted, MD_untilted, DF1sorted, DF2sorted, DFweights;

	MD_tilted.read(fntiltimage_In);
	MD_untilted.read(fnuntiltimage_In);

	DF1sorted.sort(MD_tilted,MDL_ITEM_ID,true);
	DF2sorted.sort(MD_untilted,MDL_ITEM_ID,true);

	MDIterator iter1(DF1sorted), iter2(DF2sorted);
	std::vector< Matrix1D<double> > ang1, ang2;
	Matrix1D<double> rotTiltPsi(3), z(3);
	size_t currentId;
	bool anotherImageIn2=iter2.hasNext();
	size_t id1, id2;
	bool mirror;
	Matrix2D<double> Eu, Et, R;
	double alpha, beta;
	while (anotherImageIn2)
	{
		ang1.clear();
		ang2.clear();

		// Take current id
		DF2sorted.getValue(MDL_ITEM_ID,currentId,iter2.objId);

		// Grab all the angles in DF2 associated to this id
		bool anotherIteration=false;
		do
		{
			DF2sorted.getValue(MDL_ITEM_ID,id2,iter2.objId);
			anotherIteration=false;
			if (id2==currentId)
			{
				DF2sorted.getValue(MDL_ANGLE_ROT,XX(rotTiltPsi),iter2.objId);
				DF2sorted.getValue(MDL_ANGLE_TILT,YY(rotTiltPsi),iter2.objId);
				DF2sorted.getValue(MDL_ANGLE_PSI,ZZ(rotTiltPsi),iter2.objId);
				DF2sorted.getValue(MDL_FLIP,mirror,iter2.objId);
				std::cout << "From DF2:" << XX(rotTiltPsi) << " " << YY(rotTiltPsi) << " " << ZZ(rotTiltPsi) << " " << mirror << std::endl;
				//LINEA ANTERIOR ORIGINAL
				if (mirror)
				{
					double rotp, tiltp, psip;
					Euler_mirrorX(XX(rotTiltPsi),YY(rotTiltPsi),ZZ(rotTiltPsi), rotp, tiltp, psip);
					XX(rotTiltPsi)=rotp;
					YY(rotTiltPsi)=tiltp;
					ZZ(rotTiltPsi)=psip;
				}
				ang2.push_back(rotTiltPsi);
				iter2.moveNext();
				if (iter2.hasNext())
					anotherIteration=true;
			}
		} while (anotherIteration);

		// Advance Iter 1 to catch Iter 2
		double N=0, cumulatedDistance=0;
		size_t newObjId=0;
		if (iter1.objId>0)
		{
			DF1sorted.getValue(MDL_ITEM_ID,id1,iter1.objId);
			while (id1<currentId && iter1.hasNext())
			{
				iter1.moveNext();
				DF1sorted.getValue(MDL_ITEM_ID,id1,iter1.objId);
			}

			// If we are at the end of DF1, then we did not find id1 such that id1==currentId
			if (!iter1.hasNext())
				break;

			// Grab all the angles in DF1 associated to this id
			anotherIteration=false;
			do
			{
				DF1sorted.getValue(MDL_ITEM_ID,id1,iter1.objId);
				anotherIteration=false;
				if (id1==currentId)
				{
					DF1sorted.getValue(MDL_ANGLE_ROT,XX(rotTiltPsi),iter1.objId);
					DF1sorted.getValue(MDL_ANGLE_TILT,YY(rotTiltPsi),iter1.objId);
					DF1sorted.getValue(MDL_ANGLE_PSI,ZZ(rotTiltPsi),iter1.objId);
					DF1sorted.getValue(MDL_FLIP,mirror,iter1.objId);
					std::cout << "From DF1:" << XX(rotTiltPsi) << " " << YY(rotTiltPsi) << " " << ZZ(rotTiltPsi) << " " << mirror << std::endl;
					//LINEA ANTERIOR ORIGINAL
					if (mirror)
					{
						double rotp, tiltp, psip;
						Euler_mirrorX(XX(rotTiltPsi),YY(rotTiltPsi),ZZ(rotTiltPsi), rotp, tiltp, psip);
						XX(rotTiltPsi)=rotp;
						YY(rotTiltPsi)=tiltp;
						ZZ(rotTiltPsi)=psip;
					}
					ang1.push_back(rotTiltPsi);
					iter1.moveNext();
					if (iter1.hasNext())
						anotherIteration=true;
				}
			} while (anotherIteration);

			// Process both sets of angles
			for (size_t i=0; i<ang2.size(); ++i)
			{
				const Matrix1D<double> &anglesi=ang2[i];
				double rotu=XX(anglesi);
				double tiltu=YY(anglesi);
				double psiu=ZZ(anglesi);
				Euler_angles2matrix(rotu,tiltu,psiu,Eu,false);
				/*std::cout << "------UNTILTED MATRIX------" << std::endl;
				std::cout << Eu << std::endl;
				std::cout << "vector" << std::endl;
				std::cout << Eu(2,0) << "  " << Eu(2,1) << "  " << Eu(2,2) << std::endl;*/


				for (size_t j=0; j<ang1.size(); ++j)
				{
					const Matrix1D<double> &anglesj=ang1[j];
					double rott=XX(anglesj);
					double tiltt=YY(anglesj);
					double psit=ZZ(anglesj);
					double alpha_x, alpha_y;
					Euler_angles2matrix(rott,tiltt,psit,Et,false);
					//////////////////////////////////////////////////////////////////
					double untilt_angles[3]={rotu, tiltu, psiu}, tilt_angles[3]={rott, tiltt, psit};
					angles2tranformation(untilt_angles, tilt_angles, alpha_x, alpha_y);
					//std::cout << "alpha = " << (alpha_x*alpha_x+alpha_y*alpha_y) << std::endl;
					//////////////////////////////////////////////////////////////////
					/*std::cout << "------TILTED MATRIX------" << std::endl;
					std::cout << Et << std::endl;
					std::cout << "vector" << std::endl;
					std::cout << Et(2,0) << "  " << Et(2,1) << "  " << Et(2,2) << std::endl;
					std::cout << "---------------------------" << std::endl;
					std::cout << "---------------------------" << std::endl;*/
					R=Eu*Et.transpose();
					double rotTransf, tiltTransf, psiTransf;
					Euler_matrix2angles(R, rotTransf, tiltTransf, psiTransf);
					std::cout << "Rot_and_Tilt " << rotTransf << " " << tiltTransf << std::endl;
					//LINEA ANTERIOR ORIGINAL

				XX(z) = Eu(2,0) - Et(2,0);
				YY(z) = Eu(2,1) - Et(2,1);
				ZZ(z) = Eu(2,2) - Et(2,2);

				alpha = atan2(YY(z), XX(z));        //Expressed in rad
				beta = atan2(XX(z)/cos(alpha), ZZ(z));   //Expressed in rad
				std::cout << "alpha = " << alpha*180/PI << std::endl;
				std::cout << "beta = " << beta*180/PI << std::endl;
				}
			}
		}
		else
			N=0;

		if (N>0)
		{
			double meanDistance=cumulatedDistance/ang2.size();
			DFweights.setValue(MDL_ANGLE_DIFF,meanDistance,newObjId);
		}
		else
			if (newObjId>0)
				DFweights.setValue(MDL_ANGLE_DIFF,-1.0,newObjId);
		anotherImageIn2=iter2.hasNext();
	}

	std::complex<double> qu[4], qt[4], M[4], Inv_qu[4], test[4], P1[4], P2[4], Inv_quu[4];
	double rotu=34*PI/180, tiltu=10*PI/180, psiu=5*PI/180;
	double rott=25*PI/180, tiltt=15*PI/180, psit=40*PI/180;


    quaternion2Paulibasis(rotu, tiltu, psiu, qu);
    /*std::cout << "quaternion2Pauli" << std::endl;
    std::cout << "Untilted " << qu[0] << " " << qu[1] << " " << qu[2] << " " << qu[3] << std::endl;
    std::cout << "      " << std::endl;*/

    Paulibasis2matrix(qu,M);
    /*std::cout << "Pauli2matrix" << std::endl;
    std::cout << "Matriz   " << M[0] << " " << M[1] << " " << M[2] << " " << M[3] << std::endl;
    std::cout << "      " << std::endl;*/

    inverse_matrixSU2(M, Inv_qu);
    /*std::cout << "inverse_matrixSU(2)" << std::endl;
    std::cout << "Inversa  " << Inv_qu[0] << " " << Inv_qu[1] << " " << Inv_qu[2] << " " << Inv_qu[3] << std::endl;
    std::cout << "      " << std::endl;*/

    quaternion2Paulibasis(rott, tiltt, psit, qt);
    /*std::cout << "quaternion2Pauli" << std::endl;
    std::cout << "Tilted " << qt[0] << " " << qt[1] << " " << qt[2] << " " << qt[3] << std::endl;
    std::cout << "      " << std::endl;*/

    InversefromPaulibasis(qu,Inv_quu);

    Pauliproduct(qt, Inv_qu, P1);
    /*std::cout << "Pauliproduct" << std::endl;
    std::cout << "quaternion qt  " << P1[0] << " " << P1[1] << " " << P1[2] << " " << P1[3] << std::endl;
    std::cout << "      " << std::endl;
    std::cout << "-----------------------------------" << std::endl;*/

    //double alpha_x, alpha_y;
    //extrarotationangles(P1, alpha_x, alpha_y);
    //std::cout << "alpha_x = " << alpha_x << " " << "alpha_y = " << alpha_y << std::endl;
}



