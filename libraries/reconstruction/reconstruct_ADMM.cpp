/***************************************************************************
 *
 * Authors:        Masih Nilchian (masih.nilchian@epfl.ch)
 *                 Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

/* TODO:
 * - Default parameters
 * - Parallelization
 * - Symmetrize HtKH and Htb
 */

#include "reconstruct_ADMM.h"
#include "symmetrize.h"
#include <data/metadata_extension.h>

#define SYMMETRIZE_PROJECTIONS

ProgReconsADMM::ProgReconsADMM(): XmippProgram()
{
	rank=0;
	Nprocs=1;
}

void ProgReconsADMM::defineParams()
{
    addUsageLine("Reconstruct with Alternative Direction Method of Multipliers");
    addParamsLine(" -i <metadata>: Input metadata with images and angles");
    addParamsLine(" [--oroot <root=reconstruct_admm>]: Rootname for output files");
    addParamsLine(" [--Htb <volume=\"\">]: Filename of the Htb volume");
    addParamsLine(" [--HtKH <volume=\"\">]: Filename of the HtKH volume");
    addParamsLine(" [--firstVolume <volume=\"\">]: Filename of the first initial guess");
    addParamsLine(" [--kernel <shape=KaiserBessel>]: Kernel shape");
    addParamsLine("    where <shape>");
    addParamsLine("      KaiserBessel: Kaiser-Bessel function as kernel");
    addParamsLine(" [--downsamplingV <Ti=1>]: Downsampling factor for volume");
    addParamsLine(" [--downsamplingI <Tp=1>]: Downsampling factor for projections");
    addParamsLine(" [--dontUseWeights]: Do not use weights if available in the input metadata");
    addParamsLine(" [--dontUseCTF]: Do not use CTF if available in the input metadata");
    addParamsLine("               : If the CTF information is used, the algorithm assumes that the input images are phase flipped");
    addParamsLine(" [--sampling <Ts=1>]: The sampling rate is only used to correct for the CTF");
    addParamsLine(" [--mu <mu=1e-4>]: Augmented Lagrange penalty");
    addParamsLine(" [--lambda <lambda=1e-5>]: Total variation parameter");
    addParamsLine(" [--lambda1 <lambda1=1e-7>]: Tikhonov parameter");
    addParamsLine(" [--cgiter <N=3>]: Conjugate Gradient iterations");
    addParamsLine(" [--admmiter <N=30>]: ADMM iterations");
    addParamsLine(" [--positivity]: Positivity constraint");
    addParamsLine(" [--sym <s=c1>]: Symmetry constraint");
    addParamsLine(" [--saveIntermediate]: Save Htb and HtKH volumes for posterior calls");

    mask.defineParams(this,INT_MASK);
}

void ProgReconsADMM::readParams()
{
	fnIn=getParam("-i");
	fnRoot=getParam("--oroot");
	fnFirst=getParam("--firstVolume");
	fnHtb=getParam("--Htb");
	fnHtKH=getParam("--HtKH");
	kernelShape=getParam("--kernel");
	if (kernelShape=="KaiserBessel")
	{
		a=4; //getDoubleParam("--kernel",1);
		alpha=19; // getDoubleParam("--kernel",2);
	}
	Ti=getDoubleParam("--downsamplingV");
	Tp=getDoubleParam("--downsamplingI");
	useWeights=!checkParam("--dontUseWeights");
	useCTF=!checkParam("--dontUseCTF");
	Ts=getDoubleParam("--sampling");
	mu=getDoubleParam("--mu");
	lambda=getDoubleParam("--lambda");
	lambda1=getDoubleParam("--lambda1");
	Ncgiter=getIntParam("--cgiter");
	Nadmmiter=getIntParam("--admmiter");
	positivity=checkParam("--positivity");
	saveIntermediate=checkParam("--saveIntermediate");

	applyMask=checkParam("--mask");
	if (applyMask)
		mask.readParams(this);
	SL.readSymmetryFile(getParam("--sym"));
}

void ProgReconsADMM::produceSideInfo()
{
	// Read input images
	mdIn.read(fnIn);

	// Symmetrize input images
#ifdef SYMMETRIZE_PROJECTIONS
	Matrix2D<double> Lsym(3,3), Rsym(3,3);
	MetaData mdSym;
	MDRow row;
	double rot, tilt, psi, newrot, newtilt, newpsi;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getRow(row,__iter.objId);
		row.getValue(MDL_ANGLE_ROT,rot);
		row.getValue(MDL_ANGLE_TILT,tilt);
		row.getValue(MDL_ANGLE_PSI,psi);
		mdSym.addRow(row);
		for (int isym = 0; isym < SL.symsNo(); isym++)
		{
			SL.getMatrices(isym, Lsym, Rsym, false);
			Euler_apply_transf(Lsym,Rsym,rot,tilt,psi,newrot,newtilt,newpsi);
			row.setValue(MDL_ANGLE_ROT,newrot);
			row.setValue(MDL_ANGLE_TILT,newtilt);
			row.setValue(MDL_ANGLE_PSI,newpsi);
			mdSym.addRow(row);
		}
	}
	mdIn=mdSym;
	mdSym.clear();
#endif

	// Prepare kernel
	if (kernelShape=="KaiserBessel")
		kernel.initializeKernel(alpha,a,0.0001);
	kernel.convolveKernelWithItself(0.005);

	// Get Htb reconstruction
	if (fnHtb=="")
	{
		constructHtb();
		if (saveIntermediate && rank==0)
			CHtb.write(fnRoot+"_Htb.vol");
	}
	else
	{
		CHtb.read(fnHtb);
		CHtb().setXmippOrigin();
	}

	// Adjust regularization weights
	//double Nimgs=mdIn.size();
	//double CHtb_norm=sqrt(CHtb().sum2()/MULTIDIM_SIZE(CHtb()));
	//COSS mu*=Nimgs/XSIZE(CHtb())*CHtb_norm;
	//COSS lambda*=Nimgs/XSIZE(CHtb())*CHtb_norm;
	//COSS lambda1*=Nimgs/XSIZE(CHtb());

	// Get first volume
	if (fnFirst=="")
		Ck().initZeros(CHtb());
	else
	{
		Ck.read(fnFirst);
		Ck().setXmippOrigin();
	}

	// Compute H'*K*H+mu*L^T*L+lambda_1 I
	Image<double> kernelV;
	if (fnHtKH=="")
	{
		computeHtKH(kernelV());
		if (saveIntermediate && rank==0)
			kernelV.write(fnRoot+"_HtKH.vol");
	}
	else
	{
		kernelV.read(fnHtKH);
		kernelV().setXmippOrigin();
	}
	FourierTransformer transformer;
	transformer.FourierTransform(kernelV(),fourierKernelV,true);

	// Add regularization in Fourier space
	addRegularizationTerms();

	// Resize u and d volumes
	ux.initZeros(CHtb());
	ux.setXmippOrigin();
	uy=ux;
	uz=ux;
	dx=ux;
	dy=ux;
	dz=ux;

	// Prepare Lt filters
	MultidimArray<double> L;
	L.resize(ux);
	kernel.computeGradient(L,'x');
	transformer.FourierTransform(L,fourierLx);
	kernel.computeGradient(L,'y');
	transformer.FourierTransform(L,fourierLy);
	kernel.computeGradient(L,'z');
	transformer.FourierTransform(L,fourierLz);

	ud=ux;
	transformerL.setReal(ud);

	// Prepare mask
	if (applyMask)
		mask.generate_mask(CHtb());
	synchronize();
}

void ProgReconsADMM::show()
{

}

void ProgReconsADMM::run()
{
	produceSideInfo();
	if (rank==0)
	{
		std::cout << "Running ADMM iterations ..." << std::endl;
		init_progress_bar(Nadmmiter);
		for (int iter=0; iter<Nadmmiter; ++iter)
		{
			applyConjugateGradient();
			doPOCSProjection();
			updateUD();
			progress_bar(iter);
		}
		progress_bar(Nadmmiter);
		produceVolume();
		Vk.write(fnRoot+".vol");
	}
	synchronize();
}

void ProgReconsADMM::constructHtb()
{
	size_t xdim, ydim, zdim, ndim;
	getImageSize(mdIn,xdim,ydim,zdim,ndim);
	CHtb().initZeros(xdim,xdim,xdim);
	CHtb().setXmippOrigin();

	Image<double> I;
	double rot, tilt, psi;
	size_t i=0;
	if (rank==0)
	{
		std::cerr << "Performing Htb ...\n";
		init_progress_bar(mdIn.size());
	}
	double weight=1.;
	ApplyGeoParams geoParams;
	geoParams.only_apply_shifts=true;
	geoParams.wrap=DONT_WRAP;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		if ((i+1)%Nprocs==rank)
		{
			I.readApplyGeo(mdIn,__iter.objId,geoParams);
			I().setXmippOrigin();
			mdIn.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
			mdIn.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
			mdIn.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
			if (useWeights && mdIn.containsLabel(MDL_WEIGHT))
				mdIn.getValue(MDL_WEIGHT,weight,__iter.objId);

			project(rot,tilt,psi,I(),true,weight);
		}
		i++;
		if (i%100==0 && rank==0)
			progress_bar(i);
	}
	if (rank==0)
		progress_bar(mdIn.size());

	// Symmetrize Htb
#ifndef SYMMETRIZE_PROJECTIONS
	MultidimArray<double> CHtbsym;
	symmetrizeVolume(SL, CHtb(), CHtbsym, false, false, true);
	CHtb=CHtbsym;
#endif
	shareVolume(CHtb());
}

void ProgReconsADMM::project(double rot, double tilt, double psi, MultidimArray<double> &P, bool adjoint, double weight)
{
	Matrix2D<double> E;
	Euler_angles2matrix(rot,tilt,psi,E,false);
	Matrix1D<double> r1(3), r2(3);
	E.getRow(0,r1);
	E.getRow(1,r2);
	project(r1,r2,P,adjoint,weight);
}

void ProgReconsADMM::project(const Matrix1D<double> &r1, const Matrix1D<double> &r2, MultidimArray<double> &P, bool adjoint, double weight)
{
	const MultidimArray<double> &mV=CHtb();
	if (!adjoint)
	{
		P.initZeros(std::ceil(XSIZE(mV)/Tp),std::ceil(YSIZE(mV)/Tp));
		P.setXmippOrigin();
	}

    // Tomographic projection
    for (int k=STARTINGZ(mV); k<=FINISHINGZ(mV); ++k) {
        // initial computation
        double rzn  = k;
        double rzn1= rzn*VEC_ELEM(r1,2);
        double rzn2= rzn*VEC_ELEM(r2,2);
        for (int i=STARTINGY(mV); i<=FINISHINGY(mV); ++i) {
            // initial computation
            double ryn  = i;
            double ryn1= ryn*VEC_ELEM(r1,1);
            double ryn2= ryn*VEC_ELEM(r2,1);
            double sx0=rzn1+ryn1;
            double sy0=rzn2+ryn2;

            for (int j=STARTINGX(mV); j<=FINISHINGX(mV); ++j) {
                double rxn  = j;
                double rxn1= rxn*VEC_ELEM(r1,0);
                double rxn2= rxn*VEC_ELEM(r2,0);

                double sx = Ti*(sx0+rxn1);
                double sy = Ti*(sy0+rxn2);

                int sxmin = std::floor(sx-kernel.supp);
                int sxmax = std::ceil(sx+kernel.supp);
                int symin = std::floor(sy-kernel.supp);
                int symax = std::ceil(sy+kernel.supp);
                if (sxmin<STARTINGX(P))
                    sxmin = STARTINGX(P);
                if (sxmax>FINISHINGX(P))
                    sxmax = FINISHINGX(P);
                if (symin<STARTINGY(P))
                    symin = STARTINGY(P);
                if (symax>FINISHINGY(P))
                    symax = FINISHINGY(P);

                if (adjoint)
					for (int ii=symin; ii<=symax; ii++) {
						double u=(ii-sy)*Tp;
						for (int jj=sxmin; jj<=sxmax; jj++) {
							double v=(jj-sx)*Tp;
							A3D_ELEM(mV,k,i,j)+=weight*A2D_ELEM(P,ii,jj)*kernel.projectionValueAt(u,v);
						}
					}
                else
					for (int ii=symin; ii<=symax; ii++) {
						double u=(ii-sy)*Tp;
						for (int jj=sxmin; jj<=sxmax; jj++) {
							double v=(jj-sx)*Tp;
							A2D_ELEM(P,ii,jj)+=weight*A3D_ELEM(mV,k,i,j)*kernel.projectionValueAt(u,v);
						}
					}
            }
        }
    }
}

void ProgReconsADMM::computeHtKH(MultidimArray<double> &kernelV)
{
	kernelV.initZeros(2*ZSIZE(CHtb())-1,2*YSIZE(CHtb())-1,2*XSIZE(CHtb())-1);
	kernelV.setXmippOrigin();

	if (rank==0)
	{
		std::cerr << "Calculating H'KH ...\n";
		init_progress_bar(mdIn.size());
	}
	size_t i=0;
	double rot, tilt, psi, weight=1;
	bool hasWeight=mdIn.containsLabel(MDL_WEIGHT) && useWeights;
	bool hasCTF=(mdIn.containsLabel(MDL_CTF_MODEL) || mdIn.containsLabel(MDL_CTF_DEFOCUSU)) && useCTF;
	CTFDescription ctf;
	MultidimArray<double> kernelAutocorr;
	if (!hasCTF)
		kernel.getKernelAutocorrelation(kernelAutocorr);
	Matrix2D<double> E;
	Matrix1D<double> r1(3), r2(3);
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		if ((i+1)%Nprocs==rank)
		{
			// COSS: Read also MIRROR
			mdIn.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
			mdIn.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
			mdIn.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
			if (hasWeight)
				mdIn.getValue(MDL_WEIGHT,weight,__iter.objId);
			if (hasCTF)
			{
				ctf.readFromMetadataRow(mdIn,__iter.objId);
				ctf.produceSideInfo();
				kernel.applyCTFToKernelAutocorrelation(ctf,Ts,kernelAutocorr);
			}
			kernelAutocorr.setXmippOrigin();

			// Update kernel
			Euler_angles2matrix(rot,tilt,psi,E,false);
			E.getRow(0,r1);
			E.getRow(1,r2);
			double iStep=1.0/kernel.autocorrStep;
			for (int k=((kernelV).zinit); k<=((kernelV).zinit + (int)(kernelV).zdim - 1); ++k)
			{
				double r1_z=k*ZZ(r1);
				double r2_z=k*ZZ(r2);
				for (int i=((kernelV).yinit); i<=((kernelV).yinit + (int)(kernelV).ydim - 1); ++i)
				{
					double r1_yz=i*YY(r1)+r1_z;
					double r2_yz=i*YY(r2)+r2_z;
					for (int j=((kernelV).xinit); j<=((kernelV).xinit + (int)(kernelV).xdim - 1); ++j)
					{
						double r1_xyz=j*XX(r1)+r1_yz;
						double r2_xyz=j*XX(r2)+r2_yz;
						A3D_ELEM(kernelV,k,i,j)+=weight*kernelAutocorr.interpolatedElement2D(r1_xyz*iStep,r2_xyz*iStep);
					}
				}
			}
		}

		i++;
		if (i%100==0 && rank==0)
			progress_bar(i);
	}
	if (rank==0)
		progress_bar(mdIn.size());

#ifndef SYMMETRIZE_PROJECTIONS
	MultidimArray<double> kernelVsym;
	symmetrizeVolume(SL, kernelV, kernelVsym, false, false, true);
	kernelV=kernelVsym;
#endif

	shareVolume(kernelV);
}

void addGradientTerm(double mu, AdmmKernel &kernel, MultidimArray<double> &L, FourierTransformer &transformer,
		MultidimArray<std::complex<double> >&fourierKernelV, char direction)
{
	kernel.computeGradient(L,direction);

	transformer.FourierTransform();
	double K=mu*MULTIDIM_SIZE(L);
	double xdim_2=(double)(XSIZE(L)/2);
	double Kargument=2*PI*xdim_2/XSIZE(L);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(fourierKernelV)
	{
		// Calculation of L^t*L
		// -L^2(r)*exp(-i*2*pi*(N/2)/N*(k+i+j)): the last term is a phase shift due to FFT
		double argument=Kargument*(k+i+j);
		double s, c;
		sincos(argument,&s,&c);
		A3D_ELEM(fourierKernelV,k,i,j)+=-A3D_ELEM(transformer.fFourier,k,i,j)*
				A3D_ELEM(transformer.fFourier,k,i,j)*std::complex<double>(K*c,K*s);
	}
}

void ProgReconsADMM::addRegularizationTerms()
{
	MultidimArray<double> L;
	L.initZeros(2*ZSIZE(CHtb())-1,2*YSIZE(CHtb())-1,2*XSIZE(CHtb())-1);
	L.setXmippOrigin();
	FourierTransformer transformer;
	transformer.setReal(L);

	if (mu>0)
	{
		addGradientTerm(mu, kernel,L,transformer,fourierKernelV, 'x');
		addGradientTerm(mu, kernel,L,transformer,fourierKernelV, 'y');
		addGradientTerm(mu, kernel,L,transformer,fourierKernelV, 'z');
	}

	if (lambda1>0)
	{
		L.initZeros();
		L(0,0,0)=lambda1;
		transformer.FourierTransform();
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fourierKernelV)
			DIRECT_MULTIDIM_ELEM(fourierKernelV,n)+=DIRECT_MULTIDIM_ELEM(transformer.fFourier,n);
	}
}

void ProgReconsADMM::applyKernel3D(MultidimArray<double> &x, MultidimArray<double> &AtAx)
{
	paddedx.initZeros(2*ZSIZE(x)-1,2*YSIZE(x)-1,2*XSIZE(x)-1);
	paddedx.setXmippOrigin();

	// Copy Vk into paddedx
	for (int k=STARTINGZ(x); k<=FINISHINGZ(x); ++k)
		for (int i=STARTINGY(x); i<=FINISHINGY(x); ++i)
			memcpy(&A3D_ELEM(paddedx,k,i,STARTINGX(x)),&A3D_ELEM(x,k,i,STARTINGX(x)),XSIZE(x)*sizeof(double));

	// Compute Fourier transform of paddedx
	transformerPaddedx.setReal(paddedx);
	transformerPaddedx.FourierTransform();

	// Apply kernel
	double K=MULTIDIM_SIZE(paddedx);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(transformerPaddedx.fFourier)
	{
		DIRECT_MULTIDIM_ELEM(transformerPaddedx.fFourier,n)*=DIRECT_MULTIDIM_ELEM(fourierKernelV,n);
		DIRECT_MULTIDIM_ELEM(transformerPaddedx.fFourier,n)*=K;
	}

	// Inverse Fourier transform
	transformerPaddedx.inverseFourierTransform();
	CenterFFT(paddedx,false);

	// Crop central region
	AtAx.resize(x);
	for (int k=STARTINGZ(x); k<=FINISHINGZ(x); ++k)
		for (int i=STARTINGY(x); i<=FINISHINGY(x); ++i)
			memcpy(&A3D_ELEM(AtAx,k,i,STARTINGX(x)),&A3D_ELEM(paddedx,k,i,STARTINGX(x)),XSIZE(x)*sizeof(double));
}

void ProgReconsADMM::applyLFilter(MultidimArray< std::complex<double> > &fourierL, bool adjoint)
{
	transformerL.FourierTransform();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fourierL)
		DIRECT_MULTIDIM_ELEM(transformerL.fFourier,n)*=DIRECT_MULTIDIM_ELEM(fourierL,n);
	transformerL.inverseFourierTransform();
	CenterFFT(ud,false);
	if (adjoint)
		ud*=-1;
}

void ProgReconsADMM::applyLtFilter(MultidimArray< std::complex<double> > &fourierL, MultidimArray<double> &u, MultidimArray<double> &d)
{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(ud)
		DIRECT_MULTIDIM_ELEM(ud,n)=DIRECT_MULTIDIM_ELEM(u,n)-DIRECT_MULTIDIM_ELEM(d,n);
	applyLFilter(fourierL,true);
}

void ProgReconsADMM::applyConjugateGradient()
{
	// Conjugate gradient is solving Ax=b, when A is symmetric (=> positive semidefinite) and we apply it
	// to the problem
	// (H^T K H+mu L^TL + lambda I) x = H^Tb + mu L^T(u-d)
	// Being H the projection operator
	//       K the CTF operator
	//       L the gradient operator

	// Compute H^tb+mu*L^t(u-d)
	MultidimArray<double> r;
	applyLtFilter(fourierLx,ux,dx);
	r=ud;
	applyLtFilter(fourierLy,uy,dy);
	r+=ud;
	applyLtFilter(fourierLz,uz,dz);
	r+=ud;
	r*=mu;

	r+=CHtb();

	// Apply A^tA to the current estimate of the reconstruction
	MultidimArray<double> AtAVk;
	applyKernel3D(Ck(),AtAVk);

	// Compute first residual. This is the negative gradient of ||Ax-b||^2
	r-=AtAVk;
	double d=r.sum2();

	// Search direction
	MultidimArray<double> p, AtAp;
	p=r;

	// Perform CG iterations
	MultidimArray<double> &mCk=Ck();

	for (int iter=0; iter<Ncgiter; ++iter)
	{
		applyKernel3D(p,AtAp);
		double alpha=d/p.dotProduct(AtAp);

		// Update residual and current estimate
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(r)
		{
			DIRECT_MULTIDIM_ELEM(r,n)  -=alpha*DIRECT_MULTIDIM_ELEM(AtAp,n);
			DIRECT_MULTIDIM_ELEM(mCk,n)+=alpha*DIRECT_MULTIDIM_ELEM(p,n);
		}
		double newd=r.sum2();
		double beta=newd/d;
		d=newd;

		// Update search direction
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p)
			DIRECT_MULTIDIM_ELEM(p,n)  = DIRECT_MULTIDIM_ELEM(r,n)+beta*DIRECT_MULTIDIM_ELEM(p,n);

//		Image<double> save;
//		save()=Vk();
//		save.write("PPPVk.vol");
//		std::cout << "Press any key" << std::endl;
//		char c; std::cin >> c;
	}
}

void ProgReconsADMM::doPOCSProjection()
{
	if (applyMask)
		mask.apply_mask(Ck(),Ck());

	if (positivity)
	{
		Ck().threshold("below",0.,0.);
#ifdef NEVERDEFINED
		MultidimArray<double> smallKernel3D(generate_mask2*a+1,2*a+1,2*a+1);
		smallKernel3D.setXmippOrigin();
		kernel.computeKernel3D(smallKernel3D);
		smallKernel3D/=sqrt(smallKernel3D.sum2());

		produceVolume();
		MultidimArray<double> &mVk=Vk();
		MultidimArray<double> &mCk=Ck();
		Ck.write("PPPCkBefore.vol");
		int k0=STARTINGZ(mVk)+a;
		int i0=STARTINGY(mVk)+a;
		int j0=STARTINGX(mVk)+a;
		int kF=FINISHINGZ(mVk)-a;
		int iF=FINISHINGY(mVk)-a;
		int jF=FINISHINGX(mVk)-a;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(mVk)
			if (k<k0 || k>kF || i<i0 || i>iF || j<j0 || j>jF)
				A3D_ELEM(mCk,k,i,j)=0.;
			else if (A3D_ELEM(mVk,k,i,j)<0)
			{
				// Calculate dot product
				double dot=0.;
				for (int kk=-a; kk<=a; ++kk)
					for (int ii=-a; ii<=a; ++ii)
						for (int jj=-a; jj<=a; ++jj)
							dot+=A3D_ELEM(mCk,k+kk,i+ii,j+jj)*A3D_ELEM(smallKernel3D,kk,ii,jj);

				// Update C
				for (int kk=-a; kk<=a; ++kk)
					for (int ii=-a; ii<=a; ++ii)
						for (int jj=-a; jj<=a; ++jj)
							A3D_ELEM(mCk,k+kk,i+ii,j+jj)-=dot*A3D_ELEM(smallKernel3D,kk,ii,jj);
			}
		Ck.write("PPPCkAfter.vol");
#endif
	}
}

void ProgReconsADMM::updateUD()
{
	double threshold;
	if (mu>0)
		threshold=lambda/mu;
	else
		threshold=0;
//	std::cout << "lambda=" << lambda << std::endl;
//	std::cout << "mu=" << mu << std::endl;
//	std::cout << "Threshold=" << threshold << std::endl;
	MultidimArray<double> LCk;

	// Update gradient and Lagrange parameters
	ud=Ck(); applyLFilter(fourierLx); LCk=ud;
	ux=LCk; ux+=dx; ux.threshold("soft",threshold);
	dx+=LCk; dx-=ux;

	ud=Ck(); applyLFilter(fourierLy); LCk=ud;
	uy=LCk; uy+=dy; uy.threshold("soft",threshold);
	dy+=LCk; dy-=uy;

	ud=Ck(); applyLFilter(fourierLz); LCk=ud;
	uz=LCk; uz+=dz; uz.threshold("soft",threshold);
	dz+=LCk; dz-=uz;
}

void ProgReconsADMM::produceVolume()
{
	MultidimArray<double> kernel3D;
	kernel3D.resize(Ck());
	kernel.computeKernel3D(kernel3D);
	convolutionFFT(Ck(),kernel3D,Vk());
	// Vk()*=kernel3D.sum();
}

void AdmmKernel::initializeKernel(double _alpha, double a, double _projectionStep)
{
	projectionStep=_projectionStep;
	supp=a;
	alpha=_alpha;
	size_t length=ceil(a/projectionStep)+1;
	projectionProfile.initZeros(length);
	double ia=1.0/a;
	double K=a/bessi2(alpha)*sqrt(2.0*PI/alpha);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(projectionProfile)
	{
		double s=i*projectionStep;
		double tmp=s*ia;
		tmp=sqrt(1.0 - tmp*tmp);
		VEC_ELEM(projectionProfile,i)=K*pow(tmp,2.5)*bessi2_5(alpha*tmp);
		//	function p = KaiserBesselProjection(m, alpha, a, s)
		//		tmp = sqrt(1 - (s/a).^2);
		//		p = a ./ besseli(m, alpha) .* sqrt(2*pi/alpha) .* tmp.^(m+0.5) .* besseli(m+0.5, alpha*tmp);
		//	end
	}
}

double AdmmKernel::projectionValueAt(double u, double v) {
    double r = sqrt(u*u+v*v);
    if (r > supp) {
        return 0.;
    } else {
        r = r/projectionStep ;
        int rmin = std::floor(r);
        int rmax = rmin+1    ;
        double p = r-rmin;
        return p * VEC_ELEM(projectionProfile,rmin) + (1-p) * VEC_ELEM(projectionProfile,rmax);
    }
}

void AdmmKernel::convolveKernelWithItself(double _autocorrStep)
{
	autocorrStep=_autocorrStep;
	size_t length=ceil(2.5*supp/autocorrStep)+1;
	projectionAutocorrWithCTF.initZeros(2*length+1,2*length+1);
	projectionAutocorrWithCTF.setXmippOrigin();
	double isupp=1.0/supp;
	double K=supp/bessi2(alpha)*sqrt(2.0*PI/alpha);
	FOR_ALL_ELEMENTS_IN_ARRAY2D(projectionAutocorrWithCTF)
	{
		double s=sqrt(i*i+j*j)*autocorrStep;
		double tmp=s*isupp;
		if (tmp<1)
		{
			tmp=sqrt(1.0 - tmp*tmp);
			A2D_ELEM(projectionAutocorrWithCTF,i,j)=K*pow(tmp,2.5)*bessi2_5(alpha*tmp);
		}
	}

	transformer.FourierTransform(projectionAutocorrWithCTF,FourierProjectionAutocorr);
	K=MULTIDIM_SIZE(projectionAutocorrWithCTF)*autocorrStep*autocorrStep;
	FOR_ALL_ELEMENTS_IN_ARRAY2D(FourierProjectionAutocorr)
		A2D_ELEM(FourierProjectionAutocorr,i,j)*=A2D_ELEM(FourierProjectionAutocorr,i,j)*K;

//  Debugging code
//	  transformer.inverseFourierTransform();
//    Image<double> save;
//    save()=projectionAutocorrWithCTF;
//    save.write("PPPautocorrelationWithoutCTF.xmp");
//    std::cout << "Press any key" << std::endl;
//    char c; std::cin>> c;
}

void AdmmKernel::getKernelAutocorrelation(MultidimArray<double> &autocorrelation)
{
	transformer.fFourier=FourierProjectionAutocorr;
	transformer.inverseFourierTransform();
	autocorrelation=projectionAutocorrWithCTF;
    CenterFFT(autocorrelation, false);
}

void AdmmKernel::applyCTFToKernelAutocorrelation(CTFDescription &ctf, double Ts, MultidimArray<double> &autocorrelationWithCTF)
{
	double dig2cont=1.0/Ts;
	double wx, wy;
	int xdim=(int)XSIZE(projectionAutocorrWithCTF);
	int xdim_2=xdim/2;
	double ixdim=1.0/xdim;
	double maxFreq=2*autocorrStep;
	double maxFreq2=maxFreq*maxFreq;
	// Initialize Fourier transform to 0
	memset(&A2D_ELEM(transformer.fFourier,0,0),0,MULTIDIM_SIZE(transformer.fFourier)*sizeof(std::complex<double>));
	for (int i=((transformer.fFourier).yinit); i<=((transformer.fFourier).yinit + (int)(transformer.fFourier).ydim - 1); ++i)
	{
		FFT_IDX2DIGFREQ_FAST(i,xdim,xdim_2,ixdim,wy);
		if (fabs(wy)>maxFreq)
			continue;
		double wy2=wy*wy;
		wy*=dig2cont;
	    for (int j=((transformer.fFourier).xinit); j<=((transformer.fFourier).xinit + (int)(transformer.fFourier).xdim - 1); ++j)
	    {
			FFT_IDX2DIGFREQ_FAST(j,xdim,xdim_2,ixdim,wx);
			if (fabs(wx)>maxFreq)
				continue;
			double wx2=wx*wx;
			if (wy2+wx2>maxFreq2)
				continue;
			wx*=dig2cont;
			ctf.precomputeValues(wx,wy);
			// A2D_ELEM(transformer.fFourier,i,j)=A2D_ELEM(FourierProjectionAutocorr,i,j)*fabs(ctf.getValueAt());
			A2D_ELEM(transformer.fFourier,i,j)=fabs(ctf.getValueAt());
	    }
	}
	transformer.inverseFourierTransform();
	autocorrelationWithCTF=projectionAutocorrWithCTF;
    CenterFFT(autocorrelationWithCTF, false);

// Debugging code
//    Image<double> save;
//    save()=autocorrelationWithCTF;
//    save.write("PPPautocorrelationWithCTF.xmp");
//    std::cout << "Press any key" << std::endl;
//    char c; std::cin>> c;
}

void AdmmKernel::computeGradient(MultidimArray<double> &gradient, char direction, bool adjoint)
{
	double supp2=supp*supp;
	double K=-1/(supp2*bessi2(alpha));
	gradient.initZeros();
	for (int k=-supp; k<=supp; ++k)
	{
		int k2=k*k;
		for (int i=-supp; i<=supp; ++i)
		{
			int i2=i*i;
			for (int j=-supp; j<=supp; ++j)
			{
				int j2=j*j;
				double r2=k2+i2+j2;
				if (r2>supp2)
					continue;
				double z=alpha*sqrt(1.0-r2/supp2);
				double value=K*z*bessi1(z);
				switch (direction)
				{
				case 'x': value*=j; break;
				case 'y': value*=i; break;
				case 'z': value*=k; break;
				}
				if (adjoint)
					value*=-1;
				A3D_ELEM(gradient,k,i,j)=value;
			}
		}
	}
}

void AdmmKernel::computeKernel3D(MultidimArray<double> &kernel)
{
	double supp2=supp*supp;
	double K=1/(bessi2(alpha)*alpha*alpha);
	kernel.initZeros();
	for (int k=-supp; k<=supp; ++k)
	{
		int k2=k*k;
		for (int i=-supp; i<=supp; ++i)
		{
			int i2=i*i;
			for (int j=-supp; j<=supp; ++j)
			{
				int j2=j*j;
				double r2=k2+i2+j2;
				if (r2>supp2)
					continue;
				double z=alpha*sqrt(1.0-r2/supp2);
				double value=K*z*z*bessi2(z);
				A3D_ELEM(kernel,k,i,j)=value;
			}
		}
	}
}
