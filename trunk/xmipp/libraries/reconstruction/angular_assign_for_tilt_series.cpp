/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2008)
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

#include "angular_assign_for_tilt_series.h"

#include <data/args.h>
#include <data/filters.h>
#include <data/fft.h>
#include <data/docfile.h>
#include <data/de_solver.h>
#include <fstream>
#include <queue>
#include <iostream.h>

/* Compute affine matrix --------------------------------------------------- */
class AffineFitness
{
public:
	Matrix1D<double> minAllowed;
	Matrix1D<double> maxAllowed;
	Matrix2D<double> I1;
	Matrix2D<double> I2;
	Matrix2D<double> Mask1;
	Matrix2D<double> Mask2;
	bool showMode;
    
	AffineFitness()
	{
		showMode=false;
	}
	
	double affine_fitness_individual(double *p)
    {
       // Check limits
       if (!showMode)
           FOR_ALL_ELEMENTS_IN_MATRIX1D(minAllowed)
              if (p[i]<minAllowed(i) || p[i]>maxAllowed(i)) return 1e20;

       // Separate solution
       Matrix2D<double> A12, A21;
       A12.initIdentity(3);
       A12(0,0)=p[0]; A12(0,1)=p[1]; A12(0,2)=p[4];
       A12(1,0)=p[2]; A12(1,1)=p[3]; A12(1,2)=p[5];

       A21=A12.inv();

       // Check it is approximately a rotation
       if (!showMode)
       {
           double A12det=A12.det();
           if (ABS(A12det)<0.9 || ABS(A12det)>1.1) return 1e20;
       }

       // Produce the transformed images
       Matrix2D<double> transformedI1, transformedI2;
       applyGeometry(transformedI1,A12,I1,IS_NOT_INV,DONT_WRAP);
       applyGeometry(transformedI2,A21,I2,IS_NOT_INV,DONT_WRAP);
       
       // Produce masks for the comparison
       Matrix2D<int> maskInTheSpaceOf1, maskInTheSpaceOf2;
       Matrix2D<double> maskAux;
       applyGeometry(maskAux,A12,Mask1,IS_NOT_INV,DONT_WRAP);
       maskInTheSpaceOf2.initZeros(YSIZE(maskAux),XSIZE(maskAux));
       maskInTheSpaceOf2.setXmippOrigin();
       FOR_ALL_ELEMENTS_IN_MATRIX2D(maskAux)
            maskInTheSpaceOf2(i,j)=ROUND(maskAux(i,j));
       maskAux.initZeros();
       applyGeometry(maskAux,A21,Mask2,IS_NOT_INV,DONT_WRAP);
       maskInTheSpaceOf1.initZeros(YSIZE(maskAux),XSIZE(maskAux));
       maskInTheSpaceOf1.setXmippOrigin();
       FOR_ALL_ELEMENTS_IN_MATRIX2D(maskAux)
            maskInTheSpaceOf1(i,j)=ROUND(maskAux(i,j));
       
       // Compare the two images
       double dist=0.5*(1-correlation_index(transformedI1,I2,&maskInTheSpaceOf2)+
                        1-correlation_index(transformedI2,I1,&maskInTheSpaceOf1));

       if (showMode) {
          ImageXmipp save;
          save()=I1; save.write("PPPimg1.xmp");
          save()=I2; save.write("PPPimg2.xmp");
          save()=transformedI1; save.write("PPPTransformedImg1.xmp");
          save()=transformedI2; save.write("PPPTransformedImg2.xmp");
          typeCast(maskInTheSpaceOf1,save()); save.write("PPPmaskInTheSpaceOf1.xmp");
          typeCast(maskInTheSpaceOf2,save()); save.write("PPPmaskInTheSpaceOf2.xmp");
          save()=I1-transformedI2; save.write("PPPDiffImg1.xmp");
          save()=I2-transformedI1; save.write("PPPDiffImg2.xmp");
          std::cout << "A12=\n" << A12 << "A21=\n" << A21 << std::endl;
	      std::cout << "dist=" << dist << std::endl;
          std::cout << "Press any key\n";
	      char c;
	      std::cin >> c;
       }

       return dist;
    }

    static double Powell_affine_fitness_individual(double *p, void *prm)
    {
       return ((AffineFitness*)prm)->affine_fitness_individual(p+1);
    }
};

class AffineSolver: public DESolver {
   AffineFitness *fitness;
public:
   AffineSolver(AffineFitness *newFitness, int dim,int pop) : DESolver(dim,pop), count(0), fitness(newFitness) {;}
   double EnergyFunction(double trial[], bool &bAtSolution) {
      double result=fitness->affine_fitness_individual(trial);
      #ifdef DEBUG
		if (count++ % nPop == 0)
		{
			std::cout << "EvaluationsAffine= " << count/nPop
						<< " energyAffine= " << Energy()
				<< std::endl;
		}
	  #endif
	  bAtSolution=false;
      return(result);
   }
private:
   int count;
};
#undef DEBUG

double computeAffineTransformation(const Matrix2D<double> &I1,
    const Matrix2D<double> &I2, int maxShift, int maxIterDE,
    const FileName &fn_affine, 
    Matrix2D<double> &A12, Matrix2D<double> &A21, bool show,
    double thresholdAffine, bool localAffine, bool isMirror)
{
    AffineFitness fitness;

    // Set images
    fitness.I1=I1;
    fitness.I1.setXmippOrigin();
	fitness.I2=I2;
    fitness.I2.setXmippOrigin();
   
    fitness.Mask1=fitness.I1;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I1)
    	if (I1(i,j)!=0) fitness.Mask1(i,j)=1;

    fitness.Mask2=fitness.I2;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I2)
    	if (I2(i,j)!=0) fitness.Mask2(i,j)=1;

    // Set limits for the affine matrices
    // Order: 1->2: 4 affine params+2 translations
    // Order: 2->1: 4 affine params+2 translations
    fitness.minAllowed.resize(6);
    fitness.maxAllowed.resize(6);

    // Scale factors
    fitness.minAllowed(0)=fitness.minAllowed(3)=0.5;
    fitness.maxAllowed(0)=fitness.maxAllowed(3)=1.5;
    if (isMirror)
    {
        fitness.minAllowed(3)=-1.5;
        fitness.maxAllowed(3)=-0.5;
    }

    // Rotation factors
    fitness.minAllowed(1)=fitness.minAllowed(2)=-0.5;
    fitness.maxAllowed(1)=fitness.maxAllowed(2)= 0.5;

    // Shifts
    fitness.minAllowed(4)=fitness.minAllowed(5)=-maxShift;
    fitness.maxAllowed(4)=fitness.maxAllowed(5)= maxShift;

    std::ifstream fh_in;
    fh_in.open(fn_affine.c_str());
	
	// Return result
	double cost;
	
    if (fh_in) {
    	A12.resize(3,3);
	    A21.resize(3,3);
     	fh_in >> A12 >> A21 >> cost;
	    fh_in.close();
    } else {
        // Optimize with differential evolution
        Matrix1D<double> A(6);
        if (!localAffine)
        {
            double bestEnergy=2, energy;
            int n=0;
            do
            {
                AffineSolver solver(&fitness,6,6*10);
                solver.Setup(MULTIDIM_ARRAY(fitness.minAllowed),
                             MULTIDIM_ARRAY(fitness.maxAllowed),
		             stBest2Bin, 0.5, 0.8);
                solver.Solve(maxIterDE);
                energy=solver.Energy();
                if (n==0 || bestEnergy>energy)
                {
                    FOR_ALL_ELEMENTS_IN_MATRIX1D(A)
                       A(i)=solver.Solution()[i];
                    bestEnergy=energy;
                }
                n++;
            } while ((n<3 || (n>=3 && n<10 && bestEnergy>1-thresholdAffine))
                     && bestEnergy>=0.07);
        }
        else
        {
            A(0)=A(3)=1;
            if (isMirror)
                A(3)*=-1;
            double tx, ty;
            if (!isMirror) best_shift(I1,I2,tx,ty);
            else
            {
                Matrix2D<double> auxI2=I2;
                auxI2.selfReverseY();
                best_shift(I1,auxI2,tx,ty);
                ty=-ty;
            }
            A(4)=-tx;
            A(5)=-ty;
        }

        // Optimize with Powell
        Matrix1D<double> steps(A);
        steps.initConstant(1);
        int iter;
        powellOptimizer(A, 1, XSIZE(A),
            AffineFitness::Powell_affine_fitness_individual, &fitness, 0.005,
            cost, iter, steps, false);

        // Separate solution
        A12.initIdentity(3);
        A12(0,0)=A(0); A12(0,1)=A(1); A12(0,2)=A(4);
        A12(1,0)=A(2); A12(1,1)=A(3); A12(1,2)=A(5);

        A21=A12.inv();

        std::ofstream fh_out;
        fh_out.open(fn_affine.c_str());
        fh_out << A12 << std::endl << A21 << std::endl << cost;
        fh_out.close();
    }   
    if (show)
    {
    	fitness.showMode=true;
        Matrix1D<double> p(6);
        p(0)=A12(0,0); p(1)=A12(0,1); p(4)=A12(0,2);
        p(2)=A12(1,0); p(3)=A12(1,1); p(5)=A12(1,2);
	    fitness.affine_fitness_individual(MULTIDIM_ARRAY(p));
    	fitness.showMode=false;
    }
	
	return cost;
}

/* Parameters -------------------------------------------------------------- */
void Prog_tomograph_alignment::read(int argc, char **argv) {
   fnSel=getParameter(argc,argv,"-i");
   fnSelOrig=getParameter(argc,argv,"-iorig","");
   fnRoot=getParameter(argc,argv,"-oroot","");
   if (fnRoot=="")
      fnRoot=fnSel.without_extension();
   localAffine=checkParameter(argc,argv,"-localAffine");
   seqLength=textToInteger(getParameter(argc,argv,"-seqLength","5"));
   maxStep=textToInteger(getParameter(argc,argv,"-maxStep","4"));
   gridSamples=textToInteger(getParameter(argc,argv,"-gridSamples","40"));
   psiMax=textToFloat(getParameter(argc,argv,"-psiMax","-1"));
   deltaRot=textToFloat(getParameter(argc,argv,"-deltaRot","5"));
   localSize=textToFloat(getParameter(argc,argv,"-localSize","0.04"));
   optimizeTiltAngle=checkParameter(argc,argv,"-optimizeTiltAngle");
   isCapillar=checkParameter(argc,argv,"-isCapillar");
   corrThreshold=textToFloat(getParameter(argc,argv,"-threshold","0.9"));
   maxShiftPercentage=textToFloat(getParameter(argc,argv,"-maxShiftPercentage","0.2"));
   maxIterDE=textToInteger(getParameter(argc,argv,"-maxIterDE","30"));
   showAffine=checkParameter(argc,argv,"-showAffine");
   thresholdAffine=textToFloat(getParameter(argc,argv,"-thresholdAffine","0.85"));
   numThreads = textToInteger(getParameter(argc, argv, "-thr", "1"));
   if (numThreads<1) numThreads = 1;
}

void Prog_tomograph_alignment::show() {
   std::cout << "Input images:       " << fnSel              << std::endl
             << "Original images:    " << fnSelOrig          << std::endl
             << "Output rootname:    " << fnRoot             << std::endl
             << "Local affine:       " << localAffine        << std::endl
             << "SeqLength:          " << seqLength          << std::endl
             << "MaxStep:            " << maxStep            << std::endl
             << "Grid samples:       " << gridSamples        << std::endl
             << "Maximum psi:        " << psiMax             << std::endl
             << "Delta rot:          " << deltaRot           << std::endl
             << "Local size:         " << localSize          << std::endl
             << "Optimize tilt angle:" << optimizeTiltAngle  << std::endl
             << "isCapillar:         " << isCapillar         << std::endl
             << "Threshold:          " << corrThreshold      << std::endl
             << "MaxShift Percentage:" << maxShiftPercentage << std::endl
             << "MaxIterDE:          " << maxIterDE          << std::endl
             << "Show Affine:        " << showAffine         << std::endl
             << "Threshold Affine:   " << thresholdAffine    << std::endl
             << "Threads to use:     " << numThreads         << std::endl
   ;
}

void Prog_tomograph_alignment::usage() const {
   std::cerr << "tomograph_alignment\n"
             << "   -i <selfile>                   : Input images\n"
             << "  [-iorig <selfile>]              : Selfile with images at original scale\n"
             << "  [-oroot <fn_out>]               : Output alignment\n"
             << "  [-localAffine]                  : Look for affine transformations close to I\n"
             << "  [-seqLength <n=5>]              : Sequence length\n"
             << "  [-maxStep <step=4>]             : Maximum step for chain refinement\n"
             << "  [-gridSamples <n=40>]           : Total number of samples=n*n\n"
             << "  [-psiMax <psi=-1>]              : Maximum psi in absolute value (degrees)\n"
             << "                                    -1: do not optimize for psi\n"
             << "  [-deltaRot <rot=5>]             : In degrees. For the first optimization stage\n"
             << "  [-localSize <size=0.04>]        : In percentage\n"
             << "  [-optimizeTiltAngle]            : Optimize tilt angle\n"
             << "  [-isCapillar]                   : Set this flag if the tilt series is of a capillar\n"
             << "  [-threshold <th=0.9>]           : threshold\n"
             << "  [-maxShiftPercentage <p=0.2>]   : Maximum shift as percentage of image size\n"
             << "  [-maxIterDE <n=30>]             : Maximum number of iteration in Differential Evolution\n"
             << "  [-showAffine]                   : Show affine transformations as PPP*\n"
             << "  [-thresholdAffine <th=0.85>]    : Threshold affine\n"
             << "  [-thr <num=1>]                  : Parallel processing using \"num\" threads\n"
   ;
}

/* Produce side info ------------------------------------------------------- */
void Prog_tomograph_alignment::produceSideInfo() {
   bestPreviousAlignment=new Alignment(this);
   // Read input data
   SF.read(fnSel);
   Nimg=SF.ImgNo();
   if (Nimg!=0)
   {
       // Clear the list of images if not empty
       if (img.size()!=0) {
          for (int i=0; i<img.size(); i++) delete img[i];
          img.clear();
       }

       std::cerr << "Reading input data\n";
       init_progress_bar(Nimg);
       int n=0;
       SF.go_first_ACTIVE();
       while (!SF.eof()) {
          FileName fn=SF.NextImg();
          if (fn=="") break;
          ImageXmipp imgaux(fn);
          Matrix2D<double>* img_i=new Matrix2D<double>;
          *img_i=imgaux();
          img_i->setXmippOrigin();
          img.push_back(img_i);
          tiltList.push_back(imgaux.tilt());
          name_list.push_back(imgaux.name());

          progress_bar(n++);
       }
       progress_bar(Nimg);
   }
   
    // Read images at original scale
    if (fnSelOrig!="")
    {
        SForig.read(fnSelOrig);
        if (SForig.ImgNo()!=SF.ImgNo())
            REPORT_ERROR(1,"The number of images in both selfiles (-i and -iorig) is different");
    }

    // Fill the affine transformations with empty matrices
    std::vector< Matrix2D<double> > emptyRow;
    for (int i=0; i<Nimg; i++)
    {
        Matrix2D<double> A;
        emptyRow.push_back(A);
    }
    for (int i=0; i<Nimg; i++)
    	affineTransformations.push_back(emptyRow);

    pthread_t * th_ids = new pthread_t[numThreads];
    ThreadParams * th_args= new ThreadParams[numThreads];
	
    for( int nt = 0 ; nt < numThreads ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        pthread_create( (th_ids+nt) , NULL, threadComputeTransform, (void *)(th_args+nt) );
    }
		
    // Waiting for threads to finish
    for( int nt = 0 ; nt < numThreads ; nt ++ )
        pthread_join(*(th_ids+nt), NULL);

    // Threads structures are not needed any more
    delete( th_ids );
    delete( th_args );

    // Do not show refinement
    showRefinement=false;
}

void * threadComputeTransform( void * args )
{
    ThreadParams * master = (ThreadParams *) args;

    Prog_tomograph_alignment * parent = master->parent;
    int thread_id = master->myThreadID;
    int localnumThreads = parent->numThreads;
    bool isCapillar = parent->isCapillar;
    int Nimg = parent->Nimg;
    double maxShiftPercentage = parent->maxShiftPercentage;
    int maxIterDE = parent->maxIterDE;
    bool showAffine = parent->showAffine;
    double thresholdAffine = parent->thresholdAffine;
    std::vector < Matrix2D<double> *> img = parent->img;
    std::vector< std::vector< Matrix2D<double> > > affineTransformations = parent->affineTransformations;
    double localAffine = parent->localAffine;

    int maxShift=FLOOR(XSIZE(*img[0])*maxShiftPercentage);
    int initjj=1;
    if (isCapillar) initjj=0;
	
    initjj += thread_id;

    double cost;
    for (int jj=initjj; jj<Nimg; jj+= localnumThreads)
    {
        int jj_1;
        if (isCapillar) jj_1=intWRAP(jj-1,0,Nimg-1);
        else            jj_1=jj-1;
        Matrix2D<double>& img_i=*img[jj_1];
        Matrix2D<double>& img_j=*img[jj];
        bool isMirror=(jj==0) && (jj_1==Nimg-1);
        Matrix2D<double> Aij, Aji;
        cost = computeAffineTransformation(img_i, img_j, maxShift,
            maxIterDE,
            (std::string)"affine_"+integerToString(jj_1,3)+
            "_"+integerToString(jj,3)+".txt", Aij, Aji,
            showAffine, thresholdAffine, localAffine,
            isMirror);
        
		pthread_mutex_lock( &printingMutex );
		cout << "Cost for [" << jj_1 << "] - [" << jj << "] = " << cost << endl;
		pthread_mutex_unlock( &printingMutex );
		affineTransformations[jj_1][jj]=Aij;
		affineTransformations[jj][jj_1]=Aji;
    }
	
    return NULL;
}
 
/* Generate landmark set --------------------------------------------------- */
#define DEBUG
void Prog_tomograph_alignment::generateLandmarkSet() {
    if (!exists(fnRoot+"_landmarks.txt"))
    {
        int deltaShift=FLOOR(XSIZE(*img[0])/gridSamples);
        std::vector<LandmarkChain> chainList;
        Matrix1D<double> rii(3), rjj(3);
        ZZ(rii)=1;
        ZZ(rjj)=1;
        init_progress_bar(gridSamples);
        int totalPoints=0;
        int includedPoints=0;
        for (int nx=0; nx<gridSamples; ++nx)
        {
            XX(rii)=STARTINGX(*img[0])+ROUND(deltaShift*(0.5+nx));
            for (int ny=0; ny<gridSamples; ++ny)
            {
                YY(rii)=STARTINGY(*img[0])+ROUND(deltaShift*(0.5+ny));
                for (int ii=0; ii<=Nimg-1; ++ii)
                {
                    LandmarkChain chain;
                    chain.clear();
                    Landmark l;
                    l.x=XX(rii);
                    l.y=YY(rii);
                    l.imgIdx=ii;
                    chain.push_back(l);
                    
                    // Follow this landmark backwards
                    bool acceptLandmark=true;
                    int jjleft=ii, jj;
                    if (isCapillar) jj=intWRAP(ii-1,0,Nimg-1);
                    else            jj=ii-1;
                    Matrix2D<double> Aij, Aji;
                    Matrix1D<double> rcurrent=rii;
                    while (jj>=0 && acceptLandmark && jj!=ii)
                    {
                        // Compute the affine transformation between ii and jj
                        int jj_1;
                        if (isCapillar) jj_1=intWRAP(jj+1,0,Nimg-1);
                        else            jj_1=jj+1;
                        Aij=affineTransformations[jj][jj_1];
	                Aji=affineTransformations[jj_1][jj];
                        rjj=Aji*rcurrent;
                        acceptLandmark=refineLandmark(jj_1,jj,rcurrent,rjj);
                        if (acceptLandmark)
                        {
                            l.x=XX(rjj);
                            l.y=YY(rjj);
                            l.imgIdx=jj;
                            chain.push_back(l);
                            jjleft=jj;
                            if (isCapillar) jj=intWRAP(jj-1,0,Nimg-1);
                            else            jj=jj-1;
                            rcurrent=rjj;
                        }
                    }

                    // Follow this landmark forward
                    acceptLandmark=true;
                    if (isCapillar) jj=intWRAP(ii+1,0,Nimg-1);
                    else jj=ii+1;
                    rcurrent=rii;
                    int jjright=ii;
                    while (jj<Nimg && acceptLandmark && jj!=ii)
                    {
                        // Compute the affine transformation between ii and jj
                        int jj_1;
                        if (isCapillar) jj_1=intWRAP(jj-1,0,Nimg-1);
                        else            jj_1=jj-1;
                        Aij=affineTransformations[jj_1][jj];
	                Aji=affineTransformations[jj][jj_1];
                        rjj=Aij*rcurrent;
                        acceptLandmark=refineLandmark(jj_1,jj,rcurrent,rjj);
                        if (acceptLandmark)
                        {
                            l.x=XX(rjj);
                            l.y=YY(rjj);
                            l.imgIdx=jj;
                            chain.push_back(l);
                            jjright=jj;
                            if (isCapillar) jj=intWRAP(jj+1,0,Nimg-1);
                            else            jj=jj+1;
                            rcurrent=rjj;
                        }
                    }
                    
                    std::cout << "img=" << ii << " chain length="
                              << chain.size() << " [" << jjleft 
                              << " - " << jjright << "]";
                    if (chain.size()>seqLength)
                    {
                        bool accepted=refineChain(chain);
                        if (accepted)
                        {
                            std::cout << " Accepted with length= "
                                      << chain.size();
                            chainList.push_back(chain);
                            includedPoints+=chain.size();
                        }
                    }
                    std::cout << std::endl;
                }
                #ifdef DEBUG
                    std::cout << "Point nx=" << nx << " ny=" << ny
                              << " Number of points="
                              << includedPoints
                              << " Number of chains=" << chainList.size()
                              << " ( " << ((double) includedPoints)/chainList.size() << " )\n";
                #endif
            }
            progress_bar(nx);
        }
        progress_bar(gridSamples);

        // Generate the landmark "matrix"
        allLandmarksX.resize(chainList.size(),Nimg);
        allLandmarksY.resize(chainList.size(),Nimg);
        allLandmarksX.initConstant(XSIZE(*img[0]));
        allLandmarksY.initConstant(YSIZE(*img[0]));
        for (int i=0; i<chainList.size(); i++) {
            for (int j=0; j<chainList[i].size(); j++)
            {
                int idx=chainList[i][j].imgIdx;
                allLandmarksX(i,idx)=chainList[i][j].x;
                allLandmarksY(i,idx)=chainList[i][j].y;
            }
        }
        
        // Write landmarks
        writeLandmarkSet(fnRoot+"_landmarks.txt");
    }
    else
    {
        readLandmarkSet(fnRoot+"_landmarks.txt");
    }
}
#undef DEBUG

/* Produce information from landmarks -------------------------------------- */
void Prog_tomograph_alignment::produceInformationFromLandmarks()
{
    // Produce V sets
    std::vector<int> emptyVector;

    Vseti.clear();
    for (int i=0; i<XSIZE(allLandmarksX); i++)
       Vseti.push_back(emptyVector);

    Vsetj.clear();
    for (int j=0; j<YSIZE(allLandmarksX); j++)
       Vsetj.push_back(emptyVector);

    for (int j=0; j<YSIZE(allLandmarksX); j++)
        for (int i=0; i<XSIZE(allLandmarksX); i++)
            if (allLandmarksX(j,i)!=XSIZE(*img[0]))
            {
                Vseti[i].push_back(j);
                Vsetj[j].push_back(i);
            }
    
    // Count the number of landmarks per image and average projection
    // of all landmarks in a given image
    ni.initZeros(Nimg);
    barpi.clear();
    for (int i=0; i<Nimg; i++)
    {
        ni(i)=static_cast<int>(Vseti[i].size());
        Matrix1D<double> pi(2);
        for (int jj=0; jj<ni(i); jj++)
        {
            int j=Vseti[i][jj];
            XX(pi)+=allLandmarksX(j,i);
            YY(pi)+=allLandmarksY(j,i);
        }
        pi/=ni(i);
        barpi.push_back(pi);
    }
}

/* Refine landmark --------------------------------------------------------- */
bool Prog_tomograph_alignment::refineLandmark(int ii, int jj,
    const Matrix1D<double> &rii, Matrix1D<double> &rjj) const
{
    int halfSize=ROUND(localSize*XSIZE(*img[ii]))/2;
    if (XX(rii)<STARTINGX(*img[ii])+2.6*halfSize || 
        XX(rii)>FINISHINGX(*img[ii])-2.6*halfSize ||
        YY(rii)<STARTINGY(*img[ii])+2.6*halfSize || 
        YY(rii)>FINISHINGY(*img[ii])-2.6*halfSize ||
        XX(rjj)<STARTINGX(*img[jj])+2.6*halfSize || 
        XX(rjj)>FINISHINGX(*img[jj])-2.6*halfSize ||
        YY(rjj)<STARTINGY(*img[jj])+2.6*halfSize || 
        YY(rjj)>FINISHINGY(*img[jj])-2.6*halfSize)
        return 0;

    // Select piece in image ii, compute its statistics and normalize
    Matrix2D<double> pieceii(2*halfSize+1,2*halfSize+1),
        piecejj(2*halfSize+1,2*halfSize+1);
    pieceii.setXmippOrigin();
    piecejj.setXmippOrigin();
    double mean_ii=0, stddev_ii=0;
    const Matrix2D<double> &Iii=(*img[ii]);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(pieceii)
    {
        MAT_ELEM(pieceii,i,j)=MAT_ELEM(Iii,
           (int)(YY(rii)+i),(int)(XX(rii)+j));
        mean_ii+=MAT_ELEM(pieceii,i,j);
        stddev_ii+=MAT_ELEM(pieceii,i,j)*MAT_ELEM(pieceii,i,j);
    }
    mean_ii/=MULTIDIM_SIZE(pieceii);
    stddev_ii = stddev_ii / MULTIDIM_SIZE(pieceii) - mean_ii * mean_ii;
    stddev_ii *= MULTIDIM_SIZE(pieceii) / (MULTIDIM_SIZE(pieceii) - 1);
    stddev_ii = sqrt(static_cast<double>((ABS(stddev_ii))));
    if (stddev_ii>XMIPP_EQUAL_ACCURACY)
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pieceii)
           DIRECT_MULTIDIM_ELEM(pieceii,n)=
              (DIRECT_MULTIDIM_ELEM(pieceii,n)-mean_ii)/stddev_ii;

    // Try all possible shifts
    Matrix2D<double> corr((int)(1.5*(2*halfSize+1)),(int)(1.5*(2*halfSize+1)));
    corr.setXmippOrigin();
    corr.initConstant(-1.1);
    bool accept=false;
    if (stddev_ii>XMIPP_EQUAL_ACCURACY)
    {
        double maxval=0;
        int imax=0, jmax=0;
        std::queue< Matrix1D<double> > Q;
        Q.push(vectorR2(0,0));
        while (!Q.empty())
        {
            // Get the first position to evaluate
            int shifty=(int)YY(Q.front());
            int shiftx=(int)XX(Q.front());
            Q.pop();

            // Select piece in image jj, compute its statistics and normalize
            double mean_jj=0, stddev_jj=0;
            const Matrix2D<double> &Ijj=(*img[jj]);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(piecejj)
            {
                MAT_ELEM(piecejj,i,j)=MAT_ELEM(Ijj,
                    (int)(YY(rjj)+shifty+i),
                    (int)(XX(rjj)+shiftx+j));
                mean_jj+=MAT_ELEM(piecejj,i,j);
                stddev_jj+=MAT_ELEM(piecejj,i,j)*MAT_ELEM(piecejj,i,j);
            }
            mean_jj/=MULTIDIM_SIZE(piecejj);
            stddev_jj = stddev_jj / MULTIDIM_SIZE(piecejj) - mean_jj * mean_jj;
            stddev_jj *= MULTIDIM_SIZE(piecejj) / (MULTIDIM_SIZE(piecejj) - 1);
            stddev_jj = sqrt(static_cast<double>((ABS(stddev_jj))));

            // Compute the correlation
            corr(shifty,shiftx)=0;
            double &corrRef=corr(shifty,shiftx);
            if (stddev_jj>XMIPP_EQUAL_ACCURACY)
            {
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(piecejj)
                   corrRef+=
                      DIRECT_MULTIDIM_ELEM(pieceii,n)*
                      (DIRECT_MULTIDIM_ELEM(piecejj,n)-mean_jj)/stddev_jj;
                corrRef/=MULTIDIM_SIZE(piecejj);
            }
            
            if (corrRef>maxval)
            {
                maxval=corrRef;
                imax=shifty;
                jmax=shiftx;
                for (int stepy=-1; stepy<=1; stepy++)
                    for (int stepx=-1; stepx<=1; stepx++)
                    {
                        int newshifty=shifty+stepy;
                        int newshiftx=shiftx+stepx;
                        if (newshifty>=STARTINGY(corr) &&
                            newshifty<=FINISHINGY(corr) &&
                            newshiftx>=STARTINGX(corr) &&
                            newshiftx<=FINISHINGX(corr))
                            if (corr(newshifty,newshiftx)<-1)
                                Q.push(vectorR2(newshifty,newshiftx));
                    }
            }
        }

        if (maxval>corrThreshold)
        {
            YY(rjj)+=imax;
            XX(rjj)+=jmax;
            accept=true;

            if (showRefinement)
            {
                ImageXmipp save;
                save()=pieceii; save.write("PPPpieceii.xmp");
                FOR_ALL_ELEMENTS_IN_MATRIX2D(piecejj)
                    piecejj(i,j)=(*img[jj])((int)(YY(rjj)+i),(int)(XX(rjj)+j));
                save()=piecejj; save.write("PPPpiecejj.xmp");
                save()=corr; save.write("PPPcorr.xmp");
                std::cout << "ii=" << ii << " jj=" << jj << std::endl;
                std::cout << "rii=" << rii.transpose() << std::endl;
                std::cout << "rjj=" << rjj.transpose() << std::endl;
                std::cout << "imax=" << imax << " jmax=" << jmax << std::endl;
                std::cout << "maxval=" << maxval << std::endl;
                std::cout << "Press any key\n";
                char c;
                std::cin >> c;
            }
        }
    }
    return (accept);
}

/* Refine chain ------------------------------------------------------------ */
//#define DEBUG
bool Prog_tomograph_alignment::refineChain(LandmarkChain &chain)
{
    for (int K=0; K<2 && chain.size()>seqLength; K++)
    {
        sort(chain.begin(), chain.end());
        int chainLength=chain.size();
        Matrix1D<double> rii(2), rjj(2), newrjj(2);

        // Refine every step
        for (int step=2; step<=maxStep; step++)
        {
            // Refine forwards every step
            for (int i=0; i<chainLength-step; i++)
            {
                int ii=chain[i].imgIdx;
                VECTOR_R2(rii,chain[i].x,chain[i].y);
                int jj=chain[i+step].imgIdx;
                VECTOR_R2(rjj,chain[i+step].x,chain[i+step].y);
                newrjj=rjj;
                bool accepted=refineLandmark(ii,jj,rii,newrjj);
                if ((newrjj-rjj).module()<4 && accepted)
                {
                    chain[i+step].x=XX(newrjj);
                    chain[i+step].y=YY(newrjj);
                }
            }

            #ifdef DEBUG
                showRefinement=(step==maxStep);
            #endif
            // Refine backwards all images
            for (int i=chainLength-1; i>=1; i--)
            {
                int ii=chain[i].imgIdx;
                VECTOR_R2(rii,chain[i].x,chain[i].y);
                int jj=chain[i-1].imgIdx;
                VECTOR_R2(rjj,chain[i-1].x,chain[i-1].y);
                newrjj=rjj;
                bool accepted=refineLandmark(ii,jj,rii,newrjj);
                if ((newrjj-rjj).module()<4 && accepted)
                {
                    chain[i-1].x=XX(newrjj);
                    chain[i-1].y=YY(newrjj);
                }
            }
            #ifdef DEBUG
                showRefinement=false;
            #endif
        }

        // Check that all images in the stack correlate well with the
        // average of all pieces
        int halfSize=ROUND(localSize*XSIZE(*img[0]))/2;
        Matrix2D<double> avgPiece(2*halfSize+1,2*halfSize+1);
        avgPiece.setXmippOrigin();
        // compute average piece
        for (int n=0; n<chainLength; n++)
        {
            int ii=chain[n].imgIdx;
            VECTOR_R2(rii,chain[n].x,chain[n].y);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(avgPiece)
                avgPiece(i,j)+=(*img[ii])((int)(YY(rii)+i),(int)(XX(rii)+j));
        }
        avgPiece/=chainLength;

        // Compute correlation with the average
        Matrix2D<double> pieceii(2*halfSize+1,2*halfSize+1);
        pieceii.setXmippOrigin();
        double avgCorr=0;
        LandmarkChain refinedChain;
        refinedChain.clear();
        for (int n=0; n<chainLength; n++)
        {
            int ii=chain[n].imgIdx;
            VECTOR_R2(rii,chain[n].x,chain[n].y);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(pieceii)
                pieceii(i,j)=(*img[ii])((int)(YY(rii)+i),(int)(XX(rii)+j));
            double corr=correlation_index(pieceii,avgPiece);
            if (corr>corrThreshold)
            {
               avgCorr+=corr;
               refinedChain.push_back(chain[n]);
            }
            #ifdef DEBUG
                ImageXmipp save;
                save()=pieceii; save.write("PPPpieceiiRefine.xmp");
                save()=avgPiece; save.write("PPPpieceAvg.xmp");
                std::cout << "Corr=" << corr;
                if (corr<=corrThreshold) std::cout << "  Removed!!";
                std::cout << "\nAvgCorr=" << avgCorr/refinedChain.size() << std::endl
                          << "RefinedLength=" << refinedChain.size() << std::endl
                          << "K=" << K << std::endl
                          << "Press any key\n";
                char c;
                std::cin >> c;
            #endif
        }
        avgCorr/=refinedChain.size();
        chain.clear();
        chain=refinedChain;
    }
    
    double tilt0=tiltList[chain[0].imgIdx];
    double tiltF=tiltList[chain[chain.size()-1].imgIdx];
    double lengthThreshold=FLOOR(seqLength*cos(DEG2RAD(0.5*(tilt0+tiltF))));
    return chain.size()>lengthThreshold;
}
#undef DEBUG

/* Read/Write landmark set ------------------------------------------------- */
void Prog_tomograph_alignment::writeLandmarkSet(const FileName &fnLandmark) const
{
    std::ofstream fhOut;
    fhOut.open(fnLandmark.c_str());
    if (!fhOut)
       REPORT_ERROR(1,(std::string)"Cannot open "+fnLandmark+" for output");
    fhOut << "Point     x       y       slice   color "
          << YSIZE(allLandmarksX) << " " << XSIZE(allLandmarksX) << std::endl;
    for (int i=0; i<XSIZE(allLandmarksX); i++)
    {
        int counter=0;
        for (int j=0; j<YSIZE(allLandmarksX); j++)
            if (allLandmarksX(j,i)!=XSIZE(*img[0]))
            {
                fhOut << counter << " \t"
                      << ROUND(allLandmarksX(j,i)-STARTINGX(*img[0])) << " \t"
                      << ROUND(allLandmarksY(j,i)-STARTINGY(*img[0])) << " \t"
                      << i+1 << " \t" << j << std::endl;
                counter++;
            }
    }
    fhOut.close();
}

void Prog_tomograph_alignment::readLandmarkSet(const FileName &fnLandmark)
{
    std::ifstream fhIn;
    fhIn.open(fnLandmark.c_str());
    if (!fhIn)
       REPORT_ERROR(1,(std::string)"Cannot open "+fnLandmark+" for input");
    std::string dummyStr;
    int Nlandmark;
    fhIn >> dummyStr >> dummyStr >> dummyStr >> dummyStr >> dummyStr
         >> Nlandmark >> Nimg;
    if (Nlandmark<=0)
        REPORT_ERROR(1,(std::string)"No landmarks are found in "+fnLandmark);
    allLandmarksX.resize(Nlandmark,Nimg);
    allLandmarksY.resize(Nlandmark,Nimg);
    allLandmarksX.initConstant(XSIZE(*img[0]));
    allLandmarksY.initConstant(XSIZE(*img[0]));
    while (!fhIn.eof())
    {
        int dummyInt, x, y, i, j;
        fhIn >> dummyInt >> x >> y >> i >> j;
        i=i-1;
        allLandmarksX(j,i)=x+STARTINGX(*img[0]);
        allLandmarksY(j,i)=y+STARTINGY(*img[0]);
    }
    fhIn.close();
    std::cout << "The file " << fnLandmark << " has been read for the landmarks\n"
              << Nlandmark << " landmarks are read\n";
}

/* Align images ------------------------------------------------------------ */
void Prog_tomograph_alignment::alignImages(const Alignment &alignment)
{
   DocFile DF;
   if (fnSelOrig!="")
      SForig.go_first_ACTIVE();
   DF.append_comment("in-plane rotation    Xshift      Yshift");
   DF.append_comment("First shift by -(Xshift,Yshift)");
   DF.append_comment("Then, rotate by in-plane rotation");
   for (int i=0;i<Nimg; i++) {
        // Align the normal image
	ImageXmipp I;
	I.read(name_list[i]);
        Matrix2D<double> mask;
        mask.initZeros(I());
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
            if (I(i,j)!=0) mask(i,j)=1;
	I().selfTranslate(-alignment.di[i],DONT_WRAP);
	I().selfRotate(90-alignment.rot+alignment.psi(i),DONT_WRAP);
	mask.selfTranslate(-alignment.di[i],DONT_WRAP);
	mask.selfRotate(90-alignment.rot+alignment.psi(i),DONT_WRAP);
	mask.binarize(0.5);
	Matrix2D<int> iMask;
	typeCast(mask,iMask);
        double minval, maxval, avg, stddev;
        computeStats_within_binary_mask(iMask,I(),minval, maxval, avg, stddev);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(iMask)
            if (iMask(i,j)==0) I(i,j)=0;
            else I(i,j)=(I(i,j)-avg)/stddev;
        double rot=0;
        double tilt=tiltList[i];
        double psi=0;
	I.set_eulerAngles(rot, tilt, psi);
	FileName fn_corrected=fnRoot+"_corrected_"+integerToString(i,3)+".xmp";
	I.write(fn_corrected);
	
        // Align the original image
        if (fnSelOrig!="")
        {
	    ImageXmipp Iorig;
	    Iorig.read(SForig.NextImg());
            mask.initZeros(Iorig());
            FOR_ALL_ELEMENTS_IN_MATRIX2D(Iorig())
                if (Iorig(i,j)!=0) mask(i,j)=1;
	    Iorig().selfTranslate(-alignment.di[i]*XSIZE(Iorig())/XSIZE(I()),
                DONT_WRAP);
	    Iorig().selfRotate(90-alignment.rot+alignment.psi(i),DONT_WRAP);
	    mask.selfTranslate(-alignment.di[i]*XSIZE(Iorig())/XSIZE(I()),
                DONT_WRAP);
	    mask.selfRotate(90-alignment.rot+alignment.psi(i),DONT_WRAP);
	    mask.binarize(0.5);
	    typeCast(mask,iMask);
            computeStats_within_binary_mask(iMask,Iorig(),minval, maxval,
                avg, stddev);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(iMask)
                if (iMask(i,j)==0) Iorig(i,j)=0;
                else Iorig(i,j)=(Iorig(i,j)-avg)/stddev;
	    Iorig.set_eulerAngles(rot, tilt, psi);
	    fn_corrected=fnRoot+"_corrected_originalsize_"+integerToString(i,3)+".xmp";
	    I.write(fn_corrected);
        }

        // Prepare data for the docfile
        Matrix1D<double> params(3);
        params(0)=90-alignment.rot+alignment.psi(i);
        params(1)=XX(alignment.di[i]);
        params(2)=YY(alignment.di[i]);
        DF.append_comment(fn_corrected);
        DF.append_data_line(params);
   }
   DF.write(fnRoot+"_correction_parameters.txt");
}

/* Remove outliers ---------------------------------------------------------*/
void Prog_tomograph_alignment::removeOutlierLandmarks(
    const Alignment &alignment)
{
    std::cout << "Removing outliers ...\n";

    // Compute threshold for outliers
    histogram1D hist;
    compute_hist(alignment.errorLandmark, hist, 100);
    double threshold0=hist.percentil(10);
    double thresholdF=hist.percentil(90);

    // Identify outliers
    std::vector<bool> validLandmark;
    int invalidCounter=0;
    int Nlandmark=YSIZE(allLandmarksX);
    for (int j=0; j<Nlandmark; j++)
        if (alignment.errorLandmark(j)<threshold0 ||
            alignment.errorLandmark(j)>thresholdF)
        {
            validLandmark.push_back(false);
            invalidCounter++;
        }
        else
            validLandmark.push_back(true);

    // Remove outliers
    Matrix2D<double> newAllLandmarksX(Nlandmark-invalidCounter,Nimg);
    Matrix2D<double> newAllLandmarksY(Nlandmark-invalidCounter,Nimg);
    int jj=0;
    for (int j=0; j<Nlandmark; j++)
    {
        if (!validLandmark[j]) continue;
        for (int i=0; i<Nimg; i++)
        {
            newAllLandmarksX(jj,i)=allLandmarksX(j,i);
            newAllLandmarksY(jj,i)=allLandmarksY(j,i);
        }
        jj++;
    }
    allLandmarksX=newAllLandmarksX;
    allLandmarksY=newAllLandmarksY;
    
    std::cout << invalidCounter << " out of " << Nlandmark
              << " landmarks have been removed\n";
}

/* Run --------------------------------------------------------------------- */
namespace TomographAlignment {
    const Prog_tomograph_alignment* global_prm;
}

double wrapperError(double *p, void *prm)
{
    Alignment alignment(TomographAlignment::global_prm);
    alignment=*(TomographAlignment::global_prm->bestPreviousAlignment);
    alignment.rot=p[1];
    alignment.tilt=p[2];
    return alignment.optimizeGivenAxisDirection();
}

void Prog_tomograph_alignment::run() {
    generateLandmarkSet();
    produceInformationFromLandmarks();
    Alignment *alignment=new Alignment(this);

    // Exhaustive search for rot
    double bestError=0, bestRot=-1;
    // for (double rot=100; rot<=100; rot+=deltaRot)
    for (double rot=0; rot<=180-deltaRot; rot+=deltaRot)
    {
        alignment->clear();
        alignment->rot=rot;
        double error=alignment->optimizeGivenAxisDirection();
        std::cout << "rot= " << rot
                  << " error= " << error << std::endl;
        if (bestRot<0 || bestError>error)
        {
            bestRot=rot;
            bestError=error;
            *bestPreviousAlignment=*alignment;
        }
    }
    std::cout << "Final bestRot=" << bestRot
              << " Final bestError=" << bestError << std::endl;
    
    // Continuous optimization for the axis direction
    Matrix1D<double> axisAngles(2), steps(2);
    axisAngles(0)=bestRot;
    axisAngles(1)=90;
    steps.initConstant(1);
    if (!optimizeTiltAngle) steps(1)=0;
    double fitness;
    int iter;
    TomographAlignment::global_prm=this;
    powellOptimizer(axisAngles,1,2,&wrapperError, NULL,
        0.01,fitness,iter,steps,true);

    // Outlier removal
    for (int i=0; i<3; i++)
    {
        // Compute the best alignment
        bestPreviousAlignment->rot=axisAngles(0);
        bestPreviousAlignment->tilt=axisAngles(1);
        fitness=bestPreviousAlignment->optimizeGivenAxisDirection();
        bestPreviousAlignment->computeErrorForLandmarks();

        // Remove landmarks that are outliers in the current model
        removeOutlierLandmarks(*bestPreviousAlignment);
        produceInformationFromLandmarks();
        delete bestPreviousAlignment;
        bestPreviousAlignment=new Alignment(this);
        bestPreviousAlignment->rot=axisAngles(0);
        bestPreviousAlignment->tilt=axisAngles(1);
        fitness=bestPreviousAlignment->optimizeGivenAxisDirection();

        // Optimize again
        powellOptimizer(axisAngles,1,2,&wrapperError,NULL,
            0.01,fitness,iter,steps,true);
    }
    bestPreviousAlignment->rot=axisAngles(0);
    bestPreviousAlignment->tilt=axisAngles(1);

    // Save the alignment
    std::ofstream fh_out;
    fh_out.open((fnRoot+"_alignment.txt").c_str());
    if (!fh_out)
        REPORT_ERROR(1,(std::string)"Cannot open "+fnRoot+"_alignment.txt for output");
    fh_out << *bestPreviousAlignment;
    fh_out.close();
    
    // Correct the input images
    alignImages(*bestPreviousAlignment);
}

/* Optimize for rot -------------------------------------------------------- */
#define DEBUG
double Alignment::optimizeGivenAxisDirection()
{
    double bestError;
    bool firstIteration=true, finish=false;
    int Niterations=0;
    do {
        computeGeometryDependentOfAxis();
        computeGeometryDependentOfRotation();
        double error=computeError();
        #ifdef DEBUG
            std::cout << "it=" << Niterations << " Error= " << error << std::endl;
        #endif
        updateModel();
        Niterations++;
        if (firstIteration)
        {
            bestError=error;
            firstIteration=false;
        }
        else
        {
            finish=((error>bestError) || (Niterations>1000) ||
               (ABS(error-bestError)/bestError<0.001)) && Niterations>20;
            if (error<bestError) bestError=error;
        }
    } while (!finish);
    return bestError;
}
#undef DEBUG

/* Compute the geometry part corresponding to axis ------------------------- */
void Alignment::computeGeometryDependentOfAxis()
{
    Matrix1D<double> axis;
    Euler_direction(rot, tilt, 0, axis);

    // Compute Aip, Aipt
    Matrix2D<double> Raxis;
    for (int i=0; i<Nimg; i++)
    {
        Raxis=rotation3DMatrix(prm->tiltList[i],axis);
        Aipt[i](0,0)=Aip[i](0,0)=Raxis(0,0);
        Aipt[i](1,0)=Aip[i](0,1)=Raxis(0,1);
        Aipt[i](2,0)=Aip[i](0,2)=Raxis(0,2); 
        Aipt[i](0,1)=Aip[i](1,0)=Raxis(1,0);
        Aipt[i](1,1)=Aip[i](1,1)=Raxis(1,1);
        Aipt[i](2,1)=Aip[i](1,2)=Raxis(1,2);
    }
}

/* Compute the geometry part corresponding to axis ------------------------- */
void Alignment::computeGeometryDependentOfRotation()
{
    // Compute Ai, Ait
    Matrix2D<double> Rinplane;
    for (int i=0; i<Nimg; i++)
    {
        Rinplane=rotation3DMatrix(-psi(i),'Z');
        Rinplane.resize(2,2);
        Ai[i]=Rinplane*Aip[i];
        Ait[i]=Ai[i].transpose();
    }
}

/* Compute error ----------------------------------------------------------- */
double Alignment::computeError() const
{
    Matrix1D<double> pijp;
    double error=0;
    double N=0;
    for (int i=0; i<Nimg; i++)
    {
        int jjmax=prm->Vseti[i].size();
        for (int jj=0; jj<jjmax; jj++)
        {
            int j=prm->Vseti[i][jj];
            pijp=Ai[i]*rj[j]+di[i];
            DIRECT_MAT_ELEM(allLandmarksPredictedX,j,i)=XX(pijp);
            DIRECT_MAT_ELEM(allLandmarksPredictedY,j,i)=YY(pijp);
            double diffx=DIRECT_MAT_ELEM(prm->allLandmarksX,j,i)-XX(pijp);
            double diffy=DIRECT_MAT_ELEM(prm->allLandmarksY,j,i)-YY(pijp);
            error+=diffx*diffx+diffy*diffy;
            N++;
        }
    }
    return sqrt(error/N);
}

void Alignment::computeErrorForLandmarks()
{
    Nlandmark=YSIZE(prm->allLandmarksX);
    errorLandmark.initZeros(Nlandmark);
    for (int j=0; j<Nlandmark; j++)
    {
        int counterj=0;
        for (int i=0; i<Nimg; i++)
        {
            if (prm->allLandmarksX(j,i)!=XSIZE(*(prm->img[0])))
            {
                double diffx=DIRECT_MAT_ELEM(prm->allLandmarksX,j,i)-
                    DIRECT_MAT_ELEM(allLandmarksPredictedX,j,i);
                double diffy=DIRECT_MAT_ELEM(prm->allLandmarksY,j,i)-
                    DIRECT_MAT_ELEM(allLandmarksPredictedY,j,i);
                errorLandmark(j)+=sqrt(diffx*diffx+diffy*diffy);
                counterj++;
            }
        }
        errorLandmark(j)/=counterj;
    }
}

/* Update model ------------------------------------------------------------ */
//#define DEBUG
void Alignment::updateModel()
{
    Matrix1D<double> pij(2), piN(2);
    Matrix2D<double> A(3,3), AitAi;
    Matrix1D<double> b(3);

    // Average of all landmarks from 1 to N-1 as reconstructed from the last
    // landmark
    int iimax=prm->Vsetj[Nlandmark-1].size();
    // double error0=computeError();
/*
    Matrix1D<double> rNlandmark(3);
    for (int ii=0; ii<iimax; ii++)
    {
        int i=prm->Vsetj[Nlandmark-1][ii];
        AitAi=Ait[i]*Ai[i];
        for (int j=0; j<Nlandmark-1; j++)
            rNlandmark+=AitAi*rj[j];
    }
    // std::cout << "rNlandmark=" << rNlandmark << std::endl;
*/
    // Update the 3D positions of the landmarks
//    Matrix1D<double> sumLandmarks(3);
//    for (int j=0; j<Nlandmark-1; j++)
    for (int j=0; j<Nlandmark; j++)
    {
        A.initZeros();
        b.initZeros();
        // Compute the part correponding to Vj
        iimax=prm->Vsetj[j].size();
        for (int ii=0; ii<iimax; ii++)
        {
            int i=prm->Vsetj[j][ii];
            XX(pij)=DIRECT_MAT_ELEM(prm->allLandmarksX,j,i);
            YY(pij)=DIRECT_MAT_ELEM(prm->allLandmarksY,j,i);
            A+=Ait[i]*Ai[i];
            b+=Ait[i]*(pij-di[i]);
        }

        // Compute the part corresponding to V_N
/*
        // std::cout << "A before\n" << A << "b before: " << b.transpose() << std::endl;
        iimax=prm->Vsetj[Nlandmark-1].size();
        for (int ii=0; ii<iimax; ii++)
        {
            int i=prm->Vsetj[Nlandmark-1][ii];
            XX(piN)=DIRECT_MAT_ELEM(prm->allLandmarksX,Nlandmark-1,i);
            YY(piN)=DIRECT_MAT_ELEM(prm->allLandmarksY,Nlandmark-1,i);
            AitAi=Ait[i]*Ai[i];
            A+=AitAi;
            b-=Ait[i]*(piN-di[i]);
            for (int jp=0; jp<Nlandmark-1; jp++)
            {
               if (jp!=j) b-=AitAi*rj[jp];
            }
            // b+=AitAi*rj[j];
        }
        // b-=rNlandmark;
        // std::cout << "A after\n" << A << "b after: " << b.transpose() << std::endl;
*/
        // Update rj[j]
        rj[j]=A.inv()*b;
//        sumLandmarks-=rj[j];
    }
//    rj[Nlandmark-1]=sumLandmarks;
    // double error1=computeError();

    // Compute the average landmarks seen in each image
    for (int i=0; i<Nimg; i++)
    {
        barri[i].initZeros();
        for (int jj=0; jj<prm->ni(i); jj++)
        {
            int j=prm->Vseti[i][jj];
            barri[i]+=rj[j];
        }
        barri[i]/=prm->ni(i);
    }

    // Update shifts
    for (int i=0; i<Nimg; i++)
        di[i]=prm->barpi[i]-Ai[i]*barri[i];
    // double error2=computeError();
    
    // Update rotations
    if (prm->psiMax>0)
    {
        Matrix2D<double> tmp1(2,2), tmp2(2,2), Ri(2,2);
        Matrix2D<double> Aiprj(2,1), Aiprjt(1,2), dim(2,1), Pij(2,1);
        for (int i=0; i<Nimg; i++)
        {
            tmp1.initZeros();
            tmp2.initZeros();
            dim.fromVector(di[i]);
            for (int jj=0; jj<prm->ni(i); jj++)
            {
                int j=prm->Vseti[i][jj];
                Aiprj.fromVector(Aip[i]*rj[j]);
                Aiprjt=Aiprj.transpose();
                tmp1+=Aiprj*Aiprjt;

                DIRECT_MAT_ELEM(Pij,0,0)=DIRECT_MAT_ELEM(prm->allLandmarksX,j,i);
                DIRECT_MAT_ELEM(Pij,1,0)=DIRECT_MAT_ELEM(prm->allLandmarksY,j,i);
                //tmp2+=Pij*Aiprjt+Aiprj*dit;
                tmp2+=(dim-Pij)*Aiprjt;
            }
            Ri=tmp1.inv()*tmp2;
            Ri=tmp2;
            psi(i)=CLIP(RAD2DEG(atan((Ri(0,1)-Ri(1,0))/(Ri(0,0)+Ri(1,1)))),
               -(prm->psiMax),prm->psiMax);
            Matrix2D<double> Rinplane;
            Rinplane=rotation3DMatrix(-psi(i),'Z');
            Rinplane.resize(2,2);
            // std::cout << "image i=" << i << " psi=" << psi(i) << " Ri\n"
            //          << Ri << "Rpsi\n" << Rinplane << std::endl;
        }
    }

    // Update the rotation dependent part    
    computeGeometryDependentOfRotation();

    // double error3=computeError();
    // std::cout << "Error optimization: " << error0 << " "
    //           << error1 << " " << error2 << " " << error3 << " " << std::endl;
}
#undef DEBUG

/* Print ------------------------------------------------------------------- */
std::ostream& operator << (std::ostream &out, Alignment &alignment)
{
    out << "Alignment parameters ===========================================\n"
        << "rot=" << alignment.rot << " tilt=" << alignment.tilt << std::endl;
    out << "Images ---------------------------------------------------------\n";
    for (int i=0; i<alignment.Nimg; i++)
        out << "Image " << i << " psi= " << alignment.psi(i)
            << " di= " << alignment.di[i].transpose() 
            << " (" << alignment.prm->ni(i) << ")\n";
    out << "Landmarks ------------------------------------------------------\n";
    alignment.computeErrorForLandmarks();
    for (int j=0; j<alignment.Nlandmark; j++)
        out << "Landmark " << j << " rj= " << alignment.rj[j].transpose()
            << " " << alignment.errorLandmark(j) << std::endl;
    return out;
}
