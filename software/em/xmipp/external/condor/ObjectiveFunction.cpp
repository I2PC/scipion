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

#ifdef WIN32
//#include <windows.h>
#else
#include <unistd.h>
#endif

#include <string.h>
#include "ObjectiveFunction.h"
#include "tools.h"

void ObjectiveFunction::saveStats(char *resultsFile, Vector vG, Matrix mH, Vector vLambda)
{
    FILE *ff=fopen(resultsFile,"w");
    fprintf(ff,";dimension of search-space, total NFE, NFE before best point found, value OF at solution\n"
               "%i\t%i\t(%i)\t%e\n"
               ";Solution vector\n", dim(), nfe, nfe2, valueBest);
    xBest.save(ff,2);
    fprintf(ff,";Hessian matrix at the solution\n");
    mH.save(ff,3);
    fprintf(ff,";Gradient vector at the solution (should be zero if no active constraints)\n");
    vG.save(ff,2);
    if (isConstrained)
    {
        fprintf(ff,";Lagrangian Vector at the solution (lower,upper,linear,non-linear)\n");
        vLambda.save(ff,2);
    }
    fprintf(ff,"\n");
    fclose(ff);
}

char ObjectiveFunction::isFeasible(Vector vx, double *d)
{
    if (vx==Vector::emptyVector) return 1;
    if (!isConstrained) return 1;
    int i, dim=vx.sz(), nerror; 
    char feasible=1;
    double *bbl=bl, *bbu=bu, *x=vx, t;

    if (d) *d=0.0;
    initTolLC(vx);

    for (i=0; i<dim; i++)
        if ((t=bbl[i]-x[i])>tolLC) 
            { if (d) *d=mmax(*d,t); else return 0; feasible=0; }

    for (i=0; i<dim; i++)
        if ((t=x[i]-bbu[i])>tolLC) 
            { if (d) *d=mmax(*d,t); else return 0; feasible=0; }

    for (i=0; i<A.nLine(); i++)
        if ((t=b[i]-A.scalarProduct(i,vx))>tolLC) 
            { if (d) *d=mmax(*d,t); else return 0; feasible=0; }

    for (i=0; i<nNLConstraints; i++)
        if ((t=-evalNLConstraint(i,vx,&nerror))>tolNLC ) 
            { if (d) *d=mmax(*d,t); else return 0; feasible=0; }

//    printf("");
    return feasible;
}

void ObjectiveFunction::endInit()
 // init linear tolerances and init variable "isConstrained"
{
    int i,mdim=dim();
    double *bbl=bl, *bbu=bu;

    isConstrained=0;
    for (i=0; i<mdim; i++)
    {
        if (bbl[i]>-INF) {isConstrained=1; maxNormLC=mmax(maxNormLC, condorAbs(bbl[i])); }
        if (bbu[i]< INF) {isConstrained=1; maxNormLC=mmax(maxNormLC, condorAbs(bbu[i])); }
    }
    if (b.sz()) { isConstrained=1; maxNormLC=mmax(maxNormLC,b.LnftyNorm()); }
    tolLC=(1.0+maxNormLC)*tolRelFeasibilityForLC*(mdim*2+A.nLine());

    if (nNLConstraints) isConstrained=1;
}

void ObjectiveFunction::initTolLC(Vector vX)
{
    if (!isConstrained) return;
    int i;
    double *ofb=b;

    for (i=0; i<A.nLine(); i++)
        maxNormLC=mmax(maxNormLC, condorAbs(ofb[i]-A.scalarProduct(i,vX)));

    tolLC=(1.0+maxNormLC)*tolRelFeasibilityForLC*(dim()*2+A.nLine());
}

void ObjectiveFunction::initTolNLC(Vector c, double delta)
{
    int i;

    for (i=0; i<nNLConstraints; i++) maxNormNLC=mmax(maxNormNLC,condorAbs(c[i]));
    if (delta<INF) maxNormNLC=mmax(maxNormNLC,delta*delta);

    tolNLC=(1.0+maxNormNLC)*tolRelFeasibilityForNLC*nNLConstraints;
}

void ObjectiveFunction::updateCounter(double df, Vector vX, int nerror)
{
    nfe++;
    if ((dfold==INF)&&(nerror==0)) { dfref=(1+condorAbs(df))*1e-8; dfold=df; nfe2=nfe; return; }
    if (dfold-df<dfref) return;
    if (!isFeasible(vX)) return;
    
    if (nerror==0) 
    {
        nfe2=nfe;  
        dfold=df;
    }
}

void ObjectiveFunction::setSaveFile(char *s)
{
    char buffer[300];
    if (saveFileName) free(saveFileName);
    if (s==NULL)
    {
        strcpy(buffer,name); strcat(buffer,".dat"); s=buffer;
    }
    saveFileName=(char*)malloc(strlen(s)+1);
    strcpy(saveFileName,s);
}

void ObjectiveFunction::setName(char *s)
{
    char *p=s+strlen(s)-1;
    while ((*p!='.')&&(p!=s)) p--;
    if (p==s) { strncpy(name,s, 8); name[8]=0; return; }
    *p='\0';
    while ((*p!='\\')&&(*p!='/')&&(p!=s)) p--;
    if (p==s) { strncpy(name,s, 8); name[8]=0; return; }
    p++;
    strncpy(name,p, 8); name[8]=0;
}

void ObjectiveFunction::printStats(char cc)
{
    printf("\n\nProblem Name: %s\n",name);
    printf("Dimension of the search space: %i\n",dim());
    printf("best (lowest) value found: %e\n", valueBest+objectiveConst);
    printf("Number of function Evaluation: %i (%i)\n",nfe,nfe2);
    if (xOptimal.sz())
    {
        printf("Lnfty distance to the optimum: %e\n", xBest.LnftyDistance(xOptimal));
   //   printf("Euclidian distance to the optimum: %e\n", xBest.euclidianDistance(xOptimal));
    }
    int idim=xBest.sz(),j=0;
    if (idim<20)
    {
        double *dd=xBest;
        printf("Solution Vector is : \n[%e",dd[0]);
        for (j=1; j<idim; j++) printf(", %e",dd[j]);
        printf("]\n"); j=0;
    }
    if ((cc==0)||(!isConstrained)) return;
    double *dbl=bl,*dbu=bu;
    while (idim--)
    {
        if (*(dbl++)>-INF) j++;
        if (*(dbu++)< INF) j++;
    }
    printf("number of        box constraints:%i\n"
           "number of     linear constraints:%i\n"
           "number of non-linear constraints:%i\n",j,A.nLine(),nNLConstraints);

}

void ObjectiveFunction::saveValue(Vector tmp,double valueOF, int nerror)
{
    int nl=data.nLine(), mdim=tmp.sz();
    data.setNColumn(mdim+2);
    data.append(tmp);
    ((double**)data)[nl][mdim]=valueOF;
    ((double**)data)[nl][mdim+1]=nerror;
    if (saveFileName) data.updateSave(saveFileName);
}

int ObjectiveFunction::dim() 
{
    int n=xStart.sz();
    if (n>0) return n;
    return data.nColumn()-2;
}

#ifdef NO_OPTIMIZER
void projectionIntoFeasibleSpace(Vector vFrom, Vector vBase, ObjectiveFunction *of) { vBase=vFrom.clone(); }
#else
void projectionIntoFeasibleSpace(Vector vFrom, Vector vBase, ObjectiveFunction *of);
#endif

void ObjectiveFunction::addClosestFeasiblePointInData(Vector vX)
{
    double v;
    int nerror=0;
    // closest feasible point from vX
    if (isFeasible(vX))
    {
        v=eval(vX,&nerror);
        if (nerror)
        {
            printf("Evaluation of the Obj. Funct. at the starting point as failed.\n");
            exit(255);
        }
        saveValue(vX, v, 0);
        data.swapLines(0,data.nLine()-1);
        return;
    }

    double best;
    Vector b(vX.sz());
    projectionIntoFeasibleSpace(vX,b,this);
    if (!isFeasible(b,&best))
    {
        printf("unable to start (violation=%e).\n",best);
    }
    v=eval(b,&nerror);
    if (nerror)
    {
        printf("Unable to start.\n"
            "Evaluation of the Obj. Funct. at the feasible starting point as failed.\n"
            "Feasible starting point is:\n");
        b.print();
        exit(255);
    }
    saveValue(b,v,0);
    data.swapLines(0,data.nLine()-1);
}

void ObjectiveFunction::initData()
{
    if (data.nLine()==0)
    {
        initTolLC(xStart);
        addClosestFeasiblePointInData(xStart);
        return;
    }

    if (startPointIsGiven) 
    {
        int i=data.lineIndex(xStart);
        if (i!=-1)
        {
            data.swapLines(0,i);
            return;
        }
        initTolLC(xStart);
        addClosestFeasiblePointInData(xStart);
        return;
    }

    // find THE best point in the datas.
    int mdim=dim();
    Vector r(mdim);
    data.getLine(0,r,mdim);
    initTolLC(r);

    int k=-1,j,i=data.nLine();
    double v,best=INF,best2=INF;
    while (i--)
    {
        if (((double**)data)[i][mdim+1]) continue;
        v=((double**)data)[i][mdim];
        if (v<best2) { j=i; best2=v; }
        if (!isConstrained) continue;
        data.getLine(i,r,mdim);
        if (isFeasible(r)&&(v<best)) { k=i; best=v; }
    }

    if (!isConstrained) 
    {
        data.swapLines(0,j);
        return;
    }

    if (k!=-1)
    {
        data.swapLines(0,k);
        return;
    }

    data.getLine(j,r,mdim);
    addClosestFeasiblePointInData(r);
}

void ObjectiveFunction::initBounds()
{
    int dim=this->dim();
    bl.setSize(dim);
    bu.setSize(dim);
    double *dbl=bl,*dbu=bu;
    while (dim--)
    {
        *(dbl++)=-INF;
        *(dbu++)=INF;
    }
}

Vector ObjectiveFunction::evalGradNLConstraint(int j, Vector v, int *nerror)
{
    Vector R(dim());
    evalGradNLConstraint(j, v, R, nerror);
    return R;
}


void CorrectScaleOF::saveValue(Vector X,double valueOF, int nerror)
{
    int i=dim();
    double *x=X, *xr=xTemp, *re=rescaling;
    while (i--) xr[i]=re[i]*x[i];
    of->saveValue(xTemp,valueOF, nerror);
    ObjectiveFunction::saveValue(X,valueOF, nerror);
}

double CorrectScaleOF::eval(Vector X, int *nerror)
{
    int i=dim();
    double *x=X, *xr=xTemp, *re=rescaling;
    while (i--) xr[i]=re[i]*x[i];
    double r=of->eval(xTemp,nerror);
    updateCounter(r,X,*nerror);
    return r;
}

double CorrectScaleOF::evalNLConstraint(int j, Vector X, int *nerror)
{
    int i=dim();
    double *x=X, *xr=xTemp, *re=rescaling;
    while (i--) xr[i]=re[i]*x[i];
    return of->evalNLConstraint(j,xTemp,nerror);
}

void CorrectScaleOF::evalGradNLConstraint(int j, Vector X, Vector result, int *nerror)
{
    int i=dim();
    double *x=X, *xr=xTemp, *re=rescaling;
    while (i--) xr[i]=re[i]*x[i];
    of->evalGradNLConstraint(j,xTemp,result,nerror);
}

CorrectScaleOF::CorrectScaleOF(int _t, ObjectiveFunction *_of):
    of(_of)
{
    t=_t;

    int i=of->dim();
    rescaling.setSize(i);
    double *xs=of->xStart,*r=rescaling;
    while (i--) r[i]=condorAbs(xs[i])+1.0;

    if (of->isConstrained)
    {
        double *bl=of->bl, *bu=of->bu;
        r=rescaling; i=of->dim();
        while (i--)
        {
            if ((bl[i]>-INF)&&(bu[i]<INF)) { r[i]=bu[i]-bl[i]; continue; }
            if ((r[i]==0.0 )&&(bu[i]<INF))  { r[i]=bu[i];       continue; }
            if (r[i]==0.0) r[i]=1.0;
        }
    }
    init();
}

CorrectScaleOF::CorrectScaleOF(int _t, ObjectiveFunction *_of, Vector _rescaling):
    rescaling(_rescaling), of(_of)
{
    t=_t;
    if ((int)_rescaling.sz()!=_of->dim())
    {
        printf("Error in rescaling vector: dimension do not agree.\n");
        exit(254);
    }
    init();
}

void CorrectScaleOF::init()
{
    double *xos=of->xOptimal, *xss=of->xStart, *xod, *xsd, 
           *r=rescaling, **datas=of->data, **datad,
           *bls=of->bl, *bus=of->bu, *bld, *bud, **as=of->A, **ad;
    int n=of->dim(), i=n,j;
    strcpy(name,"SCALING");

    xTemp.setSize(n);
    xOptimal.setSize(n);                xod=xOptimal; 
    xStart.setSize(n);                  xsd=xStart;
    data.setSize(of->data.nLine(),n+2); datad=data;

    while(i--)
    {
        if (xos) xod[i]=xos[i]/r[i];
        xsd[i]=xss[i]/r[i];
        j=data.nLine();
        while (j--) datad[j][i]=datas[j][i]/r[i];
    }
    j=data.nLine();
    while (j--) { datad[j][n]=datas[j][n]; datad[j][n+1]=datas[j][n+1]; }

    startPointIsGiven=of->startPointIsGiven;
    valueOptimal=of->valueOptimal;
    noiseAbsolute=of->noiseAbsolute;
    noiseRelative=of->noiseRelative;
    objectiveConst=of->objectiveConst;

    if (of->isConstrained==0) { isConstrained=0; return; }

    // there are (box&linear) constraints: scale them !
    isConstrained=of->isConstrained;
    nNLConstraints=of->nNLConstraints;
    bl.setSize(n);                      bld=bl;
    bu.setSize(n);                      bud=bu;
    A.setSize(of->A.nLine(),n);         ad=A;
    b=of->b;

    i=n;
    while(i--)
    {
        bld[i]=bls[i]/r[i];
        bud[i]=bus[i]/r[i];
        j=A.nLine();
        while (j--) ad[j][i]=as[j][i]*r[i];
    }
}

void CorrectScaleOF::finalize(Vector vG, Matrix mH, Vector vLambda)
{
    of->xBest.copyFrom(xBest);
    of->xBest.oneByOneMutiply(rescaling);

    of->valueBest=valueBest;

    // rescale vG,mH,vLambda
    rescaling.oneByOneInvert();

    vG.oneByOneMutiply(rescaling);
    mH.multiplyByDiagonalMatrix(rescaling);
    rescaling.diagonalizeAndMultiply(mH);
    vLambda.oneByOneMutiply(rescaling);
    of->finalize(vG,mH,vLambda);
}

/*
#ifdef __INCLUDE_SIF__

#include "sif/SIFFunction.h"

// extern elfunType elfunPARKCH_; extern groupType groupPARKCH_; 
extern elfunType elfunAkiva_;    extern groupType groupAkiva_;
extern elfunType elfunRosen_;    extern groupType groupRosen_;
extern elfunType elfunALLINITU_; extern groupType groupALLINITU_;
extern elfunType elfunSTRATEC_;  extern groupType groupSTRATEC_;
extern elfunType elfunTOINTGOR_; extern groupType groupTOINTGOR_;
extern elfunType elfunTOINTPSP_; extern groupType groupTOINTPSP_;
extern elfunType elfun3PK_;      extern groupType group3PK_;
extern elfunType elfunBIGGS6_;   extern groupType groupBIGGS6_;
extern elfunType elfunBROWNDEN_; extern groupType groupBROWNDEN_;
extern elfunType elfunDECONVU_;  extern groupType groupDECONVU_;
extern elfunType elfunHEART_;    extern groupType groupHEART_;
extern elfunType elfunOSBORNEB_; extern groupType groupOSBORNEB_;
extern elfunType elfunVIBRBEAM_; extern groupType groupVIBRBEAM_;
extern elfunType elfunKOWOSB_;   extern groupType groupKOWOSB_;
extern elfunType elfunHELIX_;    extern groupType groupHELIX_;

extern elfunType elfunCRAGGLVY_; extern groupType groupCRAGGLVY_;
extern elfunType elfunEIGENALS_; extern groupType groupEIGENALS_;
extern elfunType elfunHAIRY_;    extern groupType groupHAIRY_;
extern elfunType elfunPFIT1LS_;  extern groupType groupPFIT1LS_;
extern elfunType elfunVARDIM_;   extern groupType groupVARDIM_;
extern elfunType elfunMANCINO_;  extern groupType groupMANCINO_;
extern elfunType elfunPOWER_;    extern groupType groupPOWER_;
extern elfunType elfunHATFLDE_;  extern groupType groupHATFLDE_;
extern elfunType elfunWATSON_;   extern groupType groupWATSON_;
extern elfunType elfunFMINSURF_; extern groupType groupFMINSURF_;
extern elfunType elfunDIXMAANK_; extern groupType groupDIXMAANK_;
extern elfunType elfunMOREBV_;   extern groupType groupMOREBV_;
extern elfunType elfunBRYBND_;   extern groupType groupBRYBND_;
extern elfunType elfunSCHMVETT_; extern groupType groupSCHMVETT_;
extern elfunType elfunHEART6LS_; extern groupType groupHEART6LS_;
extern elfunType elfunBROWNAL_;  extern groupType groupBROWNAL_;
extern elfunType elfunDQDRTIC_;  extern groupType groupDQDRTIC_;
extern elfunType elfunGROWTHLS_; extern groupType groupGROWTHLS_;
extern elfunType elfunSISSER_;   extern groupType groupSISSER_;
extern elfunType elfunCLIFF_;    extern groupType groupCLIFF_;
extern elfunType elfunGULF_;     extern groupType groupGULF_;
extern elfunType elfunSNAIL_;    extern groupType groupSNAIL_;
extern elfunType elfunHART6_;    extern groupType groupHART6_;

#endif

#ifdef __INCLUDE_AMPL__
#include "ampl/AMPLof.h"
#endif

#include "simpleObjFunctions.h"

ObjectiveFunction *getObjectiveFunction(int i, double *rho)
{
    int n=2;
    ObjectiveFunction *of=NULL;
    double rhoEnd=-1;

    switch (i)
    {
    // first choice: internally coded functions:
    case  1: of=new Rosenbrock(i); break; // n=2;
    case  2: of=new BADScaleRosenbrock(i); break; // n=2;
    case  3: of=new FletcherTest(i); break; // n=2;
    case  4: of=new SuperSimpleConstrainedObjectiveFunction(i); break; //n=2;
    case  5: of=new FletcherTest2(i); break; // n=3;
    case  6: of=new NoisyRosenbrock(i); break; //n=2;
    case  7: of=new NoisyQuadratic(i); break; //n=2;
    case  8: of=new SimpleQuadratic(i); break; //n=2;
// second choice: create new random objective function
    case  20: of=new RandomOF(i+1,n); ((RandomOF *)of)->save("test.dat"); break; 

    // third choice : reload from disk previous random objective function
    case  21: of=new RandomOF(i,"test.dat"); break; 

#ifdef __INCLUDE_SIF__

    // fourth choice: use SIF file
    case 104: of= new SIFFunction(i,"sif/examples/akiva.d"     ,elfunAkiva_    ,groupAkiva_   ); break; // 2
    case 105: of= new SIFFunction(i,"sif/examples/allinitu.d"  ,elfunALLINITU_ ,groupALLINITU_); break; // 4
    case 106: of= new SIFFunction(i,"sif/examples/stratec.d"   ,elfunSTRATEC_  ,groupSTRATEC_ ); break; // 10 
    case 107: of= new SIFFunction(i,"sif/examples/heart.d"     ,elfunHEART_    ,groupHEART_   ); break; // 8
    case 108: of= new SIFFunction(i,"sif/examples/osborneb.d"  ,elfunOSBORNEB_ ,groupOSBORNEB_); break; // 11
    case 109: of= new SIFFunction(i,"sif/examples/vibrbeam.d"  ,elfunVIBRBEAM_ ,groupVIBRBEAM_); break; // 8
    case 110: of= new SIFFunction(i,"sif/examples/kowosb.d"    ,elfunKOWOSB_   ,groupKOWOSB_  ); break; // 4
    case 111: of= new SIFFunction(i,"sif/examples/helix.d"     ,elfunHELIX_    ,groupHELIX_   ); break; // 3

    case 112: of= new SIFFunction(i,"sif/examples/rosenbrock.d",elfunRosen_    ,groupRosen_   ); rhoEnd= 5e-3; break; // 2
    case 114: of= new SIFFunction(i,"sif/examples/sisser.d"    ,elfunSISSER_   ,groupSISSER_  ); rhoEnd= 1e-2; break; // 2
    case 115: of= new SIFFunction(i,"sif/examples/cliff.d"     ,elfunCLIFF_    ,groupCLIFF_   ); rhoEnd= 1e-3; break; // 2
    case 116: of= new SIFFunction(i,"sif/examples/hairy.d"     ,elfunHAIRY_    ,groupHAIRY_   ); rhoEnd= 2e-3; break; // 2
    case 117: of= new SIFFunction(i,"sif/examples/pfit1ls.d"   ,elfunPFIT1LS_  ,groupPFIT1LS_ ); rhoEnd= 1e-2; break; // 3
    case 118: of= new SIFFunction(i,"sif/examples/hatflde.d"   ,elfunHATFLDE_  ,groupHATFLDE_ ); rhoEnd=12e-3; break; // 3
    case 119: of= new SIFFunction(i,"sif/examples/schmvett.d"  ,elfunSCHMVETT_ ,groupSCHMVETT_); rhoEnd= 1e-2; break; // 3
    case 120: of= new SIFFunction(i,"sif/examples/growthls.d"  ,elfunGROWTHLS_ ,groupGROWTHLS_); rhoEnd= 5e-3; break; // 3
    case 121: of= new SIFFunction(i,"sif/examples/gulf.d"      ,elfunGULF_     ,groupGULF_    ); rhoEnd= 5e-2; break; // 3
    case 122: of= new SIFFunction(i,"sif/examples/brownden.d"  ,elfunBROWNDEN_ ,groupBROWNDEN_); rhoEnd=57e-2; break; // 4
    case 123: of= new SIFFunction(i,"sif/examples/eigenals.d"  ,elfunEIGENALS_ ,groupEIGENALS_); rhoEnd= 1e-2; break; // 6
    case 124: of= new SIFFunction(i,"sif/examples/heart6ls.d"  ,elfunHEART6LS_ ,groupHEART6LS_); rhoEnd= 5e-2; break; // 6
    case 125: of= new SIFFunction(i,"sif/examples/biggs6.d"    ,elfunBIGGS6_   ,groupBIGGS6_  ); rhoEnd= 6e-2; break; // 6
    case 126: of= new SIFFunction(i,"sif/examples/hart6.d"     ,elfunHART6_    ,groupHART6_   ); rhoEnd= 2e-1; break; // 6
    case 127: of= new SIFFunction(i,"sif/examples/cragglvy.d"  ,elfunCRAGGLVY_ ,groupCRAGGLVY_); rhoEnd= 6e-2; break; // 10
    case 128: of= new SIFFunction(i,"sif/examples/vardim.d"    ,elfunVARDIM_   ,groupVARDIM_  ); rhoEnd= 1e-3; break; // 10
    case 129: of= new SIFFunction(i,"sif/examples/mancino.d"   ,elfunMANCINO_  ,groupMANCINO_ ); rhoEnd= 1e-6; break; // 10
    case 130: of= new SIFFunction(i,"sif/examples/power.d"     ,elfunPOWER_    ,groupPOWER_   ); rhoEnd= 2e-2; break; // 10
    case 131: of= new SIFFunction(i,"sif/examples/morebv.d"    ,elfunMOREBV_   ,groupMOREBV_  ); rhoEnd= 1e-1; break; // 10
    case 132: of= new SIFFunction(i,"sif/examples/brybnd.d"    ,elfunBRYBND_   ,groupBRYBND_  ); rhoEnd= 6e-3; break; // 10
    case 133: of= new SIFFunction(i,"sif/examples/brownal.d"   ,elfunBROWNAL_  ,groupBROWNAL_ ); rhoEnd= 8e-3; break; // 10
    case 134: of= new SIFFunction(i,"sif/examples/dqdrtic.d"   ,elfunDQDRTIC_  ,groupDQDRTIC_ ); rhoEnd= 1e-3; break; // 10
    case 135: of= new SIFFunction(i,"sif/examples/watson.d"    ,elfunWATSON_   ,groupWATSON_  ); rhoEnd= 4e-2; break; // 12
    case 137: of= new SIFFunction(i,"sif/examples/fminsurf.d"  ,elfunFMINSURF_ ,groupFMINSURF_); rhoEnd= 1e-1; break; // 16

    case 138: of= new SIFFunction(i,"sif/examples/tointgor.d"  ,elfunTOINTGOR_ ,groupTOINTGOR_); break; // 50
    case 139: of= new SIFFunction(i,"sif/examples/tointpsp.d"  ,elfunTOINTPSP_ ,groupTOINTPSP_); break; // 50
    case 140: of= new SIFFunction(i,"sif/examples/3pk.d"       ,elfun3PK_      ,group3PK_     ); break; // 30
    case 141: of= new SIFFunction(i,"sif/examples/deconvu.d"   ,elfunDECONVU_  ,groupDECONVU_ ); break; // 61
//  case 142: of= new SIFFunction(i,"sif/examples/parkch.d"    ,elfunPARKCH_   ,groupPARKCH_  ); break; // 15 

#ifdef WIN32
    case 113: of= new SIFFunction(i,"sif/examples/snail.d"     ,elfunSNAIL_    ,groupSNAIL_   ); rhoEnd= 2e-4; break; // 2
    case 136: of= new SIFFunction(i,"sif/examples/dixmaank.d"  ,elfunDIXMAANK_ ,groupDIXMAANK_); rhoEnd= 3e-1; break; // 15
#else
    case 113: of= new SIFFunction(i,"sif/examples/snail.d"     ,elfunSNAIL_    ,groupSNAIL_   ); rhoEnd= 7e-4; break; // 2
    case 136: of= new SIFFunction(i,"sif/examples/dixmaank.d"  ,elfunDIXMAANK_ ,groupDIXMAANK_); rhoEnd= 4e-1; break; // 15
#endif

#endif
#ifdef __INCLUDE_AMPL__

    case 200: of= new AMPLObjectiveFunction(i,"ampl/examples/hs022.nl",1.0); break;
    case 201: of= new AMPLObjectiveFunction(i,"ampl/examples/hs023.nl"); break;
    case 202: of= new AMPLObjectiveFunction(i,"ampl/examples/hs026.nl"); break;
    case 203: of= new AMPLObjectiveFunction(i,"ampl/examples/hs034.nl"); break;
    case 204: of= new AMPLObjectiveFunction(i,"ampl/examples/hs038.nl"); break;
    case 205: of= new AMPLObjectiveFunction(i,"ampl/examples/hs044.nl"); break;
    case 206: of= new AMPLObjectiveFunction(i,"ampl/examples/hs065.nl"); break;
    case 207: of= new AMPLObjectiveFunction(i,"ampl/examples/hs076.nl"); break;
    case 208: of= new AMPLObjectiveFunction(i,"ampl/examples/hs100.nl"); break;
    case 209: of= new AMPLObjectiveFunction(i,"ampl/examples/hs106.nl"); break;
    case 210: of= new AMPLObjectiveFunction(i,"ampl/examples/hs108.nl"); rhoEnd= 1e-5; break;
    case 211: of= new AMPLObjectiveFunction(i,"ampl/examples/hs116.nl"); break;
    case 212: of= new AMPLObjectiveFunction(i,"ampl/examples/hs268.nl"); break;
    
    case 250: of= new AMPLObjectiveFunction(i,"ampl/examples/rosenbr.nl" ); rhoEnd= 5e-3; break; // 2
    case 251: of= new AMPLObjectiveFunction(i,"ampl/examples/sisser.nl"  ); rhoEnd= 1e-2; break; // 2
    case 252: of= new AMPLObjectiveFunction(i,"ampl/examples/cliff.nl"   ); rhoEnd= 1e-3; break; // 2
    case 253: of= new AMPLObjectiveFunction(i,"ampl/examples/hairy.nl"   ); rhoEnd= 2e-2; break; // 2
    case 254: of= new AMPLObjectiveFunction(i,"ampl/examples/pfit1ls.nl" ); rhoEnd= 1e-2; break; // 3
    case 255: of= new AMPLObjectiveFunction(i,"ampl/examples/hatflde.nl" ); rhoEnd=12e-3; break; // 3
    case 256: of= new AMPLObjectiveFunction(i,"ampl/examples/growthls.nl"); rhoEnd= 5e-3; break; // 3
    case 257: of= new AMPLObjectiveFunction(i,"ampl/examples/gulf.nl"    ); rhoEnd= 5e-2; break; // 3
    case 258: of= new AMPLObjectiveFunction(i,"ampl/examples/brownden.nl"); rhoEnd=57e-2; break; // 4
    case 259: of= new AMPLObjectiveFunction(i,"ampl/examples/eigenals.nl"); rhoEnd= 1e-2; break; // 6
    case 260: of= new AMPLObjectiveFunction(i,"ampl/examples/heart6ls.nl"); rhoEnd= 1e-2; break; // 6
    case 261: of= new AMPLObjectiveFunction(i,"ampl/examples/biggs6.nl"  ); rhoEnd= 6e-2; break; // 6
    case 262: of= new AMPLObjectiveFunction(i,"ampl/examples/hart6.nl"   ); rhoEnd= 2e-1; break; // 6
    case 263: of= new AMPLObjectiveFunction(i,"ampl/examples/cragglvy.nl"); rhoEnd= 1e-2; break; // 10
    case 264: of= new AMPLObjectiveFunction(i,"ampl/examples/vardim.nl"  ); rhoEnd= 1e-3; break; // 10
    case 265: of= new AMPLObjectiveFunction(i,"ampl/examples/mancino.nl" ); rhoEnd= 1e-6; break; // 10
    case 266: of= new AMPLObjectiveFunction(i,"ampl/examples/power.nl"   ); rhoEnd= 2e-2; break; // 10
    case 267: of= new AMPLObjectiveFunction(i,"ampl/examples/morebv.nl"  ); rhoEnd= 1e-1; break; // 10
    case 268: of= new AMPLObjectiveFunction(i,"ampl/examples/brybnd.nl"  ); rhoEnd= 3e-3; break; // 10
    case 269: of= new AMPLObjectiveFunction(i,"ampl/examples/brownal.nl" ); rhoEnd= 8e-3; break; // 10
    case 270: of= new AMPLObjectiveFunction(i,"ampl/examples/dqdrtic.nl" ); rhoEnd= 1e-3; break; // 10
    case 271: of= new AMPLObjectiveFunction(i,"ampl/examples/watson.nl"  ); rhoEnd= 4e-2; break; // 12
    case 272: of= new AMPLObjectiveFunction(i,"ampl/examples/dixmaank.nl"); rhoEnd= 3e-1; break; // 15
    case 273: of= new AMPLObjectiveFunction(i,"ampl/examples/fminsurf.nl"); rhoEnd= 1e-1; break; // 16

#endif

    }
    if ((i>=250)&&(i<=273)) of->isConstrained=0;

    if ((rho)&&(rhoEnd!=-1)) *rho=rhoEnd;
    if (of==NULL)
    {
        printf("No objective function defined under index %i\n",i);
        exit(255);
    }
    return of;
}
*/


