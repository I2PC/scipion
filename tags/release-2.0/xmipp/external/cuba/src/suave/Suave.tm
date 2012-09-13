:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Suave::usage =
	"Suave[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: NNew::usage = "NNew is an option of Suave.
	It specifies the number of new integrand evaluations in each subdivision."

:Evaluate: Flatness::usage = "Flatness is an option of Suave.
	It determines how prominently individual samples with a large fluctuation figure in the total fluctuation, which in turn determines how a region is split up.
	Explicitly, if F[i] is the individual fluctuation of sample i, the total fluctuation is computed as Sum[(1 + F[i])^p, {i, nsamples}]^(2/3/p), i.e. as the p-norm of the fluctuation vector to the power 2/3, where p is the number given by Flatness.
	Thus with increasing p, the fluctuation becomes more and more dominated by outliers, i.e. points with a large fluctuation.
	As suggested by the name Flatness, p should be chosen large for `flat' integrands and small for `volatile' integrands with high peaks.
	Note that since p appears in the exponent, one should not use too large values (say, no more than a few hundred) lest terms be truncated internally to prevent overflow."

:Evaluate: MinPoints::usage = "MinPoints is an option of Suave.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Suave.
	It can take the values Last or All which determine whether only the last (largest) or all sets of samples collected on a subregion over the iterations contribute to the final result."

:Evaluate: PseudoRandom::usage = "PseudoRandom is an option of Suave.
	If set to True, pseudo-random numbers are used instead of Sobol quasi-random numbers."

:Evaluate: Regions::usage = "Regions is an option of Suave.
	It specifies whether the regions into which the integration region has been cut are returned together with the integration results."

:Evaluate: Region::usage = "Region[ll, ur, res, df] describes a subregion:
	ll and ur are multidimensional equivalents of the region's lower left and upper right corner.
	res gives the integration results for the region in a list with entries of the form {integral, error, chi-square} for each component of the integrand.
	df is the number of degrees of freedom corresponding to the chi-square values in res."

:Evaluate: Begin["`Suave`"]

:Begin:
:Function: Suave
:Pattern: MLSuave[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, mineval_, maxeval_,
  nnew_, flatness_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, mineval, maxeval,
  nnew, flatness}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer, Integer,
  Integer, Real}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Suave] = {HoldFirst}

:Evaluate: Options[Suave] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000, NNew -> 1000, Flatness -> 50,
	Verbose -> 1, Final -> Last, PseudoRandom -> False,
	Regions -> False, Compiled -> True}

:Evaluate: Suave[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}],
	tags, vars, lower, range, jac, tmp, defs, integrand,
	rel, abs, mineval, maxeval, nnew, flatness,
	verbose, final, pseudo, regions, compiled},
	  Message[Suave::optx, #, Suave]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Suave]];
	  {rel, abs, mineval, maxeval, nnew, flatness,
	    verbose, final, pseudo, regions, compiled} =
	    tags /. {opt} /. Options[Suave];
	  {vars, lower, range} = Transpose[{v}];
	  jac = Simplify[Times@@ (range -= lower)];
	  tmp = Array[tmpvar, ndim];
	  defs = Simplify[lower + range tmp];
	  Block[{Set}, define[compiled, tmp, vars, Thread[vars = defs], jac]];
	  integrand = fun[f];
	  MLSuave[ndim, ncomp[f], 10.^-rel, 10.^-abs,
	    Min[Max[verbose, 0], 3] +
	      If[final === Last, 4, 0] +
	      If[TrueQ[pseudo], 8, 0] +
	      If[TrueQ[regions], 256, 0],
            mineval, maxeval, nnew, flatness]
	]

:Evaluate: tmpvar[n_] := ToExpression["Cuba`Suave`t" <> ToString[n]]

:Evaluate: Attributes[ncomp] = Attributes[fun] = {HoldAll}

:Evaluate: ncomp[f_List] := Length[f]

:Evaluate: _ncomp = 1

:Evaluate: define[True, tmp_, vars_, {defs__}, jac_] :=
	fun[f_] := Compile[tmp, Block[vars, defs; check[vars, Chop[f jac]//N]]]

:Evaluate: define[False, tmp_, vars_, {defs__}, jac_] :=
	fun[f_] := Function[tmp, Block[vars, defs; check[vars, Chop[f jac]//N]]]

:Evaluate: check[_, f_Real] = {f}

:Evaluate: check[_, f:{__Real}] = f

:Evaluate: check[x_, _] := (Message[Suave::badsample, ff, x]; {})

:Evaluate: sample[x_] :=
	Check[Apply[integrand, Partition[x, ndim], 1]//Flatten, {}]

:Evaluate: region[ll_, ur_, r___] :=
	Region[lower + range ll, lower + range ur, r]

:Evaluate: Suave::badsample = "`` is not a real-valued function at ``."

:Evaluate: Suave::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Suave::badcomp = "Cannot integrate `` components."

:Evaluate: Suave::accuracy =
	"Desired accuracy was not reached within `` function evaluations on `` subregions."

:Evaluate: Suave::success = "Needed `` function evaluations on `` subregions."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Suave.tm
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 1 Mar 06 th
*/


#include <setjmp.h>
#include "mathlink.h"
#include "util.c"

jmp_buf abort_;

/*********************************************************************/

static void Status(MLCONST char *msg, cint n1, cint n2)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 3);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "Suave");
  MLPutString(stdlink, msg);
  MLPutInteger(stdlink, n1);
  MLPutInteger(stdlink, n2);
}

/*********************************************************************/

static void Print(MLCONST char *s)
{
  int pkt;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Print", 1);
  MLPutString(stdlink, s);
  MLEndPacket(stdlink);

  do {
    pkt = MLNextPacket(stdlink);
    MLNewPacket(stdlink);
  } while( pkt != RETURNPKT );
}

/*********************************************************************/

static void DoSample(cnumber n, real *x, real *f)
{
  int pkt;
  real *mma_f;
  long mma_n;

  if( MLAbort ) goto abort;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Cuba`Suave`sample", 1);
  MLPutRealList(stdlink, x, n*ndim_);
  MLEndPacket(stdlink);

  while( (pkt = MLNextPacket(stdlink)) && (pkt != RETURNPKT) )
    MLNewPacket(stdlink);

  if( !MLGetRealList(stdlink, &mma_f, &mma_n) ) {
    MLClearError(stdlink);
    MLNewPacket(stdlink);
abort:
    MLPutFunction(stdlink, "Abort", 0);
    longjmp(abort_, 1);
  }

  if( mma_n != n*ncomp_ ) {
    MLDisownRealList(stdlink, mma_f, mma_n);
    MLPutSymbol(stdlink, "$Failed");
    longjmp(abort_, 1);
  }

  Copy(f, mma_f, n*ncomp_);
  MLDisownRealList(stdlink, mma_f, mma_n);

  neval_ += n;
}

/*********************************************************************/

#include "common.c"

void Suave(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) ) {
    Status("badcomp", ncomp, 0);
    MLPutSymbol(stdlink, "$Failed");
  }
  else if( BadDimension(ndim, flags) ) {
    Status("baddim", ndim, 0);
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    real integral[NCOMP], error[NCOMP], prob[NCOMP];
    count comp;
    int fail;

    neval_ = 0;

    fail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, nnew, flatness,
      integral, error, prob);

    Status(fail ? "accuracy" : "success", neval_, nregions_);

    MLPutFunction(stdlink, "List", ncomp);
    for( comp = 0; comp < ncomp; ++comp ) {
      real res[] = {integral[comp], error[comp], prob[comp]};
      MLPutRealList(stdlink, res, Elements(res));
    }
  }

  MLEndPacket(stdlink);
}

/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

