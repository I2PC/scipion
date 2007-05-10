:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Cuhre::usage =
	"Cuhre[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: Key::usage = "Key is an option of Cuhre.
	It specifies the basic integration rule:\n
	7 = use a degree-7 rule,\n
	9 = use a degree-9 rule,\n
	11 = use a degree-11 rule (available only in 3 dimensions),\n
	13 = use a degree-13 rule (available only in 2 dimensions),\n
	otherwise the default rule is used: the degree-13 rule in 2 dimensions, the degree-11 rule in 3 dimensions, else the degree-9 rule."

:Evaluate: MinPoints::usage = "MinPoints is an option of Cuhre.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Cuhre.
	It can take the values Last or All which determine whether only the last (largest) or all sets of samples collected on a subregion over the iterations contribute to the final result."

:Evaluate: Regions::usage = "Regions is an option of Cuhre.
	It specifies whether the regions into which the integration region has been cut are returned together with the integration results."

:Evaluate: Region::usage = "Region[ll, ur, res] describes a subregion:
	ll and ur are multidimensional equivalents of the region's lower left and upper right corner.
	res gives the integration results for the region in a list with entries of the form {integral, error} for each component of the integrand."

:Evaluate: Begin["`Cuhre`"]

:Begin:
:Function: Cuhre
:Pattern: MLCuhre[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, mineval_, maxeval_,
  key_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, mineval, maxeval,
  key}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer, Integer,
  Integer}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Cuhre] = {HoldFirst}

:Evaluate: Options[Cuhre] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000, Key -> 0,
	Verbose -> 1, Final -> Last, Regions -> False, Compiled -> True}

:Evaluate: Cuhre[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}],
	tags, vars, lower, range, jac, tmp, defs, integrand,
	rel, abs, mineval, maxeval, key, verbose, final, regions, compiled},
	  Message[Cuhre::optx, #, Cuhre]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Cuhre]];
	  {rel, abs, mineval, maxeval, key,
	    verbose, final, regions, compiled} =
	    tags /. {opt} /. Options[Cuhre];
	  {vars, lower, range} = Transpose[{v}];
	  jac = Simplify[Times@@ (range -= lower)];
	  tmp = Array[tmpvar, ndim];
	  defs = Simplify[lower + range tmp];
	  Block[{Set}, define[compiled, tmp, vars, Thread[vars = defs], jac]];
	  integrand = fun[f];
	  MLCuhre[ndim, ncomp[f], 10.^-rel, 10.^-abs,
	    Min[Max[verbose, 0], 3] +
	      If[final === Last, 4, 0] +
	      If[TrueQ[regions], 256, 0],
            mineval, maxeval, key]
	]

:Evaluate: tmpvar[n_] := ToExpression["Cuba`Cuhre`t" <> ToString[n]]

:Evaluate: Attributes[ncomp] = Attributes[fun] = {HoldAll}

:Evaluate: ncomp[f_List] := Length[f]

:Evaluate: _ncomp = 1

:Evaluate: define[True, tmp_, vars_, {defs__}, jac_] :=
	fun[f_] := Compile[tmp, Block[vars, defs; check[vars, Chop[f jac]//N]]]

:Evaluate: define[False, tmp_, vars_, {defs__}, jac_] :=
	fun[f_] := Function[tmp, Block[vars, defs; check[vars, Chop[f jac]//N]]]

:Evaluate: check[_, f_Real] = {f}

:Evaluate: check[_, f:{__Real}] = f

:Evaluate: check[x_, _] := (Message[Cuhre::badsample, ff, x]; {})

:Evaluate: sample[x_] :=
	Check[Apply[integrand, Partition[x, ndim], 1]//Flatten, {}]

:Evaluate: region[ll_, ur_, r___] :=
	Region[lower + range ll, lower + range ur, r]

:Evaluate: Cuhre::badsample = "`` is not a real-valued function at ``."

:Evaluate: Cuhre::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Cuhre::badcomp = "Cannot integrate `` components."

:Evaluate: Cuhre::accuracy =
	"Desired accuracy was not reached within `` function evaluations on `` subregions."

:Evaluate: Cuhre::success = "Needed `` function evaluations on `` subregions."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Cuhre.tm
		Adaptive integration using cubature rules
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
  MLPutSymbol(stdlink, "Cuhre");
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
  MLPutFunction(stdlink, "Cuba`Cuhre`sample", 1);
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

void Cuhre(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cint key)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) ) {
    Status("badcomp", ncomp, 0);
    MLPutSymbol(stdlink, "$Failed");
  }
  else if( BadDimension(ndim) ) {
    Status("baddim", ndim, 0);
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    real integral[NCOMP], error[NCOMP], prob[NCOMP];
    count comp;
    int fail;

    neval_ = 0;

    fail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key,
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

