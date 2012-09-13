:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Divonne::usage =
	"Divonne[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: Key1::usage = "Key1 is an option of Divonne.
	It determines sampling in the partitioning phase.\n
	Special cases:\n
	  Key1 = 7: use a degree-7 cubature rule,\n
	  Key1 = 9: use a degree-9 cubature rule,\n
	  Key1 = 11: use a degree-11 cubature rule (available only in 3 dimensions),\n
	  Key1 = 13: use a degree-13 cubature rule (available only in 2 dimensions),\n
	otherwise a quasi-random sample of n1 = Abs[Key1] points is used, where the sign of Key1 determines the type of sample:\n
	  Key1 > 0: use a Korobov quasi-random sample,\n
	  Key1 < 0: use a Sobol quasi-random sample."

:Evaluate: Key2::usage = "Key2 is an option of Divonne.
	It determines sampling in the final integration phase.\n
	Special cases:\n
	  Key2 = 7: use a degree-7 cubature rule,\n
	  Key2 = 9: use a degree-9 cubature rule,\n
	  Key2 = 11: use a degree-11 cubature rule (available only in 3 dimensions),\n
	  Key2 = 13: use a degree-13 cubature rule (available only in 2 dimensions),\n
	otherwise a quasi-random sample is used, where the sign of Key2 determines the type of sample:\n
	  Key2 > 0: use a Korobov quasi-random sample,\n
	  Key2 < 0: use a Sobol quasi-random sample,\n
	and n2 = Abs[Key2] determines the number of points:\n
	  n2 >= 40: sample n2 points,\n
	  n2 < 40: sample n2*nneed points, where nneed is the number of points needed to reach the prescribed accuracy, as estimated by Divonne from the results of the partitioning phase."

:Evaluate: Key3::usage = "Key3 is an option of Divonne.
	It sets the strategy for the refinement phase:\n
	  Key3 = 0: do not further treat the subregion,\n
	  Key3 = 1: split the subregion up once more,\n
	for other values the region is sampled a third time:\n
	  Key3 = 7: use a degree-7 cubature rule,\n
	  Key3 = 9: use a degree-9 cubature rule,\n
	  Key3 = 11: use a degree-11 cubature rule (available only in 3 dimensions),\n
	  Key3 = 13: use a degree-13 cubature rule (available only in 2 dimensions),\n
	otherwise a quasi-random sample is used, where the sign of Key3 determines the type of sample:\n
	  Key3 > 0: use a Korobov quasi-random sample,\n
	  Key3 < 0: use a Sobol quasi-random sample,\n
	and n3 = Abs[Key3] determines the number of points:\n
	  n3 >= 40: sample n3 points,\n
	  n3 < 40: sample n3*nneed points, where nneed is the number of points needed to reach the prescribed accuracy, as estimated by Divonne from the results of the partitioning phase."

:Evaluate: MaxPass::usage = "MaxPass is an option of Divonne.
	It controls the partitioning termination.
	The partitioning phase is terminated when the estimated total number of integrand evaluations (partitioning plus final integration) does not decrease for MaxPass successive iterations."

:Evaluate: Border::usage = "Border is an option of Divonne.
	It specifies the width of the border of the integration region.
	Points falling into this border region are not sampled directly, but are extrapolated from two samples from the interior.
	The border width always refers to the unit hypercube, i.e. it is not rescaled if the integration region is not the unit hypercube."

:Evaluate: MaxChisq::usage = "MaxChisq is an option of Divonne.
	It specifies the maximum chi-square value a single subregion is allowed to have in the final integration phase.
	Regions which fail this chi-square test and whose sample averages differ by more than MinDeviation move on to the refinement phase."

:Evaluate: MinDeviation::usage = "MinDeviation is an option of Divonne.
	Regions which fail the chi-square test are not treated further if their sample averages differ by less than MinDeviation.
	MinDeviation is specified as the fraction of the requested error of the entire integral."

:Evaluate: Given::usage = "Given is an option of Divonne.
	It provides a list of points where the integrand might have peaks.
	Divonne will consider these points when partitioning the integration region."

:Evaluate: NExtra::usage = "NExtra is an option of Divonne.
	It specifies the maximum number of points that will be considered in the output of the PeakFinder function."

:Evaluate: PeakFinder::usage = "PeakFinder is an option of Divonne.
	It specifies the peak-finder function.
	This function is called whenever a region is up for subdivision and is supposed to point out possible peaks lying in the region, thus acting as the dynamic counterpart of the static list of points supplied with Given.
	It is invoked with two arguments, the multidimensional equivalents of the lower left and upper right corners of the region being investigated, and must return a (possibly empty) list of points."

:Evaluate: MinPoints::usage = "MinPoints is an option of Divonne.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Divonne.
	It can take the values Last or All which determine whether only the last (largest) or all sets of samples collected on a subregion over the integration phases contribute to the final result."

:Evaluate: PseudoRandom::usage = "PseudoRandom is an option of Divonne.
	If set to True, pseudo-random numbers are used instead of Sobol quasi-random numbers."

:Evaluate: Regions::usage = "Regions is an option of Divonne.
	It specifies whether the regions into which the integration region has been cut are returned together with the integration results."

:Evaluate: Region::usage = "Region[ll, ur, res, df] describes a subregion:
	ll and ur are multidimensional equivalents of the region's lower left and upper right corner.
	res gives the integration results for the region in a list with entries of the form {integral, error, chi-square} for each component of the integrand.
	df is the number of degrees of freedom corresponding to the chi-square values in res."

:Evaluate: Begin["`Divonne`"]

:Begin:
:Function: Divonne
:Pattern: MLDivonne[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, mineval_, maxeval_,
  key1_, key2_, key3_, maxpass_,
  border_, maxchisq_, mindeviation_,
  xgiven_, fgiven_, nextra_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, mineval, maxeval,
  key1, key2, key3, maxpass,
  border, maxchisq, mindeviation,
  xgiven, fgiven, nextra}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer, Integer,
  Integer, Integer, Integer, Integer,
  Real, Real, Real,
  RealList, RealList, Integer}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Divonne] = {HoldFirst}

:Evaluate: Options[Divonne] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000,
	Key1 -> 47, Key2 -> 1, Key3 -> 1, MaxPass -> 5,
	Border -> 0, MaxChisq -> 10, MinDeviation -> .25,
	Given -> {}, NExtra -> 0, PeakFinder -> ({}&),
	Verbose -> 1, Final -> All, PseudoRandom -> True,
	Regions -> False, Compiled -> True}

:Evaluate: Divonne[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}],
	tags, vars, lower, range, jac, tmp, defs, integrand,
	rel, abs, mineval, maxeval, key1, key2, key3, maxpass, border,
	maxchisq, mindeviation,	given, nextra, peakfinder,
	final, verbose, pseudo, regions, compiled},
	  Message[Divonne::optx, #, Divonne]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Divonne]];
	  {rel, abs, mineval, maxeval, key1, key2, key3, maxpass, border,
	    maxchisq, mindeviation, given, nextra, peakfinder,
	    verbose, final, pseudo, regions, compiled} =
	    tags /. {opt} /. Options[Divonne];
	  {vars, lower, range} = Transpose[{v}];
	  jac = Simplify[Times@@ (range -= lower)];
	  tmp = Array[tmpvar, ndim];
	  defs = Simplify[lower + range tmp];
	  Block[{Set}, define[compiled, tmp, vars, Thread[vars = defs], jac]];
	  integrand = fun[f];
	  given = Flatten[N[(# - lower)/range]&/@ given];
	  MLDivonne[ndim, ncomp[f], 10.^-rel, 10.^-abs,
	    Min[Max[verbose, 0], 3] +
	      If[final === Last, 4, 0] +
	      If[TrueQ[pseudo], 8, 0] +
	      If[TrueQ[regions], 256, 0],
	    mineval, maxeval, key1, key2, key3, maxpass,
	    N[border], N[maxchisq], N[mindeviation],
	    given, sample[given], nextra]
	]

:Evaluate: tmpvar[n_] := ToExpression["Cuba`Divonne`t" <> ToString[n]]

:Evaluate: Attributes[ncomp] = Attributes[fun] = {HoldAll}

:Evaluate: ncomp[f_List] := Length[f]

:Evaluate: _ncomp = 1

:Evaluate: define[True, tmp_, vars_, {defs__}, jac_] :=
	fun[f_] := Compile[tmp, Block[vars, defs; check[vars, Chop[f jac]//N]]]

:Evaluate: define[False, tmp_, vars_, {defs__}, jac_] :=
	fun[f_] := Function[tmp, Block[vars, defs; check[vars, Chop[f jac]//N]]]

:Evaluate: check[_, f_Real] = {f}

:Evaluate: check[_, f:{__Real}] = f

:Evaluate: check[x_, _] := (Message[Divonne::badsample, ff, x]; {})

:Evaluate: sample[x_] :=
	Check[Apply[integrand, Partition[x, ndim], 1]//Flatten, {}]

:Evaluate: findpeak[b_] := Check[Join[#, sample[#]]& @
	N[Flatten[peakfinder@@ Transpose[lower + range Partition[b, 2]]]], {}]

:Evaluate: region[ll_, ur_, r___] :=
	Region[lower + range ll, lower + range ur, r]

:Evaluate: Divonne::badsample = "`` is not a real-valued function at ``."

:Evaluate: Divonne::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Divonne::badcomp = "Cannot integrate `` components."

:Evaluate: Divonne::accuracy =
	"Desired accuracy was not reached within `` integrand evaluations on `` subregions.
	Estimate that MaxPoints needs to be increased by `` for this accuracy."

:Evaluate: Divonne::success = "Needed `` integrand evaluations on `` subregions."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Divonne.tm
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 1 Mar 06 th
*/


#include <setjmp.h>
#include "mathlink.h"
#include "util.c"

jmp_buf abort_;

/*********************************************************************/

static void Status(MLCONST char *msg, cint n1, cint n2, cint n3)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 4);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "Divonne");
  MLPutString(stdlink, msg);
  MLPutInteger(stdlink, n1);
  MLPutInteger(stdlink, n2);
  MLPutInteger(stdlink, n3);
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

static void DoSample(cnumber n, ccount ldx, real *x, real *f)
{
  int pkt;
  real *mma_f;
  long mma_n;

  if( MLAbort ) goto abort;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Cuba`Divonne`sample", 1);
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

  neval_ += n;

  Copy(f, mma_f, n*ncomp_);
  MLDisownRealList(stdlink, mma_f, mma_n);
}

/*********************************************************************/

static count SampleExtra(cBounds *b)
{
  int pkt;
  count n, nget;
  real *mma_f;
  long mma_n;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Cuba`Divonne`findpeak", 1);
  MLPutRealList(stdlink, (real *)b, 2*ndim_);
  MLEndPacket(stdlink);

  while( (pkt = MLNextPacket(stdlink)) && (pkt != RETURNPKT) )
    MLNewPacket(stdlink);

  if( !MLGetRealList(stdlink, &mma_f, &mma_n) ) {
    MLClearError(stdlink);
    MLNewPacket(stdlink);
    MLPutFunction(stdlink, "Abort", 0);
    longjmp(abort_, 1);
  }

  neval_ += nget = mma_n/(ndim_ + ncomp_);

  n = IMin(nget, nextra_);
  if( n ) {
    Copy(xextra_, mma_f, n*ndim_);
    Copy(fextra_, mma_f + nget*ndim_, n*ncomp_);
  }

  MLDisownRealList(stdlink, mma_f, mma_n);

  return n;
}

/*********************************************************************/

#include "common.c"

void Divonne(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, cint maxpass,
  creal border, creal maxchisq, creal mindeviation,
  real *xgiven, clong nxgiven, real *fgiven, clong nfgiven,
  cnumber nextra)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) ) {
    Status("badcomp", ncomp, 0, 0);
    MLPutSymbol(stdlink, "$Failed");
  }
  else if( BadDimension(ndim, flags, key1) ||
           BadDimension(ndim, flags, key2) ||
           ((key3 & -2) && BadDimension(ndim, flags, key3)) ) {
    Status("baddim", ndim, 0, 0);
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    real integral[NCOMP], error[NCOMP], prob[NCOMP];
    int fail;
    ccount nx = nxgiven + nextra*ndim;
    ccount nf = nfgiven + nextra*ncomp;

    neval_ = ngiven_ = nxgiven/ndim;
    neval_opt_ = neval_cut_ = 0;
    nextra_ = nextra;
    ldxgiven_ = ndim;

    Alloc(xgiven_, nx + nf);
    xextra_ = xgiven_ + nxgiven;
    fgiven_ = xgiven_ + nx;
    fextra_ = fgiven_ + nfgiven;

    Copy(xgiven_, xgiven, nxgiven);
    Copy(fgiven_, fgiven, nfgiven);

    border_.lower = border;
    border_.upper = 1 - border_.lower;

    fail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key1, key2, key3, maxpass,
      maxchisq, mindeviation,
      integral, error, prob);

    if( fail >= 0 ) {
      count comp;

      Status(fail ? "accuracy" : "success", neval_, nregions_, fail);

      MLPutFunction(stdlink, "List", ncomp);
      for( comp = 0; comp < ncomp; ++comp ) {
        real res[] = {integral[comp], error[comp], prob[comp]};
        MLPutRealList(stdlink, res, Elements(res));
      }
    }

    free(xgiven_);
  }

  MLEndPacket(stdlink);
}

/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

