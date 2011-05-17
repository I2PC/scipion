#ifndef XMIPP_AUX
#define XMIPP_AUX
#include <data/ctf.h>
#include <data/multidim_array.h>

void CTFProfile(CTFDescription &ctfmodel, double angle, double FMAX, int samples, MultidimArray<double> &profiles);
void CTFAverageProfile(CTFDescription &ctfmodel, double FMAX, int samples, MultidimArray<double> &profiles);

#endif
