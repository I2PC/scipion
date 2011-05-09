#include <iostream>
#include "xmipp_Aux.h"

#define NOISE 0
#define NOISE_ERROR2 1
#define NOISE_CTF2 2
#define NOISE_PURE 3

void CTFProfile(CTFDescription &ctfmodel, double angle, double FMAX, int samples, MultidimArray<double> &profiles) {
	double step = FMAX / samples;

	profiles.resizeNoCopy(samples, 4);
	double F = 0;
	double sinus = sin(angle);
	double cosinus = cos(angle);

	for(int i = 0; i < YSIZE(profiles); i++) {
		double fx = F * cosinus;
		double fy = F * sinus;

		// Compute current frequencies.
		ctfmodel.precomputeValues(fx, fy);

		// Store values.
		double noise_at = ctfmodel.CTFnoise_at();
		double pure_at = ctfmodel.CTFpure_at();
		double E = ctfmodel.CTFdamping_at();
		double ctf = pure_at;

		profiles(i, NOISE) = noise_at;
		profiles(i, NOISE_ERROR2) = noise_at + E * E;
		profiles(i, NOISE_CTF2) = noise_at + ctf * ctf;
		profiles(i, NOISE_PURE) = pure_at;

		F += step;
	}
}
