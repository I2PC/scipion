#include <iostream>
#include "xmipp_Aux.h"

#define BGNOISE 0
#define ENVELOPE 1
#define PSD 2
#define CTF 3

void CTFProfile(CTFDescription &ctfmodel, double angle, double FMAX,
		int samples, MultidimArray<double> &profiles) {
	double step = FMAX / samples;

	profiles.resizeNoCopy(samples, 4);
	double sinus = sin(angle);
	double cosinus = cos(angle);
	double F;
	int i;

	for (i = 0, F = 0; i < YSIZE(profiles); i++, F += step) {
		double fx = F * cosinus;
		double fy = F * sinus;

		// Compute current frequencies.
		ctfmodel.precomputeValues(fx, fy);

		// Store values.
		double bgNoise = ctfmodel.CTFnoise_at();
		double ctf = ctfmodel.CTFpure_at();
		double E = ctfmodel.CTFdamping_at();

		A2D_ELEM(profiles, i, BGNOISE) = bgNoise;
		A2D_ELEM(profiles, i, ENVELOPE) = bgNoise + E * E;
		A2D_ELEM(profiles, i, PSD) = bgNoise + ctf * ctf;
		A2D_ELEM(profiles, i, CTF) = ctf;
	}
}

void CTFAverageProfile(CTFDescription &ctfmodel, double FMAX, int samples,
		MultidimArray<double> &profiles) {
	double step = FMAX / samples;

	profiles.resizeNoCopy(samples, 4);

	for (double angle = 0.0; angle < 360; angle++) {
		double sinus = sin(angle);
		double cosinus = cos(angle);
		double F;
		int i;

		for (i = 0, F = 0; i < YSIZE(profiles); i++, F += step) {
			double fx = F * cosinus;
			double fy = F * sinus;

			// Compute current frequencies.
			ctfmodel.precomputeValues(fx, fy);

			// Store values.
			double bgNoise = ctfmodel.CTFnoise_at();
			double ctf = ctfmodel.CTFpure_at();
			double E = ctfmodel.CTFdamping_at();

			profiles(i, BGNOISE) += bgNoise;
			profiles(i, ENVELOPE) += bgNoise + E * E;
			profiles(i, PSD) += bgNoise + ctf * ctf;
			profiles(i, CTF) += ctf;
		}
	}

	for (int i = 0; i < YSIZE(profiles); i++) {
		profiles(i, BGNOISE) /= 360;
		profiles(i, ENVELOPE) /= 360;
		profiles(i, PSD) /= 360;
		profiles(i, CTF) /= 360;
	}
}
