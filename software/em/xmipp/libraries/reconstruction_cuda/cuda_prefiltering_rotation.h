
template<class floatN>
__device__ floatN InitialCausalCoefficient(
	floatN* c,			// coefficients
	uint DataLength,	// number of coefficients
	int step)			// element interleave in bytes
{
	const uint Horizon = UMIN(12, DataLength);

	// this initialization corresponds to clamping boundaries
	// accelerated loop
	float zn = Pole;
	float Sum = *c;
	for (uint n = 0; n < Horizon; n++) {
		Sum += zn * *c;
		zn *= Pole;
		c = (floatN*)((uchar*)c + step);
	}
	return(Sum);
}


