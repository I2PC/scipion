% this function is a projection operator. It projects on positive region.

function x = positiveProj(x)

x(x<=0) = 0;