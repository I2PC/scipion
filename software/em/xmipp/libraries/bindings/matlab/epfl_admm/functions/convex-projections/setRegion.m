% this function is a projection operator. It projects on specific region [a,b].

function x = setRegion(x,a,b)

x(x<=a) = a;
x(x>=b) = b;