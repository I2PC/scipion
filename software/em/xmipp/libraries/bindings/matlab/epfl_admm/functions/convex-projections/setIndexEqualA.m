% this function is a projection operator. It projects on specific region to
% be zero

function x = setIndexEqualA(x,x1,x2,a)

x(1:x1,:)   = a;
x(x2:end,:) = a;
