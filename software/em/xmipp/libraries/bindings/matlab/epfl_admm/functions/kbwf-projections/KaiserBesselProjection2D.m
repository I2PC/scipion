
% Compute the 2D projection of a KBWF basis function
function X_KBWF = KaiserBesselProjection2D(param,  x_value)

% Set kbwf parameters
a        =   param.a;
alpha    =   param.alpha;
m        =   param.m;

% Initialize
[X,Y]             =   meshgrid(x_value,x_value);
X_KBWF            =   zeros(length(x_value),length(x_value));

% Compute distance to the center
s                 =   sqrt(X.^2+Y.^2);

% Fill the 2D projection inside the support
temp              =  alpha*sqrt(1-(abs(s)/a).^2);
X_KBWF(abs(s)<a)  =  (a/besseli(m,alpha))*sqrt(2*pi/alpha)*...
                        ((temp(abs(s)<a)/alpha).^(m+1/2)).*...
                            besseli(m+1/2,temp(abs(s)<a));
end