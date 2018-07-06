%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct on Fine Grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y   =   expandScaledKBWFcoeffstoFineGrid_3D(Ck,scale,param,output_size)

% The way that we fill the coarse grid to fine grid is as follows:
%    fix the boundaries, and refine the gaps in between

disp('-- Reexpanding reconstruction in fine grid.');

% Set KBWF parameters 
a     = param.a; 
alpha = param.alpha; 
m     = param.m; 


%% Basis Function in Fine Grid 

% Compute fine grid (2 times finer that needed)
x   =  linspace(0,a,a*scale+1);
x   =  x(:);
x   =  [x(end:-1:2);x(:)];

% Kaiser bessel function
KB  =  KaiserBesselWindowFunction3D(x,alpha,a,m);

% Center of the basis function 
L   =  (length(x)+1)/2    ;


%% Represent Coarse Grid in Finer Grid Dimension 

%  Coordinates for the coarse grid (its center is the center of the coefficients)
x1  = ((1:size(Ck,1))-(size(Ck,1)+1)/2)*scale;
x2  = ((1:size(Ck,2))-(size(Ck,2)+1)/2)*scale;
x3  = ((1:size(Ck,3))-(size(Ck,3)+1)/2)*scale;

% Fine grid coordinate + padding
padd_size = 2*L;
y1  = (x1(1)-scale-padd_size):(x1(end)+scale+padd_size);
y2  = (x2(1)-scale-padd_size):(x2(end)+scale+padd_size);
y3  = (x3(1)-scale-padd_size):(x3(end)+scale+padd_size);

% Filling the rest by zeros 
Aux =  zeros([length(y1),length(y2),length(y3)]) ;

% Place the coefficent at the right spots
% The rest is zeros 
Aux((scale+padd_size+1):scale:(end-scale-padd_size),(scale+padd_size+1):scale:(end-scale-padd_size),(scale+padd_size+1):scale:(end-scale-padd_size)) = Ck;
Ck                     =  Aux;


%% Zero-padding the kernel 
kernel = zeros(size(Ck)) ;
kernel((ceil((size(kernel,1)+1)/2)-L+1):(ceil((size(kernel,1)+1)/2)+L-1),...
    (ceil((size(kernel,2)+1)/2)-L+1):(ceil((size(kernel,2)+1)/2)+L-1),...
    (ceil((size(kernel,3)+1)/2)-L+1):(ceil((size(kernel,3)+1)/2)+L-1)) = KB;


%% Convolution with the kernel with FFT
Y = real(ifftn(real(fftn(ifftshift(kernel))).*fftn(Ck)));


%% Output coordinate
z1  = (1:output_size(1))-(output_size(1)+1)/2;
z2  = (1:output_size(2))-(output_size(2)+1)/2;
z3  = (1:output_size(3))-(output_size(3)+1)/2;

[Z1,Z2,Z3] = ndgrid(z1,z2,z3);

%% Interpolate on the output grid
Y    = interp3(y1,y2,y3,Y,Z2,Z1,Z3);