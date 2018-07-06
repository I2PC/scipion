% =====================        Function        ========================== %
% Description y = l1Norm(x)
%This function computes l1 norm of x. it is sum(abs(x))
%
% =====================         INPUT          ========================== %
% x: any dimension real value
% ======================== OPTIONAL INPUT PARAMETERS ======================
%
% =====================         OUTPUT         ========================== %
% y: a scalar
%
% =====================         EXAMPLE        ========================== %
% x=rand(100,100); y=l1Norm(x)

function y = l1l2Norm3D(x)
validateattributes(x,{'numeric'},{'real','nonempty','finite'},'l1l2Norm',...
    'x',1);
x = sqrt(sum(x.^2,4));
y = sum(abs(x(:)));
