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

function y = l1Norm(x)
validateattributes(x,{'numeric'},{'real','nonempty','finite'},'l1Norm',...
    'x',1);
x = x(:);
y = sum(abs(x));
