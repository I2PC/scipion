% =====================        Function        ========================== %
% Description y = l2Norm(x)
%This function computes l2 norm of x. it is sum((x.^2))
%
% =====================         INPUT          ========================== %
% x: any dimension real value
% ======================== OPTIONAL INPUT PARAMETERS ======================
%
% =====================         OUTPUT         ========================== %
% y: a scalar
%
% =====================         EXAMPLE        ========================== %
% x=rand(100,100); y=l2Norm(x)

function y = l2Norm(x)
validateattributes(x,{'numeric'},{'real','nonempty','finite'},'l1Norm',...
    'x',1);
x = x(:);
y = sqrt(sum(x.^2));
