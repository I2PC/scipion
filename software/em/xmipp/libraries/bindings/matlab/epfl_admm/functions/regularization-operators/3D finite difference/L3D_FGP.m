% =====================        Function        ========================== %
% Description [p,q,r] = L3D_FGP(X)
% this function computes the finite difference in 3D
% p(i,j,k) = x(i,j,k)-x(i+1,j,k)
% q(i,j,k) = x(i,j,k)-x(i,j+1,k)
% r(i,j,k) = x(i,j,k)-x(i,j,k+1)
%
% =====================         INPUT          ========================== %
% x : real value three dimensional data
% ======================== OPTIONAL INPUT PARAMETERS ======================
%
% =====================         OUTPUT         ========================== %
% p,q,r : finite difference along different dimension
%
% =====================         EXAMPLE        ========================== %

function [p,q,r] = L3D_FGP(X)

% this function is the transpose of the L_FGP function
% it is related to the paper of Amir Beck and Teboulle

p  =  X - circshift(X,[-1,0,0]);
p  =  p(1:(size(X,1)-1),:,:)   ;
p(size(X,1),:,:) = 0;

q  =  X - circshift(X,[0,-1,0]);
q  =  q(:,1:(size(X,2)-1),:)   ;
q(:,size(X,2),:)=0;

r  =  X - circshift(X,[0,0,-1]) ;
r  =  r(:,:,1:(size(X,3)-1))   ;
r(:,:,size(X,3))=0;