% =====================        Function        ========================== %
% Description y = LT3D_FGP(p,q,r)
% this function computes the adjoint of the finite difference in 3D
% y = (p1jk pijk-pi-1jk -pm-1jk)+ similar equations for q and r
%
% =====================         INPUT          ========================== %
% p,q,r : finite difference along different dimension
% ======================== OPTIONAL INPUT PARAMETERS ======================
%
% =====================         OUTPUT         ========================== %
% y : real value three dimensional data
% =====================         EXAMPLE        ========================== %

function   y = LT3D_FGP(p,q,r)

% this function is the linear operator defined in 
%     Fast gradient based algorithm 
%             written by 
%    Amir beck and Marc Teboulle
% L(p,q)_i,j=p(i,j)+q(i,j)-p(i-1,j)-q(i,j-1)

p1 = p - circshift(p,[1,0,0])                       ;
p  = cat(1,p(1,:,:),p1(2:(end-1),:,:),-p(end-1,:,:));


q1 = q - circshift(q,[0,1,0])                       ;
q  = cat(2,q(:,1,:),q1(:,2:(end-1),:),-q(:,end-1,:));

r1 = r - circshift(r,[0,0,1])                       ;
r  = cat(3,r(:,:,1),r1(:,:,2:(end-1)),-r(:,:,end-1));

y = p+q+r;
