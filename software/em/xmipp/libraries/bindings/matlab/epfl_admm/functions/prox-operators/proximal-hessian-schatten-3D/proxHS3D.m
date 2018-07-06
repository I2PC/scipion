function [x,P,iter,L]=proxHS3D(y,lambda,varargin)

%3D Image denoising with a Hessian Schatten-Norm regularizer under convex
%constraints
%
%  argmin  0.5*||y-x||^2+ lambda*|| Hx||_(1,p) + i_C(x), 
%   
% i_C: indicator function of the closed convex set C.
%  

% ========================== INPUT PARAMETERS (required) ==================
% Parameters    Values description
% =========================================================================
% y             Noisy image.
% lambda        Regularization penalty parameter.
% ======================== OPTIONAL INPUT PARAMETERS ======================
% Parameters    Values' description
%
% img           Original Image. (For the compution of the ISNR improvement)
% maxiter       Number of iterations (Default: 100)
% tol           Stopping threshold for denoising (Default:1e-4)
% optim         The type of gradient-based method used {'fgp'|'gp'}. 
%               (Fast Gradient Projection or Gradient Projection) 
%               (Default: 'fgp')
% P             Initialization of the dual variables. (Default: zeros([size(y) 6])).
% verbose       If verbose is set on then info for each iteration are
%               printed on screen. (Default: false)
% project       A function handle for the projection onto the convex set C.
%               (Default: project=@(x)x, which means that there is no
%               constrain on x.)
% bc            Boundary conditions for the differential operators.
%               {'reflexive'|'circular'|'zero'} (Default: 'reflexive')
% snorm         Specifies the type of the Hessian Schatten norm.
%               {'spectral'|'nuclear'|'frobenius'|'Sp'}. (Default:
%               'frobenius'). If snorm is set to Sp then the order of
%               the norm has also to be specified.
% order         The order of the Sp-norm. 1<order<inf. For order=1 set
%               snorm='nuclear', for order=2 set snorm='frobenius' and 
%               for order=inf set snorm='spectral'.
% =========================================================================
% ========================== OUTPUT PARAMETERS ============================
% x             Denoised image.
% P             The solution of the dual problem.
% iter          Number of iterations until convergence of the algorithm.
% L             Lipschitz constant of the dual objective.
% =========================================================================
%
% Author: stamatis.lefkimmiatis@epfl.ch
%
% =========================================================================

[maxiter,L,tol,optim,verbose,img,project,P,bc,snorm,order]=process_options(...
  varargin,'maxiter',100,'L',144,'tol',1e-4,'optim','fgp','verbose',...
  false,'img',[],'project',@(x)x,'P',zeros([size(y) 6]),'bc',...
  'reflexive','snorm','frobenius','order',[]);

if isempty(L)
  L=Lipschitz(y)/1.25;%Lipschitz constant
end

if isequal(snorm,'Sp') && isempty(order)
  error('proxHS3D: The order of the Sp-norm must be specified.');
end

if isequal(snorm,'Sp') && isinf(order)
  error('proxHS3D: Try spectral norm for the type-norm instead.');
end

if isequal(snorm,'Sp') && order==1
  error('proxHS3D: Try nuclear norm for the type-norm instead.');
end

if isequal(snorm,'Sp') && order < 1
  error('proxHS3D: The order of the Sp-norm should be greater or equal to 1.');
end


count=0;
flag=false;

if verbose
  fprintf('**********************************************\n');
  fprintf('** Denoising with Hessian-Based Regularizer **\n');
  fprintf('**********************************************\n');
  fprintf('#iter     relative-dif   \t fun_val         Duality Gap    \t   ISNR\n')
  fprintf('====================================================================\n');
end
switch optim
  case 'fgp'
    t=1;
    F=P;
    for i=1:maxiter
      K=y-lambda*AdjHessianOp3D(F,bc);
      %Pnew=F-step*lambda*HessianOp(lambda*AdjHessianOp(F)-y);
      %step=1/(L*lambda^2)==>
      %Pnew=F-(1/(L*lambda))*HessianOp(lambda*AdjHessianOp(F)-y);
      Pnew=F+(1/(L*lambda))*HessianOp3D(project(K),bc);
      Pnew=projectLB(Pnew,snorm,order);
      
      re=norm(Pnew(:)-P(:))/norm(Pnew(:));%relative error
      if (re<tol)
        count=count+1;
      else
        count=0;
      end
      
      tnew=(1+sqrt(1+4*t^2))/2;
      F=Pnew+(t-1)/tnew*(Pnew-P);
      P=Pnew;
      t=tnew;
      
      if verbose
        if ~isempty(img)
          k=y-lambda*AdjHessianOp3D(P,bc);
          x=project(k);
          fun_val=cost(y,x,lambda,bc,snorm,order);
          dual_fun_val=dualcost(y,k,project);
          dual_gap=(fun_val-dual_fun_val);
          ISNR=20*log10(norm(y(:)-img(:))/norm(x(:)-img(:)));
          % printing the information of the current iteration
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f \t %2.8f\n',i,re,fun_val,dual_gap,ISNR);
        else
          k=y-lambda*AdjHessianOp3D(P,bc);
          x=project(k);
          fun_val=cost(y,x,lambda,bc,snorm,order);
          dual_fun_val=dualcost(y,k,project);
          dual_gap=(fun_val-dual_fun_val);
          % printing the information of the current iteration
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f\n',i,re,fun_val,dual_gap);
        end
      end
      
      if count >=5
        flag=true;
        iter=i;
        break;
      end
    end
    
  case 'gp'
    
    for i=1:maxiter
      K=y-lambda*AdjHessianOp3D(P,bc);
      %Pnew=P-step*lambda*HessianOp(lambda*AdjHessianOp(P)-y);
      %step=1/(L*lambda^2)==>
      %Pnew=P-(1/(L*lambda))*HessianOp(lambda*AdjHessianOp(P)-y);
      Pnew=P+(1/(L*lambda))*HessianOp3D(project(K),bc);
      Pnew=projectLB(Pnew,snorm,order);
      
      re=norm(Pnew(:)-P(:))/norm(Pnew(:));%relative error
      if (re<tol)
        count=count+1;
      else
        count=0;
      end
      
      P=Pnew;
      
      if verbose
        if ~isempty(img)
          k=y-lambda*AdjHessianOp3D(P,bc);
          x=project(k);
          fun_val=cost(y,x,lambda,bc,snorm,order);
          dual_fun_val=dualcost(y,k,project);
          dual_gap=(fun_val-dual_fun_val);
          ISNR=20*log10(norm(y(:)-img(:))/norm(x(:)-img(:)));
          % printing the information of the current iteration
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f \t %2.8f\n',i,re,fun_val,dual_gap,ISNR);
        else
          k=y-lambda*AdjHessianOp3D(P,bc);
          x=project(k);
          fun_val=cost(y,x,lambda,bc,snorm,order);
          dual_fun_val=dualcost(y,k,project);
          dual_gap=(fun_val-dual_fun_val);
          % printing the information of the current iteration
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f\n',i,re,fun_val,dual_gap);
        end
      end
      
      if count >=5
        flag=true;
        iter=i;
        break;
      end
    end
end

if ~flag
  iter=maxiter;
end

x=project(y-lambda*AdjHessianOp3D(P,bc));


function Hf=HessianOp3D(f,bc)

[x,y,z]=size(f);

fxx=(f-2*shift(f,[-1,0,0],bc)+shift(f,[-2,0,0],bc));
fyy=(f-2*shift(f,[0,-1,0],bc)+shift(f,[0,-2,0],bc));
fzz=(f-2*shift(f,[0,0,-1],bc)+shift(f,[0,0,-2],bc));

fxy=(f-shift(f,[0,-1,0],bc)-shift(f,[-1,0,0],bc)+shift(f,[-1,-1,0],bc));
fxz=(f-shift(f,[0,0,-1],bc)-shift(f,[-1,0,0],bc)+shift(f,[-1,0,-1],bc));
fyz=(f-shift(f,[0,-1,0],bc)-shift(f,[0,0,-1],bc)+shift(f,[0,-1,-1],bc));

%Compute Hf (Apply to image f the Hessian Operator)
%Hf will be a cube (3D image)
Hf=zeros(x,y,z,6);
Hf(:,:,:,1)=fxx;
Hf(:,:,:,2)=fxy;
Hf(:,:,:,3)=fxz;
Hf(:,:,:,4)=fyy;
Hf(:,:,:,5)=fyz;
Hf(:,:,:,6)=fzz;

function Hf=HessianOp3D_(f,bc)

[x,y,z]=size(f);

fxx=(f-2*shift(f,[-1,0,0],bc)+shift(f,[-2,0,0],bc));
fyy=(f-2*shift(f,[0,-1,0],bc)+shift(f,[0,-2,0],bc));
fzz=(f-2*shift(f,[0,0,-1],bc)+shift(f,[0,0,-2],bc));

fxy=(f-shift(f,[0,-1,0],bc)-shift(f,[-1,0,0],bc)+shift(f,[-1,-1,0],bc));
fxz=(f-shift(f,[0,0,-1],bc)-shift(f,[-1,0,0],bc)+shift(f,[-1,0,-1],bc));
fyz=(f-shift(f,[0,-1,0],bc)-shift(f,[0,0,-1],bc)+shift(f,[0,-1,-1],bc));

%Compute Hf (Apply to image f the Hessian Operator)
%Hf will be a cube (3D image)
Hf=zeros(6,x,y,z);
Hf(1,:,:,:)=fxx;
Hf(2,:,:,:)=fxy;
Hf(3,:,:,:)=fxz;
Hf(4,:,:,:)=fyy;
Hf(5,:,:,:)=fyz;
Hf(6,:,:,:)=fzz;


function HaA=AdjHessianOp3D(A,bc)
%A is 3x3 symmetric so we dont have to store the redundant values for A.
%Thus, we store only 6 out of the 9 values.

Axx=A(:,:,:,1);
Axx=(Axx-2*shiftAdj(Axx,[-1,0,0],bc)+shiftAdj(Axx,[-2,0,0],bc));
Ayy=A(:,:,:,4);
Ayy=(Ayy-2*shiftAdj(Ayy,[0,-1,0],bc)+shiftAdj(Ayy,[0,-2,0],bc));
Azz=A(:,:,:,6);
Azz=(Azz-2*shiftAdj(Azz,[0,0,-1],bc)+shiftAdj(Azz,[0,0,-2],bc));

A2xy=2*A(:,:,:,2);
A2xy=(A2xy-shiftAdj(A2xy,[0,-1,0],bc)-shiftAdj(A2xy,[-1,0,0],bc)+...
  shiftAdj(A2xy,[-1,-1,0],bc));
A2xz=2*A(:,:,:,3);
A2xz=(A2xz-shiftAdj(A2xz,[0,0,-1],bc)-shiftAdj(A2xz,[-1,0,0],bc)+...
  shiftAdj(A2xz,[-1,0,-1],bc));
A2yz=2*A(:,:,:,5);
A2yz=(A2yz-shiftAdj(A2yz,[0,-1,0],bc)-shiftAdj(A2yz,[0,0,-1],bc)+...
  shiftAdj(A2yz,[0,-1,-1],bc));

%Compute H*A (Apply to cube A the adjoint of the Hessian Operator)
%H*A will be an image
HaA=Axx+Ayy+Azz+A2xy+A2xz+A2yz;

function Ap=projectLB(A,snorm,order)

if nargin < 3
  order=[];
end

switch snorm
  case 'spectral'
    Ap=projectSpMat3x3(A,1,1);
        
  case 'frobenius'
    Ap=projectSpMat3x3(A,2,1);
        
  case 'nuclear'
    Ap=projectSpMat3x3(A,inf,1);
        
  case 'Sp'
    Ap=projectSpMat3x3(A,order/(order-1),1);
        
  otherwise
    error('proxHS3D::Unknown type of norm.');
end


% function Px=project(x,bounds,range)
% lb=bounds(1);%lower box bound
% ub=bounds(2);%upper box bound
% 
% x(range)=0;
% 
% if isequal(lb,-Inf) && isequal(ub,Inf)
%   Px=x;
% elseif isequal(lb,-Inf) && isfinite(ub)
%   x(x>ub)=ub;
%   Px=x;
% elseif isequal(ub,Inf) && isfinite(lb)
%   x(x<lb)=lb;
%   Px=x;
% else
%   x(x<lb)=lb;
%   x(x>ub)=ub;
%   Px=x;
% end


function [Q,Hnorm]=cost(y,f,lambda,bc,snorm,p)

if ~isequal(snorm,'frobenius')
  e=eigenRoots(HessianOp3D_(f,bc));
end

switch snorm
  case 'spectral'
    %Sum of the Hessian spectral radius
    k=max(abs(e),[],1);
    Hnorm=sum(k(:));
  case 'frobenius'
    Hf=HessianOp3D(f,bc);
    Hnorm=sqrt(Hf(:,:,:,1).^2+Hf(:,:,:,4).^2+Hf(:,:,:,6).^2+...
      2*(Hf(:,:,:,2).^2+Hf(:,:,:,3).^2+Hf(:,:,:,5).^2));
    Hnorm=sum(Hnorm(:));
  case 'nuclear'
    k=sum(abs(e),1);
    Hnorm=sum(k(:));
  case 'Sp'
    k=sum(abs(e).^p,1).^(1/p);
    Hnorm=sum(k(:));
  otherwise
    error('proxHS3D::Unknown type of specified norm.');
end

Q=0.5*norm(y(:)-f(:))^2+lambda*Hnorm;


% function [Q,Hnorm]=cost(y,f,lambda,bc,snorm,order)
% 
% if nargin < 6
%   order=[];
% end
% 
% [e1,e2,e3]=eig3x3(HessianOp3D(f,bc));
% 
% switch snorm
%   case 'spectral'
%     %Sum of the Hessian spectral radius
%     Hnorm=sum(max(max(abs(e1(:)),abs(e2(:))),abs(e3(:))));
%   case 'frobenius'
%     Hnorm=sum(sqrt(e1(:).^2+e2(:).^2+e3(:).^2));
%   case 'nuclear'
%     Hnorm=sum(abs(e1(:))+abs(e2(:))+abs(e3(:)));
%   case 'Sp'
%     Hnorm=sum((abs(e1(:)).^order+abs(e2(:)).^order+abs(e3(:)).^order).^(1/order));
%   otherwise
%     error('denoiseHS3D::Unknown type of specified norm.');
% end
% 
% Q=0.5*norm(y(:)-f(:))^2+lambda*Hnorm;

function Q=dualcost(y,f,project)
r=f-project(f);
Q=0.5*(sum(r(:).^2)+sum(y(:).^2)-sum(f(:).^2));

function L=Lipschitz(f)

%Finds the Lipschitz constant for the function f(A)=0.5*||H*(A)-f||^2,
%where H* is the adjoint Hessian operator and A is a 3D image, e.g. H(A),
%where H is the Hessian operator.

[x,y,z]=size(f);

%The Lipschitz constant for f(A) is equal to ||HH*||=||H*H||=||H||^2=r(H*H)
%where r is the spectral radius of the operator.
%H*H=Dxx*Dxx+2Dxy*Dxy+Dyy*Dyy, which is a circulant operator since each one
%of the suboperators are circulant. The addition and multiplication of
%circulant operators results to a circulant operator.

hxx=[1 -2 1 0 0]';% Dxx Operator
hxx=padarray(hxx,[0 2 2]);
hyy=[1 -2 1 0 0]; %Dyy Operator
hyy=padarray(hyy,[2 0 2]);
hzz=[1 -2 1 0 0]; %Dzz Operator
hzz=reshape(hzz,[1 1 5]);
hzz=padarray(hzz,[2 2 0]);
center1=[3 3 3];

hxy=[1 -1 0;-1 1 0;0 0 0];%Dxy Operator
hxy=padarray(hxy,[0 0 1]);
hxz=zeros(3,3,3);%Dxz Operator
hxz(2,2,2)=1;hxz(2,2,1)=-1;hxz(1,2,2)=-1;hxz(1,2,1)=1;
hyz=zeros(3,3,3);%Dyz Operator
hyz(2,2,2)=1;hyz(2,2,1)=-1;hyz(2,1,2)=-1;hyz(2,1,1)=1;
center2=[2 2 2];

hxx(x,y,z)=0;hyy(x,y,z)=0;hzz(x,y,z)=0;
hxy(x,y,z)=0;hxz(x,y,z)=0;hyz(x,y,z)=0;

%Circular shift of the operator T1 so as to have consistent results by
%applying either one of the following 2 operations on an image x (T1zp is
%the zero padded operator)
%1)imfilter(x,T1,'conv','circular')
%2)ifft2(fft2(x).*fft2(T1zp));
hxx=circshift(hxx,1-center1);
hyy=circshift(hyy,1-center1);
hzz=circshift(hzz,1-center1);
hxy=circshift(hxy,1-center2);
hxz=circshift(hxz,1-center2);
hyz=circshift(hyz,1-center2);

%Operator eigenvalues;
Op_eig=abs(fftn(hxx)).^2+abs(fftn(hyy)).^2+abs(fftn(hzz)).^2+...
  2*abs(fftn(hxy)).^2+2*abs(fftn(hxz)).^2+2*abs(fftn(hyz)).^2;
L=max(Op_eig(:));

% function [e1,e2,e3]=eig3x3(A)
% 
% %Compute the eigenvalues of the symmetric 3x3 Hessian matrix for each pixel
% %of the image stack.
% m=(A(:,:,:,1)+A(:,:,:,4)+A(:,:,:,6))/3;
% 
% q=m.^3-1/6*(A(:,:,:,1).^2.*(A(:,:,:,4)+A(:,:,:,6))+A(:,:,:,4).^2.*...
%   (A(:,:,:,1)+A(:,:,:,6))+A(:,:,:,6).^2.*(A(:,:,:,1)+A(:,:,:,4)))+...
%   A(:,:,:,1)/6.*(A(:,:,:,2).^2+A(:,:,:,3).^2-2*A(:,:,:,5).^2)+...
%   A(:,:,:,4)/6.*(A(:,:,:,2).^2+A(:,:,:,5).^2-2*A(:,:,:,3).^2)+...
%   A(:,:,:,6)/6.*(A(:,:,:,3).^2+A(:,:,:,5).^2-2*A(:,:,:,2).^2)+...
%   A(:,:,:,2).*A(:,:,:,3).*A(:,:,:,5);
% 
% p=m.^2+(A(:,:,:,2).^2+A(:,:,:,3).^2+A(:,:,:,5).^2)/3-...
%   (A(:,:,:,1).*A(:,:,:,4)+A(:,:,:,1).*A(:,:,:,6)+A(:,:,:,4).*A(:,:,:,6))/3;
% 
% phi = 1/3*acos(q./p.^(3/2));
% 
% % NOTE: the following formula assumes accurate computation and therefore
% % q/p^(3/2) should be in the range of [-1,1]. In practice though, due to
% % finite computation accuracy there might be errors, so it must be checked.
% % Therefore, in case abs(q) >= abs(p^(3/2)), we have to set phi = 0;
% 
% phi(abs(q) >= abs(p.^(3/2)))=0;
% phi(phi<0)=phi(phi<0)+pi/3;
% 
% e1=m+2*sqrt(p).*cos(phi);
% e2=m-sqrt(p).*(cos(phi)+sqrt(3)*sin(phi));
% e3=m-sqrt(p).*(cos(phi)-sqrt(3)*sin(phi));


