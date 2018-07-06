%_________________________________________________________________________%
%_______________________________  Apply_ADMM _____________________________%
%
% This method apply Alternating direction method of multipliers to solve
%
%     argmin { 1/2|| AX-b ||^2 + lambda/2 ||x||^2 + alpha \Psi(Lx)}
%
% By using variable spliting the minimization problem change to
% argmin  1/2|| AX-b ||^2 + lambda/2 ||x||^2 + alpha \psi(u) s.t.  LX=u
% By using Augmented lagragian method solving this equivalent problem,
% converge to the same result.
%
%     argmin { 1/2|| AX-b ||^2 + lambda/2 ||x||^2 + alpha \psi(u) + mu/2
%            ||LX-u+d||^2 }
%
% u0 and d0 are (N*N)*2 matrix, where N is the size of the input of
% operator A and ATA and LTL are the operators that their input can be a
% vector or any dimensional matrix
% and their output has the same shape as the input.
% L is the difference operator, its input
% is a any dimensional matrix and the output is also a matrix where its number
% of columns is equal to the dimension of its input, and its number of rows
% is equal to the number of elements in input matrix. LT is the adjoint of
% finite difference, the input is two or three
% vectors depending on the dimension of input and the output is a matrix.
% You can choose symetric or asymetric TV regularization term you are
% interested in. It could be 'symetric' or 'asymetric'.
%
%
%-------------------------------   INPUTS   ------------------------------%
% ATA : AT(A())
% ATb : AT(b)

% ======================== OPTIONAL INPUT PARAMETERS ======================
% A : Forward Operator
% AT: adjoint of Forward Operator
% L : Discrete derivative operator
% LT: Adjoint of the discrete derivative operator
% ATb : AT(b)
% lambda, alpha and mu are the parameters of Augmented function
% N and M : Image size
% ========================
% b : Measurement vector (DEFAULT : 0)
% M_Pre : Preconditioner (DEFAULT : IDENTITY)
% bound : size of the zero padding area (DEFAULT : 0)
% u0 and d0: initial value for the algorithm (DEFAULT : 0)
% MaxOutIter : the number of outer iterations (DEFAULT : 30)
% Max_InIter : the number of inner CG itterations (DEFAULT : 3)
% Oracle     : oracle image (DEFAULT : 0)
% tol   : tolerance as stopping criteria. The stopping criteria consider (DEFAULT : 1e-20)
% how slow the cost function decreases in comparison with the previous one
% BsplineDegree : Since the model is based on the B-spline functions, this
% shows the degree of the B-spline which has been used (DEFAULT : 3)
% TV_type : symetric or asymetric (DEFAULT : symetric)
% verbose : to type the variation of SNR (DEFAULT : false)
% drawflag : To plot the variation of SNR with respect to the iterations (DEFAULT: false)
%
%
%-------------------------------- OUTPUTS
%-----------------------------
% y    : the soltion of the inverse problem which is the outcome of CG
% technique
% resvec:the behaviour of the inner CG step
% Cost_Fun Total Cost function
% SNR   : SNR variation with respect to the number of iteration
%__________________________________________________________________________
% =========================================================================
%
% Author: masih.nilchian@epfl.ch
%
% =========================================================================
function [y,varargout] = applyADMM(ATb,varargin)

validateattributes(ATb,{'numeric'},{'real','nonempty','finite'},'deblurMFISTA',...
    'ATb',1);

if (isvector(ATb))
    error('input should be two dimension or three dimension matrix');
else
    dimx = length(size(ATb));
end
% set finite difference as default value

if (dimx == 2)
    LFGP = @(x) L2D_OneOutput(x);
    LTFGP= @(p,q) LT2D_FGP(p,q) ;
elseif(dimx==3)
    LFGP = @(x) L3DOneOutput_FGP(x);
    LTFGP= @(p,q,r) LT3D_FGP(p,q,r);
else
    error('input should be two or three dimension vector')
end


% set the necessary parameters
Identity = @(X) X;
proxSoftThresholding = @(x,lambda) softThresholding(x,lambda);

[ATA,L,LT,prox,potentialFun,lambda,alpha,mu,M,bound,x0,u0,d0,maxInIter,maxOutIter,tol,...
    verbose,drawflag,doProjection]=process_options(varargin,...
    'ATA',Identity,'L',LFGP,'LT',LTFGP,'prox',proxSoftThresholding,'potentialFun',@l1Norm,...
    'lambda',0,'alpha',1e-6*norm(ATb(:),2),'mu',1,'M',Identity,'bound',0,...
    'x0',zeros(size(ATb)),'u0',zeros(numel(ATb),dimx),'d0',zeros(numel(ATb),dimx),...
    'maxInIter',3,'maxOutIter',30,'tol',1e-10,...
    'verbose',false,'drawflag',false,'doProjection',Identity);

% start verbose
if verbose
    
    fprintf('******************************************\n');
    fprintf('**  ADMM with Conjugate Gradient method **\n');
    fprintf('******************************************\n');
    fprintf('#iter   cost function   \n')
    fprintf('==========================================\n');
end

% initialize
u_k     =  u0;
d_k     =  d0;
costFun = zeros(maxOutIter,1)   ;


% compute necessary operators
ATAmuLTLlambda = @(X) ATA(X)+mu*LTL(X,L,LT)+lambda*X;


for k  =  1  :  maxOutIter
    
    % start computing cost
    if (verbose || ~(nargout==1) || drawflag)
        
        ATAlambdaIx0    =  ATA(x0)+lambda*x0;
        costFun(k)      =  sum(x0(:).*(.5*ATAlambdaIx0(:)-ATb(:)))+...
            alpha*potentialFun(L(x0));
        
        if verbose
            fprintf('%3d \t %10.5f  \n',k,costFun(k));
        end
        if drawflag
            figure(20);plot(1:k,costFun(1:k),'-*r'); drawnow;
            figure(30); imshow(squeeze(x0(:,:,27)),[]); drawnow;
        end
    end
    
    % Solving the quadratic part by using Conjugate gradient algorithm
    % argmin  1/2|| AX-b ||^2 + mu/2 ||LX-u+d||^2 + lambda/2 ||x||^2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%     PCG     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (dimx == 2)
        
        u1    =  reshape(squeeze(u_k(:,1)),size(ATb));
        u2    =  reshape(squeeze(u_k(:,2)),size(ATb));
        
        d1    =  reshape(squeeze(d_k(:,1)),size(ATb));
        d2    =  reshape(squeeze(d_k(:,2)),size(ATb));
        
        % apply preconditioned CG
        x_k  = applyPCG(ATb+mu*LT(u1-d1,u2-d2),'ATA',ATAmuLTLlambda...
            ,'x0',x0,'tol',tol,'numberIteration',maxInIter,'M',M);
        
    elseif (dimx == 3)
        u1    =  reshape(squeeze(u_k(:,1)),size(ATb));
        u2    =  reshape(squeeze(u_k(:,2)),size(ATb));
        u3    =  reshape(squeeze(u_k(:,3)),size(ATb));
        
        d1    =  reshape(squeeze(d_k(:,1)),size(ATb));
        d2    =  reshape(squeeze(d_k(:,2)),size(ATb));
        d3    =  reshape(squeeze(d_k(:,3)),size(ATb));
        
        % apply preconditioned CG
        
        x_k  = applyPCG(ATb+mu*LT(u1-d1,u2-d2,u3-d3),'ATA',ATAmuLTLlambda...
            ,'x0',x0,'tol',tol,'numberIteration',maxInIter,'M',M);
    end
    % do convex projection
    x_k = doProjection(x_k);
    
    % warm initializing
    x0      =  x_k         ;
    
 
    
    % compute prox
    u_k     =  prox(L(x_k)+d_k,alpha/mu);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update Lagrange parameter
    d_k     =  d_k + L(x_k) - u_k       ;
    
   
    
end
if nargout
    varargout{1} = costFun ;
end

y   =  x_k;
end

