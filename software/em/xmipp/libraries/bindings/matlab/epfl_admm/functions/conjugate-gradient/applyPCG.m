%_________________________________________________________________________%
%_______________________________  Apply_PCG   ____________________________%
%
% This method apply conjugate gradient technique to solve the inverse
% problem Ax=b. To that end, it is required to have ATA and ATb to apply
% this technique.
% The cost function which is minimized is .5*||Ax-b||^2 or (ATAx-ATb)
%
%-------------------------------   INPUTS   ------------------------------%
% The main Input is ATb.
% ATA and ATb are the required inputs. AT is the adjoint of the forward
% operator A. Therefore, ATA is the sequentional application of the
% forward operator and its adjoint. ATA is a function handle.
% The rests are not necessary but they put more options for the user to add
% which are presented in varargin which are as followings:
%
% ======================== OPTIONAL INPUT PARAMETERS ======================
% ATA  : application of AT and A. (DEFAULT : I)
% x0   : initual value
% tol  : a stopping criteria
% numberIteration: maximum number of iteration
% M    : preconditioner if exists
% verbose: to show the output
% drawflag : plot the cost function
%-------------------------------- OUTPUTS    -----------------------------%
% X    : the soltion of the inverse problem which is the outcome of CG
% technique
% costFun : the value of the Cost Function .5||Ax-b||^2 with respect to
% the iteration number
%__________________________________________________________________________
% =========================================================================
%
% Author: masih.nilchianm@epfl.ch
%
% =========================================================================
function  [X,varargout] = applyPCG(ATb,varargin)

Identity = @(X) X;
[ATA,x0,tol,numberIteration,M,verbose,drawflag]=...
    process_options(varargin,'ATA',Identity,'x0',zeros(size(ATb)),'tol',1e-10,...
    'numberIteration',30,'M',Identity,'verbose',false,...
    'drawflag',false);

if verbose
    
    fprintf('******************************************\n');
    fprintf('**  Conjugate Gradient Method  **\n');
    fprintf('******************************************\n');
    fprintf('#iter   cost function   \t  SNR\n')
    fprintf('==========================================\n');
end


ATAx0  =  ATA(x0)           ;
r0     =  ATb-ATAx0         ;      % initial gradient
z0     =  M(r0)             ;      % apply preconditioner
p0     =  z0                ;      % initial search direction
d0     =  real(sum(r0(:).*z0(:)));

costFun=  zeros(1,numberIteration);


if (verbose || (~(nargout==1)) || drawflag)
    
    costFun(1) =  sum(x0(:).*(.5*ATAx0(:)-ATb(:)))  ;
    if verbose
        fprintf('%3d \t %10.5f \n',1,costFun(1));
    end
    if drawflag
        figure(10);set(gcf,'Color',[1 1 1]);
        subplot(2,1,2)
        plot(1,costFun(1),'-*r'); title('CG cost function value');xlabel('iteration number'); drawnow;
        xlim auto
    end
end
for i = 2 : numberIteration
    
    h0    =  ATA(p0)                  ;
    alpha =  d0/real(sum(p0(:).*h0(:)));
    r1    =  r0 - alpha*h0      ;
    
%     if (norm(r1(:),2)>norm(r0(:),2))
%         break;
%     end
    x1    =  x0 + alpha*p0      ;
    z1    =  M(r1)              ;
    d1    =  real(sum(r1(:).*z1(:)))  ;
    beta  =  d1/d0             ;
    p1    =  z1 + beta*p0      ;
    
    
    p0  =  p1                ;
    d0  =  d1                ;
    r0  =  r1                ;
    x0  =  x1                ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%              Compute Cost function                %%%%%%%%%%%%%
    if (verbose || (~(nargout==1)) || drawflag)
        ATAx0  =  ATA(x0)           ;
        costFun(i) =  sum(x0(:).*(.5*ATAx0(:)-ATb(:)))  ;
        if (i>6)
            if ((abs(costFun(i)-costFun(i-1))<=tol*abs(costFun(i-1)-costFun(i-2)))&&...
                    (abs(costFun(i-1)-costFun(i-2))<=tol*abs(costFun(i-2)-costFun(i-3))))
                break;
            end
        end
        if drawflag
            figure(10);set(gcf,'Color',[1 1 1]);
            subplot(2,1,2);plot(1:i,costFun(1:i),'-r');title('CG cost function value');xlabel('iteration number'); drawnow;
            xlim auto
            subplot(2,1,1);imshow(x0,[]); title('Reconstructed phantom');drawnow;
        end
    end
    if verbose
        fprintf('%3d \t %10.5f  \n',i,costFun(i));
    end
end
if ~(nargout==1)
    varargout{1}=costFun;
end
X   =  x0                    ;