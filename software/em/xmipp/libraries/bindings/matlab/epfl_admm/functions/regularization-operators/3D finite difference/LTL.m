% this function do LTL using L and LT in 2 dimensional

function y = LTL(x,L,LT)

if (isvector(x))
    error('input should be two dimension or three dimension matrix');
else
    dimx = length(size(x));
end

if (dimx == 2)
    
    % compute Lx
    Lx = L(x);
    
    % compute LTLx
    y  = LT(reshape(Lx(:,1),size(x)),reshape(Lx(:,2),size(x)));
    
    
elseif (dimx == 3)
    
    % compute Lx
    Lx = L(x);
    
    % compute LTLx
    y  = LT(reshape(Lx(:,1),size(x)),reshape(Lx(:,2),size(x)),...
        reshape(Lx(:,3),size(x)));
end