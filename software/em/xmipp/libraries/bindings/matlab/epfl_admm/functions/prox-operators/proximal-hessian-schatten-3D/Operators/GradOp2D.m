function Df=GradOp2D(f,bc) %Gradient operator with boundary conditions bc

[r,c]=size(f);
Df=zeros(r,c,2);
Df(:,:,1)=shift(f,[-1,0],bc)-f; %f(i+1,j)-f(i,j)
Df(:,:,2)=shift(f,[0,-1],bc)-f; %f(i,j+1)-f(i,j)