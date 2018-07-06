function F=AdjVHessian2D(Hf,bc)
%Hf denotes the Hessian of a vector field. This should be of size MxNx6. 
%The first two indices indicate the spatial coordinates. The last dimension 
%holds the Hessian components for a specific coordinate. 
%The convention is that the last dimension holds the elements of the 
%Hessian in a column-wise fashion, i.e., 
%Hf(m,n,:)=[d^2F1/dxx;d^2F1/dxy;d^2F1/dyy;d^2F2/dxx;d^2F2/dxy;d^2F2/dyy];

[r,c,v]=size(Hf);

if ~isequal(v,6)
    error('AdjVHessian2D: The first input argument is not a valid Hessian of a 2D vector field');
end

F=zeros(r,c,2); %A 2D vector field. For every spatial coordinate we store 
%the two vector components of the vector field F.

for i=1:2
  Pxx=Hf(:,:,(i-1)*3+1);
  Pxx=(Pxx-2*shiftAdj(Pxx,[-1,0],bc)+shiftAdj(Pxx,[-2,0],bc));
  Pxy=Hf(:,:,(i-1)*3+2);
  Pxy=(Pxy-shiftAdj(Pxy,[0,-1],bc)-shiftAdj(Pxy,[-1,0],bc)+...
  shiftAdj(Pxy,[-1,-1],bc));
  Pyy=Hf(:,:,(i-1)*3+3);
  Pyy=(Pyy-2*shiftAdj(Pyy,[0,-1],bc)+shiftAdj(Pyy,[0,-2],bc));
  F(:,:,i)=Pxx+2*Pxy+Pyy;
end