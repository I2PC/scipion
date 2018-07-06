function Hf=VHessian2D(F,bc)
%F denotes the vector field. This should be of size MxNx2. The first two
%indices indicate the spatial coordinates. The last dimension holds the
%vector components for a specific coordinate.

[r,c,v]=size(F);

if ~isequal(v,2)
    error('VHessian2D: The first input argument is not a valid 2D vector field');
end

Hf=zeros(r,c,6); %Hessian of the vector field. For every spatial 
%coordinate we store the Vector-valued Hessian which is a 4x2 matrix, i.e.,
%Hf=[Hf1;Hf2]; Since Hfk is symmetric we can just save the upper triangular
%part of every Hfk. 
%The convention is that the last dimension holds the elements of the 
%Hessian in a column-wise fashion, i.e., 
%Hf(m,n,:)=[d^2F1/dxx;d^2F1/dxy;d^2F1/dyy;d^2F2/dxx;d^2F2/dxy;d^2F2/dyy];

for k=1:2
  %d^2F_k/dxx
  Hf(:,:,(k-1)*3+1)=(F(:,:,k)-2*shift(F(:,:,k),[-1,0],bc)+shift(F(:,:,k),[-2,0],bc));
  %d^2F_k/dxy
  Hf(:,:,(k-1)*3+2)=(F(:,:,k)-shift(F(:,:,k),[0,-1],bc)-shift(F(:,:,k),[-1,0],bc)+shift(F(:,:,k),[-1,-1],bc));
  %d^2F_k/dyy
  Hf(:,:,(k-1)*3+3)=(F(:,:,k)-2*shift(F(:,:,k),[0,-1],bc)+shift(F(:,:,k),[0,-2],bc));  
end




