function Hf=VHessian3D(F,bc)
%F denotes the vector field. This should be of size MxNxKx3. The first
%three indices indicate the spatial coordinates. The last dimension holds the
%vector components for a specific coordinate.

[r,c,s,v]=size(F);

if ~isequal(v,3)
    error('VHessian3D: The first input argument is not a valid 3D vector field');
end

Hf=zeros(r,c,s,18); %Hessian of the vector field. For every spatial 
%coordinate we store the vector-valued Hessian which is a 9x3 matrix, i.e.,
%Hf=[Hf1;Hf2;Hf3]; Since Hfk is symmetric we can just save the upper triangular
%part of every Hfk. 
%The convention is that the last dimension holds the elements of the 
%Hessian in a column-wise fashion, i.e., 
%Hf(m,n,:)=[d^2F1/dxx;d^2F1/dxy;d^2F1/dxz;d^2F1/dyy;d^2F1/dyz;d^2F1/dzz];
%          [d^2F2/dxx;d^2F2/dxy;d^2F2/dxz;d^2F2/dyy;d^2F2/dyz;d^2F2/dzz];
%          [d^2F3/dxx;d^2F3/dxy;d^2F3/dxz;d^2F3/dyy;d^2F3/dyz;d^2F3/dzz];

for k=1:3
  %d^2F_k/dxx
  Hf(:,:,:,(k-1)*6+1)=(F(:,:,:,k)-2*shift(F(:,:,:,k),[-1,0,0],bc)+shift(F(:,:,:,k),[-2,0,0],bc));
  %d^2F_k/dxy
  Hf(:,:,:,(k-1)*6+2)=(F(:,:,:,k)-shift(F(:,:,:,k),[0,-1,0],bc)-shift(F(:,:,:,k),[-1,0,0],bc)+shift(F(:,:,:,k),[-1,-1,0],bc));
  %d^2F_k/dxz
  Hf(:,:,:,(k-1)*6+3)=(F(:,:,:,k)-shift(F(:,:,:,k),[0,0,-1],bc)-shift(F(:,:,:,k),[-1,0,0],bc)+shift(F(:,:,:,k),[-1,0,-1],bc));
  %d^2F_k/dyy
  Hf(:,:,:,(k-1)*6+4)=(F(:,:,:,k)-2*shift(F(:,:,:,k),[0,-1,0],bc)+shift(F(:,:,:,k),[0,-2,0],bc));  
  %d^2F_k/dyz
  Hf(:,:,:,(k-1)*6+5)=(F(:,:,:,k)-shift(F(:,:,:,k),[0,0,-1],bc)-shift(F(:,:,:,k),[0,-1,0],bc)+shift(F(:,:,:,k),[0,-1,-1],bc));
  %d^2F_k/dzz
  Hf(:,:,:,(k-1)*6+6)=(F(:,:,:,k)-2*shift(F(:,:,:,k),[0,0,-1],bc)+shift(F(:,:,:,k),[0,0,-2],bc));  
end




