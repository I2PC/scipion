function Jf=Jacobian3D(F,bc) %Jacobian operator for a 3D Vector Field with boundary conditions

%F denotes the vector field. This should be of size MxNxKx9. The first
%three indices indicate the spatial coordinates. The last dimension holds 
%the vector components for a specific coordinate (voxel).

[r,c,s,v]=size(F);

if ~isequal(v,3)
    error('Jacobian3D: The first input argument is not a valid 3D vector field');
end

Jf=zeros(r,c,s,9); %Jacobian of the vector field. For every spatial 
%coordinate we store the Jacobian which is a 3x3 matrix. The convention is
%that the last dimension holds the elements of the Jacobian in a
%column-wise fashion, i.e., 
%Jf(m,n,k,:)=[dF1/dx;dF2/dx;dF3/dx;dF1/dy;dF2/dy;dF3/dy;dF1/dz;dF2/dz;dF3/dz];

Jf(:,:,:,1)=shift(F(:,:,:,1),[-1, 0, 0],bc)-F(:,:,:,1);%dF1/dx;
Jf(:,:,:,2)=shift(F(:,:,:,2),[-1, 0, 0],bc)-F(:,:,:,2);%dF2/dx;
Jf(:,:,:,3)=shift(F(:,:,:,3),[-1, 0, 0],bc)-F(:,:,:,3);%dF3/dx;
Jf(:,:,:,4)=shift(F(:,:,:,1),[0, -1, 0],bc)-F(:,:,:,1);%dF1/dy;
Jf(:,:,:,5)=shift(F(:,:,:,2),[0, -1, 0],bc)-F(:,:,:,2);%dF2/dy;
Jf(:,:,:,6)=shift(F(:,:,:,3),[0, -1, 0],bc)-F(:,:,:,3);%dF3/dy;
Jf(:,:,:,7)=shift(F(:,:,:,1),[0, 0, -1],bc)-F(:,:,:,1);%dF1/dz;
Jf(:,:,:,8)=shift(F(:,:,:,2),[0, 0, -1],bc)-F(:,:,:,2);%dF2/dz;
Jf(:,:,:,9)=shift(F(:,:,:,3),[0, 0, -1],bc)-F(:,:,:,3);%dF3/dz;

