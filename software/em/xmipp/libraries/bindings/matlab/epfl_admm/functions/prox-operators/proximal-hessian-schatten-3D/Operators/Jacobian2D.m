function Jf=Jacobian2D(F,bc) %Jacobian operator for a 2D Vector Field with boundary conditions

%F denotes the vector field. This should be of size MxNx2. The first two
%indices indicate the spatial coordinates. The last dimension holds the
%vector components for a specific coordinate.

[r,c,v]=size(F);

if ~isequal(v,2)
    error('Jacobian2D: The first input argument is not a valid 2D vector field');
end

Jf=zeros(r,c,4); %Jacobian of the vector field. For every spatial 
%coordinate we store the Jacobian which is a 2x2 matrix. The convention is
%that the last dimension holds the elements of the Jacobian in a
%column-wise fashion, i.e., Jf(m,n,:)=[dF1/dx;dF2/dx;dF1/dy;dF2/dy];

Jf(:,:,1)=shift(F(:,:,1),[-1,0],bc)-F(:,:,1);%dF1/dx
Jf(:,:,2)=shift(F(:,:,2),[-1,0],bc)-F(:,:,2);%dF2/dx
Jf(:,:,3)=shift(F(:,:,1),[0,-1],bc)-F(:,:,1);%dF1/dy
Jf(:,:,4)=shift(F(:,:,2),[0,-1],bc)-F(:,:,2);%dF2/dy


