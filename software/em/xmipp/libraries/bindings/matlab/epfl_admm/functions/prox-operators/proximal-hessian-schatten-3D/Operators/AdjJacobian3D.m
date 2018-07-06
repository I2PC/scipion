function F=AdjJacobian3D(Jf,bc) %Adjoint Jacobian operator for a 2D Vector Field with boundary conditions

%Jf denotes the Jacobian of a vector field. This should be of size MxNxKx9. 
%The first three indices indicate the spatial coordinates. The last dimension 
%holds the Jacobian components for a specific coordinate. The convention is
%that the last dimension holds the elements of the Jacobian in a
%column-wise fashion, i.e., 
%Jf(m,n,k,:)=[dF1/dx;dF2/dx;dF3/dx;dF1/dy;dF2/dy;dF3/dy;dF1/dz;dF2/dz;dF3/dz];

[r,c,s,v]=size(Jf);

if ~isequal(v,9)
    error('AdjJacobian3D: The first input argument is not a valid Jacobian of a 3D vector field');
end

F=zeros(r,c,s,3); %A 3D vector field. For every spatial coordinate we store 
%the three vector components of the vector field F.


P1x=Jf(:,:,:,1);
P1y=Jf(:,:,:,4);
P1z=Jf(:,:,:,7);

P2x=Jf(:,:,:,2);
P2y=Jf(:,:,:,5);
P2z=Jf(:,:,:,8);

P3x=Jf(:,:,:,3);
P3y=Jf(:,:,:,6);
P3z=Jf(:,:,:,9);


F(:,:,:,1)=(shiftAdj(P1x,[-1,0,0],bc)-P1x)+(shiftAdj(P1y,[0,-1,0],bc)-P1y)...
    +(shiftAdj(P1z,[0 0 -1],bc)-P1z);
F(:,:,:,2)=(shiftAdj(P2x,[-1,0,0],bc)-P2x)+(shiftAdj(P2y,[0,-1,0],bc)-P2y)...
    +(shiftAdj(P2z,[0 0 -1],bc)-P2z);
F(:,:,:,3)=(shiftAdj(P3x,[-1,0,0],bc)-P3x)+(shiftAdj(P3y,[0,-1,0],bc)-P3y)...
    +(shiftAdj(P3z,[0 0 -1],bc)-P3z);

