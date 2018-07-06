function F=AdjJacobian2D(Jf,bc) %Adjoint Jacobian operator for a 2D Vector Field with boundary conditions

%Jf denotes the Jacobian of a vector field. This should be of size MxNx4. 
%The first two indices indicate the spatial coordinates. The last dimension 
%holds the Jacobian components for a specific coordinate. The convention is
%that the last dimension holds the elements of the Jacobian in a
%column-wise fashion, i.e., Jf(m,n,:)=[dF1/dx;dF2/dx;dF1/dy;dF2/dy];

[r,c,v]=size(Jf);

if ~isequal(v,4)
    error('AdjJacobian2D: The first input argument is not a valid Jacobian of a 2D vector field');
end

F=zeros(r,c,2); %A 2D vector field. For every spatial coordinate we store 
%the two vector components of the vector field F.


P1x=Jf(:,:,1);
P1y=Jf(:,:,3);

P2x=Jf(:,:,2);
P2y=Jf(:,:,4);


F(:,:,1)=(shiftAdj(P1x,[-1,0],bc)-P1x)+(shiftAdj(P1y,[0,-1],bc)-P1y);
F(:,:,2)=(shiftAdj(P2x,[-1,0],bc)-P2x)+(shiftAdj(P2y,[0,-1],bc)-P2y);

