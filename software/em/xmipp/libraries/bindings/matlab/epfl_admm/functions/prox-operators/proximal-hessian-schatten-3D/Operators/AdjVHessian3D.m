function F=AdjVHessian3D(Hf,bc)
%Hf denotes the Hessian of a 3D vector field. This should be of size MxNxKx18. 
%The first three indices indicate the spatial coordinates. The last dimension 
%holds the Hessian components for a specific coordinate. 
%The convention is that the last dimension holds the elements of the 
%Hessian in a column-wise fashion, i.e., 
%Hf(m,n,:)=[d^2F1/dxx;d^2F1/dxy;d^2F1/dxz;d^2F1/dyy;d^2F1/dyz;d^2F1/dzz];
%          [d^2F2/dxx;d^2F2/dxy;d^2F2/dxz;d^2F2/dyy;d^2F2/dyz;d^2F2/dzz];
%          [d^2F3/dxx;d^2F3/dxy;d^2F3/dxz;d^2F3/dyy;d^2F3/dyz;d^2F3/dzz];

[r,c,s,v]=size(Hf);

if ~isequal(v,18)
    error('AdjVHessian2D: The first input argument is not a valid Hessian of a 2D vector field');
end

F=zeros(r,c,s,3); %A 3D vector field. For every spatial coordinate we store 
%the three vector components of the vector field F.

for i=1:3
  Pxx=Hf(:,:,:,(i-1)*6+1);
  Pxx=(Pxx-2*shiftAdj(Pxx,[-1,0,0],bc)+shiftAdj(Pxx,[-2,0,0],bc));
  Pxy=Hf(:,:,:,(i-1)*6+2);
  Pxy=(Pxy-shiftAdj(Pxy,[0,-1,0],bc)-shiftAdj(Pxy,[-1,0,0],bc)+shiftAdj(Pxy,[-1,-1,0],bc));
  Pxz=Hf(:,:,:,(i-1)*6+3);
  Pxz=(Pxz-shiftAdj(Pxz,[0,0,-1],bc)-shiftAdj(Pxz,[-1,0,0],bc)+shiftAdj(Pxz,[-1,0,-1],bc));
  Pyy=Hf(:,:,:,(i-1)*6+4);
  Pyy=(Pyy-2*shiftAdj(Pyy,[0,-1,0],bc)+shiftAdj(Pyy,[0,-2,0],bc));
  Pyz=Hf(:,:,:,(i-1)*6+5);
  Pyz=(Pyz-shiftAdj(Pyz,[0,0,-1],bc)-shiftAdj(Pyz,[0,-1,0],bc)+shiftAdj(Pyz,[0,-1,-1],bc));
  Pzz=Hf(:,:,:,(i-1)*6+6);
  Pzz=(Pzz-2*shiftAdj(Pzz,[0,0,-1],bc)+shiftAdj(Pzz,[0,0,-2],bc));
  F(:,:,:,i)=Pxx+Pyy+Pzz+2*(Pxy+Pxz+Pyz);
end