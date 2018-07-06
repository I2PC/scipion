function [aCF1,aCF2,aCF3] = curlAdj(F1,F2,F3,bc)

% compute appropriate shifts
as2F1 = shiftAdj(F1,[0,1,0],bc);
as3F1 = shiftAdj(F1,[0,0,1],bc);

as1F2 = shiftAdj(F2,[1,0,0],bc);
as3F2 = shiftAdj(F2,[0,0,1],bc);

as1F3 = shiftAdj(F3,[1,0,0],bc);
as2F3 = shiftAdj(F3,[0,1,0],bc);

% compute appropriate adjoint finite differences
ad3F2 = F2 - as3F2;
ad2F3 = F3 - as2F3;

ad1F3 = F3 - as1F3;
ad3F1 = F1 - as3F1;

ad2F1 = F1 - as2F1;
ad1F2 = F2 - as1F2;

% compute the adjoint curl components
aCF1 = -ad3F2+ad2F3;
aCF2 = -ad1F3+ad3F1;
aCF3 = -ad2F1+ad1F2;

end