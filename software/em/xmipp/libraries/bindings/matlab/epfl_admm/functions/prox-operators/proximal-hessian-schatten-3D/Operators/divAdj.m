function [gradF1,gradF2,gradF3] = divAdj(F,bc)

% compute appropriate shifts
Fs1 = shiftAdj(F,[1,0,0],bc);
Fs2 = shiftAdj(F,[0,1,0],bc);
Fs3 = shiftAdj(F,[0,0,1],bc);

% compute adjoint finite differences
gradF1 = F - Fs1;
gradF2 = F - Fs2;
gradF3 = F - Fs3;

end