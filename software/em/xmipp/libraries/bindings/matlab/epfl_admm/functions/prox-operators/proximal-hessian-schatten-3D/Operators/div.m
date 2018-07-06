function divF = div(F1,F2,F3,bc)

% compute appropriate shifts
F1s = shift(F1,[1,0,0],bc);
F2s = shift(F2,[0,1,0],bc);
F3s = shift(F3,[0,0,1],bc);

% compute finite differences
d1F1 = F1-F1s;
d2F2 = F2-F2s;
d3F3 = F3-F3s;

% compute the divergence
divF = d1F1+d2F2+d3F3;

end