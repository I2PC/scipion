function [CF1,CF2,CF3] = curl(F1,F2,F3,bc)

% compute appropriate shifts
s2F1 = shift(F1,[0,1,0],bc);
s3F1 = shift(F1,[0,0,1],bc);

s1F2 = shift(F2,[1,0,0],bc);
s3F2 = shift(F2,[0,0,1],bc);

s1F3 = shift(F3,[1,0,0],bc);
s2F3 = shift(F3,[0,1,0],bc);

% compute appropriate finite differences
d3F2 = F2 - s3F2;
d2F3 = F3 - s2F3;

d1F3 = F3 - s1F3;
d3F1 = F1 - s3F1;

d2F1 = F1 - s2F1;
d1F2 = F2 - s1F2;

% compute the curl components
CF1 = d3F2-d2F3;
CF2 = d1F3-d3F1;
CF3 = d2F1-d1F2;

end

