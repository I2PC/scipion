function xmipp_show_structure_factor(rundir)
%
%   xmipp_show_structure_factor(rundir)

figure()
[f2,logF]=xmipp_read_structure_factor(rundir);
plot(f2,logF)
xlabel('Frequency (1/A^2)')
ylabel('Log(StructureFactor)')
hold on

disp('Identify a LEFT position to fit damping factor')
[x1,y1]=ginput(1);
plot([x1 x1],[min(logF) max(logF)],'g','LineWidth',2)

disp('Identify a RIGHT position to fit damping factor')
[x2,y2]=ginput(1);
plot([x2 x2],[min(logF) max(logF)],'g','LineWidth',2)

% Compute the regression
idx=find(f2>x1 & f2<x2);
[P,S] = polyfit(f2(idx),logF(idx),1);
logFp=polyval(P,f2(idx));

plot(f2(idx),logFp,'r')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
f0=f2(idx(1));
fF=f2(idx(end));
disp(['Model for the interval f^2 (1/A^2)=[' num2str(f0) ',' num2str(fF) ']  f(A)=['...
    num2str(1/sqrt(fF)) ',' num2str(1/sqrt(f0)) ']'])
disp(['B2=' num2str(P(1))])
disp(['Fitting line : logF=' num2str(P(1)) '*f^2+(' num2str(P(2))+')'])
R2=1-var(logF(idx)-logFp)/var(logF(idx));
p=length(P);
R2adj=R2-(1-R2)*p/(length(idx)-p-1);
disp(['R^2 adjusted = ' num2str(R2adj)])
