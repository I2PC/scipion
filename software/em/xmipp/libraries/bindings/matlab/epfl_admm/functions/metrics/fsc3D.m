%------------------------- Function    -----------------------------------%
%         y = fsc(x1,x2)
% This function computes Fourier Shell Correlation. This is a measure to
% estimate the resolution of the two maps, x1 and x2. The refference is
% (FSC; Harauz & van Heel, 1986).
%-------------------------- INPUT     ------------------------------------%
% x1 and x2 : input maps
%-------------------------- OUTPUT    ------------------------------------%
% y : Fourier shell corelation between x1 and x2
% the value of nyquist frequency is assumed to be 0.5


function y = fsc3D(x1,x2,varargin)

[step,doPlot]  =process_options(varargin,'step',.01,'doPlot',false);

% compute fft
fftx1 = fftshift(fftn(x1)) ;
fftx2 = fftshift(fftn(x2)) ;

% set coordinates
[X,Y,Z]  = ndgrid(linspace(-.5,.5,size(x1,1)),linspace(-.5,.5,size(x1,2)),linspace(-.5,.5,size(x1,3)));

% set frequencies interval
freqInt= 0:step:.5;

if freqInt(end)~= .5
    freqInt = [freqInt,.5];
end

% set radius
W = X.^2+Y.^2+Z.^2;

y = zeros((length(freqInt)-1),1);
for i = 1 : (length(freqInt)-1)
    indx = (W>=((freqInt(i))^2) & W<((freqInt(i+1))^2));
    fftx1Indx = fftx1(indx);
    fftx2Indx = fftx2(indx);
    
    y(i) = abs(sum(fftx1Indx(:).*conj(fftx2Indx(:)))/...
        sqrt(sum(abs(fftx1Indx(:)).^2)*sum(abs(fftx2Indx(:)).^2)));
end

if (doPlot)
    figure
    plot(freqInt(2:end),y,'-r','LineWidth',2)
    xlim([0,.5])
    ylim([min(y),1])
end
