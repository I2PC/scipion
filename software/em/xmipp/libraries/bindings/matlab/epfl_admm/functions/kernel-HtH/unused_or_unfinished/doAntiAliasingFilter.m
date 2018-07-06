%======================================================================
function y = doAntiAliasingFilter(x, d)

% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the antialiasing filter to use on the projections


order = max(64,2^nextpow2(2*length(x)));



% First create a ramp filter - go up to the next highest
% power of 2.

filt = ones(1,(order/2)+1);
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

filt(w>pi*d) = 0;                    % Crop the frequency response
filt = [filt';  filt(end-1:-1:2)'];  % Symmetry of the filter


y    = ifft(fft(x,order).*filt)   ;
y    = y(1:length(x))             ;