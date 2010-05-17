function psd = tom_calc_periodogram(image,sz)
%TOM_CALC_PERIODOGRAM calculates a periodogram to be used in the xmipp ctf
%determination functions
%
% psd = tom_calc_periodogram(image,sz)
%
%PARAMETERS
%
%  INPUT
%   image           input 2D image
%   sz              (optional) output size of the averaged periodogram (default 512)
% 
%  OUTPUT
%   psd             periodogram
%
%EXAMPLE
%   psd = tom_calc_periodogram(image,256);
%
%REFERENCES
%
%
%
%SEE ALSO
%   
%   TOM_XMIPP_ADJUST_CTF, TOM_XMIPP_ENHANCE_PSD
%
%   created by AK & COS 27/08/07
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if nargin < 2
    sz = 512;
end

image = single(image);

imsz = size(image);
shift = sz/2;

psd = zeros(sz,'single');

counter = 0;

for i=1:shift:imsz(1)-sz;
    for j=1:shift:imsz(2)-sz;
        psd = psd + abs(fft2(image(i:i+sz-1,j:j+sz-1))).^2;
        counter = counter + 1;
    end
end

psd = (psd./counter)./sz.^4;
