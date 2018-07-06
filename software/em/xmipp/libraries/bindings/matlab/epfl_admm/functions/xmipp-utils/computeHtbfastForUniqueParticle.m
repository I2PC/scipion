
% laurene.donati@epfl.ch

function Htbp = computeHtbfastForUniqueParticle(sizeHtb, proj, tau, padding, fkbwf, r1p, r2p, K1, K2, K3)
    
    % Compute measurement on refined grid
    proj      = refine2Dimage(proj, tau);
    proj      = padarray(proj,round(size(proj)*padding),0,'both');
    
    % Compute the convolution as point-wise multiplication in Fourier
    conv_p      = ifftshift(ifft2(fkbwf.*fft2(proj)))*(tau^2);
    
    % Set meshgrid in projection domain 
    centerConv  = (size(conv_p,1)-1)/2;
    xValConv    = tau*(-centerConv:centerConv);
    [Y1,Y2]     = meshgrid(xValConv, xValConv);
    
    % Fill Htb
    kInner1     = K1*r1p(1)+K2*r1p(2)+K3*r1p(3); % inner product <K,r1> - relates to y1
    kInner2     = K1*r2p(1)+K2*r2p(2)+K3*r2p(3); % inner product <K,r2> - relates to y2
    contrib_p   = interp2(Y1, Y2, conv_p, kInner1(:), kInner2(:), 'spline');
    Htbp        = reshape(contrib_p,sizeHtb);
    
end