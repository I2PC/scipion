 % laurene.donati@epfl.ch
 
 
TESTLPF       = false;
if TESTLPF
    % Test: Low pass after
    [a,b]    = butter(3, (1/3));
    lpf_rec  = filter(a,b,rec,[],1);
    lpf_rec  = filter(a,b,lpf_rec,[],2);
    lpf_rec  = filter(a,b,lpf_rec,[],3);
end

    
    if TESTLPF
        figure, imagesc(lpf_rec(:,:,sizeIndex/2)), colormap(gray), colorbar, title('lpf_rec');
        figure, plot(fsc3D(vol,lpf_rec)), title('FSC lpf_rec');
    end