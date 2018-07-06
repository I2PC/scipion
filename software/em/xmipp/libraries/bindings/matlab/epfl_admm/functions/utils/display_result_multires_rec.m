 
% laurene.donati@epfl.ch

function display_result_multires_rec(vol,rec,sizeIndex)

slice1 = floor(sizeIndex/2); 
slice2 = floor(sizeIndex/6); 
slice3 = floor(sizeIndex/4); 

compareOrthoSlicesRecGT(vol,rec,slice1); 
compareOrthoSlicesRecGT(vol,rec,slice2); 
compareOrthoSlicesRecGT(vol,rec,slice3); 
    
figure, subplot(2,1,1), plot(fsc3D(vol,rec)), grid on; title('FSC rec'); ...
        subplot(2,1,2), plot(vol(:,slice1,slice1)); hold on;
    plot(rec(:,slice1,slice1));title('Profile lines'); legend('gt', 'rec'); grid on;
    
end 