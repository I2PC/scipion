%% display slices of 3D volumes 

function display3Dslices(vol, name)

    sliceYZ = round(size(vol,1)/2);
    sliceXZ = round(size(vol,2)/2);
    sliceXY = round(size(vol,3)/2);
    
   figure('pos',[10 10 1000 400]) 
    subplot 131, imshow(squeeze(vol(sliceYZ,:,:)),[]),colorbar, title([name ' in YZ (slice ' num2str(sliceYZ) ')'])
    subplot 132, imshow(squeeze(vol(:,sliceXZ,:)),[]),colorbar, title([name ' in XZ (slice ' num2str(sliceXZ) ')'])
    subplot 133, imshow(squeeze(vol(:,:,sliceXY)),[]),colorbar, title([name ' in XY (slice ' num2str(sliceXY) ')'])

end 
