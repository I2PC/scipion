
% laurene.donati@epfl.ch

function compareOrthoSlicesRecGT(vol,rec,slice)

figure, subplot(2,3,1), imagesc(vol(:,:,slice)), colormap(gray), colorbar, title(['gt - xy' ' - slice ' num2str(slice)]);
    subplot(2,3,2), imagesc(squeeze(vol(:,slice,:))), colormap(gray), colorbar, title(['gt - xz' ' - slice ' num2str(slice)]);
    subplot(2,3,3), imagesc(squeeze(vol(slice,:,:))), colormap(gray), colorbar, title(['gt - yz' ' - slice ' num2str(slice)]);
    subplot(2,3,4), imagesc(rec(:,:,slice)), colormap(gray), colorbar, title(['rec -xy' ' - slice ' num2str(slice)]);
    subplot(2,3,5), imagesc(squeeze(rec(:,slice,:))), colormap(gray), colorbar, title(['rec -xy' ' - slice ' num2str(slice) ]);
    subplot(2,3,6), imagesc(squeeze(rec(slice,:,:))), colormap(gray), colorbar, title(['rec -xy' ' - slice ' num2str(slice)]);
    
end 