
function mat2tiffstack(filename, vol)

for ii = 1 : size(vol, 3)
    imwrite(vol(:,:,ii) , filename , 'WriteMode' , 'append') ;
end

end