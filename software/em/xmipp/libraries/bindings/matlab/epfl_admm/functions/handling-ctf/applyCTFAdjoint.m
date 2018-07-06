
% laurene
% applies the adoint of the ctf to the measurement
% for version of Htb with mex file
% NOTE: ONLY WORKS FOR ISOTROPIC CTF RIGHT NOW

function projs = applyCTFAdjoint(projs,ctf)

if ctf.use == true
    image = ctf.image;
    parfor i = 1: size(projs,3)
        projs(:,:,i)    =  real(ifft2(fft2(squeeze(projs(:,:,i))).*ifftshift(image)));
    end
end


end