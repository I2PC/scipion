function hf = helicity(F,bc)

[cf1,cf2,cf3] = curl(F(:,:,:,1),F(:,:,:,2),F(:,:,:,3),bc);
hf = F(:,:,:,1).*cf1 + F(:,:,:,2).*cf2 + F(:,:,:,3).*cf3; 

end