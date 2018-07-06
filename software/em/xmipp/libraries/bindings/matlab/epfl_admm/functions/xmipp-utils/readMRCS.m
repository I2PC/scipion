fh=fopen('projections.mrcs');
fseek(fh,1024,'bof');

X=fread(fh,100*71*71,'float32','ieee-le');
I=squeeze(Y(:,:,50));
imshow(I',[])
fclose(fh)

fh=fopen('1BRD.vol');
fseek(fh,1136,'bof');
X=fread(fh,71*71*71,'float32','ieee-le');
X=reshape(X,[71,71,71]);
fclose(fh)
imshow(squeeze(X(:,:,35)),[])

vol = zeros(size(X));
for i = 1: size(vol,3)
    xx = squeeze(X(:,:,i));
    vol(:,:,i)= xx;
end


fh=fopen('mask.vol');
fseek(fh,1136,'bof');
X=fread(fh,71*71*71,'float32','ieee-le');
X=reshape(X,[71,71,71]);
fclose(fh)
