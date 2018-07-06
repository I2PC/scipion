function X = readVolume(fnIn,headerSize,imageSize)

fh  =fopen(fnIn)       ;
fseek(fh,headerSize,'bof')         ;
X   =fread(fh,imageSize^3,'float32','ieee-le');
X   =reshape(X,[imageSize,imageSize,imageSize]);