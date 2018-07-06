fh=fopen('BRDnoisyRec.raw','w'); fwrite(fh,X(:),'float32'); fclose(fh);

fh=fopen('TVrec_Coef.raw','w'); X=fread(fh,141*141*141,'float32'); fclose(fh);