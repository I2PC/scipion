function s = readmetadata(fnmetadata)
%   Copyright 2015 Jose Luis Vilas, jlvilas@cnb.csic.es.
%   $Date: 2015/08/12$
%clear all, close all, clc;
% clear all
% fnmetadata = 'candida.xmd';
% fnmetadata = 'candida.xmd';

% fnmetadata = 'candida.xmd';

% parts=regexp('a   0 0 0 0 3',' ','split');
% m=length(parts);
% field=parts{1};
% if m>1
%     value=parts(2:m);
%     value=cellfun(@(x) str2double(x),value);
%     ind=~isnan(value);
%     if any(ind)
%        value=value(ind); 
%     end
% end
% fnmetadata = 'Untilt_phantom.doc'

fid = fopen(fnmetadata);

tline = fgets(fid);

count = 0;
count2 = 1;
count3 = 0;
flagg=0;
while ischar(tline)
    tline2 =tline;
    tline2(ismember(tline2,' ')) = [];
    tline = fgets(fid);
    if (tline2(1)~='_')
        count3 = count3 + 1;
        if flagg==1
            break
        end
    else
        count = count + 1;
        labels{count} = strcat(tline2(2:end));
        flagg=1;
    end
end
fclose(fid);


fid = fopen(fnmetadata);

metadatastruct = struct();

for i=1:(count+count3)
    tline = fgets(fid);
end
tempfile = fopen('temp.txt','w');
while ischar(tline)
    tline2 =tline;
    if (tline2(1)~='#')|(tline2(1)~='_')
        fprintf(tempfile,'%s',tline);
        count2 = count2+1;
    end
    tline = fgets(fid);
end
fclose(tempfile);
fclose(fid);

A = textread('temp.txt', '%s', -1);
M = length(labels);
LA = length(A);

for k=1:length(labels)
    for j=0:M:(LA-M)
       %s.(labels{k}){ceil(j/M)}=A{j+k-1};
       s.(labels{k}){(j/M)+1}=A{k+j};
    end
end


for k=1:M
    if isnan(str2double(s.(labels{k})));
    	continue
    else
        s.(labels{k}) = str2double(s.(labels{k}));
    end
end

delete('temp.txt')

end






