function [newim, strain, localrot] = xmipp_calculate_local_deformation ( fn1, fn2, fnmask, fnroot)

    refim = xmipp_read(fn1);
    im = xmipp_read(fn2);
    mask = xmipp_read(fnmask);
    fh=fopen([fnroot '_final.raw'],'w'); fwrite(fh,refim(:),'float'); fclose(fh);
    fh=fopen([fnroot '_initial.raw'],'w'); fwrite(fh,im(:),'float'); fclose(fh);

    main.similarity='ssd';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS  
    main.subdivide=3;       % use 3 hierarchical levels
    main.okno=5;           % mesh window size
    main.lambda = 0.005;     % transformation regularization weight, 0 for none
    main.single=1;
    main.alpha=1;           
    main.ro=0.6;

    % Optimization settings
    optim.maxsteps = 200;    % maximum number of iterations at each hierarchical level
    optim.fundif = 1e-6;     % tolerance (stopping criterion)
    optim.gamma = 1;         % initial optimization step size 
    optim.anneal=0.8;        % annealing rate on the optimization step    

    m=min([min(im(:)) min(refim(:))]);
    M=max([max(im(:)) max(refim(:))]);
    im=(im-m)/(M-m);
    refim=(refim-m)/(M-m);

    [res, newim]=mirt3D_register(refim, im, main, optim);

    % res is a structure of resulting transformation parameters
    % newim is a deformed image 
    %
    % you can also apply the resulting transformation directly as
    % newim=mirt3D_transform(im, res);

    %figure,imagesc(refim(:,:,50)); title('Reference (fixed) image slice');
    %figure,imagesc(im(:,:,50));    title('Source (float) image slice');
    %figure,imagesc(newim(:,:,50)); title('Registered (deformed) image slice');

    % Prepare some variables to estimate the strain
    Xknots = res.Xgrid(1,:,1,1);
    Yknots = res.Xgrid(:,1,1,2);
    Yknots = Yknots';
    Zbknots = res.Xgrid(1,1,:,3);
    Zknots=reshape(Zbknots,[1 size(Zbknots,3)]);
    h = res.okno;
    c = res.X;

    dim=size(im);
    strain=zeros(dim);
    localrot=zeros(dim);
    for x=1:dim(1)
        for y=1:dim(2)
            for z=1:dim(3)
                if mask(x,y,z)>0 
                    [ux_x, ux_y,ux_z, uy_x, uy_y, uy_z, uz_x, uz_y, uz_z] = get3DMIRTDeformationDiff(z,y,x);
                    U=[ux_x ux_y ux_z; uy_x uy_y uy_z ; uz_x uz_y uz_z];
                    D=1/2*(U+U');
                    H=1/2*(U-U');
                    strain(x,y,z)=abs(det(D));            
                    localrot(x,y,z)=max(imag(eigs(H)));
                end
             end
        end
    end

    newim(isnan(newim))=0;
    fh=fopen([fnroot '_initialDeformedToFinal.raw'],'w'); fwrite(fh,newim(:),'float'); fclose(fh);
    fh=fopen([fnroot '_strain.raw'],'w'); fwrite(fh,strain(:),'float'); fclose(fh);
    fh=fopen([fnroot '_localrot.raw'],'w'); fwrite(fh,localrot(:),'float'); fclose(fh);
    disp('Registration finished')

    function [ux_x, ux_y,ux_z, uy_x, uy_y, uy_z, uz_x, uz_y, uz_z] = get3DMIRTDeformationDiff(z,y,x)
        idxX=Xknots > x - 2*h & Xknots < x + 2*h;
        idxY=Yknots > y - 2*h & Yknots < y + 2*h;
        idxZ=Zknots > z - 2*h & Zknots < z + 2*h;

        Kx = Xknots(idxX);
        Ky = Yknots(idxY);
        Kz = Zknots(idxZ);

        cx=c(idxX,idxY,idxZ,1);
        cy=c(idxX,idxY,idxZ,2);
        cz=c(idxX,idxY,idxZ,3);

        normalizedx=(x-Kx)/h;
        normalizedy=(y-Ky)/h;
        normalizedz=(z-Kz)/h;
        BKxDiff = Bspline03Diff1(normalizedx);
        BKx = Bspline03(normalizedx);
        BKyDiff = Bspline03Diff1(normalizedy);
        BKy = Bspline03(normalizedy);
        BKzDiff = Bspline03Diff1(normalizedz);
        BKz = Bspline03(normalizedz);
        ux_x=0;
        ux_y=0;
        ux_z=0;
        uy_x=0;
        uy_y=0;
        uy_z=0;
        uz_x=0;
        uz_y=0;
        uz_z=0;

        for k=1:length(BKx)
            for l=1:length(BKy)
                BKxDiffBKy=BKxDiff(k)*BKy(l);
                BKxBKyDiff=BKx(k)*BKyDiff(l);
                BKxBKy=BKx(k)*BKy(l);
                for m=1:length(BKz)
                    ux_x=ux_x+cx(k,l,m)*BKxDiffBKy*BKz(m);
                    ux_y=ux_y+cx(k,l,m)*BKxBKyDiff*BKz(m);
                    ux_z=ux_z+cx(k,l,m)*BKxBKy*BKzDiff(m);
                    uy_x=uy_x+cy(k,l,m)*BKxDiffBKy*BKz(m);
                    uy_y=uy_y+cy(k,l,m)*BKxBKyDiff*BKz(m);
                    uy_z=uy_z+cy(k,l,m)*BKxBKy*BKzDiff(m);
                    uz_x=uz_x+cz(k,l,m)*BKxDiffBKy*BKz(m);
                    uz_y=uz_y+cz(k,l,m)*BKxBKyDiff*BKz(m);
                    uz_z=uz_z+cz(k,l,m)*BKxBKy*BKzDiff(m);
                end
            end
        end    

        ux_x=ux_x/h-1;
        ux_y=ux_y/h;
        ux_z=ux_z/h;
        uy_x=uy_x/h;
        uy_y=uy_y/h-1;
        uy_z=uy_z/h;
        uz_x=uz_x/h;
        uz_y=uz_y/h;
        uz_z=uz_z/h-1;
    end
end

function [y] = Bspline03(x)
    x = abs(x);
    y=zeros(size(x));

    idx=x<1;
    xidx=x(idx);
    y(idx)=xidx.*xidx.*(xidx-2)*0.5+2/3;

    idx=x>=1 & x<2;
    xidx=x(idx)-2;
    y(idx)=-xidx.*xidx.*xidx/6;
end

function [y] = Bspline03Diff1(x)
    xorig=x;
	x = abs(x);
    y=zeros(size(x));
    
    idx=x<1;
    xidx= x(idx) ;
    y(idx) = xidx .* (xidx * 3/2 - 2);
    
    idx=x>=1 & x<2;
    xidx= 2-x(idx) ;
    y(idx) = xidx .* xidx * (-1)/(2) ;
    
    idx=xorig<0;
    y(idx)=-y(idx);
end
