% laurene.donati@epfl.ch

function  generateDataSetForTestingHtb(testVersion, sizeIndex, nberProjs, tau, kbwfProj, kbwfRec,scale, filename, indComp)

switch testVersion
    case 1
        
        % Construct spike image
        gt   = zeros(sizeIndex,sizeIndex,sizeIndex);
        center = floor(sizeIndex/2);
        gt(center, center, center) = 1;
        sizeVol  = [sizeIndex; sizeIndex; sizeIndex];
        
        % Compute orientations
        [rot,tilt,psi]  = generateEquidistributedAngles(nberProjs);
        [r1,r2]  = setProjectionPlanes(rot,tilt,psi);
        
        % Compute projections
        ltProj  = KaiserBesselProjection1D(kbwfProj, tau);
        projs   = projection3DMT(gt, r1, r2, [1,1,1], [1,1], ltProj, kbwfProj.a);
        
        % Set reconstruction size (based on scale)
        sizeRec  = floor((sizeVol-1)/scale+1);
        
        % Compute LT for rec
        ltRec1D   = KaiserBesselProjection1D(kbwfRec, tau);
        
    case 2
        % Construct gt image
        switch indComp
            case 1 % Personal computer (Mac)
                mainPath = '/Users/laurenedonati/Dropbox/multires-cryo';
            case 2 % Personal laptop (Mac)
                mainPath = '/Users/laurene/Dropbox/multires-cryo';
            case 3 % Masih's laptop
                mainPath = '/Users/masih.nilchian/Dropbox/multires-cryo';
        end
        addpath(genpath(mainPath));
        volName = 'betagal'; 
        cd([mainPath '/functions']);
        coeffPath   = setKBWFCoeffPath(mainPath,volName,sizeIndex);
        gt      = struct2array(load(coeffPath));
        sizeVol  = [sizeIndex; sizeIndex; sizeIndex];
        
        % Compute orientations
        [rot,tilt,psi]  = generateEquidistributedAngles(nberProjs);
        [r1,r2]  = setProjectionPlanes(rot,tilt,psi);
        
        % Compute projections
        ltProj  = KaiserBesselProjection1D(kbwfProj, tau);
        projs   = projection3DMT(gt, r1, r2, [1,1,1], [1,1], ltProj, kbwfProj.a);
        
        % Set reconstruction size (based on scale)
        sizeRec  = floor((sizeVol-1)/scale+1);
        
        % Compute LT for rec
        ltRec1D   = KaiserBesselProjection1D(kbwfRec, tau);
end
% save data
save(filename, 'gt', 'projs', 'r1', 'r2', 'sizeRec', 'ltRec1D');
end