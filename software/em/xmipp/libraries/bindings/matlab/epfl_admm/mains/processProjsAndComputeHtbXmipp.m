% laurene.donati@epfl.ch

function [Htb, r1, r2, sizeVol, sizeRec] = processProjsAndComputeHtbXmipp(pathXmd, kbwfRec, scale, symLR)

disp('-- Processing projs and computing Htb. ')
tic

%% INITIALIZATIONS

% Get main size information
s         = xmipp_read_metadata(pathXmd);
nberProjs  = length(s.itemId);
sizeIndex = size(xmipp_read(s.image{1}),1);

% Set reconstruction size (based on scale)
sizeVol   = [sizeIndex; sizeIndex; sizeIndex];
sizeRec   = floor((sizeVol-1)/scale+1);

% Refinement paramater for look-up tables (TO OPTIMIZE)
tauHtb    = scale*0.25;

% Padding ratio for computing FFR in Htb (TO OPTIMIZE)
paddRatio =  0.05;

% Initialize Htb and its meshgrid
k1          = sizeRec(1);
Htb         = zeros(k1,k1,k1);
sizeHtb     = size(Htb);
centerHtb   = (k1-1)/2;
xValHtb     = scale*(-centerHtb:centerHtb);
[K1,K2,K3]  = meshgrid(xValHtb, xValHtb, xValHtb);

% Read symmetry information and load rotation matrices
% No symmetry <-> symLR is empty text file
[Rsymm, nSymm] = readSymmetryMatricesXmipp(symLR);


%% PRE-COMPUTATIONS

% Compute FT of refined 2D KBWF projection
kbwf   = KaiserBesselProjection2D_RefinedScaled(kbwfRec,tauHtb,scale,[sizeIndex,sizeIndex]);
kbwf   = padarray(kbwf,round(size(kbwf)*paddRatio),0,'both');
fkbwf  = fft2(kbwf);

% Generate r1 and r2
[r1, r2]  = generateProjPlanesForSymmetryXmipp(s,nberProjs,nSymm,Rsymm);

%% LOOP OVER ALL PROJECTIONS (WITH SYMMETRY)
% IT USES PARFOR IF nberProjs > 3000 

if nberProjs<10000
    % NOTE (IF YOU MODIFY SOMETHING, MODIFY IN PARFOR TOO)
    
    for p = 1:nberProjs
        
        % Load information for a particle (if enabled)
        proj   = xmipp_read(s.image{p});
        shiftX = s.shiftX(p);
        shiftY = s.shiftY(p);
        flip   = s.flip(p);
        
        % Correct shift (for a particle)
        proj   = correctShiftAndMirrorForUniqueParticle(proj,shiftX,shiftY,flip);
        
        % Generate symmetrized orientations
        numP    = ((p-1)*nSymm)+1;
        counter = 0;
        for i=1:nSymm
            ind     = numP+counter;
            % Compute Htb (for a particle) and stack
            Htbp    = computeHtbfastForUniqueParticle(sizeHtb, proj, tauHtb, paddRatio, fkbwf, r1(:,ind), r2(:,ind), K1, K2, K3);
            Htb     = Htb + Htbp;
            counter = counter + 1;
            
        end
    end
    
else
    parfor p = 1:nberProjs
        
        % Load information for a particle (if enabled)
        proj   = xmipp_read(s.image{p});
        shiftX = s.shiftX(p);
        shiftY = s.shiftY(p);
        flip   = s.flip(p);
        
        % Correct shift (for a particle)
        proj   = correctShiftAndMirrorForUniqueParticle(proj,shiftX,shiftY,flip);
        
        % Generate symmetrized orientations
        numP    = ((p-1)*nSymm)+1;
        counter = 0;
        for i=1:nSymm
            ind     = numP+counter;
            % Compute Htb (for a particle) and stack
            Htbp    = computeHtbfastForUniqueParticle(sizeHtb, proj, tauHtb, paddRatio, fkbwf, r1(:,ind), r2(:,ind), K1, K2, K3);
            Htb     = Htb + Htbp;
            counter = counter + 1;
            
        end
    end
    
    
end

disp(['-- Time to process projs and compute Htb: ' num2str(toc) ' sec.'])


end