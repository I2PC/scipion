
% Compare the fast and regular ways of computing Htb
% laurene.donati@epfl.ch

clear all; clc; close all;

% Set default paths
indComp = 2;
switch indComp
    case 1 % Personal computer (Mac)
        mainPath = '/Users/laurenedonati/Dropbox/multires-cryo';
    case 2 % Personal laptop (Mac)
        mainPath = '/Users/laurene/Dropbox/multires-cryo';
    case 3 % Masih's laptop
        mainPath = '/Users/masih.nilchian/Dropbox/multires-cryo';
end
addpath(genpath(mainPath));
cd([mainPath '/functions']);

%% PARAMETERS
sizeIndex     = 128;
nberProjs     = 500;
scale         = 4;
tau           = scale*0.25;
paddRatio     = 0; 

GEN_DATASET   = true;
testVersion   = 2;

kbwfProj    = struct('a', 2, 'alpha', 10.8, 'm',2);  % proj
kbwfRec     = struct('a', 4, 'alpha', 19, 'm', 2);   % rec


%% LOAD/GENERATE DATA
filename = [mainPath '/functions/tests/test-datasets/testFastHtb_v' num2str(testVersion) '.mat'];
if ~GEN_DATASET
    load(filename);
else
    disp('Generating Measurements')
    generateDataSetForTestingHtb(testVersion, sizeIndex, nberProjs, tau, kbwfProj, kbwfRec,scale, filename, indComp);
    load(filename);
end


%% TEST REGULAR HTB
tic
disp('==============================================')
disp('Starting Htb (MEX FILE).')
Htb_reg   = projectionAdjoint3DMT(projs, sizeRec, r1, r2, scale*[1,1,1], [1,1], ltRec1D, kbwfRec.a);
disp(['Time to compute Htb (MEX FILE): ' num2str(toc) ' sec'])


%% TEST FAST HTB
tic
disp('==============================================')
disp('Starting Htb (fast way, VERSION 1).')
Htb_fast = computeHtbfast(projs, kbwfRec, tau, r1, r2, scale, sizeRec, paddRatio);
disp(['Time to compute Htb (fast way, VERSION 1): ' num2str(toc) ' sec'])

%% TEST FAST HTB VERSION 2
tic
disp('==============================================')
disp('Starting Htb (fast way, VERSION 2).')
Htb_fast_V2 = computeHtbfast_v2(projs, kbwfRec, tau, r1, r2, scale, sizeRec, paddRatio);
disp(['Time to compute Htb (fast way, VERSION 2): ' num2str(toc) ' sec'])


%% COMPARE RESULTS
sliceNum = round(sizeRec(1)/2);
figure, subplot(1,3,1), imagesc(Htb_reg(:,:,sliceNum)), colormap(gray), colorbar, title('Htb_mex'), ...
    subplot(1,3,2), imagesc(Htb_fast(:,:,sliceNum)), colormap(gray), colorbar,title('Htb_fast_V1');
subplot(1,3,3), imagesc(Htb_fast_V2(:,:,sliceNum)), colormap(gray), colorbar, title('Htb_fast_V2');

figure, plot(Htb_reg(:,sliceNum,sliceNum)); hold on; ...
    plot(Htb_fast(:,sliceNum,sliceNum)); hold on; ...
    plot(Htb_fast_V2(:,sliceNum,sliceNum)); ...
    title('Profile lines'); legend('Htb_mex', 'Htb_fast_V1', 'Htb_fast_V2'); grid on;
