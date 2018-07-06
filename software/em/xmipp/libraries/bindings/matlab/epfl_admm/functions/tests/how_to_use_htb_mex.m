%     % Compute/load look-up table of 1D KBWF projection (for reconstruction)
%     ltRec1D   = KaiserBesselProjection1D(kbwfRec, tauLtRec1D);
%     % Compute Htb (slow version for now)
%     projs = applyCTFAdjoint(projs,ctf);
%     Htb   = projectionAdjoint3DMT(projs, sizeRec, r1, r2, ...
%         scale*[1,1,1], [1,1], ltRec1D, kbwfRec.a);

% % This was the way I used in when testing it in the xmipp reconstruction
% function. 
% NOTE: Maybe first retry with the original Htb to check that everything is working well. 
% NOTE: TESTED, everything is working well !
% tauLtRec1D = 0.01; 
% disp('Computing Htb with mex file')
% ltRec1D   = KaiserBesselProjection1D(kbwfRec, tauLtRec1D);
% Htb_mex   = projectionAdjoint3DMT(projs, sizeRec, r1, r2, scale*[1,1,1],
% [1,1], ltRec1D, kbwfRec.a); 
% figure, imagesc(Htb_mex(:,:,floor(sizeRec(1)/2))), title('Htb Mex');
% figure, plot(Htb(:,floor(sizeRec(1)/2),floor(sizeRec(1)/2))); hold on;
% plot(Htb_mex(:,floor(sizeRec(1)/2),floor(sizeRec(1)/2))), grid on; title('Profile lines');
% legend('Htb fast', 'Htb mex')