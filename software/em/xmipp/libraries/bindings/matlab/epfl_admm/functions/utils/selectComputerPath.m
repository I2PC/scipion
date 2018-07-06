
% laurene.donati@epfl.ch 

function [mainPath] = selectComputerPath(indComp) 

switch indComp
    case 1 % Personal computer (Mac)
        mainPath = '/Users/laurenedonati/Dropbox/multires-cryo';
    case 2 % Personal laptop (Mac)
        mainPath = '/Users/laurene/Dropbox/multires-cryo';
    case 3 % Masih's laptop
        mainPath = '/Users/masih.nilchian/Dropbox/multires-cryo';
        
end 
end