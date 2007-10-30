function tom_mex_xmippbindings(xmipp_path)

if nargin < 1
    if strcmp(computer,'PCWIN')
        xmipp_path = 'y:\xmipp';
    else
        xmipp_path = '/usr/local/apps/xmipp';
    end
        
end



includes = {[xmipp_path '/lib']};
includes2 = {[xmipp_path '/libraries/'], ...
             xmipp_path, ...
             [xmipp_path '/libraries/data/'], ...
             [xmipp_path '/libraries/interface/'], ...
             [xmipp_path '/libraries/classification/'], ...
             [xmipp_path '/libraries/reconstruction/'], ...
             [xmipp_path '/libraries/legacy/']};
libs = {'-lXmippClassif', ...
        '-lOldXmipp', ...
        '-lXmippRecons', ...
        '-lXmippData', ...
        '-lXmippInterface', ...
        '-lXmippExternal'};
binding_files = {'tom_xmipp_normalize_wrapper', ...
                 'tom_xmipp_psd_enhance_wrapper', ...
                 'tom_xmipp_adjust_ctf_wrapper', ...
                 'tom_xmipp_scale_pyramid_wrapper', ...
                 'tom_xmipp_scale_wrapper', ...
                 'tom_xmipp_rotate_wrapper', ...
                 'tom_xmipp_align2d_wrapper', ...
                 'tom_xmipp_mirror_wrapper', ...
                 'tom_xmipp_morphology_wrapper', ...
                 'tom_xmipp_mask_wrapper', ...
                 'tom_xmipp_resolution_wrapper', ...
                 'tom_xmipp_ctf_correct_phase_wrapper'};

%--------------------------------------------------------------------------

if strcmp(computer,'PCWIN')
    includes = {'z:\tom_dev\xmipp_bindings'};
    includes2 = strrep(includes2,'/','\'); 
end


string = [];
for i = includes
    string = [string ' -L' i{1}];
end

for i=includes2
    string = [string ' -I' i{1}];
end

if ~strcmp(computer,'PCWIN')
    for i=libs
        string = [string ' ' i{1}];
    end
end

for file = binding_files
    disp(['Compiling ' file{1}]);
    %unix(['rm ' file{1} '.mex*']);
    if strcmp(computer,'PCWIN')
        eval(['mex -f c:\gnumex\mexopts.bat -O -D_CYGWIN ' string ' ' file{1} '.cpp']);
    else
        eval(['mex -O -D_Linux ' string ' ' file{1} '.cpp']);
    end
end

rehash toolbox;
