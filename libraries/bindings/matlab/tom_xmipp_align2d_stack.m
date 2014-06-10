function [success psi shiftx shifty ref_out ref stack] = tom_xmipp_align2d_stack(stack,ref,numiter,mode,waitbarflag)

if nargin < 5
    waitbarflag = true;
end

if nargin < 4
    mode = 'complete';
end
if nargin < 3
    numiter = 1;
end

if strcmp(mode,'complete') == true
    mode = true;
else 
    mode = false;
end

stacksz = size(stack);
refsz = size(ref);

if ~(refsz(1) == stacksz(1) && refsz(2) == stacksz(2))
    error('reference and particles must be of the same size.');
end

success = true(stacksz(3),numiter);
psi = zeros(stacksz(3),numiter,'single');
shiftx = zeros(stacksz(3),numiter,'single');
shifty = zeros(stacksz(3),numiter,'single');
ref_out = zeros(refsz(1),refsz(2),numiter,'single');

tforms = cell(stacksz(3),numiter);

workload = stacksz(3).*numiter;
fraction = ceil(workload./100);

worklauf = 1;

for iter = 1:numiter
    for particlenum = 1:stacksz(3)
    
        if mode == true
            st = tom_xmipp_align2d(stack(:,:,particlenum),ref,'complete');
            psi(particlenum,iter) = st.Psi;
            shiftx(particlenum,iter) = st.Xoff;
            shifty(particlenum,iter) = st.Yoff;
            success(particlenum,iter) = st.success;
            tforms{particlenum,iter} = st.Tform; 
        else
            st = tom_xmipp_align2d(stack(:,:,particlenum),ref,mode,'trans');
            shiftx(particlenum,iter) = st.Xoff;
            shifty(particlenum,iter) = st.Yoff;
            success(particlenum,iter) = st.success;
            if st.success == 1
                st = tom_xmipp_align2d(stack(:,:,particlenum),ref,mode,'rot');
                psi(particlenum,iter) = st.Psi;
                tforms{particlenum,iter} = st.Tform;
            end
        end
        if waitbarflag == true && mod(particlenum.*numiter,fraction) == 0
            percent = worklauf./workload;
            tom_HT_waitbar(percent, ['Iteration ' num2str(iter) ' of ' num2str(numiter)]);
        end
        worklauf = worklauf + 1;
    end
    
    %build reference
    for particlenum = 1:stacksz(3)
        stack(:,:,particlenum) = tom_xmipp_rotate(stack(:,:,particlenum),tforms{particlenum,iter});
    end
    ref_out(:,:,iter) = sum(stack,3);
end

tom_HT_waitbar(1);