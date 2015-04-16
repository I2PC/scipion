function xmipp_nma_selection_tool(rundir)
% xmipp_nma_selection_tool(rundir)

    % Set path
    if exist('xmipp_nma_read_alignment','file')==0
        addpath([getenv('XMIPP_HOME') '/libraries/bindings/matlab']);
    end
    
    xmipp_nma_selection_tool_gui('rundir',rundir)
end