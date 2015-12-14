function ctfFilter=xmipp_ctf_for_metadata_row(rowNumber,Xdim,Tm)
    % Usage: xmipp_ctf_for_metadata_row(rowNumber,Xdim,Tm)
    %
    % To avoid copying md all the time, it is assumed that md is a global
    % variable
    %
    % Example of use
    % global md
    % md=xmipp_read_metadata('mymetadata.xmd')
    % ctfFilter=xmipp_ctf_for_metadata_row(1,512,1.5)
    
    
    global md

    ctf_K=1;
    ctf_defocus_angle=0;
    ctf_ca=0;
    ctf_energy_loss=0;
    ctf_lens_stability=0;
    ctf_convergence_cone=0;
    ctf_longitudinal_displacement=0;
    ctf_transversal_displacement=0;
    ctf_Q0=0.1;
    ctf_cs=0.0;
    
    if ~isfield(md,'ctfVoltage')
        error('Cannot find ctfVoltage in the metadata')
    end
    ctf_voltage=md.ctfVoltage(rowNumber);

    if ~isfield(md,'ctfDefocusU')
        error('Cannot find ctfDefocusU in the metadata')
    end
    ctf_defocusU=md.ctfDefocusU(rowNumber);

    if ~isfield(md,'ctfDefocusV')
        ctf_defocusV=ctf_defocusU;
    else
        ctf_defocusV=md.ctfDefocusV(rowNumber);
    end
    
    if isfield(md,'ctfDefocusAngle')
        ctf_defocus_angle=md.ctfDefocusAngle(rowNumber);
    end

    if isfield(md,'ctfQ0')
        ctf_Q0=md.ctfQ0(rowNumber);
    end

    if isfield(md,'ctfSphericalAberration')
        ctf_cs=md.ctfVoltage(rowNumber);
    end
    
    ctfFilter = xmipp_ctf_generate_filter(Xdim,Tm,ctf_K,ctf_voltage,ctf_defocusU,ctf_defocusV,...
        ctf_defocus_angle,ctf_cs,ctf_ca,ctf_energy_loss,ctf_lens_stability,...
        ctf_convergence_cone,ctf_longitudinal_displacement,ctf_transversal_displacement,...
        ctf_Q0);
end