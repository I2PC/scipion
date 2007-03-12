bin_SCRIPTS = mlf_refine3d

mlf_refine3d: 
	echo "#!$(SH_PATH)" > aux_script
	cat batch_mlf_refine3d >> aux_script
	mv aux_script mlf_refine3d

EXTRA_DIST = batch_mlf_refine3d
