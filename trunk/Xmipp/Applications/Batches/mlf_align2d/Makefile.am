bin_SCRIPTS = mlf_align2d

mlf_align2d: 
	echo "#!$(SH_PATH)" > aux_script
	cat batch_mlf_align2d >> aux_script
	mv aux_script mlf_align2d

EXTRA_DIST = batch_mlf_align2d
