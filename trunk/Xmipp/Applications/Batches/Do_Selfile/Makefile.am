bin_SCRIPTS = do_selfile

do_selfile: 
	echo "#!$(SH_PATH)" > aux_script
	cat batch_do_selfile >> aux_script
	mv aux_script do_selfile

EXTRA_DIST = batch_do_selfile
