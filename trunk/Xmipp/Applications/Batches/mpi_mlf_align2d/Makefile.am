bin_SCRIPTS = mpi_mlf_align2d

mpi_mlf_align2d: 
	echo "#!$(SH_PATH)" > aux_script
	cat batch_mpi_mlf_align2d >> aux_script
	mv aux_script mpi_mlf_align2d

EXTRA_DIST = batch_mpi_mlf_align2d
