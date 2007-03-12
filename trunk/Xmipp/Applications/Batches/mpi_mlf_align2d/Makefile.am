bin_SCRIPTS = mpi_mlf_refine3d

mpi_mlf_refine3d: 
	echo "#!$(SH_PATH)" > aux_script
	cat batch_mpi_mlf_refine3d >> aux_script
	mv aux_script mpi_mlf_refine3d

EXTRA_DIST = batch_mpi_mlf_refine3d
