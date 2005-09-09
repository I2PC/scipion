bin_SCRIPTS = PDB2descr

PDB2descr: batch_PDB2descr.pl
	echo "#!$(PERL_PATH)" > aux_script
	cat batch_PDB2descr.pl >> aux_script
	mv aux_script PDB2descr

EXTRA_DIST = batch_PDB2descr.pl
