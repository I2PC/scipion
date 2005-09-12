bin_SCRIPTS = PDBsidechain

PDBsidechain: batch_PDBsidechain.pl
	echo "#!$(PERL_PATH)" > aux_script
	cat batch_PDBsidechain.pl >> aux_script
	mv aux_script PDBsidechain

EXTRA_DIST = batch_PDBsidechain.pl
