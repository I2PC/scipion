OPTIONAL_PROGS =

if BUILD_MPI
OPTIONAL_PROGS += \
	mpi_angular_discrete_assign \
	mpi_angular_continuous_assign \
	mpi_ml_align2d \
	mpi_ml_refine3d \
	mpi_angular_projection_matching \
	mpi_reconstruct_art \
	mpi_reconstruct_wbp \
	mpi_run
endif

if BUILD_TIFF
OPTIONAL_PROGS += convert_tiff2raw
endif

if BUILD_QT
OPTIONAL_PROGS += \
	ctf_view \
	edit \
	micrograph_mark \
	show \
	mask_design
endif

SUBDIRS = \
	add_noise \
	adjust_volume_grey_levels \
	align_volumes \
	align2d \
	angular_distribution_show \
	angular_distance \
	angular_discrete_assign \
	angular_continuous_assign \
	angular_assign_for_tomogram \
	header_apply \
	docfile_append \
	reconstruct_art \
	ctf_estimate_from_micrograph \
	classify_batch_som \
	break_symmetry \
	align_tilt_pairs \
	selfile_compare \
	ctf_correct_phase \
	correlation \
	selfile_copy \
	phantom_create \
	ctf_profile \
	convert_data2img \
	convert_data2vol \
	crystal_shear \
	denoise \
	micrograph_denoise \
	docfile_histogram \
	micrograph_downsample \
	crystal_euler2mrc \
	psd_enhance \
	foms_evaluate \
	classify_fcmeans \
	find_center3d \
	classify_fkcn \
	mirror \
	fourier_filter \
	classify_fsom \
	header_reset \
	header_assign \
	header_extract \
	histogram \
	ctf_correct_idr \
	convert_img2data \
	header_print \
	statistics \
	classify_kcmeans \
	classify_kerdensom \
	crystal_lattice_vectors \
	make_spectra \
	mask \
	try_symmetry \
	phantom_simulate_microscope \
	ml_align2d \
	ml_align2d_combine \
	ml_refine3d \
	morphology \
	selfile_move \
	normalize \
	operate \
	classify_pca \
	convert_pdb2vol \
	phantom_transform \
	project \
	angular_neighbourhood \
	angular_projection_matching \
	classify_project_pca \
	scale_pyramid \
	range_adjust \
	convert_raw22spi \
	phantom_test_reconstruction \
	reverse_endian \
	resolution_fsc \
	selfile_delete \
	rotate \
	classify_sammon_projection \
	scale \
	micrograph_scissor \
	volume_segment \
	selfile_select \
	selfile_statistics \
	separate_objects \
	shift \
	crystal_skew \
	classify_som \
	sort_by_statistics \
	selfile_split \
	convert_spi22em \
	convert_spi22ccp4 \
	crystal_aph2img \
	resolution_ssnr \
	average \
	symmetrize \
	threshold \
	convert_data2pdb \
	convert_data2descr \
	crystal_unbend \
	fourier_transform \
	convert_vol2data \
	convert_vol2pseudo \
	convert_voxels22blobs \
	reconstruct_wbp \
	window \
	micrograph_window \
	$(OPTIONAL_PROGS)
