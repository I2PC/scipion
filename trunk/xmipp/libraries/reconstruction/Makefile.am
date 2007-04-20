lib_LTLIBRARIES = \
	libXmippRecons.la \
	libXmippRecons_Interface.la \
	libXmippReconsINRIA.la

libXmippRecons_la_SOURCES = \
	adjust_ctf.cpp \
	adjust_volume_grey_levels.cpp \
	align2d.cpp \
	align_tilt_pairs.cpp \
	angular_assign_for_tomogram.cpp \
	angular_continuous_assign.cpp \
	angular_discrete_assign.cpp \
	angular_distance.cpp \
	angular_neighbourhood.cpp \
	angular_projection_matching.cpp \
	art_crystal.cpp\
	basic_art.cpp \
	basis.cpp \
	blobs.cpp \
	break_symmetry.cpp \
	convert_pdb2vol.cpp \
	crystal_create_surface.cpp \
	ctf_correct_idr.cpp \
	ctf_correct_phase.cpp \
	ctf.cpp \
	ctf_estimate_from_micrograph.cpp \
	directions.cpp \
	draw_surface.cpp \
	euler.cpp \
	evaluate_fscs.cpp \
	foms_evaluate.cpp \
	fourier_filter.cpp \
	grids.cpp \
	ml_align2d.cpp \
	ml_refine3d.cpp \
	phantom.cpp \
	phantom_create_random.cpp \
	phantom_simulate_microscope.cpp \
        precompute_sampling.cpp \
	project.cpp \
	project_crystal.cpp \
	projection.cpp \
	psd_enhance.cpp \
	radon.cpp \
	recons_misc.cpp \
	reconstruct_art.cpp \
	reconstruct_wbp.cpp \
	refinement.cpp \
	resolution_ssnr.cpp \
        sampling.cpp \
	spar.cpp \
	sparma.cpp \
	splines.cpp \
	surface.cpp \
	symmetries.cpp \
	volume_foms.cpp \
	volume_labelling.cpp

libXmippRecons_Interface_la_SOURCES = \
	crystal_aph2img.cpp \
	crystal_unbend.cpp\
	im_unbend.cpp \
	phantom_test_reconstruction.cpp \
	recons_spider.cpp

libXmippReconsINRIA_la_SOURCES = \
	adjust_surface.cpp

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external

EXTRA_DIST = \
	adjust_ctf.h \
	adjust_surface.h \
	adjust_volume_grey_levels.h \
	align2d.h \
	align_tilt_pairs.h \
	angular_assign_for_tomogram.h \
	angular_continuous_assign.h \
	angular_discrete_assign.h \
	angular_distance.h \
	angular_neighbourhood.h \
	angular_projection_matching.h \
	art_crystal.h \
	basic_art.h \
	basis.h \
	blobs.h \
	break_symmetry.h \
	convert_pdb2vol.h \
	crystal_aph2img.h \
	crystal_create_surface.h \
	crystal_unbend.h \
	ctf_correct_idr.h \
	ctf_correct_phase.h \
	ctf_estimate_from_micrograph.h \
	ctf.h \
	directions.h \
	draw_surface.h \
	euler.h \
	evaluate_fscs.h \
	foms_evaluate.h \
	fourier_filter.h \
	gibbs.h \
	grids.h \
	im_unbend.h \
	ml_align2d.h \
	ml_refine3d.h \
	phantom_create_random.h \
	phantom.h \
	phantom_simulate_microscope.h \
	phantom_test_reconstruction.h \
        precompute_sampling.h \
	project_crystal.h \
	project.h \
	projection.h \
	projection_matching_crystal.h \
	psd_enhance.h \
	radon.h \
	recons_misc.h \
	recons_spider.h \
	reconstruct_art.h \
	reconstruct_wbp.h \
	refinement.h \
	resolution_ssnr.h \
        sampling.h \
	spar.h \
	sparma.h \
	splines.h \
	spots2realspace3d.h \
	surface.h \
	symmetries.h \
	symmetrize.h \
	volume_foms.h \
	volume_labelling.h

docs:
	@echo "Creating documentation in Reconstruction ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
