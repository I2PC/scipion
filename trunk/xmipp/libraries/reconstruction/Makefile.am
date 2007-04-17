lib_LTLIBRARIES = \
	libXmippRecons.la \
	libXmippRecons_Interface.la \
	libXmippReconsINRIA.la

libXmippRecons_la_SOURCES = \
	grids.cpp \
	blobs.cpp \
	splines.cpp \
	symmetries.cpp \
	phantom.cpp \
	projection.cpp \
	volume_foms.cpp \
	volume_labelling.cpp \
	radon.cpp \
	directions.cpp \
	basis.cpp \
	project.cpp \
	foms_evaluate.cpp \
	evaluate_fscs.cpp \
	recons_misc.cpp \
	project_crystal.cpp \
	surface.cpp \
	euler.cpp \
	draw_surface.cpp \
	symmetrize.cpp \
	angular_distance.cpp \
	align2d.cpp \
	angular_neighbourhood.cpp \
	align_tilt_pairs.cpp \
	ctf.cpp \
	refinement.cpp \
	fourier_filter.cpp \
	ctf_correct_idr.cpp \
	ctf_correct_phase.cpp \
	spar.cpp \
	art_crystal.cpp\
	sparma.cpp \
	phantom_create_random.cpp \
	phantom_simulate_microscope.cpp \
	adjust_ctf.cpp \
	adjust_volume_grey_levels.cpp \
	ctf_estimate_from_micrograph.cpp \
	reconstruct_art.cpp \
	basic_art.cpp \
	resolution_ssnr.cpp \
	angular_discrete_assign.cpp \
	angular_continuous_assign.cpp \
	angular_assign_for_tomogram.cpp \
	ml_align2d.cpp \
	ml_refine3d.cpp \
	break_symmetry.cpp \
	reconstruct_wbp.cpp \
	angular_projection_matching.cpp \
	psd_enhance.cpp \
	convert_pdb2vol.cpp \
	crystal_create_surface.cpp

libXmippRecons_Interface_la_SOURCES = \
	crystal_aph2img.cpp \
	recons_spider.cpp \
	crystal_unbend.cpp\
	im_unbend.cpp \
	phantom_test_reconstruction.cpp

libXmippReconsINRIA_la_SOURCES = \
	adjust_surface.cpp

INCLUDES = -I$(top_srcdir)/libraries

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
