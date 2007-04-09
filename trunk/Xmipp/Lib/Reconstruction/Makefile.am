lib_LTLIBRARIES = \
	libXmippRecons_@host_os@.la \
	libXmippRecons_Interface_@host_os@.la \
	libXmippReconsINRIA_@host_os@.la

libXmippRecons_@host_os@_la_SOURCES = \
	Src/grids.cc \
	Src/blobs.cc \
	Src/splines.cc \
	Src/symmetries.cc \
	Src/phantom.cc \
	Src/projection.cc \
	Src/volume_FOMs.cc \
	Src/volume_labelling.cc \
	Src/radon.cc \
	Src/directions.cc \
	Src/basis.cc \
	Programs/Src/Prog_project.cc \
	Programs/Src/Prog_evaluate.cc \
	Programs/Src/Prog_evaluate_FSCs.cc \
	Programs/Src/recons_misc.cc \
	Programs/Src/Prog_project_crystal.cc \
	Programs/Src/Prog_surface.cc \
	Programs/Src/Prog_euler.cc \
	Programs/Src/Prog_draw_surface.cc \
	Programs/Src/Prog_symmetrize.cc \
	Programs/Src/Prog_Angular_Distance.cc \
	Programs/Src/Prog_Align2d.cc \
	Programs/Src/Prog_Projection_Neighbourhood.cc \
	Programs/Src/Prog_Centilt.cc \
	Src/CTF.cc \
	Src/refinement.cc \
	Programs/Src/Prog_FourierFilter.cc \
	Programs/Src/Prog_IDR_art.cc \
	Programs/Src/Prog_CorrectPhase.cc \
	Programs/Src/Prog_SpAR.cc \
	Programs/Src/Prog_art_crystal.cc\
	Programs/Src/Prog_SpARMA.cc \
	Programs/Src/Prog_random_phantom.cc \
	Programs/Src/Prog_microscope.cc \
	Programs/Src/Prog_adjust_CTF.cc \
	Programs/Src/Prog_adjust_volume.cc \
	Programs/Src/Prog_assign_CTF.cc \
	Programs/Src/Prog_art.cc \
	Programs/Src/Basic_art.cc \
	Programs/Src/Prog_SSNR.cc \
	Programs/Src/Prog_Angular_Predict.cc \
	Programs/Src/Prog_Angular_Predict_continuous.cc \
	Programs/Src/Prog_Angular_Predict_Tomography.cc \
	Programs/Src/Prog_MLalign2D.cc \
	Programs/Src/Prog_Refine3d.cc \
	Programs/Src/Prog_break_sym.cc \
	Programs/Src/Prog_WBP.cc \
	Programs/Src/Prog_projection_matching.cc \
	Programs/Src/Prog_Enhance_PSD.cc \
	Programs/Src/Prog_PDBPhantom.cc

libXmippRecons_Interface_@host_os@_la_SOURCES = \
	Programs/Src/Prog_Spots2RealSpace2D.cc \
	Src/recons_spider.cc \
	Programs/Src/Prog_Unbend.cc\
	Src/ImUmbend.cc \
	Programs/Src/Prog_recons_test.cc

libXmippReconsINRIA_@host_os@_la_SOURCES = \
	Programs/Src/Prog_adjust_surface.cc

INCLUDES = -I$(top_srcdir)/Lib

EXTRA_DIST = \
	blobs.hh \
	CTF.hh \
	directions.hh \
	basis.hh \
	splines.hh \
	ImUmbend.hh \
	Gibbs.hh \
	grids.hh \
	phantom.hh \
	Programs/Basic_art.hh \
	Programs/Prog_adjust_CTF.hh \
	Programs/Prog_adjust_surface.hh \
	Programs/Prog_adjust_volume.hh \
	Programs/Prog_Align2d.hh \
	Programs/Prog_Angular_Distance.hh \
	Programs/Prog_Angular_Predict_continuous.hh \
	Programs/Prog_Angular_Predict_Tomography.hh \
	Programs/Prog_Angular_Predict.hh \
	Programs/Prog_art_crystal.hh \
	Programs/Prog_art.hh \
	Programs/Prog_assign_CTF.hh \
	Programs/Prog_Centilt.hh \
	Programs/Prog_CorrectPhase.hh \
	Programs/Prog_draw_surface.hh \
	Programs/Prog_euler.hh \
	Programs/Prog_evaluate_FSCs.hh \
	Programs/Prog_evaluate.hh \
	Programs/Prog_FourierFilter.hh \
	Programs/Prog_IDR_art.hh \
	Programs/Prog_microscope.hh \
	Programs/Prog_MLalign2D.hh \
	Programs/Prog_Refine3d.hh \
	Programs/Prog_break_sym.hh \
	Programs/Prog_project_crystal.hh \
	Programs/Prog_project.hh \
	Programs/Prog_projection_matching.hh \
	Programs/Prog_Enhance_PSD.hh \
	Programs/Prog_Projection_Neighbourhood.hh \
	Programs/Prog_random_phantom.hh \
	Programs/Prog_recons_test.hh \
	Programs/Prog_SpAR.hh \
	Programs/Prog_SpARMA.hh \
	Programs/Prog_Spots2RealSpace2D.hh \
	Programs/Prog_Spots2RealSpace3D.hh \
	Programs/Prog_Unbend.hh \
	Programs/Prog_SSNR.hh \
	Programs/Prog_surface.hh \
	Programs/Prog_symmetrize.hh \
	Programs/Prog_WBP.hh \
	Programs/Prog_PDBPhantom.hh \
	Programs/recons_misc.hh \
	Programs/Src/Basic_art.inc \
	projection.hh \
	radon.hh \
	recons_spider.hh \
	refinement.hh \
	symmetries.hh \
	volume_FOMs.hh \
	volume_labelling.hh

docs:
	@echo "Creating documentation in Reconstruction ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
