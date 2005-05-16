lib_LTLIBRARIES = libXmippRecons_@host_os@.la \
	libXmippRecons_Interface_@host_os@.la \
	libXmippReconsINRIA_@host_os@.la

libXmippRecons_@host_os@_la_SOURCES = Src/grids.cc \
	Src/blobs.cc \
	Src/symmetries.cc \
	Src/phantom.cc \
	Src/projection.cc \
	Src/volume_FOMs.cc \
	Src/radon.cc \
	Src/Gibbs.cc \
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
	Programs/Src/Prog_assign_CTF.cc \
	Programs/Src/Prog_art.cc \
	Programs/Src/Basic_art.cc \
	Programs/Src/Prog_SSNR.cc \
	Programs/Src/Prog_Angular_Predict.cc \
	Programs/Src/Prog_Angular_Predict_continuous.cc \
	Programs/Src/Prog_MLalign2D.cc \
	Programs/Src/Prog_WBP.cc \
	Programs/Src/Prog_projection_matching.cc

libXmippRecons_Interface_@host_os@_la_SOURCES = Programs/Src/Prog_Spots2RealSpace2D.cc \
	Src/recons_spider.cc \
	Programs/Src/Prog_recons_test.cc

libXmippReconsINRIA_@host_os@_la_SOURCES = Programs/Src/Prog_adjust_surface.cc


INCLUDES = -I$(top_srcdir)/Lib
