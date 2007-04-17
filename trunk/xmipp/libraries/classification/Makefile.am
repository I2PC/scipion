lib_LTLIBRARIES = libXmippClassif.la

libXmippClassif_la_SOURCES = \
	batch_som.cpp \
	code_book.cpp \
	distance.cpp \
	fcmeans.cpp \
	fkcn.cpp \
	fuzzy_code_book.cpp \
	fuzzy_som.cpp \
	gaussian_kerdensom.cpp \
	kerdensom.cpp \
	map.cpp \
	pca.cpp \
	plane.cpp \
	rbf.cpp \
	rm.cpp \
	sammon.cpp \
	som.cpp \
	training_vector.cpp \
	tstudent_kerdensom.cpp \
	umatrix.cpp

INCLUDES = -I$(top_srcdir)/libraries

EXTRA_DIST = \
	base_algorithm.h \
	batch_som.h \
	code_book.h \
	data_set.h \
	data_types.h \
	distance.h \
	fcmeans.h \
	fkcn.h \
	fuzzy_code_book.h \
	fuzzy_som.h \
	gaussian_kerdensom.h \
	kerdensom.h \
	map.h \
	pca.h \
	plane.h \
	rbf.h \
	rm.h \
	sammon.h \
	som.h \
	training_set.h \
	training_vector.h \
	tstudent_kerdensom.h \
	umatrix.h \
	uniform.h \
	vector_ops.h

docs:
	@echo "Creating documentation in Classification ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
