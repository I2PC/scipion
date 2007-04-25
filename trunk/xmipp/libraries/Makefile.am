OPTIONAL_LIBS =

if BUILD_QT
OPTIONAL_LIBS += graphics
endif

if BUILD_MPI
OPTIONAL_LIBS += reconstruction_mpi
endif

SUBDIRS = \
	data \
	reconstruction \
	classification \
	interface \
	legacy \
	$(OPTIONAL_LIBS)
