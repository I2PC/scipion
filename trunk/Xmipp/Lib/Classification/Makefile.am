lib_LTLIBRARIES = libXmippClassif_@host_os@.la

libXmippClassif_@host_os@_la_SOURCES = \
	Src/xmippCTVectors.cc \
	Src/xmippCodeBook.cc \
	Src/xmippFuzzyCodeBook.cc \
	Src/xmippDistances.cc \
	Src/xmippSammon.cc \
	Src/xmippFCMeans.cc \
	Src/xmippMap.cc \
	Src/xmippSOM.cc \
	Src/xmippFKCN.cc \
	Src/xmippFuzzySOM.cc \
	Src/xmippKerDenSOM.cc \
	Src/xmippBatchSOM.cc \
	Src/xmippUmatrix.cc \
	Src/xmippPlanes.cc \
	Src/xmippTStudentKerDenSOM.cc \
	Src/xmippGaussianKerDenSOM.cc \
	Src/xmippPC.cc \
	Src/xmippRBF.cc \
	Src/xmippRM.cc

INCLUDES = -I$(top_srcdir)/Lib

EXTRA_DIST = \
	xmippBaseAlgo.hh \
	xmippFKCN.hh \
	xmippRM.hh \
	xmippBatchSOM.hh \
	xmippFuzzyCodeBook.hh \
	xmippSammon.hh \
	xmippCDataTypes.hh \
	xmippFuzzySOM.hh \
	xmippSOM.hh \
	xmippCDSet.hh \
	xmippGaussianKerDenSOM.hh \
	xmippTStudentKerDenSOM.hh \
	xmippCodeBook.hh \
	xmippKerDenSOM.hh \
	xmippUmatrix.hh \
	xmippCTSet.hh \
	xmippMap.hh \
	xmippUniform.hh \
	xmippCTVectors.hh \
	xmippPC.hh \
	xmippVectorOps.hh \
	xmippDistances.hh \
	xmippPlanes.hh \
	xmippFCMeans.hh \
	xmippRBF.hh

docs:
	@echo "Creating documentation in Classification ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
