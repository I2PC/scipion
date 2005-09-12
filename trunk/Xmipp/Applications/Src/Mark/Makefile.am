bin_PROGRAMS = mark
mark_SOURCES = mark.cc QtColorLabel.cc QtColorLabel_moc.cc \
	QtDialogFamilies.cc QtDialogFamilies_moc.cc \
	QtDialogProperties.cc QtDialogProperties_moc.cc \
	QtFileMenu.cc QtFileMenu_moc.cc \
        QtFilter.cc QtFilter_moc.cc \
        QtFilterMenu.cc QtFilterMenu_moc.cc \
        QtFiltersController.cc QtFiltersController_moc.cc \
        QtImage.cc QtImage_moc.cc \
        QtImageConverter.cc \
        QtImageMicrograph.cc QtImageMicrograph_moc.cc \
        QtImageOverviewMicrograph.cc QtImageOverviewMicrograph_moc.cc \
        QtMainWidgetMark.cc QtMainWidgetMark_moc.cc \
        QtAutoMenu.cc QtAutoMenu_moc.cc \
        QtPopupMenuMark_moc.cc \
        QtWidgetMicrograph.cc QtWidgetMicrograph_moc.cc \
        QtWidgetPSD.cc QtWidgetPSD_moc.cc \
        ../DownSample/xvsmooth.c

mark_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippGraphics/libXmippGraphics_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
        $(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	

mark_LDFLAGS = @QT_LIBS@ 

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cc

QtAutoMenu_moc.cc : QtAutoMenu.hh
	@QT_MOC@ QtAutoMenu.hh -o QtAutoMenu_moc.cc
QtColorLabel_moc.cc : QtColorLabel.hh
	@QT_MOC@ QtColorLabel.hh -o QtColorLabel_moc.cc
QtDialogFamilies_moc.cc : QtDialogFamilies.hh
	@QT_MOC@ QtDialogFamilies.hh -o QtDialogFamilies_moc.cc
QtDialogProperties_moc.cc : QtDialogProperties.hh
	@QT_MOC@ QtDialogProperties.hh -o QtDialogProperties_moc.cc
QtFileMenu_moc.cc : QtFileMenu.hh
	@QT_MOC@ QtFileMenu.hh -o QtFileMenu_moc.cc
QtFilter_moc.cc : QtFilter.hh
	@QT_MOC@ QtFilter.hh -o QtFilter_moc.cc
QtFilterMenu_moc.cc : QtFilterMenu.hh
	@QT_MOC@ QtFilterMenu.hh -o QtFilterMenu_moc.cc
QtFiltersController_moc.cc : QtFiltersController.hh
	@QT_MOC@ QtFiltersController.hh -o QtFiltersController_moc.cc
QtImage_moc.cc : QtImage.hh
	@QT_MOC@ QtImage.hh -o QtImage_moc.cc
QtImageMicrograph_moc.cc : QtImageMicrograph.hh
	@QT_MOC@ QtImageMicrograph.hh -o QtImageMicrograph_moc.cc
QtImageOverviewMicrograph_moc.cc : QtImageOverviewMicrograph.hh
	@QT_MOC@ QtImageOverviewMicrograph.hh -o QtImageOverviewMicrograph_moc.cc
QtMainWidgetMark_moc.cc : QtMainWidgetMark.hh
	@QT_MOC@ QtMainWidgetMark.hh -o QtMainWidgetMark_moc.cc
QtPopupMenuMark_moc.cc : QtPopupMenuMark.hh
	@QT_MOC@ QtPopupMenuMark.hh -o QtPopupMenuMark_moc.cc
QtWidgetMicrograph_moc.cc : QtWidgetMicrograph.hh
	@QT_MOC@ QtWidgetMicrograph.hh -o QtWidgetMicrograph_moc.cc
QtWidgetPSD_moc.cc : QtWidgetPSD.hh
	@QT_MOC@ QtWidgetPSD.hh -o QtWidgetPSD_moc.cc

EXTRA_DIST = QtAutoMenu.hh          QtFilter.hh             QtImageMicrograph.hh \
	QtColorLabel.hh        QtFilterMenu.hh         QtImageOverviewMicrograph.hh \
	QtDialogFamilies.hh    QtFiltersController.hh  QtMainWidgetMark.hh \
	QtDialogProperties.hh  QtImageConverter.hh     QtPopupMenuMark.hh \
	QtFileMenu.hh          QtImage.hh              QtWidgetMicrograph.hh \
        QtWidgetPSD.hh
