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
        QtPopupMenuMark_moc.cc \
        QtWidgetMicrograph.cc QtWidgetMicrograph_moc.cc \
        ../DownSample/xvsmooth.c

mark_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	

mark_LDFLAGS = @QT_LIBS@ 

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

QtColorLabel_moc.cc : 
	@QT_MOC@ QtColorLabel.hh -o QtColorLabel_moc.cc
QtDialogFamilies_moc.cc : 
	@QT_MOC@ QtDialogFamilies.hh -o QtDialogFamilies_moc.cc
QtDialogProperties_moc.cc : 
	@QT_MOC@ QtDialogProperties.hh -o QtDialogProperties_moc.cc
QtFileMenu_moc.cc : 
	@QT_MOC@ QtFileMenu.hh -o QtFileMenu_moc.cc
QtFilter_moc.cc : 
	@QT_MOC@ QtFilter.hh -o QtFilter_moc.cc
QtFilterMenu_moc.cc : 
	@QT_MOC@ QtFilterMenu.hh -o QtFilterMenu_moc.cc
QtFiltersController_moc.cc : 
	@QT_MOC@ QtFiltersController.hh -o QtFiltersController_moc.cc
QtImage_moc.cc : 
	@QT_MOC@ QtImage.hh -o QtImage_moc.cc
QtImageMicrograph_moc.cc : 
	@QT_MOC@ QtImageMicrograph.hh -o QtImageMicrograph_moc.cc
QtImageOverviewMicrograph_moc.cc : 
	@QT_MOC@ QtImageOverviewMicrograph.hh -o QtImageOverviewMicrograph_moc.cc
QtMainWidgetMark_moc.cc : 
	@QT_MOC@ QtMainWidgetMark.hh -o QtMainWidgetMark_moc.cc
QtPopupMenuMark_moc.cc : 
	@QT_MOC@ QtPopupMenuMark.hh -o QtPopupMenuMark_moc.cc
QtWidgetMicrograph_moc.cc : 
	@QT_MOC@ QtWidgetMicrograph.hh -o QtWidgetMicrograph_moc.cc
