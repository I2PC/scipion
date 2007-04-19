bin_PROGRAMS = micrograph_mark
micrograph_mark_SOURCES = \
	main.cpp \
	color_label.cpp \
	dialog_families.cpp \
	dialog_properties.cpp \
	file_menu.cpp \
	filter.cpp \
	filter_menu.cpp \
	filters_controller.cpp \
	image.cpp \
	image_converter.cpp \
	image_micrograph.cpp \
	image_overview_micrograph.cpp \
	main_widget_mark.cpp \
	auto_menu.cpp \
	widget_micrograph.cpp \
	widget_psd.cpp \
	../micrograph_downsample/xvsmooth.cpp
	## FIXME This is a bad, bad practice..

nodist_micrograph_mark_SOURCES = \
	color_label_moc.cpp \
	dialog_families_moc.cpp \
	dialog_properties_moc.cpp \
	file_menu_moc.cpp \
	filter_moc.cpp \
	filter_menu_moc.cpp \
	filters_controller_moc.cpp \
	image_moc.cpp \
	image_micrograph_moc.cpp \
	image_overview_micrograph_moc.cpp \
	main_widget_mark_moc.cpp \
	auto_menu_moc.cpp \
	popup_menu_mark_moc.cpp \
	widget_micrograph_moc.cpp \
	widget_psd_moc.cpp

micrograph_mark_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/graphics/libXmippGraphics.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

micrograph_mark_LDFLAGS = @QT_LIBS@

## FIXME Same as above
INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I../micrograph_downsample \
	@QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cpp

auto_menu_moc.cpp : auto_menu.h
	@QT_MOC@ auto_menu.h -o auto_menu_moc.cpp

color_label_moc.cpp : color_label.h
	@QT_MOC@ color_label.h -o color_label_moc.cpp

dialog_families_moc.cpp : dialog_families.h
	@QT_MOC@ dialog_families.h -o dialog_families_moc.cpp

dialog_properties_moc.cpp : dialog_properties.h
	@QT_MOC@ dialog_properties.h -o dialog_properties_moc.cpp

file_menu_moc.cpp : file_menu.h
	@QT_MOC@ file_menu.h -o file_menu_moc.cpp

filter_moc.cpp : filter.h
	@QT_MOC@ filter.h -o filter_moc.cpp

filter_menu_moc.cpp : filter_menu.h
	@QT_MOC@ filter_menu.h -o filter_menu_moc.cpp

filters_controller_moc.cpp : filters_controller.h
	@QT_MOC@ filters_controller.h -o filters_controller_moc.cpp

image_moc.cpp : image.h
	@QT_MOC@ image.h -o image_moc.cpp

image_micrograph_moc.cpp : image_micrograph.h
	@QT_MOC@ image_micrograph.h -o image_micrograph_moc.cpp

image_overview_micrograph_moc.cpp : image_overview_micrograph.h
	@QT_MOC@ image_overview_micrograph.h -o image_overview_micrograph_moc.cpp

main_widget_mark_moc.cpp : main_widget_mark.h
	@QT_MOC@ main_widget_mark.h -o main_widget_mark_moc.cpp

popup_menu_mark_moc.cpp : popup_menu_mark.h
	@QT_MOC@ popup_menu_mark.h -o popup_menu_mark_moc.cpp

widget_micrograph_moc.cpp : widget_micrograph.h
	@QT_MOC@ widget_micrograph.h -o widget_micrograph_moc.cpp

widget_psd_moc.cpp : widget_psd.h
	@QT_MOC@ widget_psd.h -o widget_psd_moc.cpp

EXTRA_DIST = \
	autoMenu.h \
	filter.h \
	image_micrograph.h \
	color_label.h \
	filter_menu.h \
	image_overview_micrograph.h \
	dialog_families.h \
	filters_controller.h \
	main_widget_mark.h \
	dialog_properties.h \
	image_converter.h \
	popup_menu_mark.h \
	file_menu.h \
	image.h \
	widget_micrograph.h \
	widget_psd.h

