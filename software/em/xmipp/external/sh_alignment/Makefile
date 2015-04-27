# Hand-made Makefile following the compile.sh script.

CFLAGS=-O3 -fPIC
LDFLAGS=-L$(SCIPION_HOME)/software/lib

DIRS := frm/src SpharmonicKit27
SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.c))
OBJS := $(patsubst %.c, %.o, $(SOURCES)) frm_wrap.o frm.o

PYTHON_PACKAGES := $(SCIPION_HOME)/software/lib/python2.7/site-packages
INCS := \
  -Ifrm/src -ISpharmonicKit27 \
  -I$(SCIPION_HOME)/software/include \
  -I$(SCIPION_HOME)/software/include/python2.7 \
  -I$(PYTHON_PACKAGES)/numpy/core/include


all: _swig_frm.so frm/swig/swig_frm.py

frm/swig/frm_wrap.c frm/swig/swig_frm.py: frm/swig/frm.i
	swig -python frm/swig/frm.i 

frm_wrap.o: frm/swig/frm_wrap.c
	gcc -c $(CFLAGS) $(INCS) frm/swig/frm_wrap.c

frm.o: frm/swig/frm.c
	gcc -c $(CFLAGS) $(INCS) -Ifrm/src frm/swig/frm.c

_swig_frm.so: $(OBJS)
	gcc -o _swig_frm.so -shared $(LDFLAGS) -lfftw3 $(OBJS)
	mkdir -p $(PYTHON_PACKAGES)/sh_alignment
	cp *.py frm/swig/swig_frm.py _swig_frm.so $(PYTHON_PACKAGES)/sh_alignment
	cp -r tompy $(PYTHON_PACKAGES)/sh_alignment


clean:
	rm -f $(OBJS) frm/swig/frm_wrap.c frm/swig/swig_frm.py

.PHONY: clean all
