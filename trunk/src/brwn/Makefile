include ../../makevars.in

LIBTARGET  = stokes_brwn

SRC        = $(notdir $(wildcard *.cc))
LIBOBJS    = $(addsuffix .o, $(basename $(SRC)))  

INSTALLINC = $(notdir $(wildcard *.h))

DESTDIR = $(TOPDIR)/install

INC = -I. -I../util -I../mob

include $(TOPDIR)/makerules.in
