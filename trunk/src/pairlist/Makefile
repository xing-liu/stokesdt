include ../../makevars.in

LIBTARGET  = stokes_pairlist

SRC        = $(notdir $(wildcard *.cc))
LIBOBJS    = $(addsuffix .o, $(basename $(SRC)))  

INSTALLINC = $(notdir $(wildcard *.h))

DESTDIR = $(TOPDIR)/install

INC = -I. -I../util

include $(TOPDIR)/makerules.in
