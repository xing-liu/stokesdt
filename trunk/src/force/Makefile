include ../../makevars.in

LIBTARGET  = stokes_force

SRC        = $(notdir $(wildcard *.cc))
LIBOBJS    = $(addsuffix .o, $(basename $(SRC)))  

INSTALLINC = $(notdir $(wildcard *.h))

DESTDIR = $(TOPDIR)/install

INC = -I. -I../util -I../pairlist

include $(TOPDIR)/makerules.in
