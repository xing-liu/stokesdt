include ../../makevars.in

LIBTARGET  = stokes_simbox

SRC        = $(notdir $(wildcard *.cc))
LIBOBJS    = $(addsuffix .o, $(basename $(SRC)))  

INSTALLINC = $(notdir $(wildcard *.h))

DESTDIR = $(TOPDIR)/install

INC = -I. -I../util -I../mob -I../brwn -I../force -I../pairlist

include $(TOPDIR)/makerules.in
