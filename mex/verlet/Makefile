include ../../makevars.in

MEXTARGET=compute_verlet.${MEXEXT}

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))
INC     = -I../install/include

LIBS    = $(TOPDIR)/install/lib/libstokes_verlet.a \
          $(TOPDIR)/install/lib/libstokes_util.a

DESTDIR = $(TOPDIR)/install/mex

include $(TOPDIR)/makerules.in
