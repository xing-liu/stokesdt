include ../../makevars.in

TARGET  = stokesdt

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))  

DESTDIR = $(TOPDIR)/install/bin

INC = -I. -I$(TOPDIR)/install/include/
LDFLAGS = -mkl -L$(TOPDIR)/install/lib/
LIBS =  -lstokes_simbox -lstokes_brwn -lstokes_force -lstokes_mob -lstokes_pairlist -lstokes_util

include $(TOPDIR)/makerules.in
