include ../../makevars.in

TARGET  = run_spme

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))  

DESTDIR = $(TOPDIR)/install

INC = -I. -I$(TOPDIR)/install/include/
LDFLAGS = -mkl -L$(TOPDIR)/install/lib/
LIBS =  -lstokes_simbox -lstokes_brwn -lstokes_force -lstokes_mob -lstokes_verlet -lstokes_util

include $(TOPDIR)/makerules.in
