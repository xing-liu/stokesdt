include ../../makevars.in

TARGET  = teststeric

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))  

INC = -I$(TOPDIR)/install/include
LDFLAGS = -mkl -L$(TOPDIR)/install/lib/
LIBS = -lstokes_force -lstokes_pairlist -lstokes_util

include $(TOPDIR)/makerules.in
