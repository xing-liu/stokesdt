include ../../makevars.in

TARGET  = testspme

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))  

INC = -I$(TOPDIR)/install/include
LDFLAGS = -mkl -L$(TOPDIR)/install/lib/
LIBS = -lstokes_mob -lstokes_verlet -lstokes_util

include $(TOPDIR)/makerules.in
