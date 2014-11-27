include ../../makevars.in

TARGET  = testrandom

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))  

INC = -I$(TOPDIR)/install/include
LDFLAGS = -mkl -L$(TOPDIR)/install/lib/
LIBS = -lstokes_mob -lstokes_brwn -lstokes_util -lstokes_pairlist

include $(TOPDIR)/makerules.in
