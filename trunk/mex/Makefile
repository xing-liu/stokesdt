include ../makevars.in

SOURCES = $(wildcard *_mex.cc)
TARGETS = $(SOURCES:.cc=.$(MEXEXT))

SRC     = $(notdir $(wildcard *.cc))
OBJS    = $(addsuffix .o, $(basename $(SRC)))
INC     = -I../install/include

LIBS    = $(TOPDIR)/install/lib/libstokes_brwn.a \
          $(TOPDIR)/install/lib/libstokes_force.a \
          $(TOPDIR)/install/lib/libstokes_mob.a \
          $(TOPDIR)/install/lib/libstokes_pairlist.a \
          $(TOPDIR)/install/lib/libstokes_simbox.a \
          $(TOPDIR)/install/lib/libstokes_util.a

all: $(TARGETS)

%.$(MEXEXT): %.cc
	$(MEX) $^ $(MEXFLAGS) ${LIBS}

clean:
	rm -f $(TARGETS)
