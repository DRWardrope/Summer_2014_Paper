SRCS   = $(wildcard *.cc)
DOIT   = $(SRCS:%.cc=%)
ROOTCFLAGS=$(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS=$(shell ${ROOTSYS}/bin/root-config --libs)
ROOTLIBS += -lTMVA

all: $(DOIT)

%: %.cc LittlePlotter.o
	g++ -o $@ $< $(ROOTCFLAGS) $(ROOTLIBS) LittlePlotter.o

LittlePlotter.o: LittlePlotter.cpp 
	g++ -g -c $(ROOTCFLAGS) $<


clean:
	rm -f TMVAClassification FinalPlotter LittlePlotter *.gch *.o

.PHONY: clean
