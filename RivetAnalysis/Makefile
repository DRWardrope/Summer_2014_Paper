SRCS   = $(wildcard *.cc)
LIBS   = $(SRCS:%.cc=Rivet%.so)
DATA   = $(wildcard *.yoda)
PLOTS  = $(DATA:%.yoda=plots/%)
MKHTML = $(DATA:%.yoda=%.yoda:Title=%)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
ROOTLIBS += -lTMVA

all: $(LIBS)

Rivet%.so: %.cc
	rivet-buildplugin $@ $< $(EXTRA_CFLAGS) $(CPPEXPFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

plots: $(PLOTS)
	cp -r plots/* ~/www/HH-rivet/

allplots: $(DATA)
	rivet-mkhtml --no-ratio -o plots/ $^
	cp -r plots/* ~/www/HH-rivet/

plots/%: %.yoda
	rivet-mkhtml --no-ratio -o plots/ $<:"Title=$*"
