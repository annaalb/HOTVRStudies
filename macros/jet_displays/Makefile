# Makefile for intro to C++ with ROOT

CPP=$(shell root-config --cxx)
LD=$(shell root-config --ld)

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
ROOTINC=$(shell root-config --incdir)

CFLAGS =
LFLAGS =

RCXX=$(CFLAGS) $(ROOTCFLAGS) -I$(ROOTINC)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC1= style.o jet_display_helpers.o hist_helpers.o make_history.o
SRC2= style.o jet_display_helpers.o jet_display_akt.o jet_display_hotvr.o make_displays.o

%.o: %.C
	$(CPP) $(RCXX) -c $<

all: make_history make_displays

make_history: $(SRC1)
	$(LD) $(SRC1) $(RLXX) -o make_history
	@echo '-> make_history created!'

make_displays: $(SRC2)
	$(LD) $(SRC2) $(RLXX) -o make_displays
	@echo '-> make_displays created!'

clean:
		@rm -f make_displays make_history
		@rm -f *~
		@rm -f *.o
		@rm -f *.x
