
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -O2 -Wall -Wextra $(ROOTCFLAGS) -I../main/include/
#-I/home/pitzl/GBL/V01-17-00/cpp/include

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pitzl/GBL/V01-17-00/cpp/lib/

scope53m: scope53m.cc
	g++ scope53m.cc $(CXXFLAGS) -o scope53m \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: scope53m'

tele: tele.cc
	g++ tele.cc $(CXXFLAGS) -fopenmp -o tele \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: tele'

rswfm: rswfm.cc
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) rswfm.cc -o rswfm $(ROOTLIBS)
	@echo 'done: rswfm'

tekwfm: tekwfm.cc
	g++ $(CXXFLAGS) tekwfm.cc -o tekwfm \
	$(ROOTLIBS)
	@echo 'done: tekwfm'

edg53: edg53.cc
	g++ edg53.cc $(CXXFLAGS) -o edg53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: edg53'

edgDiode: edgDiode.cc
	g++ edgDiode.cc $(CXXFLAGS) -o edgDiode \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: edgDiode'

edg53Diode: edg53Diode.cc
	g++ edg53Diode.cc $(CXXFLAGS) -o edg53Diode \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: edg53Diode'

ed53: ed53.cc
	g++ ed53.cc $(CXXFLAGS) -o ed53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: ed53'

evdt: evdt.cc
	g++ evdt.cc $(CXXFLAGS) -o evdt \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: evdt'

ed3d: ed3d.cc
	g++ ed3d.cc $(CXXFLAGS) -o ed3d \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: ed3d'

kink: kink.cc
	g++ kink.cc $(CXXFLAGS) -o kink \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: kink'
