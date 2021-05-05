
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -O2 -Wall -Wextra $(ROOTCFLAGS) -I/eudaq/eudaq/include/

scope53m: scope53m.cc
	g++ $(CXXFLAGS) scope53m.cc -o scope53m \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: scope53m'

scopes: scopes_2017.cc
	g++ $(CXXFLAGS) scopes_2017.cc -o scopes \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: scopes (2017 version)'

old_scopes: scopes.cc
	g++ $(CXXFLAGS) scopes.cc -o scopes \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: scopes'

rswfm: rswfm.cc
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) rswfm.cc -o rswfm $(ROOTLIBS)
	@echo 'done: rswfm'

tekwfm: tekwfm.cc
	g++ $(CXXFLAGS) tekwfm.cc -o tekwfm \
	$(ROOTLIBS)
	@echo 'done: tekwfm'

edgDiode: edgDiode.cc
	g++ edgDiode.cc $(CXXFLAGS) -o edgDiode \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: edgDiode'

edg53Diode: edg53Diode.cc
	g++ edg53Diode.cc $(CXXFLAGS) -o edg53Diode \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: edg53Diode'

edg53: edg53.cc
	g++ $(CXXFLAGS) edg53.cc -o edg53 \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: edg53'

scope53: scope53.cc
	g++ $(CXXFLAGS) scope53.cc -o scope53 \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: scope53'

tele: tele.cc
	g++ tele.cc $(CXXFLAGS) -fopenmp -o tele \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: tele'

ed53: ed53.cc
	g++ $(CXXFLAGS) ed53.cc -o ed53 \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: ed53'

evdt: evdt.cc
	g++ evdt.cc $(CXXFLAGS) -o evdt \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: evdt'

ed3d: ed3d.cc
	g++ ed3d.cc $(CXXFLAGS) -o ed3d \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: ed3d'

kink: kink.cc
	g++ kink.cc $(CXXFLAGS) -o kink \
	$(ROOTLIBS) -L/eudaq/eudaq/lib -lEUDAQ
	@echo 'done: kink'

