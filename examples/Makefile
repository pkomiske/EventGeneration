CXX = g++
ifeq ($(shell uname), Darwin)
  CXX = clang++
endif

CXXFLAGS = -O3 -Wall -std=c++14 -g

FASTJETINC = $(shell fastjet-config --cxxflags)
PYTHIA8INC = $(shell pythia8-config --cxxflags)
INCLUDE   += $(FASTJETINC) $(PYTHIA8INC) -I../include

FASTJETLIB = $(shell fastjet-config --libs)
PYTHIA8LIB = $(shell pythia8-config --libs)
LIBRARIES += $(PYTHIA8LIB) $(FASTJETLIB) -L../lib -lEventGenerator -lboost_program_options

.PHONY: all
all:
	make events

events: Events.o
	$(CXX) -o $@ $^ $(LIBRARIES) $(CXXFLAGS)

Events.o: Events.cc
	$(CXX) -o $@ $< -c $(INCLUDE) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -fv events *.o
