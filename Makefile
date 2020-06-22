CXX = g++
ifeq ($(shell uname), Darwin)
  CXX = clang++
endif

CXXFLAGS = -O3 -Wall -std=c++14 -g

FASTJETINC = $(shell fastjet-config --cxxflags)
PYTHIA8INC = $(shell pythia8-config --cxxflags)
INCLUDE   += $(FASTJETINC) $(PYTHIA8INC) -Iinclude

# ensure lib directory exists
$(shell mkdir -p lib)

.PHONY: all
all: libEventGenerator.a
	make examples

libEventGenerator.a: src/EventGenerator.o src/JetMatcher.o
	ar crus lib/$@ $^

src/EventGenerator.o: src/EventGenerator.cc include/EventGenerator.hh include/JetMatcher.hh
	$(CXX) -o $@ $< -c $(INCLUDE) $(CXXFLAGS)

src/JetMatcher.o: src/JetMatcher.cc include/JetMatcher.hh
	$(CXX) -o $@ $< -c $(INCLUDE) $(CXXFLAGS)

.PHONY: examples
examples:
	make -C examples

.PHONY: clean
clean:
	rm -fv lib/* src/*.o
	make -C examples clean
