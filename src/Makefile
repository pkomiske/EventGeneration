CXXFLAGS = -O3 -Wall -std=c++14 -g
SRCS = $(shell ls *.cc)

FASTJETINC = $(shell fastjet-config --cxxflags)
PYTHIA8INC = $(shell pythia8-config --cxxflags)
INCLUDES  += $(FASTJETINC) $(PYTHIA8INC) -I../include

# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#combine
DEPDIR = .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCLUDES) -c
%.o : %.cc
%.o : %.cc $(DEPDIR)/%.d | $(DEPDIR)
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

.PHONY: all clean
all: $(SRCS:.cc=.o)

clean:
	rm -fv *.o

$(DEPDIR): ; @mkdir -p $@
DEPFILES := $(SRCS:%.cc=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
