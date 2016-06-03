CXXFLAGS = -Wall -O2 -I. -std=c++14
LDFLAGS = -lz -lm

.PHONY: all
all: resite

resite: CXXFLAGS += -fopenmp 

%: %.cc | kseq.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
