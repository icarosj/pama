# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -g  -fopenmp
OPT_LEVEL = -O3

SUITE = print bc bfs cc pr pr_atomic sssp tc converter

.PHONY: all
#all: $(SUITE)

all: pr print

debug: CXX_FLAGS += -DDEBUG 
debug: OPT_LEVEL = -O0
debug: all



% : %.cc *.h
	$(CXX) $(CXX_FLAGS) $(OPT_LEVEL) $< -o $@$

.PHONY: clean
clean:
	rm -f $(SUITE)
