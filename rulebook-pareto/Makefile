CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -march=native -Wall -Wextra -Wno-unused-parameter

bin/rbsearch: src/main.cpp src/graph.hpp src/rulebook.hpp src/solvers.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

.PHONY: clean
clean:
	rm -rf bin
