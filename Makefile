CXX ?= clang++
CXXFLAGS ?= -O3 -std=c++11 -Wall
# OpenMP: plain -fopenmp on Linux (clang/gcc); Apple clang needs the libomp form.
ifeq ($(shell uname -s),Darwin)
OMPFLAGS ?= -Xpreprocessor -fopenmp -lomp
else
OMPFLAGS ?= -fopenmp
endif

HDRS = src/fold_simd.h src/energy.h src/t2004.h

all: simdmfe verify

simdmfe: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

# OpenMP batch build (needs libomp: brew install libomp)
simdmfe_omp: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ src/main.cpp

# NEON disabled, isolates the scalar constant factor
simdmfe_noneon: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -DDISABLE_NEON -o $@ src/main.cpp

verify: src/verify.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/verify.cpp

clean:
	rm -f simdmfe simdmfe_omp simdmfe_noneon verify

.PHONY: all clean
