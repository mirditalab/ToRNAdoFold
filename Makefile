CXX ?= clang++
CXXFLAGS ?= -O3 -std=c++17 -Wall
OMPFLAGS ?= -Xpreprocessor -fopenmp -lomp

all: simdmfe verify kbench

simdmfe: src/main.cpp src/fold_simd.h src/params.h
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

# OpenMP batch build (needs libomp: brew install libomp)
simdmfe_omp: src/main.cpp src/fold_simd.h src/params.h
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ src/main.cpp

# NEON disabled, isolates the scalar constant factor
simdmfe_noneon: src/main.cpp src/fold_simd.h src/params.h
	$(CXX) $(CXXFLAGS) -DDISABLE_NEON -o $@ src/main.cpp

verify: src/verify.cpp src/fold.h src/fold_simd.h src/params.h
	$(CXX) $(CXXFLAGS) -o $@ src/verify.cpp

kbench: src/kernel_bench.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/kernel_bench.cpp

clean:
	rm -f simdmfe simdmfe_omp simdmfe_noneon verify kbench

.PHONY: all clean
