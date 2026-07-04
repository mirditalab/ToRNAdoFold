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

EMCC ?= emcc
EMFLAGS ?= -O3 -std=c++17 -Isrc -lembind \
	-sMODULARIZE=1 -sEXPORT_NAME=createSimdMFE -sENVIRONMENT=web,worker \
	-sINITIAL_MEMORY=268435456

wasm: web/simdmfe.js web/simdmfe-simd.js
web/simdmfe.js: web/fold_wasm.cpp $(HDRS)
	$(EMCC) $(EMFLAGS) -o $@ web/fold_wasm.cpp
web/simdmfe-simd.js: web/fold_wasm.cpp $(HDRS)
	$(EMCC) $(EMFLAGS) -msimd128 -o $@ web/fold_wasm.cpp

clean:
	rm -f simdmfe simdmfe_omp simdmfe_noneon verify \
	  web/simdmfe.js web/simdmfe.wasm web/simdmfe-simd.js web/simdmfe-simd.wasm

.PHONY: all clean wasm
