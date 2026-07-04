# SimdMFE

Fast RNA secondary-structure prediction by minimum free energy (Turner et al, 2004)
nearest-neighbour model with dangling ends scored on both sides of every helix.

SIMD-accelerated and header-only C++11 library. It runs on the command line, as an
embedded library, and in the browser.

- **Exact** the full Turner 2004 nearest-neighbour model in exact integer arithmetic.
- **Fast** SIMD-accelerated folding (ARM NEON, x86 SSE2/AVX2/AVX-512,
  WebAssembly SIMD) with an OpenMP batch mode.
- **Portable** header-only C++11 core, no dependencies; embed via CMake.
- **In the browser** an interactive WebAssembly folding and visualization. 

## Command line

```sh
make
# prints the sequence, structure, and MFE
echo GGGGAAAACCCC | ./simdmfe
# FASTA, or one sequence per line
./simdmfe < sequences.fa
```

Fold a batch in parallel:

```sh
make simdmfe_omp
OMP_NUM_THREADS=8 ./simdmfe_omp < sequences.fa
```

## As a library

```cmake
add_subdirectory(path/to/lib EXCLUDE_FROM_ALL)
target_link_libraries(myapp PRIVATE simdmfe::simdmfe)
```

```cpp
#include "fold_simd.h"

mfe::FoldSimd f;
f.setSeq("GGGGAAAACCCC");
// minimum free energy, in 0.01 kcal/mol
int mfe = f.fold();
// dot-bracket
std::string structure = f.traceback(mfe);
```

## In the browser

```sh
make wasm
python3 -m http.server -d web
# open http://localhost:8000
```

Paste a sequence to fold and view its structure, or upload a FASTA to fold in
parallel across web workers.

## Benchmark

Full ArchiveII (3966 reference structures;
[Sloma & Mathews, 2016](https://rnajournal.cshlp.org/content/22/12/1808)) against
RNAfold ([Lorenz et al., 2011](https://link.springer.com/article/10.1186/1748-7188-6-26))
2.7.2 `-d2`, on an NVIDIA DGX Spark (20-core ARM CPU):

| tool    | 1 thread | 20 threads | sensitivity | PPV   | F1    | F1 (macro) |
|---------|---------:|-----------:|------------:|------:|------:|-----------:|
| SimdMFE |   41.5 s |      3.3 s |       0.577 | 0.506 | 0.539 |      0.577 |
| RNAfold |  112.3 s |      9.4 s |       0.576 | 0.506 | 0.539 |      0.577 |

