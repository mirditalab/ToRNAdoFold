// Regression self-test: the int16 fast path must equal the exact int32 path
// on sequences short enough to stay in 16-bit range (no fallback).
#include "fold_simd.h"
#include <cstdio>
#include <random>

int main() {
    std::mt19937 rng(42);
    const char* B = "ACGU";
    std::uniform_int_distribution<int> d(0, 3);
    int fails = 0, tot = 0;
    for (int L : {20, 35, 50, 80, 120, 200}) {
        for (int t = 0; t < 40; ++t) {
            std::string s;
            for (int i = 0; i < L; ++i) s += B[d(rng)];
            mfe::FoldSimd   f32; f32.setSeq(s); int e32 = f32.fold();
            mfe::FoldSimd16 f16; f16.setSeq(s); int e16 = f16.fold();
            tot++;
            if (e32 != e16) { if (fails < 8) printf("MISMATCH L=%d i32=%d i16=%d\n%s\n", L, e32, e16, s.c_str()); fails++; }
        }
    }
    printf("verify: %d/%d agree (%d mismatches)\n", tot - fails, tot, fails);
    return fails ? 1 : 0;
}
