// Cross-check: scalar fold.h (2D storage) vs SIMD fold_simd.h (diagonal+NEON)
// must agree on the MFE for random sequences. Both use the exact energy model.
#include "fold.h"
#include "fold_simd.h"
#include <cstdio>
#include <random>

int main() {
    std::mt19937 rng(42);
    const char* B = "ACGU";
    int fails = 0, tot = 0;
    std::uniform_int_distribution<int> d(0, 3);
    for (int L : {20, 35, 50, 80, 120, 200}) {
        for (int t = 0; t < 40; ++t) {
            std::string s;
            for (int i = 0; i < L; ++i) s += B[d(rng)];
            mfe::Fold fa; fa.setSeq(s); int ea = fa.fold();
            mfe::FoldSimd fb; fb.setSeq(s); int eb = fb.fold();
            tot++;
            if (ea != eb) {
                if (fails < 8) printf("MISMATCH L=%d scalar=%d simd=%d\n%s\n", L, ea, eb, s.c_str());
                fails++;
            }
        }
    }
    printf("verify: %d/%d agree (%d mismatches)\n", tot - fails, tot, fails);
    return fails ? 1 : 0;
}
