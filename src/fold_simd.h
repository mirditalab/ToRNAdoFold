// SimdMFE — wavefront (anti-diagonal) SIMD folder with exact Turner 2004 energy.
// Matrices stored diagonal-by-diagonal: D[s][i] holds cell (i, i+s). The O(N^3)
// multibranch bifurcation FMbif(i,j)=min_k FM(i,k)+FM1(k+1,j) becomes, per split
// span, a unit-stride vector add + min over i (NEON, 4x int32).
#pragma once
#include "energy.h"
#include <vector>
#include <string>
#include <algorithm>
#if (defined(__aarch64__) || defined(__ARM_NEON)) && !defined(DISABLE_NEON)
#include <arm_neon.h>
#define HAVE_NEON 1
#endif

namespace mfe {
using en::INF;

// Lower bound on any single internal-loop closure energy: the most stabilizing
// case is the best nearest-neighbour stack (-3.40 kcal/mol = -340 in centikcal
// units) across the Turner 2004 tables (stack/int11/int21/int22/general all
// >= -340). Used for an exact prune of the O(MAXLOOP^2) internal-loop scan.
static const int SINGLE_LOOP_LB = -340;

struct FoldSimd {
    int n = 0;
    std::vector<int> b;
    std::string seq;
    en::EM em;
    std::vector<std::vector<int32_t>> V, FM, FM1, FMbif;
    std::vector<int32_t> Erow;
    static const int PAD = 8;

    void setSeq(const std::string& str) {
        seq = str;
        n = (int)str.size();
        b.resize(n);
        for (int i = 0; i < n; ++i) {
            char c = str[i];
            b[i] = (c == 'A')   ? 0
                   : (c == 'C') ? 1
                   : (c == 'G') ? 2
                                : ((c == 'U' || c == 'T') ? 3 : 4);
        }
        em.set(seq, b.data(), n);
    }
    int pt(int i, int j) const {
        return em.pt(i, j);
    }
    int eHairpin(int i, int j) const {
        return em.eHairpin(i, j);
    }
    int eIntLoop(int i, int j, int p, int q) const {
        return em.eIntLoop(i, j, p, q);
    }
    int mlStem(int p, int q) const {
        return em.mlStem(p, q);
    }
    int extStem(int k, int j) const {
        return em.extStem(k, j);
    }
    int mlClose(int i, int j) const {
        return em.mlClose(i, j);
    }

    void bifKernel(int s) {
        int m = n - s;
        int32_t* out = FMbif[s].data();
        for (int i = 0; i < m; ++i) {
            out[i] = INF;
        }
        for (int a = 0; a < s; ++a) {
            int bsp = s - 1 - a;
            const int32_t* pa = FM[a].data();
            const int32_t* pb = FM1[bsp].data() + (a + 1);
            int i = 0;
#ifdef HAVE_NEON
            for (; i + 4 <= m; i += 4) {
                int32x4_t va = vld1q_s32(pa + i), vb = vld1q_s32(pb + i);
                int32x4_t sum = vaddq_s32(va, vb);
                sum = vminq_s32(sum, vdupq_n_s32(2 * INF));
                int32x4_t cur = vld1q_s32(out + i);
                vst1q_s32(out + i, vminq_s32(cur, sum));
            }
#endif
            // scalar tail (NEON remainder), or the whole span when NEON is off
            for (; i < m; ++i) {
                int va = pa[i], vb = pb[i];
                if (va >= INF || vb >= INF) {
                    continue;
                }
                int sm = va + vb;
                if (sm < out[i]) {
                    out[i] = sm;
                }
            }
        }
    }

    int fold() {
        V.assign(n, {});
        FM.assign(n, {});
        FM1.assign(n, {});
        FMbif.assign(n, {});
        for (int s = 0; s < n; ++s) {
            int m = n - s;
            V[s].assign(m + PAD, INF);
            FM[s].assign(m + PAD, INF);
            FM1[s].assign(m + PAD, INF);
            FMbif[s].assign(m + PAD, INF);
        }
        const int MLB = t2004::ML_BASE;
        for (int s = 1; s < n; ++s) {
            int m = n - s;
            for (int i = 0; i < m; ++i) {
                int j = i + s, p = pt(i, j), best = INF;
                if (p && s >= 4) {
                    best = eHairpin(i, j);
                    en::EM::IntPre pre = em.intPre(i, j); // (i,j)-invariant, hoisted
                    int maxp = std::min(j - 1, i + 31);
                    for (int a = i + 1; a <= maxp; ++a) {
                        int n1 = a - i - 1;
                        int minq = std::max(a + 1, j - 1 - (30 - n1));
                        for (int c = j - 1; c >= minq; --c) {
                            // minq bounds n1+n2 <= 30, so no explicit MAXLOOP check needed
                            int v2 = V[c - a][a];
                            if (v2 >= INF) {
                                continue;
                            }
                            // Exact lower-bound prune: eIntLoop >= SINGLE_LOOP_LB,
                            // so if v2 + SINGLE_LOOP_LB already >= best this enclosed
                            // pair cannot improve V(i,j); skip its energy evaluation.
                            if (v2 + SINGLE_LOOP_LB >= best) {
                                continue;
                            }
                            int e = em.eIntLoop(pre, i, j, a, c);
                            if (e < INF && e + v2 < best) {
                                best = e + v2;
                            }
                        }
                    }
                    if (s - 2 >= 1) {
                        int bif = FMbif[s - 2][i + 1];
                        if (bif < INF) {
                            int cl = mlClose(i, j);
                            if (bif + cl < best) {
                                best = bif + cl;
                            }
                        }
                    }
                }
                V[s][i] = best;
            }
            for (int i = 0; i < m; ++i) {
                int j = i + s, f1 = INF;
                if (pt(i, j) && V[s][i] < INF) {
                    f1 = V[s][i] + mlStem(i, j);
                }
                int ext = FM1[s - 1][i];
                if (ext < INF && ext + MLB < f1) {
                    f1 = ext + MLB;
                }
                FM1[s][i] = f1;
            }
            bifKernel(s);
            for (int i = 0; i < m; ++i) {
                int fm = FM1[s][i];
                int u = FM[s - 1][i + 1];
                if (u < INF && u + MLB < fm) {
                    fm = u + MLB;
                }
                if (FMbif[s][i] < fm) {
                    fm = FMbif[s][i];
                }
                FM[s][i] = fm;
            }
        }
        Erow.assign(n, 0);
        for (int j = 0; j < n; ++j) {
            int e = (j > 0) ? Erow[j - 1] : 0;
            for (int k = 0; k <= j; ++k) {
                int v = V[j - k][k];
                if (v >= INF) {
                    continue;
                }
                int pre = (k > 0) ? Erow[k - 1] : 0;
                int cand = pre + v + extStem(k, j);
                if (cand < e) {
                    e = cand;
                }
            }
            Erow[j] = e;
        }
        return n ? Erow[n - 1] : 0;
    }

    int getV(int i, int j) const {
        return V[j - i][i];
    }
    int getFM(int i, int j) const {
        return FM[j - i][i];
    }
    int getFM1(int i, int j) const {
        return FM1[j - i][i];
    }
    int getFMbif(int i, int j) const {
        return FMbif[j - i][i];
    }

    std::string traceback(int) {
        std::string dot(n, '.');
        int j = n - 1;
        while (j >= 0) {
            int e = Erow[j];
            if (j > 0 && Erow[j - 1] == e) {
                --j;
                continue;
            }
            bool found = false;
            for (int k = 0; k <= j; ++k) {
                int v = getV(k, j);
                if (v >= INF) {
                    continue;
                }
                int pre = (k > 0) ? Erow[k - 1] : 0;
                if (pre + v + extStem(k, j) == e) {
                    dot[k] = '(';
                    dot[j] = ')';
                    traceV(k, j, dot);
                    j = k - 1;
                    found = true;
                    break;
                }
            }
            if (!found) {
                --j;
            }
        }
        return dot;
    }
    void traceV(int i, int j, std::string& dot) {
        int v = getV(i, j);
        if (eHairpin(i, j) == v) {
            return;
        }
        int maxp = std::min(j - 1, i + 31);
        for (int a = i + 1; a <= maxp; ++a) {
            int n1 = a - i - 1;
            int minq = std::max(a + 1, j - 1 - (30 - n1));
            for (int c = j - 1; c >= minq; --c) {
                // minq bounds n1+n2 <= 30 (see fold())
                int v2 = getV(a, c);
                if (v2 >= INF) {
                    continue;
                }
                int e = eIntLoop(i, j, a, c);
                if (e < INF && e + v2 == v) {
                    dot[a] = '(';
                    dot[c] = ')';
                    traceV(a, c, dot);
                    return;
                }
            }
        }
        if ((j - 1) - (i + 1) >= 1) {
            int cl = mlClose(i, j);
            if (getFMbif(i + 1, j - 1) + cl == v) {
                traceFMbif(i + 1, j - 1, dot);
                return;
            }
        }
    }
    void traceFMbif(int i, int j, std::string& dot) {
        int val = getFMbif(i, j);
        for (int k = i; k < j; ++k) {
            int a = getFM(i, k), c = getFM1(k + 1, j);
            if (a >= INF || c >= INF) {
                continue;
            }
            if (a + c == val) {
                traceFM(i, k, dot);
                traceFM1(k + 1, j, dot);
                return;
            }
        }
    }
    void traceFM1(int i, int j, std::string& dot) {
        int f1 = getFM1(i, j);
        if (pt(i, j) && getV(i, j) < INF && getV(i, j) + mlStem(i, j) == f1) {
            dot[i] = '(';
            dot[j] = ')';
            traceV(i, j, dot);
            return;
        }
        if (j - 1 >= i && getFM1(i, j - 1) < INF && getFM1(i, j - 1) + t2004::ML_BASE == f1) {
            traceFM1(i, j - 1, dot);
        }
    }
    void traceFM(int i, int j, std::string& dot) {
        int fm = getFM(i, j);
        if (getFM1(i, j) == fm) {
            traceFM1(i, j, dot);
            return;
        }
        if (i + 1 <= j && getFM(i + 1, j) < INF && getFM(i + 1, j) + t2004::ML_BASE == fm) {
            traceFM(i + 1, j, dot);
            return;
        }
        traceFMbif(i, j, dot);
    }
};

} // namespace mfe
