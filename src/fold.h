// SimdMFE core folder — Turner 2004, unconditional (-d2) dangles, affine
// multiloop/exterior scoring. Scalar reference recurrences with backtrace.
// Recurrences implemented from the SimdMFE design description only.
#pragma once
#include "energy.h"
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace mfe {

using en::INF;

struct Fold {
    int n;
    std::vector<int> b;                 // base indices 0..3 (4=N)
    std::string seq;
    // DP matrices flattened row-major [i*n+j], only i<=j used.
    std::vector<int> V, FM, FM1;        // V, multiloop-interior, single-branch
    std::vector<int> Erow;              // exterior prefix E(0,j), size n
    // Backtrace helpers recomputed on demand.

    static int idx(int i, int j, int n) { return i * n + j; }
    int I(int i, int j) const { return i * n + j; }

    en::EM em;
    void setSeq(const std::string& s) {
        seq = s; n = (int)s.size();
        b.resize(n);
        for (int i = 0; i < n; ++i) { char c=s[i];
            b[i] = (c=='A')?0:(c=='C')?1:(c=='G')?2:((c=='U'||c=='T')?3:4); }
        em.set(seq, b.data(), n);
    }

    int pt(int i, int j) const { return em.pt(i, j); }
    int eHairpin(int i, int j) const { return em.eHairpin(i, j); }
    int eIntLoop(int i, int j, int p, int q) const { return em.eIntLoop(i, j, p, q); }
    int multiBranch(int p, int q) const { return em.mlStem(p, q); }
    int mlClose(int i, int j) const { return em.mlClose(i, j); }

    // Exterior branch (k,j): delegates to the exact model.
    int extBranch(int k, int j) const {
        int t = pt(k, j);
        if (!t) return INF;
        return em.extStem(k, j);
    }

    static int addsat(int a, int c) {
        if (a >= INF || c >= INF) return INF;
        return a + c;
    }

    int fold() {
        V.assign((size_t)n * n, INF);
        FM.assign((size_t)n * n, INF);
        FM1.assign((size_t)n * n, INF);
        // Fill by increasing span.
        for (int span = 1; span < n; ++span) {
            for (int i = 0; i + span < n; ++i) {
                int j = i + span;
                // ---- V ----
                int p = pt(i, j);
                int best = INF;
                if (p && span >= 4) {
                    best = eHairpin(i, j);
                    // internal/bulge/stack
                    int maxp = std::min(j - 1, i + 31);
                    for (int a = i + 1; a <= maxp; ++a) {
                        int n1 = a - i - 1;
                        int minq = std::max(a + 1, j - 1 - (30 - n1));
                        for (int c = j - 1; c >= minq; --c) {
                            if ((a - i - 1) + (j - c - 1) > 30) continue;
                            int v2 = V[I(a, c)];
                            if (v2 >= INF) continue;
                            int e = eIntLoop(i, j, a, c);
                            if (e < INF) best = std::min(best, e + v2);
                        }
                    }
                    // multiloop: closing (i,j) around FM(i+1,k)+FM1(k+1,j-1)
                    int close = mlClose(i, j);
                    for (int k = i + 2; k <= j - 2; ++k) {
                        int a = FM[I(i + 1, k)];
                        int c = FM1[I(k + 1, j - 1)];
                        if (a >= INF || c >= INF) continue;
                        int cand = a + c + close;
                        if (cand < best) best = cand;
                    }
                }
                V[I(i, j)] = best;

                // ---- FM1: exactly one branch starting at i ----
                int f1 = INF;
                if (p) {
                    int br = V[I(i, j)];
                    if (br < INF) f1 = br + multiBranch(i, j);
                }
                if (j - 1 >= i) {
                    int ext = FM1[I(i, j - 1)];
                    if (ext < INF) f1 = std::min(f1, ext + t2004::ML_BASE);
                }
                FM1[I(i, j)] = f1;

                // ---- FM: >=1 branch ----
                int fm = FM1[I(i, j)];
                if (i + 1 <= j) {
                    int u = FM[I(i + 1, j)];
                    if (u < INF) fm = std::min(fm, u + t2004::ML_BASE);
                }
                for (int k = i; k < j; ++k) {
                    int a = FM[I(i, k)];
                    int c = FM1[I(k + 1, j)];
                    if (a >= INF || c >= INF) continue;
                    int cand = a + c;
                    if (cand < fm) fm = cand;
                }
                FM[I(i, j)] = fm;
            }
        }
        // ---- Exterior: E(0,j) ----
        Erow.assign(n, 0);
        for (int j = 0; j < n; ++j) {
            int e = (j > 0) ? Erow[j - 1] : 0;      // j unpaired
            for (int k = 0; k <= j; ++k) {
                int v = V[I(k, j)];
                if (v >= INF) continue;
                int pre = (k > 0) ? Erow[k - 1] : 0;
                int cand = pre + v + extBranch(k, j);
                if (cand < e) e = cand;
            }
            Erow[j] = e;
        }
        return n ? Erow[n - 1] : 0;
    }

    // ---------- Backtrace ----------
    std::string traceback(int mfeVal) {
        std::string dot(n, '.');
        // exterior trace
        traceExt(0, n - 1, dot);
        return dot;
    }

    void traceExt(int lo, int hi, std::string& dot) {
        // Reconstruct E over [lo..hi] via the same recurrence (lo assumed 0).
        int j = hi;
        // recompute Erow already available as Erow[]; walk it.
        while (j >= lo) {
            int e = Erow[j];
            if (j > lo && Erow[j - 1] == e) { --j; continue; }     // j unpaired
            bool found = false;
            for (int k = lo; k <= j; ++k) {
                int v = V[I(k, j)];
                if (v >= INF) continue;
                int pre = (k > lo) ? Erow[k - 1] : 0;
                if (pre + v + extBranch(k, j) == e) {
                    dot[k] = '('; dot[j] = ')';
                    traceV(k, j, dot);
                    j = k - 1; found = true; break;
                }
            }
            if (!found) --j;
        }
    }

    void traceV(int i, int j, std::string& dot) {
        int v = V[I(i, j)];
        int p = pt(i, j);
        // hairpin?
        if (eHairpin(i, j) == v) return;
        // internal/bulge/stack?
        int maxp = std::min(j - 1, i + 31);
        for (int a = i + 1; a <= maxp; ++a) {
            int n1 = a - i - 1;
            int minq = std::max(a + 1, j - 1 - (30 - n1));
            for (int c = j - 1; c >= minq; --c) {
                if ((a - i - 1) + (j - c - 1) > 30) continue;
                int v2 = V[I(a, c)];
                if (v2 >= INF) continue;
                int e = eIntLoop(i, j, a, c);
                if (e < INF && e + v2 == v) {
                    dot[a] = '('; dot[c] = ')';
                    traceV(a, c, dot);
                    return;
                }
            }
        }
        // multiloop
        int close = mlClose(i, j);
        for (int k = i + 2; k <= j - 2; ++k) {
            int a = FM[I(i + 1, k)];
            int c = FM1[I(k + 1, j - 1)];
            if (a >= INF || c >= INF) continue;
            if (a + c + close == v) {
                traceFM(i + 1, k, dot);
                traceFM1(k + 1, j - 1, dot);
                return;
            }
        }
    }

    void traceFM1(int i, int j, std::string& dot) {
        int f1 = FM1[I(i, j)];
        int p = pt(i, j);
        if (p && V[I(i, j)] < INF && V[I(i, j)] + multiBranch(i, j) == f1) {
            dot[i] = '('; dot[j] = ')';
            traceV(i, j, dot);
            return;
        }
        if (j - 1 >= i && FM1[I(i, j - 1)] < INF && FM1[I(i, j - 1)] + t2004::ML_BASE == f1) {
            traceFM1(i, j - 1, dot);
        }
    }

    void traceFM(int i, int j, std::string& dot) {
        int fm = FM[I(i, j)];
        if (FM1[I(i, j)] == fm) { traceFM1(i, j, dot); return; }
        if (i + 1 <= j && FM[I(i + 1, j)] < INF && FM[I(i + 1, j)] + t2004::ML_BASE == fm) {
            traceFM(i + 1, j, dot); return;
        }
        for (int k = i; k < j; ++k) {
            int a = FM[I(i, k)];
            int c = FM1[I(k + 1, j)];
            if (a >= INF || c >= INF) continue;
            if (a + c == fm) { traceFM(i, k, dot); traceFM1(k + 1, j, dot); return; }
        }
    }
};

} // namespace mfe
