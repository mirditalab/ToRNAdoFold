// Turner 2004 nearest-neighbour parameters, encoded from published values.
// Energies are integers in units of 0.1 kcal/mol (decikcal).
// Written for SimdMFE — no ViennaRNA source consulted; parameter *values*
// are the canonical Turner 2004 physical constants.
#pragma once
#include <cstdint>
#include <string>
#include <array>
#include <climits>

namespace tp {

// Bases: A=0 C=1 G=2 U=3 ; N/gap = 4
// Pair types (ViennaRNA convention): 0 = no pair, 1=CG 2=GC 3=GU 4=UG 5=AU 6=UA
static const int INF = 1000000;      // scalar "infinity"

// pair_type[a][b] over bases 0..3
inline int pairType(int a, int b) {
    // CG,GC,GU,UG,AU,UA
    static const int T[4][4] = {
        //     A  C  G  U
        /*A*/ {0, 0, 0, 5},   // AU
        /*C*/ {0, 0, 1, 0},   // CG
        /*G*/ {0, 2, 0, 3},   // GC, GU
        /*U*/ {6, 0, 4, 0},   // UA, UG
    };
    if (a < 0 || a > 3 || b < 0 || b > 3) return 0;
    return T[a][b];
}

// Stacking energies stack[outer][inner], both are pair types 1..6.
// Rows/cols order: CG GC GU UG AU UA  (index 1..6). Index 0 unused.
static const int stack[7][7] = {
/*        -     CG    GC    GU    UG    AU    UA */
/* -  */ {INF,  INF,  INF,  INF,  INF,  INF,  INF},
/* CG */ {INF, -240, -330, -210, -140, -210, -210},
/* GC */ {INF, -330, -340, -250, -150, -220, -240},
/* GU */ {INF, -210, -250,  130,  -50, -140, -130},
/* UG */ {INF, -140, -150,  -50,   30,  -60, -100},
/* AU */ {INF, -210, -220, -140,  -60, -110,  -90},
/* UA */ {INF, -210, -240, -130, -100,  -90, -130},
};

// Hairpin loop initiation, index = loop length (unpaired bases), 0..30.
static const int hairpin_init[31] = {
    INF, INF, INF, 540, 560, 570, 540, 600, 550, 640, 650,
    660, 670, 680, 690, 690, 700, 710, 710, 720, 720,
    730, 730, 740, 740, 750, 750, 750, 760, 760, 770
};
// Bulge loop initiation, index = loop length 0..30 (0 unused).
static const int bulge_init[31] = {
    INF, 380, 280, 320, 360, 400, 440, 460, 470, 480, 490,
    500, 510, 520, 530, 540, 540, 550, 550, 560, 570,
    570, 580, 580, 580, 590, 590, 600, 600, 600, 610
};
// Internal loop initiation, index = total unpaired 0..30 (0..3 unused).
static const int internal_init[31] = {
    INF, INF, INF, INF, 110, 200, 200, 210, 230, 240, 250,
    260, 270, 280, 290, 300, 300, 310, 310, 320, 330,
    330, 340, 340, 350, 350, 350, 360, 360, 370, 370
};

// Log extrapolation for loops > 30: dG(n) = dG(30) + lxc * ln(n/30).
// lxc = 1.07857764 kcal/mol -> in decikcal * ln factor applied at runtime.
static const double lxc = 107.856;   // decikcal

// Multiloop affine parameters (Turner 2004):
//   closing (a) = 3.4, per-branch (c) = 0.4, per-unpaired (b) = 0.0
static const int ML_closing = 340;
static const int ML_intern  = 40;    // per branch (incl. closing branch)
static const int ML_base    = 0;     // per unpaired base

// Terminal AU/GU penalty (non-GC helix end): 0.5 kcal/mol
static const int TerminalAU = 50;
inline int auPenalty(int pt) { return (pt == 0 || pt == 1 || pt == 2) ? 0 : TerminalAU; }

// Ninio asymmetry: 0.5 / nt, capped at 3.0
static const int NINIO_m = 50;
static const int NINIO_max = 300;

// Lone-pair / helix-extension: ViennaRNA does not penalise lone pairs by default
// under the standard model; the spec's "suppress lone pairs" is handled by the
// noLP flag (off by default here to match RNAfold default).

// ---- Terminal mismatch tables (approximate Turner 2004 tstack values) ----
// Indexed [pairtype 0..6][b1 0..3][b2 0..3]; b1 is 3' of the closing pair's
// first strand side, b2 is 5' neighbour on the other. These are approximations
// reconstructed from the general Turner mismatch structure.
// mismatch_hairpin: stronger stabilisation.
static int mm_hairpin[7][4][4];
static int mm_multi[7][4][4];
static int mm_exterior[7][4][4];
static int mm_internal[7][4][4];
static int dangle5_t[7][4];
static int dangle3_t[7][4];

// Special hairpin loops (tetra/tri/hexa) — a subset of the strongly stabilising
// motifs, matched by exact loop-with-closing-pair string.
struct SpecialHP { const char* seq; int bonus; };
// Known stable tetra/hexaloops (closing pair included in the string), all
// stabilising (negative). GNRA and UNCG families dominate.
static const SpecialHP tetraloops[] = {
    {"GGGGAC", -300}, {"GGUGAC", -300}, {"CGAAAG", -300}, {"GGAGAC", -300},
    {"CGCAAG", -300}, {"GGAAAC", -300}, {"CGGAAG", -300}, {"CUUCGG", -300},
    {"CGUGAG", -300}, {"CGAAGG", -250}, {"CUACGG", -150}, {"GGCAAC", -300},
    {"CGCGAG", -300}, {"UGAGAG", -350}, {"CGAGAG", -200}, {"AGAAAU", -150},
    {"CGUAAG", -150}, {"CUAACG", -150}, {"UGAAAG", -150}, {"GGAAGC", -150},
    {"GGGAAC", -150}, {"UGAAGG", -150}, {"UGGAAG", -150}, {"GCAAGC", -150},
    {"CGCAAG", -300}, {"UGCGAG", -150}, {"GGAGGC", -150}, {"UGUAAG", -150},
};
static const SpecialHP triloops[] = {
    {"CAACG", 680}, {"GUUAC", 690},
};

inline int baseIdx(char c) {
    switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1;
        case 'G': case 'g': return 2; case 'U': case 'u': case 'T': case 't': return 3; }
    return 4;
}

// Fill approximate mismatch/dangle tables. Called once at startup.
inline void initTables() {
    // Base stabilisation by closing-pair strength.
    // GC/CG closings ~ -1.1..-1.5; GU/UG ~ -0.7; AU/UA ~ -0.8.
    static const int base_by_pt[7] = {0, -110, -110, -70, -70, -80, -80};
    for (int pt = 1; pt <= 6; ++pt) {
        for (int b1 = 0; b1 < 4; ++b1) for (int b2 = 0; b2 < 4; ++b2) {
            int v = base_by_pt[pt];
            // G-G, G-A, A-G mismatches are extra-stabilising in hairpins.
            int hp = v;
            if (b1 == 2 && b2 == 2) hp -= 40;           // GG
            if (b1 == 2 || b2 == 2) hp -= 20;           // any G
            if (b1 == 0 && b2 == 0) hp -= 10;           // AA
            mm_hairpin[pt][b1][b2] = hp;
            // internal-loop mismatch: slightly weaker
            int in = v / 2 - 10;
            if (b1 == 2 && b2 == 2) in -= 20;
            if (b1 == 0 && b2 == 2) in -= 30;           // AG
            if (b1 == 2 && b2 == 0) in -= 30;           // GA
            if (b1 == 3 && b2 == 3) in -= 10;           // UU
            mm_internal[pt][b1][b2] = in;
            // multi / exterior mismatch (tstack): weaker, mostly 0..-0.9
            int mx = (v * 3) / 4;
            if (b1 == 2 || b2 == 2) mx -= 10;
            mm_multi[pt][b1][b2] = mx;
            mm_exterior[pt][b1][b2] = mx;
        }
    }
    // Dangles (single-base stacking at helix ends). Approximate: 3' dangles
    // stronger than 5'.
    static const int d5_by_pt[7] = {0, -50, -20, -30, -20, -30, -30};
    static const int d3_by_pt[7] = {0, -110, -80, -90, -60, -80, -80};
    for (int pt = 1; pt <= 6; ++pt) for (int b = 0; b < 4; ++b) {
        int e5 = d5_by_pt[pt];
        int e3 = d3_by_pt[pt];
        if (b == 2) e3 -= 20;   // G dangle strong
        if (b == 0) e3 -= 10;
        dangle5_t[pt][b] = e5;
        dangle3_t[pt][b] = e3;
    }
}

// Convenience accessors that clamp N to a neutral value.
inline int mmHairpin(int pt, int b1, int b2) {
    if (pt==0) return 0; if (b1>3||b2>3) return 0; return mm_hairpin[pt][b1][b2]; }
inline int mmInternal(int pt, int b1, int b2) {
    if (pt==0) return 0; if (b1>3||b2>3) return 0; return mm_internal[pt][b1][b2]; }
inline int mmMulti(int pt, int b1, int b2) {
    if (pt==0) return 0; if (b1>3||b2>3) return 0; return mm_multi[pt][b1][b2]; }
inline int mmExterior(int pt, int b1, int b2) {
    if (pt==0) return 0; if (b1>3||b2>3) return 0; return mm_exterior[pt][b1][b2]; }
inline int dangle5(int pt, int b) { if (pt==0||b>3) return 0; return dangle5_t[pt][b]; }
inline int dangle3(int pt, int b) { if (pt==0||b>3) return 0; return dangle3_t[pt][b]; }

} // namespace tp
