// Exact Turner 2004 / -d2 energy model, using the parsed t2004 tables.
// Conventions reconstructed from the standard nearest-neighbour model:
//  bases A,C,G,U = 0..3 (N=4); table base index = base+1 (N->0).
//  pair types CG,GC,GU,UG,AU,UA,NN = 1..7; table pair index = type-1.
#pragma once
#include "t2004.h"
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace en {
using t2004::INF;

inline int bidx(int b) { return b < 4 ? b + 1 : 0; }   // N=0,A=1,C=2,G=3,U=4
inline int pidx(int pt) { return pt - 1; }             // type 1..7 -> 0..6

inline int pairType(int a, int b) {
    static const int T[4][4] = {
        {0,0,0,5},{0,0,1,0},{0,2,0,3},{6,0,4,0}
    };
    if (a<0||a>3||b<0||b>3) return 0;
    return T[a][b];
}
inline int auPen(int pt) { return (pt==1||pt==2) ? 0 : (pt==0?0:t2004::TerminalAU); }

struct EM {
    const int* s = nullptr;   // base codes 0..3 (4=N)
    int n = 0;
    std::string seq;

    void set(const std::string& str, const int* bases, int len) { seq=str; s=bases; n=len; }
    int pt(int i,int j) const { return pairType(s[i],s[j]); }

    int hpInit(int L) const {
        if (L<=30) return t2004::hairpin[L];
        return (int)std::lround(t2004::hairpin[30] + t2004::lxc*std::log(L/30.0));
    }
    int intInit(int L) const {
        if (L<=30) return t2004::internal[L];
        return (int)std::lround(t2004::internal[30] + t2004::lxc*std::log(L/30.0));
    }
    int blgInit(int L) const {
        if (L<=30) return t2004::bulge[L];
        return (int)std::lround(t2004::bulge[30] + t2004::lxc*std::log(L/30.0));
    }

    int eHairpin(int i,int j) const {
        int p = pt(i,j); if(!p) return INF;
        int L = j-i-1; if(L<3) return INF;
        // Special loops (tetra/tri/hexa): the tabulated value is the complete
        // loop free energy and replaces the generic estimate. Allocation-free.
        const char* sp = seq.data() + i;
        if (L==3){ for(int k=0;k<t2004::tri_n;++k)
                if(!std::memcmp(sp,t2004::tri[k].seq,5)) return t2004::tri[k].bonus; }
        else if (L==4){ for(int k=0;k<t2004::tetra_n;++k)
                if(!std::memcmp(sp,t2004::tetra[k].seq,6)) return t2004::tetra[k].bonus; }
        else if (L==6){ for(int k=0;k<t2004::hexa_n;++k)
                if(!std::memcmp(sp,t2004::hexa[k].seq,8)) return t2004::hexa[k].bonus; }
        int e;
        if (L==3) e = hpInit(3) + auPen(p);
        else e = hpInit(L) + t2004::mm_hairpin[pidx(p)][bidx(s[i+1])][bidx(s[j-1])];
        return e;
    }

    int eIntLoop(int i,int j,int p,int q) const {
        int t1 = pt(i,j);
        int t2f = pt(p,q);             // enclosed pair 5'->3'
        int t2r = pt(q,p);             // enclosed pair reversed (loop's view)
        if(!t1||!t2f) return INF;
        int n1=p-i-1, n2=j-q-1;
        // enclosed pair is indexed reversed (nearest-neighbour orientation)
        if (n1==0 && n2==0) return t2004::stack[pidx(t1)][pidx(t2r)];
        if (n1==0 || n2==0) {          // bulge
            int L=n1+n2, e=blgInit(L);
            if (L==1) e += t2004::stack[pidx(t1)][pidx(t2r)];
            else e += auPen(t1) + auPen(t2f);
            return e;
        }
        int si1=s[i+1], sj1=s[j-1], sp1=s[p-1], sq1=s[q+1];
        if (n1==1 && n2==1)
            return t2004::int11[pidx(t1)][pidx(t2r)][bidx(si1)][bidx(sj1)];
        if (n1==1 && n2==2)
            return t2004::int21[pidx(t1)][pidx(t2r)][bidx(si1)][bidx(sq1)][bidx(sj1)];
        if (n1==2 && n2==1)   // mirror of 1x2 (verified against RNAeval)
            return t2004::int21[pidx(t2r)][pidx(t1)][bidx(sq1)][bidx(si1)][bidx(sp1)];
        if (n1==2 && n2==2)
            return t2004::int22[pidx(t1)][pidx(t2r)][si1][sp1][sq1][sj1];
        // general
        int L=n1+n2;
        int e = intInit(L) + std::min(t2004::NINIO_max, std::abs(n1-n2)*t2004::NINIO_m);
        const int (*mm)[5][5];
        if (n1==1 || n2==1) mm = t2004::mm_internal_1n;
        else if ((n1==2&&n2==3)||(n1==3&&n2==2)) mm = t2004::mm_internal_23;
        else mm = t2004::mm_internal;
        e += mm[pidx(t1)][bidx(si1)][bidx(sj1)];
        e += mm[pidx(t2r)][bidx(sq1)][bidx(sp1)];
        return e;
    }

    // Multiloop interior branch (p,q): d2 terminal mismatch on both flanks.
    int mlStem(int p,int q) const {
        int t = pt(p,q); if(!t) return INF;
        int e = t2004::ML_intern + auPen(t);
        int mm5 = (p>0)? s[p-1] : 4;
        int mm3 = (q<n-1)? s[q+1] : 4;
        e += t2004::mm_multi[pidx(t)][bidx(mm5)][bidx(mm3)];
        return e;
    }
    // Closing pair of a multiloop, viewed from inside (reversed).
    int mlClose(int i,int j) const {
        int tr = pt(j,i); if(!tr) return INF;
        int e = t2004::ML_closing + t2004::ML_intern + auPen(tr);
        e += t2004::mm_multi[pidx(tr)][bidx(s[j-1])][bidx(s[i+1])];
        return e;
    }
    // Exterior stem (k,j): mismatch if both neighbors present, else single dangle.
    int extStem(int k,int j) const {
        int t = pt(k,j); if(!t) return INF;
        int e = auPen(t);
        int mm5 = (k>0)? s[k-1] : -1;
        int mm3 = (j<n-1)? s[j+1] : -1;
        if (mm5>=0 && mm3>=0) e += t2004::mm_exterior[pidx(t)][bidx(mm5)][bidx(mm3)];
        else if (mm5>=0) e += t2004::dangle5[pidx(t)][bidx(mm5)];
        else if (mm3>=0) e += t2004::dangle3[pidx(t)][bidx(mm3)];
        return e;
    }
    // Evaluate a dot-bracket structure the way the DP scores it (mirror of the
    // loop decomposition). Used to localise scoring bugs against RNAeval.
    int evalStructure(const std::string& db, bool verbose=false) const {
        int N=(int)db.size();
        std::vector<int> pr(N,-1), st;
        for (int i=0;i<N;++i){ char c=db[i];
            if (c=='('||c=='['||c=='{') st.push_back(i);
            else if (c==')'||c==']'||c=='}'){ if(st.empty()) return INF; pr[st.back()]=i; pr[i]=st.back(); st.pop_back(); } }
        int total=0;
        // exterior: top-level pairs
        for (int i=0;i<N;){ if(pr[i]>i){ total+=getVvalRegion(pr,i,pr[i],verbose); total+=extStem(i,pr[i]);
                if(verbose) fprintf(stderr,"ext stem (%d,%d) +%d\n",i,pr[i],extStem(i,pr[i])); i=pr[i]+1; } else ++i; }
        return total;
    }
    // energy of the loop closed by (i,j) plus everything inside it
    int getVvalRegion(const std::vector<int>& pr, int i, int j, bool verbose) const {
        // find children
        std::vector<std::pair<int,int>> ch;
        for (int k=i+1;k<j;){ if(pr[k]>k){ ch.push_back({k,pr[k]}); k=pr[k]+1; } else ++k; }
        int e=0;
        if (ch.empty()) { e=eHairpin(i,j); if(verbose) fprintf(stderr,"hairpin (%d,%d) %d\n",i,j,e); }
        else if (ch.size()==1){ auto[p,q]=ch[0]; e=eIntLoop(i,j,p,q)+getVvalRegion(pr,p,q,verbose);
            if(verbose) fprintf(stderr,"intloop (%d,%d;%d,%d) %d\n",i,j,p,q,eIntLoop(i,j,p,q)); }
        else { e=mlClose(i,j); if(verbose) fprintf(stderr,"ml close (%d,%d) +%d\n",i,j,mlClose(i,j));
            for (auto&c:ch){ e+=mlStem(c.first,c.second)+getVvalRegion(pr,c.first,c.second,verbose);
                if(verbose) fprintf(stderr,"  ml stem (%d,%d) +%d\n",c.first,c.second,mlStem(c.first,c.second)); } }
        return e;
    }
};

} // namespace en
