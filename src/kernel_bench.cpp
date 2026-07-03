// Microbenchmark: the bifurcation kernel min_a(FM[a][i]+FM1[b][i+a+1]) in
// int32 (4 NEON lanes) vs int16 (8 NEON lanes, saturating). Tests the
// "16-bit doubles SIMD throughput" claim on this M-series core.
#include <arm_neon.h>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <chrono>
#include <random>

static const int32_t INF32 = 1000000;
static const int16_t INF16 = 30000;

// Simulate a full triangular DP fill of the bifurcation kernel for size n.
template<class T> struct Bench;

template<> struct Bench<int32_t> {
    static double run(int n, int reps) {
        std::vector<std::vector<int32_t>> FM(n), FM1(n);
        std::mt19937 rng(7);
        for (int s=0;s<n;++s){FM[s].assign(n-s+8,0);FM1[s].assign(n-s+8,0);
            for(auto&x:FM[s])x=rng()%2000; for(auto&x:FM1[s])x=rng()%2000;}
        volatile int64_t sink=0;
        auto t0=std::chrono::high_resolution_clock::now();
        for(int r=0;r<reps;++r)
        for(int s=2;s<n;++s){
            int m=n-s; std::vector<int32_t> out(m+8,INF32);
            for(int a=0;a<s;++a){int bsp=s-1-a;
                const int32_t*pa=FM[a].data();const int32_t*pb=FM1[bsp].data()+(a+1);
                int i=0;for(;i+4<=m;i+=4){int32x4_t va=vld1q_s32(pa+i);int32x4_t vb=vld1q_s32(pb+i);
                    int32x4_t sm=vaddq_s32(va,vb);int32x4_t cu=vld1q_s32(out.data()+i);
                    vst1q_s32(out.data()+i,vminq_s32(cu,sm));}
                for(;i<m;++i){int v=pa[i]+pb[i];if(v<out[i])out[i]=v;}}
            sink+=out[0];
        }
        auto t1=std::chrono::high_resolution_clock::now();
        (void)sink;
        return std::chrono::duration<double,std::milli>(t1-t0).count();
    }
};

template<> struct Bench<int16_t> {
    static double run(int n, int reps) {
        std::vector<std::vector<int16_t>> FM(n), FM1(n);
        std::mt19937 rng(7);
        for (int s=0;s<n;++s){FM[s].assign(n-s+16,0);FM1[s].assign(n-s+16,0);
            for(auto&x:FM[s])x=rng()%2000; for(auto&x:FM1[s])x=rng()%2000;}
        volatile int64_t sink=0;
        auto t0=std::chrono::high_resolution_clock::now();
        for(int r=0;r<reps;++r)
        for(int s=2;s<n;++s){
            int m=n-s; std::vector<int16_t> out(m+16,INF16);
            for(int a=0;a<s;++a){int bsp=s-1-a;
                const int16_t*pa=FM[a].data();const int16_t*pb=FM1[bsp].data()+(a+1);
                int i=0;for(;i+8<=m;i+=8){int16x8_t va=vld1q_s16(pa+i);int16x8_t vb=vld1q_s16(pb+i);
                    int16x8_t sm=vqaddq_s16(va,vb);int16x8_t cu=vld1q_s16(out.data()+i);
                    vst1q_s16(out.data()+i,vminq_s16(cu,sm));}
                for(;i<m;++i){int v=pa[i]+pb[i];if(v<out[i])out[i]=v;}}
            sink+=out[0];
        }
        auto t1=std::chrono::high_resolution_clock::now();
        (void)sink;
        return std::chrono::duration<double,std::milli>(t1-t0).count();
    }
};

int main(){
    for(int n:{300,600,1000}){
        int reps = n<=300?20:(n<=600?6:2);
        double t32=Bench<int32_t>::run(n,reps)/reps;
        double t16=Bench<int16_t>::run(n,reps)/reps;
        printf("n=%-5d  int32=%7.2f ms  int16=%7.2f ms  throughput x%.2f\n",n,t32,t16,t32/t16);
    }
    return 0;
}
