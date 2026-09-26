// Real-q and joint input/q extension of the unchanged archived MPFR kernel.
// The archive is retained verbatim beside this file. No original API is changed.
// merge_mode=0: archived Hamacher final merge; merge_mode=1: source Java qROFPAA.
#include "mpfr_interval_kernel.cpp"

static I distance_real(const Pair&a,const Pair&b,const I&inverse_q) {
    I d0=absv(pow01(a[0],inverse_q)-pow01(b[0],inverse_q));
    I d1=absv(pow01(a[1],inverse_q)-pow01(b[1],inverse_q));
    return unit(root(rat(1,2)*(pown(d0,3)+pown(d1,3)),3));
}

static std::vector<I> weights_real(const std::vector<Pair>&v,
                                  const std::vector<I>&base,const I&inverse_q) {
    std::vector<I> z; I den(0);
    for(size_t i=0;i<v.size();i++) {
        I t(1);
        for(size_t j=0;j<v.size();j++) if(i!=j)
            t=t+(I(1)-distance_real(v[i],v[j],inverse_q));
        z.push_back(base[i]*t); den=den+z.back();
    }
    for(auto&x:z) x=x/den;
    return z;
}

static Pair merge_criteria(const std::vector<Pair>&v,const std::vector<I>&w,
                           int merge_mode) {
    Pair p=partition({v[0],v[2],v[4]},{w[0],w[2],w[4]},{1,2,3});
    Pair t=partition({v[1],v[3]},{w[1],w[3]},{1,2});
    if(merge_mode==1) {
        // Canonical form of Qin's Java qROFPAA across the two partitions.
        return {unit(I(1)-root(unit((I(1)-p[0])*(I(1)-t[0])),2)),
                root(unit(p[1]*t[1]),2)};
    }
    if(merge_mode!=0) throw std::runtime_error("unknown final merge mode");
    Pair bp=B(p),bt=B(t);
    return can({root(bp[0]*bt[0],2),root(bp[1]*bt[1],2)});
}

static I score_real(const std::vector<Pair>&z,int a,const I&inverse_qe,
                    const I&inverse_qc,int merge_mode) {
    std::vector<I> ew={rat(30,100),rat(22,100),rat(28,100),rat(20,100)};
    std::vector<I> cw={rat(20,100),rat(20,100),rat(15,100),rat(25,100),rat(20,100)};
    std::vector<Pair> coll;
    for(int c=0;c<5;c++) {
        std::vector<Pair> vals;
        for(int e=0;e<4;e++) vals.push_back(z[e*25+a*5+c]);
        auto w=weights_real(vals,ew,inverse_qe);
        coll.push_back(expert_average(vals,w));
    }
    auto w=weights_real(coll,cw,inverse_qc);
    Pair f=merge_criteria(coll,w,merge_mode);
    return f[0]-f[1];
}

extern "C" int kernel_scores_real(const double*lo,const double*hi,
                                   double qe_lo,double qe_hi,
                                   double qc_lo,double qc_hi,
                                   int merge_mode,double*out) {
    try {
        if(!std::isfinite(qe_lo)||!std::isfinite(qe_hi)||
           !std::isfinite(qc_lo)||!std::isfinite(qc_hi)||
           qe_lo<1||qc_lo<1||qe_lo>qe_hi||qc_lo>qc_hi)
            throw std::runtime_error("finite ordered q intervals with lower >=1 required");
        I ie=I(1)/I(qe_lo,qe_hi),ic=I(1)/I(qc_lo,qc_hi);
        auto z=state(lo,hi);
        for(int a=0;a<5;a++) result(score_real(z,a,ie,ic,merge_mode),out+2*a);
        return 0;
    } catch(const std::exception&e) { ERROR=e.what();return 1; }
}
