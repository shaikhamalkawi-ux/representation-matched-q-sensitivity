// Compact parameter t=1/q, including the continuous t=0 limit on positive inputs.
// All original interval primitives and both source-qualified H/A merges retained.
#include "continuous_q_kernel.cpp"

static I positive_power_closed(const I&a,const I&t) {
    if(mpfr_cmp_ui(a.l.v,0)<=0 || mpfr_cmp_ui(a.h.v,1)>0 ||
       mpfr_cmp_ui(t.l.v,0)<0 || mpfr_cmp(t.l.v,t.h.v)>0)
        throw std::runtime_error("positive base in (0,1] and nonnegative exponent required");
    I z;
    // x^t increases in x and decreases in t for 0<x<=1, t>=0.
    mpfr_pow(z.l.v,a.l.v,t.h.v,MPFR_RNDD);
    mpfr_pow(z.h.v,a.h.v,t.l.v,MPFR_RNDU);
    return unit(z);
}

static I distance_compact(const Pair&a,const Pair&b,const I&t) {
    I d0=absv(positive_power_closed(a[0],t)-positive_power_closed(b[0],t));
    I d1=absv(positive_power_closed(a[1],t)-positive_power_closed(b[1],t));
    return unit(root(rat(1,2)*(pown(d0,3)+pown(d1,3)),3));
}

static std::vector<I> weights_compact(const std::vector<Pair>&v,
                                     const std::vector<I>&base,const I&t) {
    std::vector<I> z; I den(0);
    for(size_t i=0;i<v.size();i++) {
        I support(1);
        for(size_t j=0;j<v.size();j++) if(i!=j)
            support=support+(I(1)-distance_compact(v[i],v[j],t));
        z.push_back(base[i]*support); den=den+z.back();
    }
    for(auto&x:z) x=x/den;
    return z;
}

static I score_compact(const std::vector<Pair>&z,int a,const I&t,
                       int merge_mode,bool limit_only=false) {
    std::vector<I> ew={rat(30,100),rat(22,100),rat(28,100),rat(20,100)};
    std::vector<I> cw={rat(20,100),rat(20,100),rat(15,100),rat(25,100),rat(20,100)};
    std::vector<Pair> coll;
    for(int c=0;c<5;c++) {
        std::vector<Pair> vals;
        for(int e=0;e<4;e++) vals.push_back(z[e*25+a*5+c]);
        auto w=limit_only?ew:weights_compact(vals,ew,t);
        coll.push_back(expert_average(vals,w));
    }
    auto w=limit_only?cw:weights_compact(coll,cw,t);
    Pair f=merge_criteria(coll,w,merge_mode);
    return f[0]-f[1];
}

extern "C" int kernel_scores_compact(const double*lo,const double*hi,
                                      double tlo,double thi,int merge_mode,double*out) {
    try {
        if(!std::isfinite(tlo)||!std::isfinite(thi)||tlo<0||thi>1||tlo>thi)
            throw std::runtime_error("0<=tlo<=thi<=1 required");
        auto z=state(lo,hi); I t(tlo,thi);
        for(int a=0;a<5;a++) result(score_compact(z,a,t,merge_mode),out+2*a);
        return 0;
    } catch(const std::exception&e) {ERROR=e.what();return 1;}
}

extern "C" int kernel_scores_limit(const double*lo,const double*hi,
                                    int merge_mode,double*out) {
    try {
        auto z=state(lo,hi);
        for(int a=0;a<5;a++) result(score_compact(z,a,I(0),merge_mode,true),out+2*a);
        return 0;
    } catch(const std::exception&e) {ERROR=e.what();return 1;}
}
