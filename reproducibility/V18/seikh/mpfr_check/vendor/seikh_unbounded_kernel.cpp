// Exploratory extension only. Does not change any V14 file.
// Original entropy weights recomputed, t=1/q in [0,1/4].
#include "inherited_mpfr_interval_kernel.cpp"
extern "C" {int mpfr_log(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);}
static I log_I(const I&x){if(mpfr_cmp_ui(x.l.v,0)<=0)throw std::runtime_error("strictly positive log input required");I z;mpfr_log(z.l.v,x.l.v,MPFR_RNDD);mpfr_log(z.h.v,x.h.v,MPFR_RNDU);return z;}
static I exp_I(const I&x){I z;mpfr_exp(z.l.v,x.l.v,MPFR_RNDD);mpfr_exp(z.h.v,x.h.v,MPFR_RNDU);return z;}
static const int SEIKH[40]={
#include "data_values.inc"
};
extern "C" int seikh_t_scores(const double*lo,const double*hi,double tl,double th,double*out,double*wout){
try{
 if(!std::isfinite(tl)||!std::isfinite(th)||tl<0||th>0.25||tl>th)throw std::runtime_error("t domain");
 I t(tl,th);std::vector<I> z;
 for(int k=0;k<40;k++){
  I x=(lo&&hi)?I(lo[k],hi[k]):pown(rat(SEIKH[k],100),4);
  if(mpfr_cmp_ui(x.l.v,0)<=0||mpfr_cmp_ui(x.h.v,1)>=0)throw std::runtime_error("interior canonical input required");
  z.push_back(x);
 }
 std::vector<I> weights;I denominator(0);
 for(int j=0;j<5;j++){
  I a(0);
  for(int i=0;i<4;i++)for(int c=0;c<2;c++){I y=t*log_I(z[(i*5+j)*2+c]);a=a+y*exp_I(y);}
  I e=I(1)+rat(1,4)*a;
  // x log x >= -1/e for x in [0,1]; e_j >= 1-2/e > 1/4.
  // This analytic range intersection is independent of observed results.
  N quarter(0.25),one(1.0);
  if(mpfr_cmp(e.l.v,quarter.v)<0)e.l=quarter;
  if(mpfr_cmp(e.h.v,one.v)>0)e.h=one;
  if(mpfr_cmp(e.l.v,e.h.v)>0)throw std::runtime_error("invalid entropy interval");
  weights.push_back(e);denominator=denominator+e;
 }
 for(int j=0;j<5;j++){weights[j]=weights[j]/denominator;result(weights[j],wout+2*j);}
 for(int i=0;i<4;i++){
  I p(1),n(1);
  for(int j=0;j<5;j++){p=p*pow01(I(1)-z[(i*5+j)*2],weights[j]);n=n*pow01(z[(i*5+j)*2+1],weights[j]);}
  result(I(1)-p-n,out+2*i);
 }
 return 0;
}catch(const std::exception&e){ERROR=e.what();return 1;}}
