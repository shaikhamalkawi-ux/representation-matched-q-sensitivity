// V17 directed MPFR adapter. The inherited kernel is copied byte-for-byte.
// Scientific improvement: scalar g(y)=y exp(y) has its only interior extremum
// at y=-1; evaluate endpoints and include -exp(-1) when the interval crosses it.
#include "vendor/seikh_unbounded_kernel.cpp"

static I yexp_range(const I& y) {
 if (mpfr_cmp_ui(y.h.v,0)>0) throw std::runtime_error("yexp requires y<=0");
 I lp(y.l,y.l),hp(y.h,y.h);
 I lv=lp*exp_I(lp),hv=hp*exp_I(hp);
 I out=lv;
 if (mpfr_cmp(hv.l.v,out.l.v)<0) out.l=hv.l;
 if (mpfr_cmp(hv.h.v,out.h.v)>0) out.h=hv.h;
 N minus_one(-1.0);
 if(mpfr_cmp(y.l.v,minus_one.v)<=0 && mpfr_cmp(y.h.v,minus_one.v)>=0) {
  I critical=I(-1)*exp_I(I(-1));
  if(mpfr_cmp(critical.l.v,out.l.v)<0) out.l=critical.l;
 }
 return out;
}

extern "C" int v17_yexp(double lo,double hi,double*out) {
 try {
  if(!out || !std::isfinite(lo)||!std::isfinite(hi)||lo>hi) throw std::runtime_error("yexp domain");
  result(yexp_range(I(lo,hi)),out);return 0;
 }catch(const std::exception&e){ERROR=e.what();return 1;}
}

extern "C" int v17_components(const double*lo,const double*hi,double tl,double th,double*eout,double*logout){
 try {
  if(!lo||!hi||!eout||!logout||!std::isfinite(tl)||!std::isfinite(th)||tl<0||th>0.25||tl>th)
   throw std::runtime_error("component input domain");
  I t(tl,th);std::vector<I>z;
  for(int k=0;k<40;k++){
   if(!std::isfinite(lo[k])||!std::isfinite(hi[k])||lo[k]<=0||hi[k]>=1||lo[k]>hi[k])
    throw std::runtime_error("canonical box domain");
   z.push_back(I(lo[k],hi[k]));
  }
  for(int j=0;j<5;j++){
   I total(0);
   for(int i=0;i<4;i++)for(int c=0;c<2;c++){
    I y=t*log_I(z[(i*5+j)*2+c]); total=total+yexp_range(y);
   }
   I e=I(1)+rat(1,4)*total;N quarter(0.25),one(1.0);
   // Global analytic bound 1-2/e > 1/4 for eight x log(x) terms.
   if(mpfr_cmp(e.l.v,quarter.v)<0)e.l=quarter;
   if(mpfr_cmp(e.h.v,one.v)>0)e.h=one;
   if(mpfr_cmp(e.l.v,e.h.v)>0)throw std::runtime_error("entropy range empty");
   result(e,eout+2*j);
  }
  for(int k=0;k<40;k++)result(log_I(k%2 ? z[k] : I(1)-z[k]),logout+2*k);
  return 0;
 }catch(const std::exception&e){ERROR=e.what();return 1;}
}

extern "C" int v17_exp(double lo,double hi,double*out){
 try{if(!out||!std::isfinite(lo)||!std::isfinite(hi)||lo>hi)throw std::runtime_error("exp input");
 result(exp_I(I(lo,hi)),out);return 0;}
 catch(const std::exception&e){ERROR=e.what();return 1;}
}
