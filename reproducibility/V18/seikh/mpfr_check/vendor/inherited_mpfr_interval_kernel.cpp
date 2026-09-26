// Directed-rounding interval kernel for the archived Qin (2019) specification.
// Uses the stable MPFR 4.x C ABI. No Python numerical package is required.
// This is an independent bounded-transform implementation; source formulas and
// high-precision replay are cross-checked separately. See MATHEMATICAL_NOTE.md.
#include <gmp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <string>
#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
extern "C" {
// Public MPFR ABI declarations (a subset of mpfr.h, to support runtime-only installs).
typedef long mpfr_prec_t; typedef long mpfr_exp_t;
typedef struct {mpfr_prec_t _mpfr_prec; int _mpfr_sign; mpfr_exp_t _mpfr_exp; mp_limb_t *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1]; typedef __mpfr_struct *mpfr_ptr; typedef const __mpfr_struct *mpfr_srcptr;
typedef enum {MPFR_RNDN=0, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD, MPFR_RNDA, MPFR_RNDF} mpfr_rnd_t;
void mpfr_init2(mpfr_ptr,mpfr_prec_t); void mpfr_clear(mpfr_ptr);
int mpfr_set(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_set_d(mpfr_ptr,double,mpfr_rnd_t);
int mpfr_set_ui(mpfr_ptr,unsigned long,mpfr_rnd_t); int mpfr_set_si(mpfr_ptr,long,mpfr_rnd_t);
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_pow(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_pow_ui(mpfr_ptr,mpfr_srcptr,unsigned long,mpfr_rnd_t);
int mpfr_rootn_ui(mpfr_ptr,mpfr_srcptr,unsigned long,mpfr_rnd_t); int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_neg(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t); int mpfr_cmp(mpfr_srcptr,mpfr_srcptr); int mpfr_cmp_ui(mpfr_srcptr,unsigned long);
int mpfr_number_p(mpfr_srcptr); double mpfr_get_d(mpfr_srcptr,mpfr_rnd_t); const char *mpfr_get_version(void);
}
#endif
static thread_local unsigned BITS=160;
struct N {mpfr_t v; N(){mpfr_init2(v,BITS);mpfr_set_ui(v,0,MPFR_RNDN);} explicit N(double d):N(){mpfr_set_d(v,d,MPFR_RNDN);} N(const N&x):N(){mpfr_set(v,x.v,MPFR_RNDN);} N&operator=(const N&x){mpfr_set(v,x.v,MPFR_RNDN);return *this;} ~N(){mpfr_clear(v);} };
struct I {N l,h; I(){} explicit I(double d):l(d),h(d){} I(double a,double b):l(a),h(b){if(a>b)throw std::runtime_error("reversed input");} I(const N&a,const N&b):l(a),h(b){} };
static I rat(unsigned long a,unsigned long b){I x((double)a),y((double)b); I z;mpfr_div(z.l.v,x.l.v,y.h.v,MPFR_RNDD);mpfr_div(z.h.v,x.h.v,y.l.v,MPFR_RNDU);return z;}
static I operator+(const I&a,const I&b){I z;mpfr_add(z.l.v,a.l.v,b.l.v,MPFR_RNDD);mpfr_add(z.h.v,a.h.v,b.h.v,MPFR_RNDU);return z;}
static I operator-(const I&a,const I&b){I z;mpfr_sub(z.l.v,a.l.v,b.h.v,MPFR_RNDD);mpfr_sub(z.h.v,a.h.v,b.l.v,MPFR_RNDU);return z;}
static I operator*(const I&a,const I&b){I z;N t;bool first=true;for(auto u:{a.l.v,a.h.v})for(auto v:{b.l.v,b.h.v}){mpfr_mul(t.v,u,v,MPFR_RNDD);if(first||mpfr_cmp(t.v,z.l.v)<0)z.l=t;mpfr_mul(t.v,u,v,MPFR_RNDU);if(first||mpfr_cmp(t.v,z.h.v)>0)z.h=t;first=false;}return z;}
static I operator/(const I&a,const I&b){if(mpfr_cmp_ui(b.l.v,0)<=0)throw std::runtime_error("nonpositive denominator");I inv;N one(1.0);mpfr_div(inv.l.v,one.v,b.h.v,MPFR_RNDD);mpfr_div(inv.h.v,one.v,b.l.v,MPFR_RNDU);return a*inv;}
static I unit(I a){if(mpfr_cmp_ui(a.l.v,0)<0)mpfr_set_ui(a.l.v,0,MPFR_RNDN);if(mpfr_cmp_ui(a.h.v,1)>0)mpfr_set_ui(a.h.v,1,MPFR_RNDN);if(mpfr_cmp(a.l.v,a.h.v)>0)throw std::runtime_error("invalid unit range");return a;}
static I pown(const I&a,unsigned n){if(mpfr_cmp_ui(a.l.v,0)<0)throw std::runtime_error("negative base");I z;mpfr_pow_ui(z.l.v,a.l.v,n,MPFR_RNDD);mpfr_pow_ui(z.h.v,a.h.v,n,MPFR_RNDU);return z;}
static I root(const I&a,unsigned n){if(!n||mpfr_cmp_ui(a.l.v,0)<0)throw std::runtime_error("root domain");I z;mpfr_rootn_ui(z.l.v,a.l.v,n,MPFR_RNDD);mpfr_rootn_ui(z.h.v,a.h.v,n,MPFR_RNDU);return z;}
static I pow01(I a,const I&c){a=unit(a);if(mpfr_cmp_ui(c.l.v,0)<=0)throw std::runtime_error("exponent domain");I z;mpfr_pow(z.l.v,a.l.v,c.h.v,MPFR_RNDD);mpfr_pow(z.h.v,a.h.v,c.l.v,MPFR_RNDU);return unit(z);}
static I absv(const I&a){N zero(0.0);I z;if(mpfr_cmp_ui(a.l.v,0)>=0)return a;if(mpfr_cmp_ui(a.h.v,0)<=0){mpfr_neg(z.l.v,a.h.v,MPFR_RNDD);mpfr_neg(z.h.v,a.l.v,MPFR_RNDU);return z;}mpfr_neg(z.h.v,a.l.v,MPFR_RNDU);if(mpfr_cmp(a.h.v,z.h.v)>0)z.h=a.h;return z;}
// Evaluate monotone rational functions at endpoint points. Dependency is avoided,
// while every primitive and endpoint conversion remains outward rounded.
static I rational_map(const I&x,int code){auto f=[&](const N&t){I z(t,t);if(code==0)return (I(1)-z)/(I(1)+I(8)*z);if(code==1)return (I(1)-z)/(I(1)+I(2)*z);if(code==2)return z/(I(3)-I(2)*z);return I(3)*z/(I(1)+I(2)*z);};bool dec=code<2;I lo=f(dec?x.h:x.l),hi=f(dec?x.l:x.h);return unit(I(lo.l,hi.h));}
static I K(const I&x){return rational_map(x,0);} static I bg(const I&u){return rational_map(u,1);} static I bf(const I&v){return rational_map(v,2);} static I invbf(const I&b){return rational_map(b,3);}
using Pair=std::array<I,2>;
static Pair B(const Pair&z){return {bg(z[0]),bf(z[1])};} static Pair can(const Pair&z){return {bg(z[0]),invbf(z[1])};}
static I distance(const Pair&a,const Pair&b,unsigned q){I d0=absv(root(a[0],q)-root(b[0],q)),d1=absv(root(a[1],q)-root(b[1],q));return unit(root(rat(1,2)*(pown(d0,3)+pown(d1,3)),3));}
static std::vector<I> weights(const std::vector<Pair>&v,const std::vector<I>&base,unsigned q){std::vector<I> z;I den(0);for(size_t i=0;i<v.size();i++){I t(1);for(size_t j=0;j<v.size();j++)if(i!=j)t=t+(I(1)-distance(v[i],v[j],q));z.push_back(base[i]*t);den=den+z.back();}for(auto&x:z)x=x/den;return z;}
static Pair expert_average(const std::vector<Pair>&v,const std::vector<I>&w){Pair acc={I(1),I(1)};for(size_t i=0;i<v.size();i++){Pair b=B(v[i]);for(int c=0;c<2;c++)acc[c]=acc[c]*pow01(b[c],w[i]);}return can(acc);}
static Pair partition(const std::vector<Pair>&v,const std::vector<I>&w,const std::vector<unsigned>&d){unsigned m=v.size(),fact=1,sum=0;for(unsigned i=1;i<=m;i++)fact*=i;for(auto x:d)sum+=x;std::vector<Pair> z;for(unsigned i=0;i<m;i++){Pair b=B(v[i]);z.push_back({K(pow01(b[0],I(m)*w[i])),K(pow01(b[1],I(m)*w[i]))});}std::vector<unsigned> p;for(unsigned i=0;i<m;i++)p.push_back(i);Pair acc={I(1),I(1)};do{Pair prod={I(1),I(1)};for(unsigned k=0;k<m;k++)for(int c=0;c<2;c++)prod[c]=prod[c]*pown(z[p[k]][c],d[k]);for(int c=0;c<2;c++)acc[c]=acc[c]*K(prod[c]);}while(std::next_permutation(p.begin(),p.end()));for(int c=0;c<2;c++)acc[c]=K(pow01(K(pow01(acc[c],rat(1,fact))),rat(1,sum)));return can(acc);}
static Pair aggregate_criteria(const std::vector<Pair>&v,const std::vector<I>&w){Pair p=B(partition({v[0],v[2],v[4]},{w[0],w[2],w[4]},{1,2,3}));Pair t=B(partition({v[1],v[3]},{w[1],w[3]},{1,2}));return can({root(p[0]*t[0],2),root(p[1]*t[1],2)});}
// Exact printed decimals: every numerator below is divided by ten.
static const int DATA[200]={7,2,8,2,4,5,7,1,9,2,8,6,7,6,4,5,5,3,7,2,6,5,5,4,5,6,8,5,8,3,7,2,6,5,6,5,6,2,6,5,6,4,7,5,5,6,7,4,7,3,7,1,7,3,3,5,6,1,8,2,8,3,7,5,2,7,9,2,6,2,9,2,8,3,5,3,8,4,9,2,6,3,7,5,4,6,8,3,7,6,7,3,8,4,3,4,7,2,6,3,8,2,6,2,3,7,8,1,7,1,8,6,7,4,4,8,7,4,6,4,7,2,8,3,2,6,9,2,8,2,4,7,9,2,2,5,5,5,9,6,6,4,7,3,3,6,7,3,8,5,9,3,6,2,3,6,7,3,7,2,6,6,7,7,2,6,6,2,6,2,7,1,9,4,4,7,8,3,9,4,7,3,5,4,4,5,9,4,6,2,8,4,8,5,2,5,7,4,8,3};
static std::vector<Pair> state(const double*lo,const double*hi){std::vector<Pair> z(100);for(int i=0;i<200;i++){int a=DATA[i];z[i/2][i%2]=(lo&&hi&&!std::isnan(lo[i]))?I(lo[i],hi[i]):rat(a*a*a,1000);}return z;}
static I score(const std::vector<Pair>&z,int a,unsigned qe,unsigned qc){std::vector<I> ew={rat(30,100),rat(22,100),rat(28,100),rat(20,100)},cw={rat(20,100),rat(20,100),rat(15,100),rat(25,100),rat(20,100)};std::vector<Pair> coll;for(int c=0;c<5;c++){std::vector<Pair> vals;for(int e=0;e<4;e++)vals.push_back(z[e*25+a*5+c]);auto w=weights(vals,ew,qe);coll.push_back(expert_average(vals,w));}auto w=weights(coll,cw,qc);Pair f=aggregate_criteria(coll,w);return f[0]-f[1];}
static void result(const I&x,double*out){if(!mpfr_number_p(x.l.v)||!mpfr_number_p(x.h.v))throw std::runtime_error("nonfinite output");out[0]=mpfr_get_d(x.l.v,MPFR_RNDD);out[1]=mpfr_get_d(x.h.v,MPFR_RNDU);}
static thread_local std::string ERROR;
extern "C" {
const char* kernel_version(){return mpfr_get_version();} const char* kernel_error(){return ERROR.c_str();}
void kernel_precision(unsigned bits){if(bits>=64&&bits<=4096)BITS=bits;}
int kernel_scores(const double*lo,const double*hi,unsigned qe,unsigned qc,double*out){try{if(qe<1||qc<1)throw std::runtime_error("integer q required");auto z=state(lo,hi);for(int a=0;a<5;a++)result(score(z,a,qe,qc),out+2*a);return 0;}catch(const std::exception&e){ERROR=e.what();return 1;}}
int kernel_margin(const double*lo,const double*hi,int j,double*out){try{if(j<1||j>4)throw std::runtime_error("challenger index");auto z=state(lo,hi);result(score(z,0,3,3)-score(z,j,3,3),out);return 0;}catch(const std::exception&e){ERROR=e.what();return 1;}}
int kernel_contract(int n,const int*idx,double*lo,double*hi,double R){try{I radius(R),rr=pown(radius,2);std::vector<I> centre;for(int i=0;i<n;i++){int a=DATA[idx[i]];centre.push_back(rat(a*a*a,1000));}for(int iter=0;iter<2;iter++){std::vector<N> sq;N tot(0.0);for(int i=0;i<n;i++){N dev(0.0);N l(lo[i]),h(hi[i]);if(mpfr_cmp(l.v,centre[i].h.v)>0)mpfr_sub(dev.v,l.v,centre[i].h.v,MPFR_RNDD);else if(mpfr_cmp(h.v,centre[i].l.v)<0)mpfr_sub(dev.v,centre[i].l.v,h.v,MPFR_RNDD);N s;mpfr_mul(s.v,dev.v,dev.v,MPFR_RNDD);sq.push_back(s);mpfr_add(tot.v,tot.v,s.v,MPFR_RNDD);}if(mpfr_cmp(tot.v,rr.h.v)>0)return 1;for(int i=0;i<n;i++){N other,rem2,rem,l,h;mpfr_sub(other.v,tot.v,sq[i].v,MPFR_RNDD);mpfr_sub(rem2.v,rr.h.v,other.v,MPFR_RNDU);if(mpfr_cmp_ui(rem2.v,0)<0)throw std::runtime_error("contract exclusion mismatch");mpfr_sqrt(rem.v,rem2.v,MPFR_RNDU);mpfr_sub(l.v,centre[i].l.v,rem.v,MPFR_RNDD);mpfr_add(h.v,centre[i].h.v,rem.v,MPFR_RNDU);lo[i]=std::max({lo[i],mpfr_get_d(l.v,MPFR_RNDD),0.0});hi[i]=std::min({hi[i],mpfr_get_d(h.v,MPFR_RNDU),1.0});if(lo[i]>hi[i])return 1;}}return 0;}catch(const std::exception&e){ERROR=e.what();return -1;}}
// Elementary test interface: +, -, *, /, integer power, integer root, [0,1]^positive interval, abs, K.
int kernel_op(int op,double al,double ah,double bl,double bh,double*out){try{I a(al,ah),b(bl,bh),z;switch(op){case 0:z=a+b;break;case 1:z=a-b;break;case 2:z=a*b;break;case 3:z=a/b;break;case 4:z=pown(a,(unsigned)bl);break;case 5:z=root(a,(unsigned)bl);break;case 6:z=pow01(a,b);break;case 7:z=absv(a);break;case 8:z=K(a);break;default:throw std::runtime_error("unknown op");}result(z,out);return 0;}catch(const std::exception&e){ERROR=e.what();return 1;}}
}
