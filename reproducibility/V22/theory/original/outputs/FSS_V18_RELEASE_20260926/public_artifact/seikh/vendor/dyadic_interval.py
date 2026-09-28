"""Outward intervals using ONLY Python integer arithmetic.

Endpoints are integers / 2**PREC. No MPFR, float transcendental or mpmath.
log: binary argument reduction + atanh Taylor series, explicit positive tail.
exp: scaling to [0,1/8] + positive Taylor series, explicit tail + squaring;
     negative arguments use reciprocal. Roots use exact integer inequalities.
See ARITHMETIC_PROOF.md. Intended domain is bounded, finite, positive when
logarithms are used. It fails closed outside supported domains.
"""
from __future__ import annotations
from fractions import Fraction
from functools import lru_cache
from math import factorial, isqrt

PREC = 224
SCALE = 1 << PREC

def ceildiv(a:int,b:int)->int:
    if b<=0: raise ValueError('positive denominator required')
    return -((-a)//b)

class IV:
    __slots__=('lo','hi')
    def __init__(self,lo:int,hi:int):
        if not isinstance(lo,int) or not isinstance(hi,int) or lo>hi:
            raise ValueError('invalid dyadic interval')
        self.lo,self.hi=lo,hi
    @classmethod
    def frac(cls,a,b=None):
        x=Fraction(a); y=x if b is None else Fraction(b)
        if x>y: raise ValueError('reversed exact input')
        return cls((x.numerator*SCALE)//x.denominator,
                   ceildiv(y.numerator*SCALE,y.denominator))
    @classmethod
    def integer(cls,n:int): return cls(n*SCALE,n*SCALE)
    def __add__(self,b):
        b=coerce(b);return IV(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return IV(-self.hi,-self.lo)
    def __sub__(self,b): return self+-coerce(b)
    def __rsub__(self,b):return coerce(b)+-self
    def __mul__(self,b):
        b=coerce(b);p=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return IV(min(p)//SCALE,ceildiv(max(p),SCALE))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=coerce(b)
        if b.lo<=0: raise ValueError('division denominator contains nonpositive value')
        low=min((a*SCALE)//c for a in (self.lo,self.hi) for c in (b.lo,b.hi))
        high=max(ceildiv(a*SCALE,c) for a in (self.lo,self.hi) for c in (b.lo,b.hi))
        return IV(low,high)
    def __rtruediv__(self,b):return coerce(b)/self
    def div_int(self,n:int):
        if n<=0:raise ValueError('nonpositive divisor')
        return IV(self.lo//n,ceildiv(self.hi,n))
    def intersect(self,lo:int,hi:int):
        return IV(max(self.lo,lo*SCALE),min(self.hi,hi*SCALE))
    def abs(self):
        if self.lo>=0:return self
        if self.hi<=0:return -self
        return IV(0,max(-self.lo,self.hi))
    def pow_int(self,n:int):
        if n<0 or self.lo<0:raise ValueError('nonnegative base and integer power required')
        if n==0:return IV.integer(1)
        den=SCALE**(n-1)
        return IV(self.lo**n//den,ceildiv(self.hi**n,den))
    def root(self,n:int):
        if self.lo<0 or n<1:raise ValueError('nonnegative root input')
        a=iroot(self.lo*SCALE**(n-1),n);b=iroot(self.hi*SCALE**(n-1),n)
        if b**n<self.hi*SCALE**(n-1):b+=1
        return IV(a,b)
    def log(self):
        if self.lo<=0:raise ValueError('log needs strictly positive lower endpoint')
        return IV(log_point(self.lo).lo,log_point(self.hi).hi)
    def exp(self):return IV(exp_point(self.lo).lo,exp_point(self.hi).hi)
    def pow_positive(self,p):
        p=coerce(p)
        if self.lo<=0:raise ValueError('positive power base required')
        return (p*self.log()).exp()
    def pow_unit(self,p):
        """Monotone power enclosure for base in [0,1], exponent >0."""
        p=coerce(p)
        if self.lo<0 or self.hi>SCALE or p.lo<=0:raise ValueError('unit power domain')
        l=0 if self.lo==0 else (IV(self.lo,self.lo).log()*IV(p.hi,p.hi)).exp().lo
        h=0 if self.hi==0 else (IV(self.hi,self.hi).log()*IV(p.lo,p.lo)).exp().hi
        return IV(max(0,l),min(SCALE,h))
    def endpoints(self):return [str(Fraction(self.lo,SCALE)),str(Fraction(self.hi,SCALE))]
    def __repr__(self):return f'IV({float(Fraction(self.lo,SCALE))}, {float(Fraction(self.hi,SCALE))})'

def coerce(b):
    if isinstance(b,IV):return b
    if isinstance(b,int):return IV.integer(b)
    if isinstance(b,Fraction):return IV.frac(b)
    raise TypeError('only intervals, integers and exact Fractions are accepted')

def iroot(a:int,n:int)->int:
    if a<0 or n<1:raise ValueError('root domain')
    if n==1 or a<=1:return a
    if n==2:return isqrt(a)
    x=1<<ceildiv(a.bit_length(),n)
    while True:
        y=((n-1)*x+a//(x**(n-1)))//n
        if y>=x:break
        x=y
    while (x+1)**n<=a:x+=1
    while x**n>a:x-=1
    return x

LOG_N=PREC//3+12
EXP_N=PREC//3+12
# tail log(1+z)/(1-z), 0<=z<=1/3, after N terms:
# <= (9/4) * 3^(-2N-1)/(2N+1).
LOG_TAIL=IV.frac(Fraction(9,4*(2*LOG_N+1)*3**(2*LOG_N+1)))
# exp y, 0<=y<=1/8; remainder after degree N:
# <= (8/7)*(1/8)^(N+1)/(N+1)!.
EXP_TAIL=IV.frac(Fraction(8,7*8**(EXP_N+1)*factorial(EXP_N+1)))

def atanh_log(m:IV)->IV:
    if m.lo<SCALE or m.hi>2*SCALE:raise ValueError('log reduction outside [1,2]')
    z=(m-1)/(m+1); z2=z*z; term=z;out=IV.integer(0)
    for k in range(LOG_N):
        out=out+term.div_int(2*k+1)
        term=term*z2
    out=2*out
    return IV(out.lo,out.hi+LOG_TAIL.hi)

LN2=atanh_log(IV.integer(2))
@lru_cache(maxsize=30000)
def log_point(a:int)->IV:
    if a<=0:raise ValueError('log domain')
    if a==SCALE:return IV.integer(0)
    k=a.bit_length()-1-PREC
    # Exact a/SCALE = m*2^k. Enclose m on the fixed lattice.
    if k>=0:m=IV(a//(1<<k),ceildiv(a,1<<k))
    else:m=IV(a<<(-k),a<<(-k))
    return atanh_log(m)+k*LN2

@lru_cache(maxsize=30000)
def exp_point(a:int)->IV:
    if a==0:return IV.integer(1)
    if abs(a)>1024*SCALE:raise ValueError('exp argument exceeds declared bounded implementation domain')
    if a<0:return 1/exp_point(-a)
    k=max(0,a.bit_length()+3-PREC)
    y=IV(a//(1<<k),ceildiv(a,1<<k))
    if y.hi>SCALE//8:raise ValueError('exp reduction failed')
    term=IV.integer(1);out=term
    for n in range(1,EXP_N+1):
        term=(term*y).div_int(n);out=out+term
    out=IV(out.lo,out.hi+EXP_TAIL.hi)
    for _ in range(k):out=out.pow_int(2)
    return out
