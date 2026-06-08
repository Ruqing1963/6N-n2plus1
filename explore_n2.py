import numpy as np, time
t0=time.time()
Nmax=3_000_000

# primes up to Nmax (for the root sieve; need q up to ~sqrt(Nmax^2+1)=Nmax)
s=np.ones(Nmax+1,bool); s[:2]=False
for p in range(2,int(Nmax**0.5)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
print("primes<=Nmax: %d  (%.1fs)"%(len(primes),time.time()-t0))

def sqrtm1(q):                       # sqrt(-1) mod q for q prime, q%4==1
    z=2
    while pow(z,(q-1)//2,q)!=q-1: z+=1
    return pow(z,(q-1)//4,q)

# root sieve: is_comp[n]=True if n^2+1 has a prime factor q<=Nmax (so composite, unless n^2+1==q)
is_comp=np.zeros(Nmax+1,bool)
is_comp[1::2]=True                   # odd n -> n^2+1 even
cnt1=0
for q in primes:
    q=int(q)
    if q==2: continue
    if q&3==1:
        r=sqrtm1(q); cnt1+=1
        is_comp[r::q]=True
        is_comp[q-r::q]=True
print("primes=1mod4 used: %d  sieve done %.1fs"%(cnt1,time.time()-t0))

# n^2+1 prime indicator; override small n where n^2+1<=Nmax could equal a sieving prime
isq=int(Nmax**0.5)+2
is_n1=~is_comp
is_n1[0]=is_n1[1]=False
for n in range(2,isq+1):
    v=n*n+1; pr=True
    for p in primes:
        if p*p>v: break
        if v%p==0: pr=False; break
    is_n1[n]=pr
n1=np.nonzero(is_n1)[0]               # all n in [2,Nmax] with n^2+1 prime
Q=len(n1)
print("\nQ(N)=#{n<=%d : n^2+1 prime} = %d"%(Nmax,Q))

# wings by n mod 6
import collections
m6=n1%6
for a in range(6):
    c=int((m6==a).sum())
    wing={0:"right 6N+1",2:"left 6N-1",4:"left 6N-1"}.get(a,"(odd->should be 0)")
    print("  n=%d mod6: %7d   %s"%(a,c,wing))
right=int((m6==0).sum()); left=int(((m6==2)|(m6==4)).sum())
print("  RIGHT(6|n)=%d  LEFT(n=2,4)=%d  right/left=%.4f  (pure-count predicts 0.5)"%(right,left,right/left))

# Bateman-Horn constant: C = prod (1-w(p)/p)/(1-1/p)
# p=2 factor =1 ; p=1mod4: (p-2)/(p-1) ; p=3mod4: p/(p-1)
pp=primes[primes>2].astype(np.float64)
fac=np.where(pp%4==1,(pp-2)/(pp-1),pp/(pp-1))
C_run=np.prod(fac)
print("\nLandau-Shanks C (product to %d) = %.5f   (known ~1.37281)"%(Nmax,C_run))

# Bateman-Horn integral  Q ~ C * int_2^N dt/log(t^2+1)
from math import log
tt=np.arange(2,Nmax+1,dtype=np.float64)
BHint=np.sum(1.0/np.log(tt*tt+1.0))
print("int_2^N dt/log(t^2+1) = %.1f"%BHint)
print("C_known * int = %.1f     measured Q = %d     ratio Q/(C*int)=%.4f"%(1.37281*BHint,Q,Q/(1.37281*BHint)))
print("=> empirical C = Q/int = %.5f"%(Q/BHint))
print("\nTOTAL %.1fs"%(time.time()-t0))
