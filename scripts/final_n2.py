import numpy as np, csv, time, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
t0=time.time()

Nmax=30_000_000

s=np.ones(Nmax+1,bool); s[:2]=False
for p in range(2,int(Nmax**0.5)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
print("primes<=Nmax: %d (%.1fs)"%(len(primes),time.time()-t0))

def sqrtm1(q):
    z=2
    while pow(z,(q-1)//2,q)!=q-1: z+=1
    return pow(z,(q-1)//4,q)

is_comp=np.zeros(Nmax+1,bool); is_comp[1::2]=True
for q in primes:
    q=int(q)
    if q==2: continue
    if q&3==1:
        r=sqrtm1(q); is_comp[r::q]=True; is_comp[q-r::q]=True
print("root sieve done %.1fs"%(time.time()-t0))

isq=int(Nmax**0.5)+2
is_n1=~is_comp; is_n1[0]=is_n1[1]=False
for n in range(2,isq+1):
    v=n*n+1; pr=True
    for p in primes:
        if p*p>v: break
        if v%p==0: pr=False; break
    is_n1[n]=pr
n1=np.nonzero(is_n1)[0]; Q=len(n1)
print("Q(N)=%d  (%.1fs)"%(Q,time.time()-t0))

# ---------- (1) wings ----------
m6=n1%6
wings=[]
for a in range(6):
    c=int((m6==a).sum())
    lab={0:"right_6N+1",2:"left_6N-1",4:"left_6N-1"}.get(a,"odd_annihilated")
    wings.append((a,c,lab))
right=int((m6==0).sum()); left=int(((m6==2)|(m6==4)).sum())
with open(os.path.join(_DATA,"n2_wings.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["n_mod6","count_n2plus1_prime","wing"])
    for a,c,lab in wings: w.writerow([a,c,lab])

# ---------- (2) per-prime local factors g(q) (character spectrum) ----------
# g(q)=(1-omega/q)/(1-1/q): q=2 ->1 ; q=1mod4 ->(q-2)/(q-1) ; q=3mod4 -> q/(q-1)
with open(os.path.join(_DATA,"n2_local_factors.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["q","q_mod4","omega","g_predicted","g_measured"])
    for q in primes[primes<=2000]:
        q=int(q)
        if q==2: om,g=1,1.0
        elif q&3==1: om,g=2,(q-2)/(q-1)
        else: om,g=0,q/(q-1)
        # measured: fraction of n in [1,Nmax] with q | n^2+1 -> count/Nmax ; g_meas=(1-frac)/(1-1/q)
        if om==0: frac=0.0
        else:
            # count n in [1,Nmax] with n^2+1==0 mod q
            if q==2: frac=(Nmax//2)/Nmax
            else:
                r=sqrtm1(q); c1=len(range(r,Nmax+1,q)); c2=len(range(q-r,Nmax+1,q)); frac=(c1+c2)/Nmax
        gmeas=(1-frac)/(1-1/q)
        w.writerow([q,q%4,om,"%.6f"%g,"%.6f"%gmeas])

# ---------- (3) Bateman-Horn cumulative & running constant ----------
# cumulative Q(N) and integral at sampled N
csum=np.cumsum(is_n1)                          # cumulative count over n
tt=np.arange(0,Nmax+1,dtype=np.float64)
with np.errstate(divide='ignore'):
    integ=np.concatenate([[0.0],np.cumsum(1.0/np.log(tt[2:]*tt[2:]+1.0))])  # int from 2
# align: integ index corresponds to n starting at 2; build sampled rows
Ns=[10**k for k in range(3,8) if 10**k<=Nmax]+[Nmax]
with open(os.path.join(_DATA,"n2_bh_cumulative.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["N","Q_cumulative","BH_integral","C_empirical_Q_over_int"])
    for N in Ns:
        QN=int(csum[N]); 
        I=float(np.sum(1.0/np.log(np.arange(2,N+1,dtype=np.float64)**2+1.0)))
        w.writerow([N,QN,"%.2f"%I,"%.5f"%(QN/I)])

# running product C_x at sampled primes
ppall=primes[primes>2].astype(np.float64)
facs=np.where(ppall%4==1,(ppall-2)/(ppall-1),ppall/(ppall-1))
runC=np.cumprod(facs)
with open(os.path.join(_DATA,"n2_running_C.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["x_prime","C_product_up_to_x"])
    idx=np.unique(np.clip((10**np.linspace(0.7,np.log10(len(ppall)-1),60)).astype(int),0,len(ppall)-1))
    for i in idx: w.writerow([int(ppall[i]),"%.6f"%runC[i]])

# ---------- (4) summary ----------
C_prod=float(runC[-1]); I_full=float(np.sum(1.0/np.log(np.arange(2,Nmax+1,dtype=np.float64)**2+1.0)))
with open(os.path.join(_DATA,"n2_summary.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["parameter","value"])
    for k,v in [("Nmax",Nmax),("n2plus1_max",Nmax*Nmax+1),("n_primes_le_Nmax",len(primes)),
        ("Q_count",Q),("right_wing_6Np1",right),("left_wing_6Nm1",left),
        ("right_over_left","%.5f"%(right/left)),("predicted_right_over_left","0.50000"),
        ("C_product_to_Nmax","%.5f"%C_prod),("C_known_Landau_Shanks","1.37281"),
        ("BH_integral","%.1f"%I_full),("C_empirical_Q_over_int","%.5f"%(Q/I_full)),
        ("BH_ratio_Q_over_Cknown_int","%.5f"%(Q/(1.37281*I_full)))]:
        w.writerow([k,v])
print("\nQ=%d right=%d left=%d r/l=%.4f"%(Q,right,left,right/left))
print("C_product=%.5f  C_empirical=%.5f  BH_ratio=%.5f"%(C_prod,Q/I_full,Q/(1.37281*I_full)))
print("TOTAL %.1fs"%(time.time()-t0))
