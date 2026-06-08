import csv, numpy as np, matplotlib, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,"figure.dpi":150,"savefig.bbox":"tight"})

def rd(n): return list(csv.DictReader(open(os.path.join(_DATA,n))))
summ={r["parameter"]:r["value"] for r in rd("n2_summary.csv")}
Cknown=1.37281

# ============ FIGURE 1 ============
fig,ax=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) wings
wg=rd("n2_wings.csv"); a=[int(r["n_mod6"]) for r in wg]; c=[int(r["count_n2plus1_prime"]) for r in wg]
cols=["#1f77b4" if r["wing"]=="right_6N+1" else ("#2ca02c" if r["wing"]=="left_6N-1" else "0.7") for r in wg]
ax[0].bar(a,c,color=cols,width=0.7)
ax[0].set_xlabel(r"$n \ \mathrm{mod}\ 6$"); ax[0].set_ylabel(r"# of primes $n^2+1$  ($n\leq3\times10^7$)")
ax[0].set_title("(A) wing decomposition: odd $n$ annihilated",fontsize=9)
ax[0].annotate("right\n$6N{+}1$",(0,c[0]),fontsize=7,ha="center",va="bottom",color="#1f77b4")
ax[0].annotate("left $6N{-}1$",(3,max(c)*0.55),fontsize=8,ha="center",color="#2ca02c")
ax[0].text(0.5,0.92,"right:left $=0.5001$ (predicted $1{:}2$)",transform=ax[0].transAxes,fontsize=7.5,ha="center")
# (B) character spectrum g(q)
lf=rd("n2_local_factors.csv")
q=np.array([int(r["q"]) for r in lf]); g=np.array([float(r["g_predicted"]) for r in lf]); m4=np.array([int(r["q_mod4"]) for r in lf])
sel=q<=160
ax[1].axhline(1.0,color="0.5",lw=0.8,ls="--")
s1=sel&(m4==1); s3=sel&(m4==3)
ax[1].plot(q[s3],g[s3],"o",ms=4,color="#d62728",label=r"$q\equiv3\,(4)$: never divides, $\frac{q}{q-1}>1$")
ax[1].plot(q[s1],g[s1],"s",ms=4,color="#1f77b4",label=r"$q\equiv1\,(4)$: two roots, $\frac{q-2}{q-1}<1$")
ax[1].plot([2],[1.0],"^",ms=5,color="0.4",label="$q=2$")
ax[1].set_xlabel(r"prime $q$"); ax[1].set_ylabel(r"local factor $g(q)=\frac{1-\omega(q)/q}{1-1/q}$")
ax[1].set_title(r"(B) the $(-1|q)$ character spectrum",fontsize=9)
ax[1].legend(fontsize=6.5,loc="upper right")
fig.suptitle(r"Primes of the form $n^2+1$ on the $6N$ skeleton",fontsize=10)
fig.savefig(os.path.join(_FIG,"p31_fig1.pdf")); print("fig1 done")

# ============ FIGURE 2 ============
fig2,ax2=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) Bateman-Horn: C_empirical(N) -> C
bh=rd("n2_bh_cumulative.csv")
N=np.array([float(r["N"]) for r in bh]); Cemp=np.array([float(r["C_empirical_Q_over_int"]) for r in bh])
ax2[0].axhline(Cknown,color="#d62728",lw=1.0,ls="--",label=r"Landau--Shanks $C=1.37281$")
ax2[0].semilogx(N,Cemp,"o-",ms=4,color="#1f77b4",label=r"$Q(N)\,/\!\int_2^N\!dt/\log(t^2{+}1)$")
ax2[0].set_xlabel("$N$"); ax2[0].set_ylabel("empirical constant")
ax2[0].set_title("(A) Bateman--Horn density verification",fontsize=9)
ax2[0].legend(fontsize=7.5,loc="lower right"); ax2[0].set_ylim(1.20,1.45)
# (B) running product C_x -> C
rc=rd("n2_running_C.csv")
x=np.array([float(r["x_prime"]) for r in rc]); Cx=np.array([float(r["C_product_up_to_x"]) for r in rc])
ax2[1].axhline(Cknown,color="#d62728",lw=1.0,ls="--",label=r"$C=1.37281$")
ax2[1].semilogx(x,Cx,"-",color="#2ca02c",lw=1.3,label=r"$\prod_{p\leq x} g(p)$")
ax2[1].axhline(float(summ["C_empirical_Q_over_int"]),color="#1f77b4",lw=0.9,ls=":",label=r"measured $C=%.4f$"%float(summ["C_empirical_Q_over_int"]))
ax2[1].set_xlabel("$x$"); ax2[1].set_ylabel(r"$\prod_{p\leq x} g(p)$")
ax2[1].set_title("(B) the character product converges to $C$",fontsize=9)
ax2[1].legend(fontsize=7.5,loc="lower right")
fig2.savefig(os.path.join(_FIG,"p31_fig2.pdf")); print("fig2 done")
