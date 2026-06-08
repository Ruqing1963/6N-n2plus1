# 6N-n2plus1

**Part XXXI — Primes of the Form n²+1 on the 6N Skeleton: the Wing Decomposition and the Character Spectrum of the Landau–Shanks Constant**

Ruqing Chen · GUT Geoservice Inc., Montreal · June 2026

Companion code and data for Part XXXI of *Arithmetic Geodynamics on the 6N Skeleton* — the first
**nonlinear** sequence in the series. **Everything here is a measured sieve result — no fitted
parameters, no fabricated numbers.**

## The question

Primes of the form `n²+1` (Landau's fourth problem). We enumerate all `n ≤ 3×10⁷` with `n²+1` prime
(**1,275,229** of them, `n²+1` as large as `9×10¹⁴`) via a **square-root sieve**: a prime `q` divides
`n²+1` iff `n ≡ ±r (mod q)` with `r² ≡ −1`, which has solutions only for `q = 2` and `q ≡ 1 (mod 4)`.

## The result, in one paragraph

1. **Wing decomposition.** `n` must be even (odd `n` → `n²+1` even → annihilated). Among even `n`:
   `6 | n` ⟹ `n²+1 ≡ 1 (mod 6)` (right wing `6N+1`); `n ≡ 2,4 (mod 6)` ⟹ `n²+1 ≡ 5 (mod 6)`
   (left wing `6N−1`). Measured right:left = **0.50015** — exactly the 1:2 of pure counting, **no
   arithmetic bias** between wings (the three even classes hold 425159 / 425348 / 424722, equal to
   0.07%).

2. **Character spectrum.** The local factor `g(q) = (1 − ω(q)/q)/(1 − 1/q)` splits the odd primes by
   the quadratic character of −1: `q ≡ 3 (mod 4)` **never** divide `n²+1` and *enhance* it by
   `q/(q−1) > 1`; `q ≡ 1 (mod 4)` divide on two residues and *suppress* by `(q−2)/(q−1) < 1`. (First
   constellation in the series whose primes help rather than hinder.)

3. **Landau–Shanks constant & Bateman–Horn.** The product `C = Π_p g(p) = 1.37281`. It controls the
   count: measured `Q(N) / ∫₂ᴺ dt/log(t²+1) = 1.37253`, matching `C` to **0.02%**, and the empirical
   constant converges onto `C` as `N` grows. The nonlinear (degree-two) local–global heuristic is
   verified to a part in 10³ at `n²+1 ~ 10¹⁵`.

**Scope (honest).** None of this is new as theory — ω(q), C, and the asymptotic are due to
Hardy–Littlewood, Bateman–Horn, and Shanks, and C has been computed before. What this adds is the 6N
wing decomposition, the character-spectrum reading, and the high-precision verification via the
square-root sieve. **No claim is made about Landau's problem** (the infinitude of such primes), which
no finite computation can reach.

## Reproducing

```bash
pip install -r requirements.txt
cd code
python3 explore_n2.py     # quick look at Nmax=3e6: wings + Bateman-Horn (console)
python3 final_n2.py       # square-root sieve to Nmax=3e7 (n²+1 up to 9e14)
                          #   -> data/n2_wings.csv, n2_local_factors.csv,
                          #      n2_running_C.csv, n2_bh_cumulative.csv, n2_summary.csv   (~15 s)
python3 makefigs_n2.py    # reads the CSVs -> figures/p31_fig1.pdf, p31_fig2.pdf
```

Paths resolve relative to the script (outputs land in `../data` and `../figures`). The sieve needs
~0.1 GB RAM; single-threaded.

## Files

```
code/    explore_n2.py   final_n2.py   makefigs_n2.py
data/    n2_wings.csv          n_mod6, count, wing
         n2_local_factors.csv  q, q_mod4, omega, g_predicted, g_measured
         n2_running_C.csv      x_prime, C_product_up_to_x
         n2_bh_cumulative.csv  N, Q_cumulative, BH_integral, C_empirical
         n2_summary.csv        parameter, value
figures/ p31_fig1.pdf  p31_fig2.pdf
paper/   paper31.tex   paper31.pdf
```

All data files are plain CSV — openable in any text editor or spreadsheet.

## Citation

See `CITATION.cff`. The paper is archived on Zenodo (DOI in the citation file once minted).

## License

MIT (see `LICENSE`).
