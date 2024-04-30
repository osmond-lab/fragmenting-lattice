# ---------------------------------------------------------------
# Figure 4 (and S1): abundance of modes vs switching rate
# ---------------------------------------------------------------

Ts = [5000] #number of iterations, >0 int
Ms = [100] #length of lattice, >0 int
ns = [4,10] #max group sizes, >0 int
bs = [0.3] #birth rate of single cell, >0 float
ass = [0,1] #strength of cooperator dependent birth, >=0 float
ds = [0.03] #group death rate, >0 float
rs = [2] #replicative advantage of cheater cells, >1 float
ss = [1] #strength of size-dependent competition, [0,1] float
ms = [0,0.1,0.2,0.3,0.4,0.5] #switching rates, >=0 float
localss= [False,True] #local dispersal, boolean
nreps = 10 #number of replicates, >0 int

outfile = 'data/fig4_{T}T_{M}M_{n}n_{b}b_{a}a_{d}d_{r}r_{s}s_{m}m_{m}w_{local}local_{rep}rep.npy'

rule fig_fours:
  input:
    expand(outfile, T=Ts, M=Ms, n=ns, b=bs, a=ass, d=ds, r=rs, s=ss, m=ms, local=localss, rep=range(nreps))

# snakemake fig_fours --profile slurm --group-components fig_fours=80 --jobs 6 -n
# they take about 15m

rule fig_four:
  input:
    'sim.py'
  output:
    outfile
  resources:
    runtime=30
  threads: 1
  group: "fig_fours" 
  run:
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    from sim import sim, non_trivial_partitions
    import numpy as np
    modes = non_trivial_partitions(int(wildcards.n)) #all modes with adult size less than or equal to n
    brs = [0] + [float(wildcards.b) + float(wildcards.a)*i/(int(wildcards.n)-1) for i in range(int(wildcards.n))] #birth rate of cells in group with i cooperators (zero indexed so start with 0)
    results = sim(modes=modes, T=int(wildcards.T), M=int(wildcards.M), N0=int(wildcards.M)**2, bs=brs, d=float(wildcards.d), r=float(wildcards.r), s=float(wildcards.s), m=float(wildcards.m), w=float(wildcards.m), local=wildcards.local)
    np.save(output[0], results)

