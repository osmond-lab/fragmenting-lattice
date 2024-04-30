# ---------------------------------------------------------------
# Figure 3: dominant modes without switching 
# ---------------------------------------------------------------

Ts = [5000] #number of iterations, >0 int
Ms = [100] #length of lattice, >0 int
ns = [4] #max group sizes, >0 int
b1s = [0.3] #birth rate of single cell, >0 float
b2s = [1,2,3,4,5] #birth rate of pair of cells relative to single cell (b2/b1), >0 float
b3s = b2s #birth rate of triplet of cells relative to single cell (b3/b1), >0 float
ds = [0.03] #group death rate, >0 float
rs = [2] #replicative advantage of cheater cells, >1 float
ss = [0,1] #strength of size-dependent competition, [0,1] float
ms = [0] #switching rates, >=0 float
localss = [False,True] #local dispersal, boolean
reps = range(10) #number of replicates, >0 int

outfile = 'data/fig3_{T}T_{M}M_{n}n_{b1}b1_{b2}b2_{b3}b3_{d}d_{r}r_{s}s_{m}m_{m}w_{local}local_{rep}rep.npy'

rule fig_threes:
  input:
    expand(outfile, T=Ts, M=Ms, n=ns, b1=b1s, b2=b2s, b3=b3s, d=ds, r=rs, s=ss, m=ms, local=localss, rep=reps)

# snakemake fig_threes --profile slurm --group-components fig_threes=80 --jobs 6 -n
# take less than 5m i think

rule fig_three:
  input:
    'sim.py'
  output:
    outfile
  resources:
    runtime=30
  threads: 1
  group: "fig_threes" 
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
    brs = [0] + [float(wildcards.b1)*i for i in [1,float(wildcards.b2),float(wildcards.b3)]] #birth rate of cells in group with i cooperators (zero indexed so start with 0)
    results = sim(modes=modes, T=int(wildcards.T), M=int(wildcards.M), N0=int(wildcards.M)**2, bs=brs, d=float(wildcards.d), r=float(wildcards.r), s=float(wildcards.s), m=float(wildcards.m), w=float(wildcards.m), local=wildcards.local)
    np.save(output[0], results)

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
localss = [False,True] #local dispersal, boolean
reps = range(10) #number of replicates, >0 int

outfile = 'data/fig4_{T}T_{M}M_{n}n_{b}b_{a}a_{d}d_{r}r_{s}s_{m}m_{m}w_{local}local_{rep}rep.npy'

rule fig_fours:
  input:
    expand(outfile, T=Ts, M=Ms, n=ns, b=bs, a=ass, d=ds, r=rs, s=ss, m=ms, local=localss, rep=reps)

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
    brs = [0] + [float(wildcards.b) + float(wildcards.a)*i/(int(wildcards.n)-1) for i in range(int(wildcards.n)-1)] #birth rate of cells in group with i cooperators (zero indexed so start with 0)
    results = sim(modes=modes, T=int(wildcards.T), M=int(wildcards.M), N0=int(wildcards.M)**2, bs=brs, d=float(wildcards.d), r=float(wildcards.r), s=float(wildcards.s), m=float(wildcards.m), w=float(wildcards.m), local=wildcards.local)
    np.save(output[0], results)

# ---------------------------------------------------------------
# Figure 5 (and S2): fraction of groups with cheathers
# ---------------------------------------------------------------

Ts = [1000] #number of iterations, >0 int
Ms = [50] #length of lattice, >0 int
bs = [0.3] #birth rate of single cell, >0 float
ass = [0] #strength of cooperator dependent birth, >=0 float
ds = [0.03] #group death rate, >0 float
rs = [2] #replicative advantage of cheater cells, >1 float
ss = [1] #strength of size-dependent competition, [0,1] float
ms = [0.2] #switching rates, >=0 float
localss = [False] #local dispersal, boolean
reps = range(1) #number of replicates, >0 int

outfile = 'data/fig5_{mode}mode_{T}T_{M}M_{n}n_{b}b_{a}a_{d}d_{r}r_{s}s_{m}m_{m}w_{local}local_{rep}rep.npy'

from sim import non_trivial_partitions
rule fig_fives:
  input:
    expand(outfile, mode=range(len(non_trivial_partitions(4))), T=Ts, M=Ms, n=[4], b=bs, a=ass, d=ds, r=rs, s=ss, m=ms, local=localss, rep=reps),
    expand(outfile, mode=range(len(non_trivial_partitions(10))), T=Ts, M=Ms, n=[10], b=bs, a=ass, d=ds, r=rs, s=ss, m=ms, local=localss, rep=reps),

# snakemake fig_fives --profile slurm --group-components fig_fives=80 --jobs 6 -n

rule fig_five:
  input:
    'sim.py'
  output:
    outfile
  resources:
    runtime=15
  threads: 1
  group: "fig_fives" 
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
    brs = [0] + [float(wildcards.b) + float(wildcards.a)*i/(int(wildcards.n)-1) for i in range(int(wildcards.n)-1)] #birth rate of cells in group with i cooperators (zero indexed so start with 0)
    results = sim(modes=[modes[int(wildcards.mode)]], T=int(wildcards.T), M=int(wildcards.M), N0=int(wildcards.M)**2, bs=brs, d=float(wildcards.d), r=float(wildcards.r), s=float(wildcards.s), m=float(wildcards.m), w=float(wildcards.m), local=wildcards.local)
    np.save(output[0], results)

