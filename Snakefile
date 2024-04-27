nreps = 10 #number of replicates
ms = [0,0.1,0.2,0.3,0.4,0.5] #[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5] #mutation rates
ns = [4,10] #max sizes
As = [0,1] #0 for size-independent birth, 1 for cooperator-dependent birth
globalls = [0,1] #0 for local dispersal, 1 for global dispersal
b = 0.3

outfile = 'data/sim_n{n}_m{m}_globall{globall}_a{a}_rep{rep}.npy'

rule sim:
  input:
    'sim.py'
  output:
    outfile
  resources:
    runtime=15
  threads: 1
  group: "sim" 
  run:
    import os
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["GOTO_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    from sim import sim
    import numpy as np
    bs = [0] + [b + float(wildcards.a)*i/(int(wildcards.n)-1) for i in range(int(wildcards.n))] #birth rate of cells in group with i cooperators
    results = sim(T=5000, M=50, N0=2000, n=int(wildcards.n), bs=bs, d=0.03, r=2, s=1, m=float(wildcards.m), w=float(wildcards.m), globall=int(wildcards.globall))
    np.save(output[0], results)

rule sims:
  input:
    expand(outfile, n=ns, m=ms, a=As, globall=globalls, rep=range(nreps))
