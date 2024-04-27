import numpy as np
from itertools import product

# determine modes

def partition(number):
    answer = set()
    answer.add((number, ))
    for x in range(1, number):
        for y in partition(number - x):
            answer.add(tuple(sorted((x, ) + y)))
    return answer

def all_partitions_below(number):
    answer = set()
    answer.add((number, ))
    for x in range(1, number):
        partition_lower = all_partitions_below(number - x)
        answer = answer.union(partition_lower)
        for y in partition_lower:
            answer.add(tuple(sorted((x, ) + y)))
    return answer

def non_trivial_partitions(n):
    all_partitions = all_partitions_below(n)
    for i in range(1, n + 1):
        all_partitions.remove((i,))
        pn = sorted(list(all_partitions), key=lambda part: np.sum(part))
    return pn

# simulation

def sim(T=5000, M=50, N0=2000, n=4, bs=[0,0.3,0.3,0.3,0.3], d=0.03, r=2, s=0, m=0, w=0, globall=1):

    # modes   
    modes = non_trivial_partitions(n) #modes
    nmodes = len(modes) #number of modes
    max_sizes = [sum(mode) for mode in modes] #max sizes
 
    # birth rates
    pbs = [1 - np.exp(-b) for b in bs] #probability a cooperator cell divides
    pbs_cancer = [1 - np.exp(-r * b) for b in bs] #probability a cancer cell divides

    # death rate
    pd = 1 - np.exp(-d) #probability a group dies

    # lattice
    group_size = np.zeros((M, M), int)  #number of cells at vertices
    mode = np.zeros((M, M), int) # fragmentation genotype for points on the lattice that are inhabited by cells
    num_cancer = np.zeros((M, M), int) #number of cancer cells
    locations_of_parents = np.zeros((M,M), int) #temporary record of which vertices fragmenting

    # -----------------------------------------------------------------------------|
    # create N0 randomly scattered individuals on an M by M grid
   
    ixs = np.random.choice(M**2, size=N0, replace=False) #N0 unique random numbers
    locations = np.array(list(product(range(M),range(M))))[ixs] #convert to unique vertices

    for i,j in locations:
        group_size[i,j] = 1 #start all individuals as single cells
        mode[i,j] = np.random.randint(1,nmodes+1) #choose mode at random (note we are saving 0 for empty vertex)
        num_cancer[i,j] = 0
    
    # -----------------------------------------------------------------------------|
    # ITERATE
    # -----------------------------------------------------------------------------|
    for t in range(T):

        if t%100==0:
            print(t, np.sum(group_size), np.unique(mode))
   
        # cell division and group death 
        for i in range(M):
            for j in range(M):
                
                gsize = group_size[i,j] #number of cells at this vertex
    
                if gsize > 0:  # if there is an individual (/group) at that spot
    
                    # info
                    ncanc = num_cancer[i,j] #number of cancer cells
                    ncoop = gsize - ncanc #number of cooperator cells
                    msize = max_sizes[mode[i,j]-1]
 
                    # births
                    coop_births = np.random.binomial(ncoop, pbs[ncoop])  #number of cooperator cells that want to divide
                    canc_births = np.random.binomial(ncanc, pbs_cancer[ncoop]) #number of cancer cells that want to divide

                    #mutations
                    coop_muts = np.random.binomial(coop_births, m)
                    canc_muts = np.random.binomial(canc_births, w)
                    new_coop = coop_births - coop_muts + canc_muts
                    new_canc = canc_births - canc_muts + coop_muts     

                    # make sure we dont overshoot max size
                    new_gsize = gsize + new_coop + new_canc
                    if new_gsize > msize:  # if the max size is exceeded
                        new_gsize = msize
                        picks = np.random.choice(new_coop + new_canc, new_gsize - gsize) #else choose randomly among cells that want to divide until max size reached
                        new_canc=0
                        for pick in picks:
                            if pick >= new_coop: 
                                new_canc += 1 
                    group_size[i,j] = new_gsize
                    num_cancer[i,j] += new_canc
            
                    # death
                    if np.random.random() > (1-pd)*(1-(ncanc+new_canc)/new_gsize):
                        group_size[i,j] = 0
                        mode[i,j] = 0
                        num_cancer[i,j] = 0

                    # fragmenters 
                    if group_size[i,j] == msize:
                        locations_of_parents[i,j] = 1
   
        # iterate through the grid and make groups fragment -----------------------|
        for i in range(M):
            for j in range(M):
    
                # if there is a parent at that spot
                if locations_of_parents[i,j] == 1:
   
                   #info 
                   mij = mode[i,j]
                   msize = max_sizes[mij-1]
                   ncanc = num_cancer[i,j] 
                   offs = modes[mij-1] #offspring sizes
                   
                   #remove parent so that offspring don't compete with it
                   group_size[i,j] = 0 
                   mode[i,j] = 0
                   num_cancer[i,j] = 0
    
                   # randomly distribute cancer cells among offspring
                   cells = [1]*ncanc + [0]*(msize-ncanc) #list of cells, 1=cancer
                   np.random.shuffle(cells) #shuffle them up
                   cs = np.hstack([[0], np.cumsum(offs)]) #indices for offspring
                   ncanc_in_offs = [sum(cells[cs[i]:cs[i+1]]) for i in range(len(cs[:-1]))] #number of cancer cells in each offspring 
                   
                   # disperse offspring 
                   for o,off in enumerate(offs):
                           
                       # global dispersal        
                       if globall == 1:
                           ival = np.random.randint(M) 
                           jval = np.random.randint(M) 
    
                       # local dispersal 
                       else: 
                           num = np.random.randint(5)  
                           if num == 0:  # place to the right of parent
                               ival = i + 1
                               jval = j
                           elif num == 1:  # place to the left of parent
                               ival = i - 1
                               jval = j
                           elif num == 2:  # place above parent
                               ival = i
                               jval = j + 1
                           elif num == 3:  # place below parent
                               ival = i
                               jval = j - 1
                           elif num == 4:  # place on the parent spot
                               ival = i
                               jval = j

                       if 0<=ival and ival<M and 0<=jval and jval<M: #if offspring in bounds 

                           # compete for spot   
                           win = False

                           # random competition
                           if s == 0:
                               if group_size[ival,jval] > 0:
                                   if np.random.random() > 0.5: #50/50 chance of replacing resident
                                      win=True

                           # size dependent competition
                           else:                 
                               if group_size[ival,jval] == off: 
                                   if np.random.random() > 0.5: #50/50 chance of replacing when same size
                                       win=True 
                               elif group_size[ival,jval] < off:
                                   win=True

                           if win:
                               group_size[ival,jval] = off
                               mode[ival,jval] = mij
                               num_cancer[ival,jval] = ncanc_in_offs[o]

                   # remove parent
                   locations_of_parents[i,j] = 0

        if np.sum(group_size) == 0:
            break

    return group_size, mode, num_cancer

