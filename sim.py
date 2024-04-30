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

def sim(modes=[(1,1)], T=1000, M=10, N0=100, bs=[0,1], d=0, r=2, s=0, m=0, w=0, local=False):

    # modes   
    nmodes = len(modes) #number of modes
    max_sizes = [sum(mode) for mode in modes] #max sizes
 
    # birth rates
    pbs = [1 - np.exp(-b) for b in bs] #probability a cooperator cell divides for given number of cooperators
    pbs_cancer = [1 - np.exp(-r * b) for b in bs] #probability a cancer cell divides for given number of cooperators

    # death rate
    pd = 1 - np.exp(-d) #probability a group dies

    # lattice
    group_size = np.zeros((M, M), int)  #number of cells at vertices
    mode = np.zeros((M, M), int) # fragmentation genotype for points on the lattice that are inhabited by cells
    num_cancer = np.zeros((M, M), int) #number of cancer cells
    parents = np.zeros((M,M), int) #temporary record of which vertices fragmenting

    # create N0 randomly scattered individuals on an M by M grid
    ixs = np.random.choice(M**2, size=N0, replace=False) #N0 unique random numbers
    locations = np.array(list(product(range(M),range(M))))[ixs] #convert to unique vertices
    for i,j in locations:
        group_size[i,j] = 1 #start all individuals as single cells
        mode[i,j] = np.random.randint(1,nmodes+1) #choose mode at random (note we are saving 0 for empty vertex)
        num_cancer[i,j] = 0
    
    # ITERATE
    for t in range(T):

        # cell division and group death 
        for i in range(M):
            for j in range(M):
                
                gsize = group_size[i,j] #number of cells at this vertex
    
                if gsize > 0:  # if there is a group at that spot
    
                    # info
                    ncanc = num_cancer[i,j] #number of cancer cells
                    ncoop = gsize - ncanc #number of cooperator cells
                    msize = max_sizes[mode[i,j]-1] #max size
 
                    # births
                    coop_births = np.random.binomial(ncoop, pbs[ncoop])  #number of cooperator cells that want to divide
                    canc_births = np.random.binomial(ncanc, pbs_cancer[ncoop]) #number of cancer cells that want to divide

                    #mutations
                    coop_muts = np.random.binomial(coop_births, m) #number of new cooperators that switch to cancer cells
                    canc_muts = np.random.binomial(canc_births, w) #number of new cancer cells that switch to cooperators
                    new_coop = coop_births - coop_muts + canc_muts #number of new cooperators
                    new_canc = canc_births - canc_muts + coop_muts #number of new cancer cells

                    # make sure we dont overshoot max size
                    new_gsize = gsize + new_coop + new_canc #size if all new cells added
                    if new_gsize > msize:  # if the max size is exceeded
                        new_gsize = msize #reach max size
                        picks = np.random.choice(new_coop + new_canc, new_gsize - gsize) #choose randomly among cells that want to divide until max size reached
                        new_canc = 0 #initialize tally of new cancer cells to be added
                        for pick in picks:
                            if pick >= new_coop: #if chosen index greater than or equal to number of cooperators
                                new_canc += 1 #then add a new cancer cell
                    group_size[i,j] = new_gsize #update group size
                    num_cancer[i,j] += new_canc #add new cancer cells
            
                    # death
                    if np.random.random() < 1 - (1-pd)*(1-(ncanc+new_canc)/new_gsize): #if dont survive both forms of mortality
                        group_size[i,j] = 0 #remove
                        mode[i,j] = 0
                        num_cancer[i,j] = 0

                    # fragmenters 
                    if group_size[i,j] == msize:
                        parents[i,j] = 1
   
        # iterate through the grid and make groups fragment
        for i in range(M):
            for j in range(M):
    
                # if there is a parent at that spot
                if parents[i,j] == 1:
   
                   #info 
                   mij = mode[i,j] #mode index
                   msize = max_sizes[mij-1] #max size
                   ncanc = num_cancer[i,j]  #number of cancer cells
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
                       if local == False or local == "False":
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
                           win = True
                           if group_size[ival,jval] > 0:
                               num = np.random.random()
                               if group_size[ival,jval] < off:
                                   if num < 1 - (1+s)/2: #lose to smaller resident with prob 1 - (1+s)/2 
                                       win=False
                               elif group_size[ival,jval] == off:
                                   if num < 0.5: #50/50 chance of losing when same size as resident
                                       win=False
                               elif group_size[ival,jval] > off:
                                   if num < (1+s)/2: #lose to larger resident with prob (1+s)/2 
                                       win=False
      
                           #replace resident
                           if win:
                               group_size[ival,jval] = off
                               mode[ival,jval] = mij
                               num_cancer[ival,jval] = ncanc_in_offs[o]

                   # remove parent
                   parents[i,j] = 0
     
        # end simulation if extinct
        if np.sum(group_size) == 0:
            break 

    return group_size, mode, num_cancer
