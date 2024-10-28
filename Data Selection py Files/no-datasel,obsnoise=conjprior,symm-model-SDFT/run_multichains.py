try:  # on server
    import mkl
    mkl.set_num_threads(1)  # this is the number of threads per chain 
except:
    # print('Number of threads not set.')
    pass

import multiprocessing as mp
import run_fitting as rf

def run(fullfile, file):
    num_chains = 1

    workdata = [] #create list to append work jobs to
    for i in range(num_chains):
        workdata.append([i+1, fullfile, file])
    workdata = tuple(workdata)

    print('Current file:', file)
    with mp.Pool(num_chains) as p:
        p.starmap(rf.run, workdata)

