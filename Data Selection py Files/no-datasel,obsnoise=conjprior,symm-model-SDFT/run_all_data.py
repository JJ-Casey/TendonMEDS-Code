import numpy as np
import os
import run_multichains as rm

# script to automaticall trim raw data to inlcude stress-strain data up to maximum stress where damage is assumed to have initiated

np.seterr(all='ignore')
dir = os.path.join(os.getcwd(), 'trimmed_data')  # filedialog.askdirectory()
for file in os.listdir(dir):
    filename = os.fsdecode(file)
    # if filename.endswith('.txt'):
    datafile = os.path.join(dir, filename)

    if __name__ == '__main__':  # this is only executed when this code is run as a script (rather than imported as a module)
        rm.run(datafile, filename)
