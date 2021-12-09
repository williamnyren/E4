# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

file = 'gaussian_results.dat';

# Read normal gaussian from results.dat
data = np.genfromtxt(file, delimiter=',')/300000;

n_bins = 100;


fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True);

# We can set the number of bins with the *bins* keyword argument.
axs.hist(data, bins=n_bins);


# %%
