#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f=np.load("energy_regen_prophist_matter.npz")
print f['_num_vec'].shape
plt.plot(f['_dists_'],f['_num_vec'][:,2])
plt.show()

