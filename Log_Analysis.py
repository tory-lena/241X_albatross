#This file reading data

import numpy as np

log = np.genfromtxt('log001.csv', delimiter=",",filling_values = 0.0, skip_header=5)
print(log.shape)

"""
skip_header=5,
"""