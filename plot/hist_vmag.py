#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab

plt.xlabel('Magnitude (V)')  
plt.ylabel('Number')  

x = np.genfromtxt('vmag.csv', dtype=float)

n, bins, patches = plt.hist(x, 20, facecolor='green', alpha=0.75)
#n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

plt.savefig('vmag.png', dpi=75)
plt.show()

