# -*- coding: utf-8 -*-
"""
Spyder Editor

Quickly plot the output RZ coordinates
Author: Laura Xin Zhang
"""
import numpy as np
import matplotlib.pyplot as plt

dirin = './output/'
prefix = 'example'

#
coordRZ = np.loadtxt(dirin + prefix + '_coords_RZ.txt', delimiter = ',')
energy = np.loadtxt(dirin + prefix + '_energy.txt')
time = np.loadtxt(dirin + prefix + '_time.txt')
print(coordRZ.shape)
plt.scatter(coordRZ[::10, 0], coordRZ[::10, 1], s = 0.3, 
            c = time[::10]*1000, cmap = 'plasma', vmax = time[-1]*1000*1.1)
cb = plt.colorbar()
cb.set_label('Time [ms]')
plt.title('Collisionless Orbit')
plt.ylabel('Z')
plt.xlabel('R')
plt.savefig('./figures/'+prefix+'.pdf', bbox_inches = 'tight')

plt.show()