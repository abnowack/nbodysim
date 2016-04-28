# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 02:51:57 2016

@author: Aaron
"""
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

import kepler as kp

def plot_perihelion_semimajor(fname, title=''):
    f = h5.File(fname, 'r')
    
    data = f["data"]
    names = f["name"]
    kbos = np.array([i[0] == 'a' for i in names])
    
    kbo_data = data[:, kbos, :]
    
    nsteps = np.size(kbo_data, 0)
    
    steps = (0, int(nsteps / 2), nsteps-1)
    colors = ['red', 'green', 'blue']
    titles = ['Initial', 'Middle', 'End']
    
    plt.figure()
    
    for q, step in enumerate(steps):
        a_ = np.zeros(np.size(kbo_data))
        w_ = np.zeros(np.size(kbo_data))
        
        for ijk in xrange(np.size(kbo_data, 1)):
            r = kbo_data[step, ijk, :3]
            rdot = kbo_data[step, ijk, 3:]
            (a, e, i, O, w, M) = kp.cart_to_kepler(r, rdot)
            
            a_[ijk] = a
            w_[ijk] = w
        
        plt.plot(a_, w_, 'o', c=colors[q], label=titles[q])
        
    plt.xlim(0, 600)
    plt.ylim(0., 360)
    plt.xlabel('Semi Major Axis')
    plt.ylabel('Argument of Perihelion')
    plt.legend()
    plt.title(title)
    plt.tight_layout()
        
    f.close()
        
if __name__ == '__main__':
    fname = r"H:\nbody\build\out.h5"
    
#    plot_perihelion_semimajor(r"H:\nbody\build\out.h5", "Solar System Dynamics, 10k Years")
#    plot_perihelion_semimajor(r"H:\nbody\build\unperturbed.h5", "KBO Dynamics Unperturbed, 10k Years")
    plot_perihelion_semimajor(r"H:\nbody\build\perturbed.h5", "KBO Dynamics Perturbed (10 Earth Mass), 10k Years")