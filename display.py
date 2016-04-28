# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:26:15 2016

@author: Aaron
"""

import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import kepler as kp

class AnimatedPaths(object):
    
    def __init__(self, file_name, step_size=1, legend=True, elements=False, select=None):
        self.file_name = file_name
        self.step = 0
        self.step_size = step_size
        self.legend = legend
        self.elements = elements

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        self.h5f = h5.File(self.file_name, mode='r')
        if select:
            self.data = self.h5f["data"][::self.step_size, select]
            self.names = self.h5f["name"][select]
        else:
            self.data = self.h5f["data"][::self.step_size]
            self.names = self.h5f["name"][:]
            
        self.suns = np.array([i[0] == 's' for i in self.names])
        self.planets = np.array([i[0] == 'p' for i in self.names])
        self.kbos = np.array([i[0] == 'a' for i in self.names])
        
        self.nframes = np.size(self.data, 0) - 2
   
    def is_finished(self):
        while (self.step == 0) or (self.step < np.size(self.data, 0)):
            yield self.step

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        self.step = 0
        
        self.colors = np.zeros((np.size(self.data, 1), 4))
        self.colors[:, 3] = 1
        self.colors[self.suns, :3] = [1.0, 0.0, 0.0]
        self.colors[self.planets, :3] = [0.0, 0.5, 0.0]
        self.colors[self.kbos, :4] = [0.0, 0.75, 0.75, 0.3]
#        self.colors = plt.cm.jet(np.linspace(0, 1, np.size(self.data, 1)))
        
        min_x, max_x = np.min(self.data[:, :, 0]), np.max(self.data[:, :, 0])
        min_y, max_y = np.min(self.data[:, :, 1]), np.max(self.data[:, :, 1])
        min_z, max_z = np.min(self.data[:, :, 2]), np.max(self.data[:, :, 2])
        min_x, max_x = (-1000, 1000)
        min_y, max_y = (-1000, 1000)
        min_z, max_z = (-1000, 1000)
        if min_x != max_x:
            self.ax.set_xlim3d(min_x, max_x)
        if min_y != max_y:
            self.ax.set_ylim3d(min_y, max_y)   
        if min_z != max_z:
            self.ax.set_zlim3d(min_z, max_z)

        plt.tight_layout()
        
        if self.elements:
            self.lines = sum([self.ax.plot([], [], [], '-', c=c) for c in self.colors], [])
            
        self.pts = sum([self.ax.plot([], [], [], 'o', c=c, label=name) for (c, name) in zip(self.colors, self.names)], [])
        
        if self.legend:
            self.ax.legend(loc='best', numpoints=1, prop={'size':12})        
        
        if self.elements:
            return self.pts, self.lines
        else:
            return self.pts, 

    def update(self, i):
        """Update the scatter plot."""

        self.step += 1
        
        if self.step % int(self.nframes / 10) == 0:
            print self.step, self.nframes
        
        for j in xrange(np.size(self.data, 1)):
            if self.elements:
                r = self.data[self.step, j, :3]
                v = self.data[self.step, j, 3:]
                t, tdot = kp.orbital_trace(r, v)
                
                self.lines[j].set_data(t[0, :], t[1, :])
                self.lines[j].set_3d_properties(t[2, :])
                
#                self.lines[j].set_data(self.data[:self.step, j, 0], self.data[:self.step, j, 1])
#                self.lines[j].set_3d_properties(self.data[:self.step, j, 2])
            
            self.pts[j].set_data(self.data[self.step, j, 0], self.data[self.step, j, 1])
            self.pts[j].set_3d_properties(self.data[self.step, j, 2])

        if self.elements:
            return self.pts, self.lines
        else:
            return self.pts, 
        

def show_trajectories(fname, select=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    h5f = h5.File(fname, mode='r')
    if select:
        data = h5f["data"][:, select]
        names = h5f["name"][select]
    else:
        data = h5f["data"][:]
        names = h5f["name"][:]

    colors = plt.cm.jet(np.linspace(0, 1, np.size(data, 1)))
    
    min_x, max_x = np.min(data[:, :, 0]), np.max(data[:, :, 0])
    min_y, max_y = np.min(data[:, :, 1]), np.max(data[:, :, 1])
    min_z, max_z = np.min(data[:, :, 2]), np.max(data[:, :, 2])
    
    for i in xrange(np.size(data, 1)):
        ax.plot(data[:, i, 0], data[:, i, 1], data[:, i, 2], c=colors[i], label=names[i])
    
    if min_x != max_x:
        ax.set_xlim3d(min_x, max_x)
    if min_y != max_y:
        ax.set_ylim3d(min_y, max_y)   
    if min_z != max_z:
        ax.set_zlim3d(min_z, max_z)
    
    ax.legend()
    plt.tight_layout()
    h5f.close()
    
def show_trajectories2d(fname, select=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    h5f = h5.File(fname, mode='r')
    if not select:
        particle_names = h5f.keys()            
    else:
        particle_names = select
    nparticles = len(particle_names)
    
    for i, pname in enumerate(particle_names):
        ax.plot(h5f[pname][:, 0], h5f[pname][:, 1], label=pname)
    
    max_x = max([np.max(h5f[p][:, 0]) for p in particle_names])
    min_x = min([np.min(h5f[p][:, 0]) for p in particle_names])
    ax.set_xlim(min_x, max_x)
    max_y = max([np.max(h5f[p][:, 1]) for p in particle_names])
    min_y = min([np.min(h5f[p][:, 1]) for p in particle_names])
    ax.set_ylim(min_y, max_y)
    
    ax.legend()
    
    h5f.close()

if __name__ == '__main__':
#    fname = r"H:\nbody\build\out.h5"
#    fname = r"H:\nbody\build\solar.h5"
    fname = r"H:\nbody\build\perturbed.h5"
    
    a = AnimatedPaths(fname, step_size=5, legend=False, elements=True)
    ani = animation.FuncAnimation(a.fig, a.update, frames=a.nframes, interval=20, 
                                  init_func=a.setup_plot, blit=False, repeat=False)
    ani.save(r'C:\Users\Aaron\Documents\perturbed.mp4', extra_args=['-vcodec', 'libx264'])
    plt.show()
    
    a.h5f.close()
    
#    show_trajectories(fname)
#    plt.show()
